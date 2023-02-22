#!/usr/bin/env python3

# Setup the environment
import lalsimulation as lalsim
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import xcp
from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict,template_amp_phase, gc
from numpy.linalg import norm
from scipy.optimize import curve_fit,minimize,fmin

'''
NOTE(s) ...
'''


# Let the user know where lalsimulation lives

#
lalsim_path = lalsim.__path__[0]
lalsuite_repo_path = lalsim_path.split('lib')[0]+'src/lalsuite/'
branch_name = bash('cd %s && git status'%lalsuite_repo_path).decode("utf-8").split('On branch ')[-1].split('\n')[0]

#
alert('We are getting our LALSimulation from:\n%s'%magenta(lalsim_path))
alert('We think that the related lalsuite source files are here:\n%s'%green(lalsuite_repo_path))
alert('Lastly, we are currently on this branch: %s'%bold(magenta(branch_name)))

#
datadir = '/Users/book/KOALA/PhenomXCP/data/version4/'
data_path = datadir+'calibration_data_dict.pickle'
calibration_data_dict = pickle.load( open( data_path, "rb" ) )

# For all pairs of l and m in the global config file
# for ll,mm in [(2,2)]:#gc.lmlist:
# for ll,mm in gc.lmlist:
for ll,mm in gc.lmlist:
        
        
    #
    files = glob( datadir+'*_l%im%i.txt'%(ll,mm) )
    files.sort()
    files = files[::-1]

    #
    file_map = []
    for j,sn_ in enumerate(metadata_dict['simname']):
        
        # Find index location of this metadata in the list of file paths
        k = [ k for k,val in enumerate(files) if sn_ in val ][0]

        # 
        file_map.append(k)
        

    # Sort file paths like matadata so that we can then sort both according to distance between physical parameters
    files = array(files)[file_map]
        
    reference_index = 0
    # theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z, M_final, chi_final
    mask = [0,3,9] # theta, eta, a1 -- the parameters that will be used to calculate distance for sorting
    coordinates = metadata_dict['array_data'][:,mask]
    coordinates[:,0] = cos(coordinates[:,0]) # use cos theta 
    index_map, sorted_coordinates = distance_sort( coordinates, reference_index, center=not False )

    #
    files = list(files[index_map])
    
    #
    if ll==3:
        files = [ f for f in files if 'q1' not in f  ]

    # #
    # files = [ f for f in files if 'q8a08t60' in f ]

    # Get number of parameters to be tuned
    scarecrow = template_amp_phase(0.5, 0.5,zeros(3),zeros(3),lm=(ll,mm))
    var_count = scarecrow.__code__.co_argcount - 1

    #
    fig,ax = subplots( len(files), 2, figsize=3*array([ 2.5*2/(0.618), 2.0*len(files) ]) )
    ax = ax.flatten()

    #
    tight_layout(pad=1,w_pad=4,h_pad=4)

    #
    foo = {}

    #
    lmlist = [ (ll,mm) ]

    #
    p = 0
    popt_array  = zeros( (len(files),var_count+1) ) # NOTE that the "+1" is for nu0; see code below
    fit_obj_list = []
    physical_param_array = zeros( (len(files), 19) )
    alert('Plotting ...')
    for j,f_ in enumerate(files):

        #
        simname = f_.split('/')[-1].split('_l%im%i.'%(ll,mm))[0]
        
        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* #
        # Load data that has been processed in issue 3c
        # ---
        (metadata,f,dphi_fd,amp_fd,xphm_dphi,dphi_fd_enforced_min,nr_dphi_lm_shift,min_xphm_dphi_l2m2) = calibration_data_dict[ll,mm][simname]
        #
        dphi_fd_floored = dphi_fd - dphi_fd_enforced_min
        # *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-* #

        # Collect params for this case 
        theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z,Mf,Xf = metadata
        chi1_vec = array([chi1_x,chi1_y,chi1_z])
        chi2_vec = array([chi2_x,chi2_y,chi2_z])
        
        #
        physical_param_array[j,:] = metadata
        
        # Generate PhenomXHM waveform generators that allow for model deviations
        # NOTE that the phase derivatives output have their min values set to zero
        action_helper_1 = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(ll,mm),option_shorthand='5-xhm-tuning',include_nu0=False)
        # Create version that includes nu0 (NOTE that varcount must be increased by 1 for this object)
        action_helper_2 = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(ll,mm),option_shorthand='5-xhm-tuning',include_nu0=True)
        # Create version that does not set min dphi to zero
        action_helper_3 = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(ll,mm),option_shorthand='5-xhm-tuning',include_nu0=True,floor_dphi=False)
        
        # ######################################################### #
        # NOTE that the waveform is tuned relative to XHM but the relative locations are tuned relative to XPHM (see issue 3c script)
        # ######################################################### #
        
        # GENERATE TEMPLATE FUNCTIONS
        # ---
        def action( p, verbose=False, output_vars=False, modeling_phase = 1, amp_only=False ):
            
            #
            amplitude,phase_derivative = action_helper_1( f, *p ) if modeling_phase==1 else action_helper_3( f, *p )
            
            #
            calibration_dphi = dphi_fd_floored if modeling_phase==1 else dphi_fd
            
            #
            calibration_dphi -= calibration_dphi[0]
            phase_derivative -= phase_derivative[0]
            
            # -- Calculate residual of phase derivative -- #
            # phase_derivative -= min(phase_derivative)
            # NOTE that we do not shift the phase derivative as suggested above because we have enforced that  min(dphi_fd) == xhm_min_phase_derivative
            residual_phase_derivative = sum((calibration_dphi - phase_derivative)**2) / std(dphi_fd_floored)
            # -- Calculate residual of amplitude --------- #
            amp_scale = f ** (7.0/6)
            inv_amp_scale = f ** (-7.0/6)
            log_scaled_amp_fd = log( amp_fd * amp_scale )
            log_scaled_amplitude = log( amplitude * amp_scale )
            residual_amplitude = sum((log_scaled_amp_fd - log_scaled_amplitude)**2) / std(log_scaled_amp_fd)
            # -- Combine residuals ----------------------- #
            '''
            I am very happy to NOTE that multiplying the amplitude residual by 100 has the
            effect of forcing the fit to always prioritize the amplitude fit, to the great
            effect of overcoming the sometimes troublesome data quality issues in the phase derivative.
            '''
            combined_residual = (1-amp_only)*residual_phase_derivative + 100*residual_amplitude
            
            #
            if output_vars:
                return (combined_residual,p)
            return combined_residual
        
        # PERFORM FIT
        # --- 
        
        # Perform fit -- PHASE 1: Do not include dphi offset parameter nu0
        # # NOTE that carrying forward previous solutions as initial guesses causes problems UNLESS the simulations are distance sorted as can be seen above
        # --- parameters for dphi shape and amp, not dphi shift
        guess_1 = zeros(var_count) #if j==0 else foo[1]
        foo_shape = minimize( lambda p: action(p,modeling_phase=1), guess_1 )
        foo_shape = (foo_shape.fun,foo_shape.x)
        
        #
        best_shape_popt = list(foo_shape[1]) + [0]
        best_shape_fit_amp,best_shape_fit_dphi = action_helper_3( f, *best_shape_popt )
        
        #
        if ll == 2:
            dev_parameter = a1 * sin(theta) 
        if ll == 3:
            dev_parameter = a1 * sin(theta) * delta 
        
        final_dphi_shidt_opt = (dphi_fd_enforced_min-min(best_shape_fit_dphi))# if ll==2 else mean(dphi_fd-best_shape_fit_dphi)
        nu0_opt_guess = final_dphi_shidt_opt / dev_parameter
        
        # # Perform fit -- PHASE 2: Optimize over nu0 ONLY
        # # --- parameter for dphi shift
        # guess_2 = nu0_opt_guess #if j==0 else foo[1]
        # foo_shift = minimize( lambda p: action(list(guess_1) + [ p ],modeling_phase=2), guess_2 )
        # foo_shift = (foo_shift.fun,foo_shift.x)
        #
        foo = list(foo_shape)
        nu0_opt = nu0_opt_guess # foo_shift[1][0]
        foo[1] = list(foo_shape[1]) + list([nu0_opt])
        
        
        # foos, boundary_par = jac_sort_minimize(action,var_count, verbose=True, initial_guess=guess)
        # foo = foos[ -1 ]
        # foo = foos[ argmin( [ f[0] for f in foos ] ) ]
        
        #
        alert('Calling best fit waveform (%f) '%final_dphi_shidt_opt)
        best_fit_amp,best_fit_dphi = action_helper_3( f, *foo[1] )
        alert('END of calling best fit waveform ')
        # best_fit_dphi += nu0_opt
        
        # Store fit params and cov 
        alert(simname,header=True)
        # print('~> ',theta,delta,a1)
        print('>> ',simname,nu0_opt,nu0_opt_guess)
        alert(foo[1],header=False)
        popt_array[j,:] = foo[1]
        fit_obj_list.append( foo )
        
        # PLOTTING
        # ---
        
        # Calculate default model amp and phase derivative
        xhm_amp,xhm_dphi = action_helper_3(f)
        
        #
        sca(ax[p]); p+=1
        plot( f, dphi_fd-dphi_fd[0], label='Calibration Data (NR)', lw=4,ls='-', alpha=0.15, color='k' )
        
        plot( f, xphm_dphi-xphm_dphi[0], label='PhenomXPHM', ls=':', lw=1, alpha=0.85, color='tab:blue', zorder=-10 )
        
        plot( f, xhm_dphi-xhm_dphi[0], label='PhenomXHM', ls='--',lw=1,alpha=0.85,color='k',zorder=-10 )
        
        plot( f, best_fit_dphi-best_fit_dphi[0], label='Best Fit', color='r', ls='-' )
        xscale('log')
        xlim(lim(f,dilate=1.1,dilate_with_multiply=True))
        ylim( limy(f, dphi_fd-dphi_fd[0],dilate=0.1) )
        title(simname,size=12,loc='left')
        legend(ncol=2,loc=1)
        ylabel(r'$\frac{d}{df}\arg(\tilde{h}_{%i%i})$'%(ll,mm))
        xlabel('$fM$')
        title(simname,loc='left',size=12)
        
        axhline(0-dphi_fd[0],ls=':',color='k')
        axhline(min_xphm_dphi_l2m2-dphi_fd[0],ls='--',color='green')
        
        #
        sca(ax[p]); p+=1
        plot( f, amp_fd, label='Calibration Data (NR)', lw=4,ls='-', alpha=0.15, color='k' )
        plot( f, xhm_amp, label='PhenomXHM', ls='--',lw=1,alpha=0.85,color='k',zorder=-10 )
        plot( f, best_fit_amp, label='Best Fit', color='r', ls='-' )
        yscale('log')
        xscale('log')
        legend(ncol=2,loc=3)
        ylim( limy(f, amp_fd,dilate=1.2) )
        xlabel('$fM$')
        ylabel(r'$|\tilde{h}_{%i%i}(f)|$'%(ll,mm))
        title(simname,loc='left',size=12)
        
        #
        print( '.',end='')
            
    #
    print( '')
    file_path = datadir+'waveform_fit_diagnostic_l%im%i.pdf'%(ll,mm)
    alert('Saving batch plot to %s'%magenta(file_path))
    savefig(file_path,pad_inches=2, bbox_inches = "tight")
    alert('Done.')

    # SAVE FIT DATA
    # --

    # Initial binary parameters
    data_path = datadir+'fit_initial_binary_parameters_l%im%i.txt'%(ll,mm)
    alert('Saving %s to %s'%( magenta('physical_param_array'), magenta(data_path)) )
    savetxt( data_path, physical_param_array, header='see "issues/3a_collect_metadata.py"; columns are theta, m1, m2, eta, delta, chi_eff, chi_p, chi1, chi2, a1, a2, chi1_L_x, chi1_L_y, chi1_L_z, chi2_L_x, chi2_L_y, chi2_L_z' )

    # Get parameter names in order to use for the next file header
    scarecrow = template_amp_phase(0.5, 0.5,zeros(3),zeros(3),lm=(ll,mm))
    parameter_names_in_order = scarecrow.__code__.co_varnames[1:scarecrow.__code__.co_argcount]

    # Fit parameters
    data_path = datadir+'fit_opt_parameters_l%im%i.txt'%(ll,mm)
    alert('Saving %s to %s'%( magenta('dphi_popt_array'), magenta(data_path)) )
    savetxt( data_path, popt_array, header='see "template_together()" in core.py; columns are %s'%(' '.join(parameter_names_in_order)) )

    #
    data_path = datadir+'fit_objects_l%im%i.pickle'%(ll,mm)
    alert('Saving fit_obj_list to %s'%magenta(data_path))
    pickle.dump( fit_obj_list, open( data_path, "wb" ) )

#
alert('Fitting complete.')
alert('Plotting complete.')
alert('Saving complete.')