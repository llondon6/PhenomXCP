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
if branch_name != 'pnrv1-ll':
    alert('We are not on the expected branch. This may cause unexpected behavior.',say=not True)


# For all pairs of l and m in the global config file
for ll,mm in gc.lmlist:

    #
    datadir = '/Users/book/KOALA/PhenomXCP/data/version2/'
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
    # theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z
    mask = [0,3,9] # theta, eta, a1 -- the parameters that will be used to calculate distance for sorting
    coordinates = metadata_dict['array_data'][:,mask]
    coordinates[:,0] = cos(coordinates[:,0]) # use cos theta 
    index_map, sorted_coordinates = distance_sort( coordinates, reference_index, center=not False )

    #
    files = list(files[index_map])
    
    #
    if ll==3:
        files = [ f for f in files if 'q1' not in f  ]
        
        
    reverse_index_map = argsort( index_map ).astype(int) # array([ int(c) for c in argsort( index_map ) ]),dtype=int

    # Get number of parameters to be tuned
    scarecrow = template_amp_phase(0.5, 0.5,zeros(3),zeros(3),ell=2)
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
    popt_array  = zeros( (len(files),var_count) )
    fit_obj_list = []
    physical_param_array = zeros( (len(files), 17) )
    alert('Plotting ...')
    for j,f_ in enumerate(files):

        #
        simname = f_.split('/')[-1].split('_l%im%i.'%(ll,mm))[0]

        # Find index location of metadata for simname 
        k = [ k for k,val in enumerate(metadata_dict['simname']) if val in simname ][0]

        # Load data for this case
        raw_data = loadtxt(f_).T
        calibration_data, dphi_lorentzian_min, f_min, f_max, f_lorentzian_min = determine_data_fitting_region( raw_data, simname=simname, lm=(ll,mm) )

        # Collect params for this case 
        metadata = metadata_dict['array_data'][k,:]

        #
        f,amp_fd,dphi_fd,alpha,beta,gamma = calibration_data.T
        theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z = metadata_dict['array_data'][k]
        
        # dphi_fd -= mean(dphi_fd)
            
        # # Compute the phase from the phase derivative
        # phi_fd = spline_antidiff( f, dphi_fd )
        # # Compute complex strain
        # h_fd = amp_fd * exp( 1j * phi_fd )
        # amp_scale = f ** (7.0/6)
        # std_h_fd_2 = std(h_fd*amp_scale)**2
        
        #
        chi1_vec = array([chi1_x,chi1_y,chi1_z])
        chi2_vec = array([chi2_x,chi2_y,chi2_z])
        
        #
        physical_param_array[j,:] = metadata_dict['array_data'][k]
        
        # GENERATE TEMPLATE FUNCTIONS
        # ---
        # Generate PhenomX (not PhenomXP) waveform generators that allow for model deviations
        action_helper = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(ll,mm),option_shorthand='4-xhm')
        def action( p, verbose=False, output_vars=False ):
            
            #
            # print('>> ',p)
            amplitude,phase_derivative = action_helper( f, *p )
            
            # # -- Compute FD waveforms, time-shift align, and then compute residual -- #
            # # Compute optimal time shift 
            # phase_shift = mean( phase_derivative - dphi_fd )
            # # Compute model phase, shifted accordingly 
            # phase = phase_derivative - phase_shift 
            # # Compute complex strain 
            # complex_strain = amplitude * exp( 1j * phase ) 
            # #
            # abs_res = abs(h_fd-complex_strain) * amp_scale
            # combined_residual = sum(abs_res**2 / std_h_fd_2)
            
            # -- Calculate residual of phase derivative -- #
            phase_derivative -= min(phase_derivative)
            residual_phase_derivative = sum((dphi_fd - phase_derivative)**2) / std(dphi_fd)
            # -- Calculate residual of amplitude --------- #
            amp_scale = f ** (7.0/6)
            inv_amp_scale = f ** (-7.0/6)
            log_scaled_amp_fd = log( amp_fd * amp_scale )
            log_scaled_amplitude = log( amplitude * amp_scale )
            residual_amplitude = sum((log_scaled_amp_fd - log_scaled_amplitude)**2) / std(log_scaled_amp_fd)
            # -- Combine residuals ----------------------- #
            combined_residual = residual_phase_derivative + residual_amplitude
            
            #
            if output_vars:
                return (combined_residual,p)
            return combined_residual
        
        # PERFORM FIT
        # --- 
        
        # Calculate default model amp and phase derivative
        mod_xhm0_amp,mod_xhm0_dphi = action_helper(f)
        
        # Perform fit 
        # # NOTE that carrying forward previous solutions as initial guesses causes problems UNLESS the simulations are distance sorted as can be seen above
        guess = zeros(var_count) #if j==0 else foo[1]
        foo = minimize( action, guess )
        foo = (foo.fun,foo.x)
        
        
        # foos, boundary_par = jac_sort_minimize(action,var_count, verbose=True, initial_guess=guess)
        # foo = foos[ -1 ]
        # foo = foos[ argmin( [ f[0] for f in foos ] ) ]
        
        #
        best_fit_amp,best_fit_dphi = action_helper( f, *foo[1] )
        
        # Store fit params and cov 
        alert(simname,header=True)
        alert(foo[1],header=False)
        popt_array[j,:] = foo[1]
        fit_obj_list.append( foo )
        
        # PLOTTING
        # ---
        
        # align by requiring zero min
        # dphi_fd -= min(dphi_fd) # NOTE that this is already done for the calibration data
        best_fit_dphi -= min(best_fit_dphi)
        mod_xhm0_dphi -= min(mod_xhm0_dphi)
        
        #
        sca(ax[p]); p+=1
        plot( f, dphi_fd, label='Calibration Data (NR)', lw=4,ls='-', alpha=0.15, color='k' )
        plot( f, mod_xhm0_dphi, label='PhenomXHM', ls='--',lw=1,alpha=0.85,color='k',zorder=-10 )
        plot( f, best_fit_dphi, label='Best Fit', color='r', ls='-' )
        xscale('log')
        xlim(lim(f,dilate=1.1,dilate_with_multiply=True))
        ylim( limy(f, dphi_fd,dilate=0.1) )
        title(simname,size=12,loc='left')
        legend(ncol=2,loc=1)
        ylabel(r'$\frac{d}{df}\arg(\tilde{h}_{22})$')
        xlabel('$fM$')
        title(simname,loc='left',size=12)
        
        #
        sca(ax[p]); p+=1
        plot( f, amp_fd, label='Calibration Data (NR)', lw=4,ls='-', alpha=0.15, color='k' )
        plot( f, mod_xhm0_amp, label='PhenomXHM', ls='--',lw=1,alpha=0.85,color='k',zorder=-10 )
        plot( f, best_fit_amp, label='Best Fit', color='r', ls='-' )
        yscale('log')
        xscale('log')
        legend(ncol=2,loc=3)
        ylim( limy(f, amp_fd,dilate=1.2) )
        xlabel('$fM$')
        ylabel(r'$|\tilde{h}_{22}(f)|$')
        title(simname,loc='left',size=12)
        
        #
        print( '.',end='')
            
    #
    print( '')
    file_path = datadir+'waveform_fit_diagnostic_l%im%i.pdf'%(ll,mm)
    alert('Saving batch plot to %s'%magenta(file_path))
    savefig(file_path,pad_inches=2, bbox_inches = "tight")
    alert('Done.')

    # #
    # physical_param_array = physical_param_array[ reverse_index_map,: ]
    # popt_array = popt_array[ reverse_index_map,: ]
    # fit_obj_list  = list(array(fit_obj_list)[ reverse_index_map ])

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