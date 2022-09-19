#!/usr/bin/env python3

# Setup the notebook's environment
import lalsimulation as lalsim
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import xcp
from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict,template_amp_phase,gc,generate_model_params

# Let the user know where lalsimulation lives

# from "/Users/book/KOALA/PhenomXCP/xcp/parameter_space_fits.py" import generate_model_params

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
    
    
#

#
from numpy.linalg import norm
from scipy.optimize import curve_fit

# For all pairs of l and m in the global config file
for ll,mm in gc.lmlist:

    # --------------------------------------- #
    # Preliminaries
    # --------------------------------------- #

    #
    alert('Loading parameter space data: (l,m) == (%i,%i)'%(ll,mm))

    # Define data location
    package_dir = parent( xcp.__path__[0] )
    datadir = package_dir + 'data/version2/'
    files = glob( datadir+'*_l%im%i.txt'%(ll,mm) )
    files.sort()
    # files = files[::-1]#
    
    # Ignore select files
    files = [ f for f in files if ('fit' not in f)  ]

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
        files = [ f for f in files if ('q1' not in f)  ]

    # Load and unpack physical parameter space
    # NOTE that the order of these parameters must be the same as that in files
    opt_parameter_range = loadtxt(datadir+'fit_opt_parameters_l%im%i.txt'%(ll,mm))
    # scarecrow = template_amp_phase(0.5, 0.5,zeros(3),zeros(3),lm=(ll,mm))
    # parameter_names_in_order = scarecrow.__code__.co_varnames[1:scarecrow.__code__.co_argcount]
    # model_range = {  parameter_names_in_order[k]:var for k,var in enumerate(opt_parameter_range.T) }

    #
    fig,ax = subplots( len(files), 2, figsize=3*array([ 2.5*2/(0.618), 1.5*len(files) ]) )
    ax = ax.flatten()

    #
    tight_layout(w_pad=4,h_pad=4)

    #
    foo = {}

    #
    lmlist = [ (ll,mm) ]

    #
    p = 0
    for j,f_ in enumerate(files):

        #
        simname = f_.split('/')[-1].split('_l%im%i.'%(ll,mm))[0]

        # Find index location of metadata for simname 
        k = [ k for k,val in enumerate(metadata_dict['simname']) if val in simname ][0]

        # Load data for this case
        raw_data = loadtxt(f_).T
        alert(simname)
        calibration_data, dphi_lorentzian_min, f_min, f_max, f_lorentzian_min = determine_data_fitting_region( raw_data, simname=simname, lm=(ll,mm) )

        # Collect params for this case 
        metadata = metadata_dict['array_data'][k,:]

        #
        f,amp_fd,dphi_fd,alpha,beta,gamma = calibration_data.T
        theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z = metadata_dict['array_data'][k]
        chi1_vec = array([chi1_x,chi1_y,chi1_z])
        chi2_vec = array([chi2_x,chi2_y,chi2_z])
        
        #
        action_helper = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(ll,mm), option_shorthand='4-xhm')
        mod_xhm0_amp,mod_xhm0_dphi = action_helper(f)
        
        #
        popt = opt_parameter_range[j,:]
        opt_amp,opt_dphi = action_helper(f,*popt)
        
        #
        (mu1_l2m2,mu2_l2m2,mu3_l2m2,nu4_l2m2,nu5_l2m2,nu6_l2m2,zeta1_l2m2,zeta2_l2m2,mu1_l3m3,mu2_l3m3,mu3_l3m3,mu4_l3m3,nu4_l3m3,nu5_l3m3,nu6_l3m3,zeta1_l3m3,zeta2_l3m3) = generate_model_params(theta,eta,a1)
        if ll == 2:
            ppy = (mu1_l2m2,mu2_l2m2,mu3_l2m2,nu4_l2m2,nu5_l2m2,nu6_l2m2,zeta1_l2m2,zeta2_l2m2)
        else:
            ppy = (mu1_l3m3,mu2_l3m3,mu3_l3m3,mu4_l3m3,nu4_l3m3,nu5_l3m3,nu6_l3m3,zeta1_l3m3,zeta2_l3m3)
        py_amp,py_dphi = action_helper(f,*ppy)
        py_dphi -= min(py_dphi)
        
        #
        mod_xhm_dict = xcp.get_phenomxphm_coprecessing_multipoles( f,lmlist, m1, m2, chi1_vec, chi2_vec, option_shorthand='2-xphm' )
        mod_xhm = mod_xhm_dict[ll,mm]
        mod_xhm_amp = abs(mod_xhm)
        mod_xhm_phi = unwrap( angle(mod_xhm) )
        mod_xhm_dphi = spline_diff(f,mod_xhm_phi)
        mod_xhm_dphi -= min( mod_xhm_dphi[ (f>0.03*ll*0.5)&(f<0.12*ll*0.5) ] )
        # mod_xhm_dphi -= mean( mod_xhm_dphi )
        
        #
        tuned_xhm_dict = xcp.get_phenomxphm_coprecessing_multipoles( f,lmlist, m1, m2, chi1_vec, chi2_vec, option_shorthand='1-pnr' )
        tuned_xhm = tuned_xhm_dict[ll,mm]
        tuned_xhm_amp = abs(tuned_xhm)
        tuned_xhm_phi = unwrap( angle(tuned_xhm) )
        tuned_xhm_dphi = spline_diff(f,tuned_xhm_phi)
        tuned_xhm_dphi -= min( tuned_xhm_dphi[ (f>0.03*ll*0.5)&(f<0.12*ll*0.5) ] )
        # tuned_xhm_dphi -= mean( tuned_xhm_dphi )

        # PLOTTING
        # ---
        
        
        # dphi_fd -= mean(dphi_fd)
        # opt_dphi -= mean(opt_dphi)
        # tuned_xhm_dphi -= mean(tuned_xhm_dphi)
        # mod_xhm0_dphi -= mean(mod_xhm0_dphi)

        #
        sca(ax[p]); p+=1
        plot( f, dphi_fd, label='Calibration Data (NR)', lw=4, alpha=0.2, color='k' )
        plot( f, opt_dphi, label='Direct Fit', ls='--',lw=2,alpha=1,color='dodgerblue' )
        plot( f, tuned_xhm_dphi, label='End Model (PNR)', ls='-',lw=2,alpha=1,color='r' )
        plot( f, py_dphi, label='End Model (PNR-py)', ls=':',lw=2,alpha=1,color='tab:orange' )
        plot( f, mod_xhm_dphi, label='PhenomXPHM', ls='-',lw=2,alpha=0.85,color='m' )
        plot( f, mod_xhm0_dphi, label='PhenomXHM', ls='--',lw=1,alpha=1,color='k',zorder=-10 )
        xscale('log')
        xlim(lim(f,dilate=1.1,dilate_with_multiply=True))
        ylim( limy(f, mod_xhm_dphi,dilate=0.1) )
        title(simname,size=12,loc='left')
        legend(ncol=2,loc=1)
        ylabel(r'$\frac{d}{df}\arg(\tilde{h}_{%i%i})$'%(ll,mm))
        xlabel('$fM$')
        title(simname,loc='left',size=12)

        sca(ax[p]); p+=1
        plot( f, amp_fd, label='Calibration Data (NR)', lw=4, alpha=0.2, color='k' )
        plot( f, opt_amp, label='Direct Fit', ls='--',lw=2,alpha=1,color='dodgerblue' )
        plot( f, tuned_xhm_amp, label='End Model (PNR)', ls='-',lw=2,alpha=1,color='r' )
        plot( f, py_amp, label='End Model (PNR-py)', ls=':',lw=2,alpha=1,color='tab:orange' )
        plot( f, mod_xhm_amp, label='PhenomXPHM', ls='-',lw=2,alpha=0.85,color='m' )
        plot( f, mod_xhm0_amp, label='PhenomXHM', ls='--',lw=1,alpha=1,color='k',zorder=-10 )
        yscale('log')
        xscale('log')
        legend(ncol=2)
        ylim( limy(f, amp_fd,dilate=1.2) )
        xlabel('$fM$')
        ylabel(r'$|\tilde{h}_{%i%i}(f)|$'%(ll,mm))
        #
        title(simname,loc='left',size=12)
            
            
    #
    file_path = datadir+'waveform_postfit_diagnostic_l%im%i.pdf'%(ll,ll)
    alert('Saving batch plot to %s'%magenta(file_path))
    savefig(file_path,pad_inches=2, bbox_inches = "tight")
    alert('Done.')