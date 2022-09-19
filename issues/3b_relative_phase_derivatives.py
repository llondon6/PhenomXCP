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
from numpy.linalg import norm
from scipy.optimize import curve_fit
from os.path import exists

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
for ll,mm in [ lm for lm in gc.lmlist if lm != (2,2)]:

    # --------------------------------------- #
    # Preliminaries
    # --------------------------------------- #

    #
    alert('Loading parameter space data: (l,m) == (%i,%i)'%(ll,mm))

    # Define data location
    package_dir = parent( xcp.__path__[0] )
    datadir = package_dir + 'data/version4/'
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

    #
    fig,ax = subplots( len(files), 1, figsize=3*array([ 2.5*1/(0.618), 1.5*len(files) ]) )
    ax = ax.flatten()

    #
    tight_layout(w_pad=4,h_pad=4)

    #
    foo = {}

    #
    lmlist = [ (ll,mm) ]

    #
    p = 0
    # c0list = 
    shift_dict = {}
    for j,f_ in enumerate(files):

        #
        simname = f_.split('/')[-1].split('_l%im%i.'%(ll,mm))[0]

        # Find index location of metadata for simname 
        k = [ k for k,val in enumerate(metadata_dict['simname']) if val in simname ][0]
        
        
        # Collect params for this case 
        metadata = metadata_dict['array_data'][k,:]
        theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z,Mf,Xf = metadata
        chi1_vec = array([chi1_x,chi1_y,chi1_z])
        chi2_vec = array([chi2_x,chi2_y,chi2_z])
        
        # ----
        
        #
        # ax_ = ax[p]; p+=1
        
        # Load QNM info
        qnmo_p = qnmobj( Mf, Xf, ll, mm,0,p=1,use_nr_convention=True,verbose=False,calc_slm=False,calc_rlm=False )
        fring  = qnmo_p.CW.real / (2*pi)
        # Load data for this case
        raw_data = loadtxt(f_).T
        # Determine data fitting region
        calibration_data, dphi_lorentzian_min, f_min, f_max, f_lorentzian_min = determine_data_fitting_region( raw_data, 
                                                                                                            fring, 
                                                                                                            lm=(ll,mm), 
                                                                                                            floor_dphi=True, 
                                                                                                            plot=not True, 
                                                                                                            simname=simname)
        
        # Load QNM info
        qnmo_p = qnmobj( Mf, Xf, 2, 2,0,p=1,use_nr_convention=True,verbose=False,calc_slm=False,calc_rlm=False )
        fring_22  = qnmo_p.CW.real / (2*pi)
        # Load data for this case
        raw_data_22 = loadtxt(f_.replace('l3m3','l2m2')).T
        # Determine data fitting region
        calibration_data_22, dphi_lorentzian_min_22, f_min_22, f_max_22, f_lorentzian_min_22 = determine_data_fitting_region( raw_data_22, 
                                                                                                            fring_22, 
                                                                                                            lm=(2,2), 
                                                                                                            floor_dphi=True, 
                                                                                                            plot=not True,
                                                                                                            simname=simname)
        # xlim(f_min_22,f_max)
        

        #
        f_22,amp_fd_22,dphi_fd_22,alpha_22,beta_22,gamma_22 = calibration_data_22.T
        f,amp_fd,dphi_fd,alpha,beta,gamma = calibration_data.T
        
        #
        action_helper_22 = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(2,2), option_shorthand='2-xphm',turn_on_relative_dphi_mode=True)
        action_helper = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(ll,mm), option_shorthand='2-xphm',turn_on_relative_dphi_mode=True)
        #
        mod_xhm0_amp_22,mod_xhm0_dphi_22,mod_xhm0_min_dphi_22 = action_helper_22(f)
        mod_xhm0_amp,mod_xhm0_dphi,mod_xhm0_min_dphi          = action_helper(f)
        
        #
        shifted_mod_xhm0_dphi = mod_xhm0_dphi - mod_xhm0_min_dphi_22
        
        #
        # plot( f, shifted_mod_xhm0_dphi, label='PhenomXHM', ls='--',lw=1,alpha=1,color='k',zorder=-10 )
        
        #
        shifted_f_lorentzian_min = dphi_lorentzian_min - dphi_lorentzian_min_22
        
        #
        shift_dict[simname] = shifted_f_lorentzian_min
        
        # ----
        
        # #
        # qnmo_p = qnmobj( Mf, Xf, ll, mm,0,p=1,use_nr_convention=True,verbose=False,calc_slm=False,calc_rlm=False )
        # fring  = qnmo_p.CW.real / (2*pi)

        # # Load data for this case
        # # --- (2,2)
        # raw_data_22 = loadtxt(f_.replace('l%im%i'%(ll,mm),'l2m2')).T
        # alert(simname)
        # calibration_data_22, dphi_lorentzian_min_22, f_min_22, f_max_22, f_lorentzian_min_22 = determine_data_fitting_region( raw_data_22, fring, lm=(2,2), floor_dphi=False, plot=False, simname=simname)
        
        # #
        # # NOTE: smooth_dphi=True is set because data must be smoothed before the fitting region is selected in order to avoid smoothing related boundary effects 
        
        # # --- (ll,mm)
        # raw_data = loadtxt(f_).T
        # calibration_data, dphi_lorentzian_min, f_min, f_max, f_lorentzian_min = determine_data_fitting_region( raw_data, fring, lm=(ll,mm), floor_dphi=False, plot=False, simname=simname)

        # #
        # f_22,amp_fd_22,dphi_fd_22_raw,alpha_22,beta_22,gamma_22 = calibration_data_22.T
        # f,amp_fd,dphi_fd_raw,alpha,beta,gamma = calibration_data.T
        
        # #
        # theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z,_,_ = metadata_dict['array_data'][k]
        # chi1_vec = array([chi1_x,chi1_y,chi1_z])
        # chi2_vec = array([chi2_x,chi2_y,chi2_z])
        
        # #
        # action_helper_22 = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(2,2), option_shorthand='4-xhm',turn_on_relative_dphi_mode=True)
        # action_helper = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(ll,mm), option_shorthand='4-xhm',turn_on_relative_dphi_mode=True)
        # #
        # mod_xhm0_amp_22,mod_xhm0_dphi_22,mod_xhm0_min_dphi_22 = action_helper_22(f)
        # mod_xhm0_amp,mod_xhm0_dphi,mod_xhm0_min_dphi          = action_helper(f)
        
        
        # # NOTE: smooth_dphi=True is set because data must be smoothed before the fitting region is selected in order to avoid smoothing related boundary effects. See code above.
        # dphi_fd_22 = dphi_fd_22_raw  # smooth(dphi_fd_22_raw,width=20).answer
        # dphi_fd    = dphi_fd_raw     # smooth(dphi_fd_raw,width=10).answer
        
        # #
        shifted_mod_xhm0_dphi = mod_xhm0_dphi - mod_xhm0_min_dphi_22
        # # shifted_opt_dphi = opt_dphi - min_dphi_22
        
        # # ------------------------------- #
        # df = f[1]-f[0]
        # N = len(f)
        # T = 1.0 / df
        shifted_dphi_fd = dphi_fd + shifted_f_lorentzian_min
        # shifted_f_lorentzian_min = dphi_lorentzian_min - dphi_lorentzian_min_22
        # print('-> ',T,dphi_lorentzian_min_22,dphi_lorentzian_min)
        # print('*> ',mod(dphi_lorentzian_min_22,T))
        # print('*> ',mod(dphi_lorentzian_min,T))
        # print('*> ',mod(dphi_lorentzian_min,T)-mod(dphi_lorentzian_min_22,T),'\n')
        # # ------------------------------- #
        
        
        
        # # ******************************** #
        # # Determine the optimal offset such that PNR overlaps with the calibration data (shifted_dphi_fd)
        # # opt_pnr_shift = mean(shifted_dphi_fd-shifted_opt_dphi)
        # # shifted_opt_dphi
        # # shifted_opt_dphi_final = shifted_opt_dphi + opt_pnr_shift
        # # ******************************** #

        # PLOTTING
        # ---

        #
        sca(ax[p]); p+=1
        plot( f, shifted_dphi_fd, label='Calibration Data (NR)', lw=4, alpha=0.2, color='k' )
        
        plot( f, shifted_mod_xhm0_dphi, label='PhenomXPHM', ls='--',lw=1,alpha=1,color='k',zorder=-10 )
        axhline(shifted_f_lorentzian_min,ls='--',color='tab:blue',label=r'Min for $(\ell,m)=(%i,%i)$'%(ll,mm))
        
        xscale('log')
        xlim( max(f)*0.7,max(f) )
        
        # ylim(-50,100)
        yl = lim(list(limy(f,shifted_dphi_fd,dilate=0.1)) + list(limy(f,shifted_mod_xhm0_dphi,dilate=0.1)))
        ylim( yl )
        axhline(0,ls=':', c='k',alpha=0.4,label=r'Min for $(\ell,m)=(%i,%i)$'%(2,2))
        title(simname,size=12,loc='left')
        legend(ncol=2,loc=1)
        ylabel(r'$\frac{d}{df}\arg(\tilde{h}_{%i%i})$'%(ll,mm))
        xlabel('$fM$')
        title(simname,loc='left',size=12)
        
        # #
        # shift_dict[simname] = shifted_f_lorentzian_min
            
            
    #
    file_path = datadir+'time_shift_diagnostic_l%im%i.pdf'%(ll,mm)
    alert('Saving batch plot to %s'%magenta(file_path))
    savefig(file_path,pad_inches=2, bbox_inches = "tight")

    #
    data_path = datadir+'dphi_shift_dict_l%im%i.pickle'%(ll,mm)
    alert('Saving dphi_shift_dict to %s'%magenta(data_path))
    print('>> ',shift_dict)
    pickle.dump( shift_dict, open( data_path, "wb" ) )
        
    alert('Done.')