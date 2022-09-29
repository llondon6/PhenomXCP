#!/usr/bin/env python3

# Setup the notebook's environment
import lalsimulation as lalsim
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import xcp
from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict,template_amp_phase,gc

# Let the user know where lalsimulation lives

#
lalsim_path = lalsim.__path__[0]
lalsuite_repo_path = lalsim_path.split('lib')[0]+'src/lalsuite/'
branch_name = bash('cd %s && git status'%lalsuite_repo_path).decode("utf-8").split('On branch ')[-1].split('\n')[0]

#
alert('We are getting our LALSimulation from:\n%s'%magenta(lalsim_path))
alert('We think that the related lalsuite source files are here:\n%s'%green(lalsuite_repo_path))
alert('Lastly, we are currently on this branch: %s'%bold(magenta(branch_name)))

# #
# if branch_name != 'pnrv1-ll':
#     alert('We are not on the expected branch. This may cause unexpected behavior.',say=True)
    
    
#

#
from numpy.linalg import norm
from scipy.optimize import curve_fit

#
calibration_data_dict = {}

#
datadir = '/Users/book/KOALA/PhenomXCP/data/version4/'
foo_path = datadir+'dphi_shift_dict_l%im%i.pickle'%(3,3)
dphi_shift_dict_l3m3 = pickle.load( open( foo_path, "rb" ) )
#
dphi_shift_dict_l2m2 = { simname:0 for simname in dphi_shift_dict_l3m3 }

# For all pairs of l and m in the global config file
for ll,mm in gc.lmlist:
    
    #
    calibration_data_dict[ll,mm] = {}

    #
    files = glob( datadir+'*_l%im%i.txt'%(ll,mm) )
    files.sort()
    
    # Ignore select files
    files = [ f for f in files if ('fit' not in f)  ]
    if ll==3:
        files = [ f for f in files if ('q1' not in f)  ]

    #
    fig,ax = subplots( len(files), 2, figsize=3*array([ 2.5*2/(0.618), (2.5)*len(files) ]) )
    ax = ax.flatten()

    #
    tight_layout(w_pad=4,h_pad=4)

    #
    foo = {}
    
    # Load time shifts
    if ll==3:
        dphi_shift_dict = dphi_shift_dict_l3m3
    else:
        dphi_shift_dict = dphi_shift_dict_l2m2
    # print('l%im%i> '%(ll,mm),dphi_shift_dict)
    

    #
    p = 0
    for f_ in files[::-1]:

        #
        simname = f_.split('/')[-1].split('_l%im%i.'%(ll,mm))[0]

        # Find index location of metadata for simname 
        k = [ k for k,val in enumerate(metadata_dict['simname']) if val in simname ][0]

        # Collect params for this case 
        theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z,Mf,Xf = metadata_dict['array_data'][k,:]
        chi1_vec = array([chi1_x,chi1_y,chi1_z])
        chi2_vec = array([chi2_x,chi2_y,chi2_z])

        # ****************************** #
        # Determine data fitting region  #
        # ****************************** #
        # Load QNM info
        qnmo_p = qnmobj( Mf, Xf, ll, mm,0,p=1,use_nr_convention=True,verbose=False,calc_slm=False,calc_rlm=False )
        fring  = qnmo_p.CW.real / (2*pi)
        # Load data for this case
        raw_data = loadtxt(f_).T
        calibration_data, dphi_lorentzian_min, f_min, f_max, f_lorentzian_min = determine_data_fitting_region( raw_data, fring, lm=(ll,mm), floor_dphi=True, plot=not True, simname=simname)
        
        #
        f,amp_fd,dphi_fd,alpha,beta,gamma = calibration_data.T
        
        # tell the XCP package to generate a PhenomXPHM waveform
        action_helper_l2m2 = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(2,2),option_shorthand='2-xphm',floor_dphi=False)
        action_helper = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(ll,mm),option_shorthand='2-xphm',floor_dphi=False)
        # action_helper = template_amp_phase(m1, m2, chi1_vec, chi2_vec,lm=(ll,mm))
        _,xphm_dphi_l2m2 = action_helper_l2m2(f)
        xphm_amp,xphm_dphi = action_helper(f)
        
        
        ''' 
        NOTE
        ---
        A choice must be made about where to set the minimum value of the NR (2,2) minimum phase derivative value. The motivations guiding this choice are:
        * The tuned coprecessing (2,2) model should as much as possible maintain the relative min phase derivative (min dphi) values that are encoded within *XPHM*, not the *XHM* value. This is key as 
        * The min dphi value of the (3,3) multipole should be tuned to NR away from the value defined within *XPHM*
        * Enforcing these to criteria will, for other multipole moments, change their dphi relationships relative to the (3,3) moment, but maintain the min dphi relationship with (2,2) up to modeling error within the XCP calibration region.
        
        THUS the XPHM min dphi relationships will be used for (2,2) and (3,3) tuning, NOT those from XHM, despite the fact that all other tuning aspect will be done relative to XPHM .... This actually implies that the best thing to do is to use XPHM as a base model, *not* XHM. However, this would have the primary complication of causing structure in deviations to the ringdown frequency that cannot be modeled with polynomials. ...
        '''
        
        #
        min_xphm_dphi_l2m2 = min(xphm_dphi_l2m2)
        
        # 
        nr_dphi_lm_shift = dphi_shift_dict[simname] if not ('q1' in simname) else 0
        dphi_fd_enforced_min = min_xphm_dphi_l2m2 + nr_dphi_lm_shift
        dphi_fd += dphi_fd_enforced_min
        
        
        # print('** ',simname)
        # print('>> ',dphi_fd_enforced_min)
        # print('>> ',dphi_shift_dict[simname])
        # print('>> ',sum(dphi_shift_dict.values()))
        # error()
        
        # EZH effective ringdown 
        mod_xhm_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, [(ll,mm)], m1, m2, chi1_vec, chi2_vec, option_shorthand='3-xphm-ezh' )
        # mod_xhm_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, [(ll,mm)], m1, m2, chi1_vec, chi2_vec, pflag= 501 if ll==2 else 500, fsflag=None )
        mod_xhm = mod_xhm_dict[ll,mm]
        mod_xhm_amp = abs(mod_xhm)
        mod_xhm_phi = unwrap( angle(mod_xhm) )
        mod_xhm_dphi = spline_diff(f,mod_xhm_phi)
        # mod_xhm_dphi -= min( mod_xhm_dphi[ (f>0.03*ll/2)&(f<0.12*ll/2) ] )

        # XAS/XHM
        xhm_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, [(ll,mm)], m1, m2, chi1_vec, chi2_vec, option_shorthand='4-xhm' )
        xhm = xhm_dict[ll,mm]
        xhm_amp = abs(xhm)
        xhm_phi = unwrap( angle(xhm) )
        xhm_dphi = spline_diff(f,xhm_phi)
        # xhm_dphi -= min( xhm_dphi[ (f>0.03*ll/2)&(f<0.12*ll/2) ] )

        # PLOTTING
        # ---

        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&% #
        #           Store all useful data for output         #
        calibration_data_dict[ll,mm][simname] = (metadata_dict['array_data'][k,:],f,dphi_fd,amp_fd,xphm_dphi,dphi_fd_enforced_min,nr_dphi_lm_shift,min_xphm_dphi_l2m2)
        # &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&% #

        #
        sca(ax[p]); p+=1
        plot( f, dphi_fd, label='Calibration Data (NR)', lw=2, alpha=1, color='k' )
        plot( f, xhm_dphi, label='PhenomXHM', ls=':',lw=2,alpha=0.85,color='k' )
        plot( f, mod_xhm_dphi, label='PhenomXPHM(501:EZH-EffRD)', ls='--',lw=2,alpha=0.85,color='r' )
        plot( f, xphm_dphi, label='PhenomXPHM', ls='--',lw=4,alpha=0.25,color='k',zorder=-10 )
        
        # xscale('log')
        #xlim(lim(f,dilate=1.1,dilate_with_multiply=True))
        xlim( max(f)*0.7,max(f) )
        
        # ylim( limy(f, mod_xhm_dphi,dilate=0.1) )
        yl = lim(list(limy(f,dphi_fd,dilate=0.1)) + list(limy(f,xhm_dphi,dilate=0.1)) + list(limy(f,mod_xhm_dphi,dilate=0.1)) + list(limy(f,xphm_dphi,dilate=0.1)))
        ylim( yl )
        
        title(simname,size=12,loc='left')
        ylabel(r'$\frac{d}{df}\arg(\tilde{h}_{%i%i})$'%(ll,mm))
        xlabel('$fM$')
        title(simname+', nr_dphi_lm_shift: %f'%nr_dphi_lm_shift,loc='left',size=12)
        # axhline(0,ls='--',color='k',alpha=0.7)
        if ll!=2: axhline(min_xphm_dphi_l2m2,ls=':',color='k',alpha=0.7,label='XPHM $(2,2)$ Min')
        #
        axhline(dphi_fd_enforced_min,c='tab:green',ls='--',label='NR Relative Value' if ll==3 else 'XPHM Relative Value')
        # alert(dphi_fd_enforced_min)
        legend(ncol=2,loc=1)

        sca(ax[p]); p+=1
        plot( f, amp_fd, label='Calibration Data (NR)', lw=2, alpha=1, color='k' )
        plot( f, xhm_amp, label='PhenomXHM', ls=':',lw=2,alpha=0.85,color='k' )
        plot( f, mod_xhm_amp, label='PhenomXPHM(501:EZH-EffRD)', ls='--',lw=2,alpha=0.85,color='r' )
        plot( f, xphm_amp, label='PhenomXPHM', ls='--',lw=4,alpha=0.25,color='k',zorder=-10 )
        yscale('log')
        xscale('log')
        legend(ncol=2)
        ylim( limy(f, amp_fd,dilate=1.2) )
        xlabel('$fM$')
        ylabel(r'$|\tilde{h}_{%i%i}(f)|$'%(ll,mm))
        #
        title(simname,loc='left',size=12)
            
            
    #
    file_path = datadir+'waveform_prefit_diagnostic_l%im%i.pdf'%(ll,mm)
    alert('Saving batch plot to %s'%magenta(file_path))
    savefig(file_path,pad_inches=2, bbox_inches = "tight")
    alert('Done.')

#
data_path = datadir+'calibration_data_dict.pickle'
alert('Saving calibration_data_dict to %s'%magenta(data_path))
pickle.dump( calibration_data_dict, open( data_path, "wb" ) )