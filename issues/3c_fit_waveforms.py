#!/usr/bin/env python3

# Setup the environment
import lalsimulation as lalsim
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import xcp
from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict,template_amp_phase

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
    alert('We are not on the expected branch. This may cause unexpected behavior.',say=True)

#
from numpy.linalg import norm
from scipy.optimize import curve_fit,minimize,fmin

#
ll = 2

#
datadir = '/Users/book/KOALA/PhenomXCP/data/version2/'
files = glob( datadir+'*_l%im%i.txt'%(ll,ll) )
files.sort()

files = files[::-1]

#
fig,ax = subplots( len(files), 2, figsize=3*array([ 2.5*2/(0.618), 2.0*len(files) ]) )
ax = ax.flatten()

#
tight_layout(1,2,4)

#
foo = {}

#
lmlist = [ (ll,ll) ]

#
p = 0
popt_array  = zeros( (len(files),8) )
fit_obj_list = []
physical_param_array = zeros( (len(files), 17) )
alert('Plotting ...')
for j,f_ in enumerate(files):

    #
    simname = f_.split('/')[-1].split('_l%im%i.'%(ll,ll))[0]

    # Find index location of metadata for simname 
    k = [ k for k,val in enumerate(metadata_dict['simname']) if val in simname ][0]

    # Load data for this case
    raw_data = loadtxt(f_).T
    calibration_data, dphi_lorentzian_min, f_min, f_max, f_lorentzian_min = determine_data_fitting_region( raw_data )

    # Collect params for this case 
    metadata = metadata_dict['array_data'][k,:]

    #
    f,amp_fd,dphi_fd,alpha,beta,gamma = calibration_data.T
    theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z = metadata_dict['array_data'][k]
    
    #
    chi1_vec = array([chi1_x,chi1_y,chi1_z])
    chi2_vec = array([chi2_x,chi2_y,chi2_z])
    
    #
    physical_param_array[j,:] = metadata_dict['array_data'][k]
    
    # GENERATE TEMPLATE FUNCTIONS
    # ---
    action_helper = template_amp_phase(m1, m2, chi1_vec, chi2_vec,ell=2)
    def action( p ):
        #
        amplitude,phase_derivative = action_helper( f, *p )
        # -- Calculate residual of phase derivative -- #
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
        return combined_residual
    
    # PERFORM FIT
    # --- 
    
    # Calculate default model amp and phase derivative
    mod_xhm0_amp,mod_xhm0_dphi = action_helper(f)
    
    # Perform fit
    foo = minimize( action,[0,0,0,0,0,0,0,0] )
    best_fit_amp,best_fit_dphi = action_helper( f, *foo.x )
    
    # Store fit params and cov 
    popt_array[j,:] = foo.x
    fit_obj_list.append( foo )
    
    # PLOTTING
    # ---
    
    #
    sca(ax[p]); p+=1
    plot( f, dphi_fd, label='Calibration Data (NR)', lw=4,ls='-', alpha=0.15, color='k' )
    plot( f, mod_xhm0_dphi, label='PhenomX(0)', ls='--',lw=1,alpha=0.85,color='k',zorder=-10 )
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
    plot( f, mod_xhm0_amp, label='PhenomX(0)', ls='--',lw=1,alpha=0.85,color='k',zorder=-10 )
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
file_path = datadir+'waveform_fit_diagnostic_l%im%i.pdf'%(ll,ll)
alert('Saving batch plot to %s'%magenta(file_path))
savefig(file_path,pad_inches=2, bbox_inches = "tight")
alert('Done.')

# SAVE FIT DATA
# --

# Initial binary parameters
data_path = datadir+'fit_initial_binary_parameters.txt'
alert('Saving %s to %s'%( magenta('physical_param_array'), magenta(data_path)) )
savetxt( data_path, physical_param_array, header='see "issues/3a_collect_metadata.py"; columns are theta, m1, m2, eta, delta, chi_eff, chi_p, chi1, chi2, a1, a2, chi1_L_x, chi1_L_y, chi1_L_z, chi2_L_x, chi2_L_y, chi2_L_z' )

# Fit parameters
data_path = datadir+'fit_opt_parameters.txt'
alert('Saving %s to %s'%( magenta('dphi_popt_array'), magenta(data_path)) )
savetxt( data_path, popt_array, header='see "template_together()" in core.py; columns are mu1, mu2, mu3, nu4, nu5, nu6, zeta1, zeta2' )

#
data_path = datadir+'fit_objects.pickle'
alert('Saving fit_obj_list to %s'%magenta(data_path))
pickle.dump( fit_obj_list, open( data_path, "wb" ) )

#
alert('Fitting complete.')
alert('Plotting complete.')
alert('Saving complete.')