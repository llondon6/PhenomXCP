

#
import pickle
from positive import alert, magenta, parent, red
from numpy import loadtxt, load
from os.path import exists
import xcp
import json

#
package_dir = parent(xcp.__path__[0])
data_dir = package_dir + 'data/'

# Detect python version
import sys
PYTHON3 = True
if sys.version_info[0] < 3:
    PYTHON3 = False

# Always load catalog list for calibration runs 
calibration_catalog_path = data_dir+'calibration_catalog.pickle'
if PYTHON3:
    calibration_catalog = pickle.load( open( calibration_catalog_path, "rb" ), encoding='latin1' )
else:
    calibration_catalog = pickle.load( open( calibration_catalog_path, "rb" ) )
alert('Catalog of calibration runs stored to %s'%magenta('"xcp.calibration_catalog"'),fname='xcp.core')

# Always load curated metadata for calibration runs 
'''NOTE that this files is created by issues/3a_collect_metadata.py'''
metadata_dict_path = data_dir+'metadata_dict.pickle'
metadata_dict = load(metadata_dict_path,allow_pickle=True)
alert('Metadata dictionary for calibration runs stored to %s'%magenta('"xcp.metadata_dict"'),fname='xcp.core')
        
#
catalog_paper_md_path = data_dir+'catalog_paper_metadata.json'
with open(catalog_paper_md_path, 'r') as f:
    catalog_paper_metadata = json.load(f)
alert('Metadata dictionary for Ed\'s catalog paper stored to %s'%magenta('"xcp.catalog_paper_metadata"'),fname='xcp.core')


#
def LALPolarizationsFD(approximant, lmlist, m1, m2, s1, s2, delta_f, phiRef=0,nu0 = 0,pflag=501, ReturnCoPrec=1):
    
    #
    import lalsimulation as lalsim
      
    lalparams = lal.CreateDict()
    
    #
    output_modes = {}
    
    #
    ModeArray = lalsim.SimInspiralCreateModeArray()
    for mode in lmlist:
        
        #
        l,m = mode
        
        #
        lalsim.SimInspiralModeArrayActivateMode(ModeArray, l,m)
        lalsim.SimInspiralWaveformParamsInsertModeArray(lalparams, ModeArray)

        #
        lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(lalparams, 0)

        #
        lalsim.SimInspiralWaveformParamsInsertPhenomXReturnCoPrec(lalparams, ReturnCoPrec)
        
        #
        if pflag:
            lalsim.SimInspiralWaveformParamsInsertPhenomXPrecVersion( lalparams, pflag )

        #
        f_min       = 10.0
        f_max       = 2048.0
        Omega       = 0.
        inclination = 0 # Chosen so that it doesn't correxpond to a spherical hamonic root
        distance_Mpc= 100.0
        distance    = distance_Mpc*1.0e6*lal.PC_SI

        Hp, Hc = lalsim.SimInspiralChooseFDWaveform(m1=lal.MSUN_SI*m1,
                                                m2=lal.MSUN_SI*m2, 
                                                S1x=s1[0], S1y=s1[1], S1z=s1[2],
                                                S2x=s2[0], S2y=s2[1], S2z=s2[2],
                                                distance=distance, 
                                                inclination=inclination, 
                                                LALpars=lalparams,
                                                phiRef=phiRef, 
                                                f_ref=f_min,
                                                deltaF=delta_f,
                                                f_min=f_min,
                                                f_max=f_max,
                                                longAscNodes=Omega,
                                                eccentricity=0.0,
                                                meanPerAno=0.0,
                                                approximant=approximant) 

        #
        freqs = np.arange(len(Hp.data.data)) * delta_f
        hp = Hp.data.data
        hc = Hc.data.data
        
        if not ( approximant in (lalsim.IMRPhenomXP,lalsim.IMRPhenomXPHM) ):

            #
            s = -2
            spherical_harmonic = sYlm(s,l,m,inclination,phiRef)
            hp /= spherical_harmonic
            hc /= spherical_harmonic
        
        #
        Mtot = m1+m2
        hp = codehf(hp,Mtot,distance_Mpc)
        hc = codehf(hc,Mtot,distance_Mpc)
        f  = codef(freqs,Mtot) 
        
        #
        output_modes[l,m] = (hp,hc,f)
    
    #
    return output_modes


#
def get_phenomxphm_coprecessing_multipoles(freqs, lmlist, m1, m2, s1, s2, phiRef=0, pflag=500, mu1=0, mu2=0, mu3=0, mu4=0, nu4=0, nu5=0, nu6=0, zeta1=0, zeta2=0 ):
    '''
    Generate dictionary of waveform arrays corresponding to input multipole list (i.e. list of [l,m] pairs ). If a single l,m pair is provided, then a single waveform array will be returned (i.e. we have opted to not have a lower-level function called "phenomxhm_multipole").
    
    USAGE
    ---
    output_modes_dict = phenomxhm_multipoles(approximant, lmlist, m1, m2, s1, s2, delta_f, phiRef, nu0=0)
    
    NOTES
    ---
    IF lmlist list of tuples, e.g. [(2,2),(2,1)], 
    THEN output_modes_dict is dictionary of [freq, hplus, hcross] in code units. 
    ELSE IF lmlist is single list of ell and m, e.g. (2,2) or [2,2], 
    THEN a single array, i.e. [freq, hplus, hcross], is returned.
    * pflag=500 is the default PhenomXPNR framework
    '''
    
    # Import usefuls 
    from numpy import ndarray, array, arange, double
    from positive.units import codef,codeh,codehf,physf
    from positive import sYlm
    import lalsimulation as lalsim
    import lal
    
    #
    
    # Validate inputs 
    # ---
    
    # Check type of lmlist 
    if not isinstance(lmlist,(list,tuple,ndarray)):
        error('lmlist input must be iterable ')
    # If a single mode is given, make lmlist a list of that mode's indices 
    single_mode_requested = False
    if len(lmlist)==2:
        l,m = lmlist
        if isinstance(l,int) and isinstance(m,int):
            single_mode_requested = True
            lmlist = [ lmlist ]
    
    #
    Mtot = 100.0
    M1 = m1 * Mtot/ ( m1 + m2 )
    M2 = m2 * Mtot / ( m1 + m2 )
    
    # Create physical frequencies as a LAL REAL* sequence
    freqs_Hz = lal.CreateREAL8Sequence( len(freqs) )
    freqs_Hz.data = physf(freqs,Mtot)

    # Main routine 
    # ---
    
    
    #
    lalparams = lal.CreateDict()

    #
    output_modes = {}
    
    #
    ModeArray = lalsim.SimInspiralCreateModeArray()
    for lm in lmlist:
        
        #
        l,m = lm
        
        #
        lalsim.SimInspiralModeArrayActivateMode(ModeArray, l,m)
        lalsim.SimInspiralWaveformParamsInsertModeArray(lalparams, ModeArray)

        # Turn off multibanding
        lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(lalparams, 0)

        # Tell the model to return the coprecessing mode -- only works on our development branches
        lalsim.SimInspiralWaveformParamsInsertPhenomXReturnCoPrec(lalparams, 1)
        
        #
        
        # Set deviations from base model based on inputs
        # mu2,mu3,mu4,nu4,nu5,nu6,zeta2 = [ double(k) for k in (mu2,mu3,mu4,nu4,nu5,nu6,zeta2) ]
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU1(lalparams, mu1)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU2(lalparams, mu2)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU3(lalparams, mu3)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU4(lalparams, mu4)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU4(lalparams, nu4)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU5(lalparams, nu5)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU6(lalparams, nu6)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPZETA1(lalparams, zeta1)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPZETA2(lalparams, zeta2)
        
        #
        if pflag:
            lalsim.SimInspiralWaveformParamsInsertPhenomXPrecVersion( lalparams, pflag )

        #
        distance_Mpc= 100.0
        distance_SI    = distance_Mpc*1.0e6*lal.PC_SI
        m1_SI = lal.MSUN_SI*M1
        m2_SI = lal.MSUN_SI*M2
        chi1x, chi1y, chi1z = s1
        chi2x, chi2y, chi2z = s2
        fRef_In = 0
        
        #
        Hp, Hc = lalsim.SimIMRPhenomXPHMFrequencySequence( freqs_Hz, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, distance_SI, 0, phiRef, fRef_In, lalparams)

        #
        hp = Hp.data.data
        hc = Hc.data.data
        
        #
        hp = codehf(hp,Mtot,distance_Mpc)
        hc = codehf(hc,Mtot,distance_Mpc)
        
        #
        output_modes[l,m] = hp + 1j * hc

    #
    return output_modes


#
def template_amp_phase(m1, m2, chi1_vec, chi2_vec, ell=2):
    
    # NOTE that mu4 is no longer to be used as it is completely degenerate with nu5 in PhenomX
    
    #
    import xcp
    from numpy import unwrap,angle 
    from positive import spline_diff
    
    #
    lmlist = [ (ell,ell) ]
    
    #
    def template_together( f, mu1=0, mu2=0, mu3=0, nu4=0, nu5=0, nu6=0, zeta1=0, zeta2=0 ):
        
        # Set phase deviations to zero
        mu4 = 0 # No longer to be used as it is completely degenerate with nu5 in PhenomX
        
        # Calculate PhenomXPHM with the input amplitude deviations
        # NOTE that pflag=0 means that we use the default setting of PhenomXPHM as a reference model
        try:
            multipole_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, lmlist, m1, m2, chi1_vec, chi2_vec, pflag=0, mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4, nu4=nu4, nu5=nu5, nu6=nu6, zeta1=zeta1, zeta2=zeta2 )
        except:
            multipole_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, lmlist, m1, m2, chi1_vec, chi2_vec, pflag=0 )
        
        # 
        complex_strain = multipole_dict[ell,ell]
            
        # Given the complex FD waveform, compute its amplitude
        amplitude = abs(complex_strain)
        # Given the complex FD waveform, compute its phase derivative
        complex_strain = multipole_dict[ell,ell]
        phase = unwrap( angle(complex_strain) )
        phase_derivative = spline_diff(f,phase)
        # Find min phase derivative
        mask = (f>0.03)&(f<0.12)
        min_phase_derivative = min( phase_derivative[ mask ] )
        # Adjust phase derivative 
        phase_derivative -= min_phase_derivative
        
        #
        return amplitude,phase_derivative
    
    # #
    # def template_amp( f, mu2=0, nu5=0 ):
        
    #     # Set phase deviations to zero
    #     # NOTE that mu4 is no longer to be used as it is completely degenerate with nu5 in PhenomX
    #     nu4=0
    #     nu6=0
    #     zeta2=0
        
    #     # Calculate PhenomXPHM with the input amplitude deviations
    #     # NOTE that pflag=0 means that we use the default setting of PhenomXPHM as a reference model
    #     multipole_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, lmlist, m1, m2, chi1_vec, chi2_vec, pflag=0, mu2=mu2, nu4=nu4, nu5=nu5, nu6=nu6, zeta2=zeta2 )
        
    #     # Given the complex FD waveform, compute its amplitude
    #     complex_strain = multipole_dict[ell,ell]
    #     amplitude = abs(complex_strain)
        
    #     #
    #     return amplitude
        
    # #
    # def template_dphi( f, nu4=0, nu5=0, nu6=0, zeta2=0 ):
        
    #     # Set amplitude deviations to zero
    #     mu2=0
    #     mu4=0
        
    #     # Calculate PhenomXPHM with the input phase deviations
    #     # NOTE that pflag=0 means that we use the default setting of PhenomXPHM as a reference model
    #     multipole_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, lmlist, m1, m2, chi1_vec, chi2_vec, pflag=0, mu2=mu2, mu4=mu4, nu4=nu4, nu5=nu5, nu6=nu6, zeta2=zeta2 )
        
    #     # Given the complex FD waveform, compute its amplitude
    #     complex_strain = multipole_dict[ell,ell]
    #     phase = unwrap( angle(complex_strain) )
    #     phase_derivative = spline_diff(f,phase)
    #     # Shift such that the min is zero
    #     phase_derivative -= min( phase_derivative[ (f>0.03)&(f<0.12) ] )
        
    #     #
    #     return phase_derivative
        
    # #
    # def make_template_dphi( mu2=0, nu5=0 ):
        
    #     #
    #     def template_dphi( f, nu6=0, nu4=0, zeta2=0 ):
            
    #         # Set amplitude deviations to zero
    #         mu4=0
            
    #         # Calculate PhenomXPHM with the input phase deviations
    #         # NOTE that pflag=0 means that we use the default setting of PhenomXPHM as a reference model
    #         multipole_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, lmlist, m1, m2, chi1_vec, chi2_vec, pflag=0, mu2=mu2, mu4=mu4, nu4=nu4, nu5=nu5, nu6=nu6, zeta2=zeta2 )
            
    #         # Given the complex FD waveform, compute its amplitude
    #         complex_strain = multipole_dict[ell,ell]
    #         phase = unwrap( angle(complex_strain) )
    #         phase_derivative = spline_diff(f,phase)
    #         # Shift such that the min is zero
    #         phase_derivative -= min( phase_derivative[ (f>0.03)&(f<0.12) ] )
            
    #         #
    #         return phase_derivative
            
    #     #
    #     return template_dphi
        
        
    #
    return template_together # template_amp, make_template_dphi


# Function to determine version2 data fitting region
def determine_data_fitting_region_legacy( data, fmin=0.03, fmax=0.12 ):
    '''
    Given version2 data array, determine fitting region dynamically.
    This function assumes clean data within the default or input values of fmin and fmax, and then uses the phase derivate to determine new fmin and fmax values that are ultimately used to define a fitting region.
    '''
    
    # Import usefuls
    from numpy import argmin,log
    from positive import smooth,find,lim
    from matplotlib.pyplot import figure,plot,show,axhline,xlim,ylim
    
    # Extract data 
    f,amp_fd,dphi_fd,alpha,beta,gamma = data

    # Use default domain bounds to determine a mask
    mask = (f>=fmin) & (f<=fmax)
    
    # Determine the minimum dphi
    # Smooth dphi using postiive's savgol filter
    y = smooth(dphi_fd[mask]).answer
    knot = argmin(y)
    y_knot = y[knot]
    data[2] = dphi_fd - smooth(dphi_fd[mask]).answer[knot] + y_knot
    
    # Determine new fmin and max using heuristic 
    f_knot = f[mask][knot]
    new_fmin = max(f_knot * 0.22,0.018) # 0.5 # 0.325
    new_fmax = f_knot + 0.020 # 0.025 
    
    #
    new_mask = (f>=new_fmin) & (f<=new_fmax)
    new_data = data.T[new_mask,:]
    
    #
    new_knot = find(f>=fmin)[0]+knot
    
    #
    return new_data,new_knot,new_fmin,new_fmax,f_knot


# Function to determine version2 data fitting region
def determine_data_fitting_region(data, threshold=0.015):
    '''
    Given version2 data array, determine fitting region dynamically.
    This function assumes clean data within the default or input values of fmin and fmax, and then uses the phase derivate to determine new fmin and fmax values that are ultimately used to define a fitting region.
    '''

    # Import usefuls
    from numpy import argmin, log, arange
    from positive import smooth, find, lim, smoothest_part_by_threshold,findpeaks
    from matplotlib.pyplot import figure, plot, show, axhline, xlim, ylim

    # DETERMINE DATA FITTING REGION
    # ---
    
    # 0. Select and unpack
    f,amp_fd,dphi_fd,alpha,beta,gamma = data
    
    # 1. Find smoothest part
    pre_mask = (f>0) #& (f<0.12)
    mask_1     = smoothest_part_by_threshold( dphi_fd[pre_mask], threshold=threshold, smooth_width=20, plot=False )
    dphi_fd_1  = smooth( dphi_fd[pre_mask][mask_1], width=30 ).answer
    f_1 = f[pre_mask][mask_1]
    
    # 2. Handle unstable end behavior and remask
    peaks,peak_locations = findpeaks( dphi_fd_1 )
    if len(peak_locations):
        mask_2    = range(0,peak_locations[-1]+1)
        dphi_fd_2 = dphi_fd_1[mask_2]
        f_2 = f_1[mask_2]
        mask_3 = f_2<0.12
        dphi_fd_3 = dphi_fd_2[ mask_3 ]
    else:
        mask_2 = mask_1
        dphi_fd_3 = dphi_fd_1
    
    # 3. Determine location of the lorentzian min
    lorentzian_mindex = mask_1[mask_2[argmin( dphi_fd_3 )]]
    dphi_lorentzian_min = min( dphi_fd_3 )
    
    # 4. Use the lorentzian_mindex to define the start and end of the fitting region
    f_lorentzian_min = f[pre_mask][ lorentzian_mindex ]
    # f_lorentzian_min = f[mask_1][mask_2][ lorentzian_mindex ]
    f_min = f_lorentzian_min * 0.2 
    # f_min = max(f_lorentzian_min * 0.22, 0.018) 
    f_max = f_lorentzian_min + 0.012# 0.018 # 0.012
    calibration_mask = arange(len(f))[(f>=f_min) & (f<=f_max)]
    
    # 5. Select region for output 
    calibration_data = data.T[ calibration_mask ]
    # 6. Smooth phase derivative data and subtract min
    calibration_data.T[2] = dphi_fd[calibration_mask] - dphi_lorentzian_min
    # calibration_data.T[2] = smooth( dphi_fd, width=10 ).answer[calibration_mask] - dphi_lorentzian_min
    # calibration_data.T[2] = smooth( calibration_data.T[2], width=30 ).answer
    

    #
    return calibration_data, dphi_lorentzian_min, f_min, f_max, f_lorentzian_min


#
def select_scenty_metadata( sceo ):
    
    '''
    Given nrutils' scentry object, collect metedata useful for generating model waveforms 
    '''
    
    #
    from numpy.linalg import norm
    from numpy import arccos,dot,pi,array
    from positive.physics import calc_chi_eff,calc_chi_p
    
    #
    a = sceo
    
    #
    l = a.L/norm(a.L)
    if (abs(a.m1-a.m2)<1e-3) and (norm(a.X1)<norm(a.X2)):
        a.X1,a.X2 = [ array(k) for k in (a.X2,a.X1) ]
        a.m1,a.m2 = [ float(k) for k in (a.m2,a.m1) ]

    #
    m1,m2 = [ k/(a.m1+a.m2) for k in (a.m1,a.m2) ] 
    eta = m1*m2/(m1+m2)
    
    #
    X1,X2,L,S = a.X1,a.X2,a.L,a.S
    
    #
    a1,a2 = norm(a.X1),norm(a.X2)
    
    #
    l = L/norm(L)
    s = S/norm(S)
    
    # NOTE that theta is in radians
    theta = arccos( dot( l, s ) ) 
    
    #
    chi1 = dot(X1,l)
    chi2 = dot(X2,l)
    
    #
    chi_p   = calc_chi_p(   m1,X1, m2,X2, L )
    chi_eff = calc_chi_eff( m1,X1, m2,X2, L )
    
    #
    delta = (m1-m2)/(m1+m2)
    
    #
    return theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2 

   
# Advanced gloss atop mvpolyfit.plot and mvrfit.plot
def advanced_gmvx_plot( fit_object ):
    
    '''
    Advanced gloss atop mvpolyfit.plot and mvrfit.plot
    '''
    
    from matplotlib.pyplot import subplots, plot, xlabel, ylabel, title, sca, gca, figaspect, tight_layout
    from numpy import cos,sin,array,around,ones_like,sort,pi,linspace
    from positive import eta2q,q2eta,eta2delta,rgb
    from glob import glob
    from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict
    
    # Load and unpack physical parameter space
    raw_domain = loadtxt(data_dir+'version2/fit_initial_binary_parameters.txt')
    theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_ = raw_domain.T


    # Define desired model domain variables and array 
    u = cos(theta)
    v = sin(theta)
    q = 1.0/eta2q(eta)

    # --------------------------------------- #
    # Plot ans save fits 
    # --------------------------------------- #

    # Collect set of unique a1 values
    a1_point = around(a1,2)
    a1_set = array(sort(list( set(a1_point) )))

    # Collect set of unique angle values
    degree_point = (around( (theta*180/pi)/10 )*10).astype(int)
    theta_point = degree_point*pi/180
    theta_set = array(sort(list( set(theta_point) )))

    # Collect set of unique mass-ratio values
    q_point = around(array([eta2q(n) for n in eta]),2)
    q_set = array(sort(list( set(q_point) )))

    # Collect set of unique eta values
    eta_point = q2eta( q_point )
    eta_set = q2eta(q_set)

    # Summary figure for internal diagnostics 
    summary_fig = fit_object.plot(size_scale=1.5)
    ax = summary_fig.axes
    
    #
    sca(ax[0])
    title(fit_object.labels['python'][0])
    
    #
    num_figs = len(a1_set)*len(theta_set)
    eta_set_figs,set_fig_ax = subplots( len(a1_set), len(theta_set), figsize=5*array([ len(theta_set),len(a1_set) ]) )
    set_fig_ax = set_fig_ax.flatten();
    tight_layout(4,4)
    ax_counter = 0
    
    #
    colors = rgb(len(a1_set),jet=True)

    #
    for k,_a1 in enumerate(a1_set):
        for _theta in theta_set:

            #
            theta_mask = (_theta==theta_point)
            a1_mask = (_a1==a1_point)
            mask = a1_mask & theta_mask

            #
            _eta = eta_point[mask]
            _u = cos(_theta) 

            #
            case_eta   = linspace( min(_eta),max(_eta),1000 ) 
            case_delta = eta2delta( case_eta )
            case_q     = 1.0/eta2q(case_eta)  
            case_theta = _theta * ones_like(case_eta)
            case_u     = cos(case_theta)
            case_a1    = _a1    * ones_like(case_eta)

            #
            case_domain = array([case_u,case_eta,case_delta,case_a1]).T
            case_range = fit_object.eval(case_domain)
            opt_range  = fit_object.eval(fit_object.domain[mask,:])

            #
            sca(ax[0])
            ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color=colors[k] )
            # ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color = 'tab:blue' if _a1==a1_set[0] else 'red' )
            
            #
            sca( set_fig_ax[ax_counter] ); ax_counter += 1
            plot( eta[mask], fit_object.range[mask] if hasattr(fit_object,'range') else fit_object.scalar_range[mask], marker='o',ls='none'  )
            plot( eta[mask], opt_range, marker='o',ms=10,mfc='none', color='r',ls='none'  )
            plot( case_eta, case_range, ls='-', color='r' )
            title( r'$a_1=%1.2f$, $\theta=%1.2f$'%(_a1,round(_theta*180.0/pi,0)) )
            xlabel(r'$\eta$')
            ylabel(r'$\%s$'%fit_object.labels['python'][0])
            

    
    #
    num_figs = len(a1_set)*len(eta_set)
    theta_set_figs,set_fig_ax = subplots( len(a1_set), len(eta_set), figsize=5*array([ len(eta_set),len(a1_set) ]) )
    set_fig_ax = set_fig_ax.flatten()
    tight_layout(4,4)
    ax_counter = 0

    #
    for k,_a1 in enumerate(a1_set):
        for _eta in eta_set:

            #
            eta_mask = (_eta==eta_point)
            a1_mask = (_a1==a1_point)
            mask = a1_mask & eta_mask

            #
            _theta = theta_point[mask]
            _u = cos(_theta) 

            #
            case_theta   = linspace( min(_theta),max(_theta),1000 ) # 
            case_u     = cos(case_theta)
            case_eta   = _eta * ones_like(case_theta)
            case_delta = eta2delta( case_eta )
            case_a1    = _a1  * ones_like(case_theta)

            #
            case_domain = array([case_u,case_eta,case_delta,case_a1]).T
            case_range = fit_object.eval(case_domain)
            opt_range  = fit_object.eval(fit_object.domain[mask,:])

            #
            sca(ax[0])
            ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color=colors[k] )
            # ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color = 'tab:blue' if _a1==a1_set[0] else 'red' )
            
            #
            sca( set_fig_ax[ax_counter] ); ax_counter += 1
            plot( cos(theta[mask]), fit_object.range[mask] if hasattr(fit_object,'range') else fit_object.scalar_range[mask], marker='o',ls='none',color='r'  )
            plot( cos(theta[mask]), opt_range, marker='o',ms=10,mfc='none', color='b',ls='none'  )
            plot( cos(case_theta), case_range, ls='-', color='b' )
            title( r'$a_1=%1.2f$, $q=%1.2f$'%(_a1,eta2q(_eta)) )
            xlabel(r'$\cos(\theta)$')
            ylabel(r'$\%s$'%fit_object.labels['python'][0])
            
    #
    return summary_fig, eta_set_figs, theta_set_figs
 

# Remove non l=m moments and symmetrize in co-precessing frame
# QUESTION: Do the angles resulting from this procedure equate to simple averages of the original angles?
def gwylmo_cpclean( gwylmo, verbose=False, safe_domain_range=None, cp_domain=None, cp_kind=None ):

    '''
    Given a gwylm object:
     * put in co-precessing frame
     * remove higher multipoles
     * reconstruct mass-quadrupoles by symmetrization
     * put back in simulation frame
    '''

    # Import usefuls
    from numpy import array

    #
    if cp_domain in (None,'td'):
        cp_domain = 'td'
        alert('Using the '+red('Time Domain')+' coprecessing frame for cleaning.')
    elif cp_domain in ('fd','FD','frequency_domain','freq'):
        cp_domain = 'fd'
        alert('Using the '+red('Frequency Domain') +
              ' coprecessing frame for cleaning.')
    else:
        error('Unknown coprecessing domain option:'+str(cp_domain))
    #
    if cp_kind is None:
        cp_kind = 'psi4'

    #
    if safe_domain_range is None:
        safe_domain_range=[0.009,0.3]

    #
    if verbose:
        alert('Now removing current-quadrupole content in coprecessing frame, symmetrizing mass-quadrupoles and then putting back in original frame.')

    #
    y = gwylmo

    # Put in frame where j init is zhat
    y0 = y.__calc_initial_j_frame__()

    # Calc coprecessing frame wf TD
    y1 = y0.__calc_coprecessing_frame__(verbose=True, safe_domain_range=safe_domain_range,transform_domain=cp_domain,kind=cp_kind)

    # Duplicate for manipulation
    y2 = y1.copy()

    # Symmetrize
    ll,mm = y.root_lm
    for kind in y2[ll,mm]:
        yy = 0.5 * (  y2[ll,mm][kind].y + ((-1)**ll)*y2[ll,-mm][kind].y.conj()  )
        wfarr22 = array( [ y2.t, yy.real, yy.imag ] ).T
        wfarr2m2 = array( [ y2.t, yy.real,-yy.imag ] ).T
        y2[ll,mm][kind].setfields(wfarr22)
        y2[ll,-mm][kind].setfields(wfarr2m2)

    # Remove higher multipole data
    for lm in y2.lm:
        l,m = lm
        if (l,abs(m)) != (ll,mm):
            for kind in y2.lm[lm]:
                wfarr = y2[lm][kind].wfarr
                wfarr[:,1:] *= 1e-10
                y2[lm][kind].setfields( wfarr )

    # Re-wrap with the original TD or FD angles
    alpha  = y1.radiation_axis_info.radiation_axis['%s_alpha'%cp_domain]
    beta   = y1.radiation_axis_info.radiation_axis['%s_beta' %cp_domain]
    gamma  = y1.radiation_axis_info.radiation_axis['%s_gamma'%cp_domain]
    angles = ( alpha,beta,gamma )
    y3     = y2.__rotate_frame_at_all_times__( angles, transform_domain=cp_domain )
    y3.frame = y.frame

    # Return answer
    ans = y3
    return ans



# Given underlying physical parameters, calculate ones useful form modeling
def parama_party( eta,theta,a1 ):
    '''
    PARAMA-PARTY:
    If L || z and m1>m2 and q=m1/m2, then 

    S2 = 0
    S1 = m1**2 a1 * exp( 1j * theta ) = Sz + 1j*Sperp
    X1 = X = S1/m1**2

    chi_eff = m1*a1*cos(theta)/(m1+m2) = a1*cos(theta)*/(1+1.0/q)

    A1 = 2 + (3*m2)/(2*m1)
    A2 = 2 + (3*m1)/(2*m2)
    B1 = A1 * a1*sin(theta)
    B2 = 0
    chi_p = max( B1,B2 ) / ( A1 * m1*m1 )
    L = L

    '''
    
    #
    from positive import eta2m1m2
    from numpy import cos,sin,maximum
    
    #
    m1,m2 = eta2m1m2(eta)
    
    #
    q = m1/m2
    chi_eff = m1*a1*cos(theta)/(m1+m2)
    
    #
    A1 = 2 + (3.0*m2)/(2.0*m1)
    Norm_S1_perp = abs(a1*sin(theta)*m1*m1)
    B1 = A1 * Norm_S1_perp 
    chi_p = maximum( B1,0 ) / ( A1 * m1*m1 )
    
    '''
    NOTE that the above is basically a1*sin(theta), but we retain the additional steps to illustrate consistency with the more general formula
    '''
    
    #
    return chi_eff, chi_p