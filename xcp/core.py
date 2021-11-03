

#
import pickle
from positive import alert, magenta, parent
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

# # Always load curated metadata for calibration runs 
# metadata_dict_path = data_dir+'metadata_dict.pickle'
# metadata_dict = load(metadata_dict_path,allow_pickle=True)
# alert('Metadata dictionary for calibration runs stored to %s'%magenta('"xcp.metadata_dict"'),fname='xcp.core')
        
#
catalog_paper_md_path = data_dir+'catalog_paper_metadata.json'
with open(catalog_paper_md_path, 'r') as f:
    catalog_paper_metadata = json.load(f)
alert('Metadata dictionary for Ed\'s catalog paper stored to %s'%magenta('"xcp.catalog_paper_metadata"'),fname='xcp.core')


#
def LALPolarizationsFD(approximant, modeList, m1, m2, s1, s2, delta_f, phiRef=0,nu0 = 0,pflag=501, ReturnCoPrec=1):
    
    #
    import lalsimulation as lalsim
      
    lalparams = lal.CreateDict()
    
    #
    output_modes = {}
    
    #
    ModeArray = lalsim.SimInspiralCreateModeArray()
    for mode in modeList:
        
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
def phenomxhm_multipole( l, m, m1, m2, s1, s2, fmin=0.005, fmax=1.0, df = 5e-5,appx=None):
    
    #
    from numpy import array
    import numpy as np
    import lal, lalsimulation as lalsim
    from positive import physf,codehf,codef,sYlm
    
    #
    lalparams = lal.CreateDict()
    ModeArray = lalsim.SimInspiralCreateModeArray()
    
    #
    if appx is None:
        appx = lalsim.IMRPhenomXPHM
    
    #
    lalsim.SimInspiralModeArrayActivateMode(ModeArray, l, m)
    lalsim.SimInspiralWaveformParamsInsertModeArray(lalparams, ModeArray)
    
    #
    lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(
        lalparams, 0)
        
    #
    M_Sol = 100.0
    distance_Mpc = 100.0
    
    #
    m1_SI       = m1 * M_Sol * lal.MSUN_SI
    m2_SI       = m2 * M_Sol * lal.MSUN_SI
    f_min       = physf(fmin,M_Sol)
    f_max       = physf(fmax,M_Sol)
    delta_f     = physf(df,M_Sol)
    Omega       = 0.
    inclination = 1.3232  # Chosen so that it doesn't correxpond to a spherical harmonic root
    distance_SI = distance_Mpc*1.0e6*lal.PC_SI

    #
    Hp, Hc = lalsim.SimInspiralChooseFDWaveform(m1=m1_SI,
                                                m2=m2_SI,
                                                S1x=s1[0], S1y=s1[1], S1z=s1[2],
                                                S2x=s2[0], S2y=s2[1], S2z=s2[2],
                                                distance=distance_SI,
                                                inclination=inclination,
                                                LALpars=lalparams,
                                                phiRef=0,
                                                f_ref=f_min,
                                                deltaF=delta_f,
                                                f_min=f_min,
                                                f_max=f_max,
                                                longAscNodes=Omega,
                                                eccentricity=0.0,
                                                meanPerAno=0.0,
                                                approximant=appx)

    #
    freqs = np.arange(len(Hp.data.data)) * delta_f

    #
    s = -2
    spherical_harmonic_lm = sYlm(s, l, m, inclination, 0)
    hp = Hp.data.data / spherical_harmonic_lm
    hc = Hc.data.data / spherical_harmonic_lm

    #
    hp = codehf(hp, M_Sol, distance_Mpc)
    hc = codehf(hc, M_Sol, distance_Mpc)
    f = codef(freqs, M_Sol)

    #
    ans = array( [f, hp, hc] ).T

    #    
    return ans

#
def phenomxhm_multipoles(approximant, modeList, m1, m2, s1, s2, delta_f, phiRef, nu0=0):
    '''
    Generate dictionary of waveform arrays corresponding to input multipole list (i.e. list of [l,m] pairs ). If a single l,m pair is provided, then a single waveform array will be returned (i.e. we have opted to not have a lower-level function called "phenomxhm_multipole").
    
    USAGE
    ---
    output_modes_dict = phenomxhm_multipoles(approximant, modeList, m1, m2, s1, s2, delta_f, phiRef, nu0=0)
    
    NOTES
    ---
    IF modeList list of tuples, e.g. [(2,2),(2,1)], 
    THEN output_modes_dict is dictionary of [freq, hplus, hcross] in code units. 
    ELSE IF modeList is single list of ell and m, e.g. (2,2) or [2,2], 
    THEN a single array, i.e. [freq, hplus, hcross], is returned.
    '''
    
    # Import usefuls 
    from numpy import ndarray 
    
    # Validate inputs 
    # ---
    
    # Check type of modeList 
    if not isinstance(modeList,(list,tuple,ndarray)):
        error('modeList input must be iterable ')
    # If a single mode is given, make modeList a list of that mode's indices 
    single_mode_requested = False
    if len(modeList)==2:
        l,m = modeList
        if isinstance(l,int) and isinstance(m,int):
            single_mode_requested = True
            modeList = [ modeList ]
    # #
    # for lm in modeList:
    #     if isinstance(lm,int)

    # Main routine 
    # ---
    
    
    #
    lalparams = lal.CreateDict()

    #
    output_modes = {}

    #
    ModeArray = lalsim.SimInspiralCreateModeArray()
    for mode in modeList:

        #
        l, m = mode

        #
        lalsim.SimInspiralModeArrayActivateMode(ModeArray, l, m)
        lalsim.SimInspiralWaveformParamsInsertModeArray(lalparams, ModeArray)

        #
        lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(
            lalparams, threshold)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPFlag(lalparams, 1)
        lalsim.SimInspiralWaveformParamsInsertPhenomXCPfRing22Deviation(
            lalparams, nu0)

        # NOTE that all of the values below ore fiducial -- we ultimately want the waveforms in code units here.
        f_min = 10.0
        f_max = 2048.0
        Omega = 0.
        inclination = 1.3232  # Chosen so that it doesn't correxpond to a spherical harmonic's root
        distance_Mpc = 100.0
        distance = distance_Mpc*1.0e6*lal.PC_SI

        # Tell LAL to generate polarizations for the current mode 
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

        # Create a frequency array
        freqs = np.arange(len(Hp.data.data)) * delta_f

        # Divide by the relate spherical harmonic function
        s = -2
        spherical_harmonic = sYlm(s, l, m, inclination, 0)
        hp = Hp.data.data / spherical_harmonic
        hc = Hc.data.data / spherical_harmonic

        # Put in code units
        Mtot = m1+m2
        hp = codehf(hp, Mtot, distance_Mpc)
        hc = codehf(hc, Mtot, distance_Mpc)
        f = codef(freqs, Mtot)

        # Store to out dictionary 
        output_modes[l, m] = (hp, hc, f)

    #
    return output_modes



# Function to determine version2 data fitting region
def determine_data_fitting_region(data, fmin=0.03, fmax=0.12):
    '''
    Given version2 data array, determine fitting region dynamically.
    This function assumes clean data within the default or input values of fmin and fmax, and then uses the phase derivate to determine new fmin and fmax values that are ultimately used to define a fitting region.
    '''

    # Import usefuls
    from numpy import argmin, log
    from positive import smooth, find, lim
    from matplotlib.pyplot import figure, plot, show, axhline, xlim, ylim

    # Extract data
    f, amp_td, amp_fd, dphi_td, dphi_fd, phi_td, phi_fd = data

    # Use default domain bounds to determine a mask
    mask = (f >= fmin) & (f <= fmax)

    # Determine the minimum dphi
    # Smooth dphi using postiive's savgol filter
    x = log(f[mask])
    y = smooth(dphi_td[mask]).answer
    knot = argmin(y)
    y_knot = y[knot]
    data[3] = dphi_td - smooth(dphi_td[mask]).answer[knot] + y_knot
    data[4] = dphi_fd - smooth(dphi_fd[mask]).answer[knot] + y_knot

    # Determine new fmin and max using heuristic
    f_knot = f[mask][knot]
    new_fmin = max(f_knot * 0.22, 0.018)  # 0.5 # 0.325
    new_fmax = f_knot + 0.020  # 0.025

    #
    new_mask = (f >= new_fmin) & (f <= new_fmax)
    new_data = data.T[new_mask, :]

    #
    new_knot = find(f >= fmin)[0]+knot

    #
    return new_data, new_knot, new_fmin, new_fmax, f_knot


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
    from positive import eta2q,q2eta,eta2delta
    from glob import glob
    from pwca import determine_data_fitting_region,pwca_catalog,metadata_dict
    
    # Load and unpuack physical parameter space
    raw_domain = loadtxt(data_dir+'version2/fit_intial_binary_parameters.txt')
    theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2 = raw_domain.T


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
    degree_point = (theta*180/pi).astype(int)
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
    for _a1 in a1_set:
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
            ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color = 'tab:blue' if _a1==a1_set[0] else 'red' )
            
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
    set_fig_ax = set_fig_ax.flatten();
    tight_layout(4,4)
    ax_counter = 0

    #
    for _a1 in a1_set:
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
            ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color = 'tab:blue' if _a1==a1_set[0] else 'red' )
            
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
 
