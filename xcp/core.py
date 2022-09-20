

#
import pickle
from positive import alert, magenta, parent, red, smart_object, blue
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

#
config_path = package_dir + 'xcp/global_config.ini'
global_config = gc = smart_object( config_path, verbose=True )
gc.lmlist = eval(','.join(gc.lmlist))

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
alert('According to the global config, the XCP package is configured to model the %s coprecessing moment multipole moment(s)'%magenta(str(gc.lmlist)),fname='xcp.core')


def get_xphm_coprec(ell, emm, Mtotal, q, chi1, chi2, pnr=False):
    
    #
    import lal 
    import lalsimulation as lalsim 
    import numpy as np

    m1 = Mtotal * q / (1.0 + q)
    m2 = Mtotal / (1.0 + q)

    # generate co-precessing waveform
    lalparams = lal.CreateDict( )
    lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(lalparams, 0)
    
    lalsim.SimInspiralWaveformParamsInsertPhenomXPHMPrecModes(lalparams, 1)
    
    if pnr:
        lalsim.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalparams, 500)
    
    #
    ModeArray = lalsim.SimInspiralCreateModeArray()
    lalsim.SimInspiralModeArrayActivateMode(ModeArray, ell,emm)
    lalsim.SimInspiralWaveformParamsInsertModeArray(lalparams, ModeArray)

    distance = 1e6*lal.PC_SI 
    phiRef = 0.0 
    inclination = 0.0

    delta_F=0.125 
    f_lower=20.
    fRef=20.
    f_max = 2**(int(np.log(3000./delta_F)/np.log(2))+1) * delta_F

    hlmpos, hlmneg = lalsim.SimIMRPhenomXPHMOneMode(    
        ell,
        emm,
        lal.MSUN_SI*m1,
        lal.MSUN_SI*m2, 
        chi1[0],chi1[1],chi1[2],
        chi2[0],chi2[1],chi2[2],                                                  
        distance,
        inclination,
        phiRef,
        delta_F,
        f_lower,
        f_max,
        fRef,
        lalparams)

    freqs = np.arange(len(hlmpos.data.data)) * delta_F

    return freqs, hlmpos.data.data, hlmneg.data.data

#
def get_phenomxphm_coprecessing_multipoles(freqs, lmlist, m1, m2, s1, s2, phiRef=0, pflag=0, fsflag=None, mu1=0, mu2=0, mu3=0, mu4=0, nu0=0, nu4=0, nu5=0, nu6=0, zeta1=0, zeta2=0,__set_XPHMThresholdMband__=True, option_shorthand=None ):
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
    
    We will want this function to do three things:
    
    1. Output PhenomXPHM for comparisons
    2. Output PhenomXHM/AS with/without deviations for calibration and comparison 
    3. Output PhenomXCP/PNR 
    
    
    
    '''
    
    # Import usefuls 
    from numpy import ndarray, array, arange, double
    from positive.units import codef,codeh,codehf,physf
    from positive import sYlm, error,warning
    import lalsimulation as lalsim
    import lal
    
    #
    #
    if len(lmlist)>1:
        error('this function only works with one mode input! otherwise there is an issue')
    #
    lm = lmlist[0]
    ll,mm = lm
            
    #
    option_shorthand_strings =  ('1-pnr','2-xphm','3-xphm-ezh','4-xhm')
    
    #
    if option_shorthand:
        
        #
        if isinstance(option_shorthand,str):
            option_shorthand = option_shorthand.lower()
        else:
            error('option shorthand key value must be string')
        #
        if not (option_shorthand in option_shorthand_strings):
            error('option shorthand key value should be in %s'%magenta(option_shorthand_strings))
            
        #
        # alert('Valid option-shorthand input found. We will now overwrite Phenom-X options according to %s'%option_shorthand)
        
        if   option_shorthand == '1-pnr':
            
            #
            pflag = 500 # Use XCP which is XAS or XHM with select deviations in model parameters
            fsflag = 0 # Use default behavior for XPHM
            
        elif option_shorthand == '2-xphm':
            
            #
            pflag = 0 # Use default PhenomXPHM behavior
            fsflag = 0 # Use default behavior for XPHM
            
        elif option_shorthand == '3-xphm-ezh':
            
            #
            pflag = 501  # Use EZH's effective ringdown frequency for dominant multipole moments
            fsflag = 0 # Use default behavior for XPHM
            # alert('(l,m) = (%i,%i), pflag = %i'%(ll,mm,pflag))
            
        elif option_shorthand == '4-xhm':
            
            # NOTE see use of fsflag in LALSimIMRPhenomX_precession.c
            fsflag = 5 # Tells XP interface to use non-spinning fit for remnant BH
            pflag = 0 # Use default PhenomXPHM behavior
            # alert('(l,m) = (%i,%i), pflag = %i, fsflag = %i'%(ll,mm,pflag,fsflag))
            
        else:
            
            #
            error('unhandled option shorthand')
            
    else:
        
        error('It is very unwise to not use this functions option_shorthand keyword input. NOTE that this can be done through template_amp_phase\'s kwargs interface')
        
    # alert('option_shorthand = %s'%option_shorthand)
    
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
    
    # #
    # if pflag is None:
    #     error('Precession version flag, or "pflag", must be input. To return PhenomX with default settings, please use pflag=0. To return PhenomXCP please use pflag=500.')
    
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
    for lm in lmlist:
        
        #
        l,m = lm
        
        # alert('(l,m) = (%i,%i), pflag = %i, fsflag = %i'%(l,m,pflag,fsflag))
        
        #
        lalparams = lal.CreateDict()
        ModeArray = lalsim.SimInspiralCreateModeArray()
        lalsim.SimInspiralModeArrayActivateMode(ModeArray, l,m)
        lalsim.SimInspiralWaveformParamsInsertModeArray(lalparams, ModeArray)

        # Turn off multibanding
        lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(lalparams, 0)
        if __set_XPHMThresholdMband__:
            lalsim.SimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(lalparams, 0)

        # Tell the model to return the coprecessing mode -- only works on our development branches
        lalsim.SimInspiralWaveformParamsInsertPhenomXReturnCoPrec(lalparams, 1)
        
        #
        
        # Set deviations from base model based on inputs
        # mu2,mu3,mu4,nu4,nu5,nu6,zeta2 = [ double(k) for k in (mu2,mu3,mu4,nu4,nu5,nu6,zeta2) ]
        # print('3>> ',*(mu1, mu2, mu3, mu4, nu4, nu5, nu6, zeta1, zeta2),'\n')
        if (l,m) == (2,2):
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU1(lalparams, mu1)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU2(lalparams, mu2)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU3(lalparams, mu3)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU4(lalparams, mu4)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU0(lalparams, nu0)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU4(lalparams, nu4)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU5(lalparams, nu5)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU6(lalparams, nu6)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPZETA1(lalparams, zeta1)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPZETA2(lalparams, zeta2)
        elif (l,m) == (3,3):
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU1l3m3(lalparams, mu1)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU2l3m3(lalparams, mu2)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU3l3m3(lalparams, mu3)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPMU4l3m3(lalparams, mu4)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU0l3m3(lalparams, nu0)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU4l3m3(lalparams, nu4)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU5l3m3(lalparams, nu5)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPNU6l3m3(lalparams, nu6)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPZETA1l3m3(lalparams, zeta1)
            lalsim.SimInspiralWaveformParamsInsertPhenomXCPZETA2l3m3(lalparams, zeta2)
        else:
            error('(%i,%i) unhandled for deviations relative to PhenomXHM'%(l,m))
        
        #
        '''
        The default value in this function is pflag=0 (see function def above). This value does not equate to True, so the true default behavior is set in LALSuite; i.e. the LAL code sets its own default value for the pflag which is 300 (see LALSimInspiralWaveformParams.c at DEFINE_ISDEFAULT_FUNC(PhenomXPrecVersion, INT4, "PrecVersion", 300). NOTE that pflag is then set to a new default value of 223 in LALSimIMRPhenomX_precession.c.   )
        '''
        if pflag:
            lalsim.SimInspiralWaveformParamsInsertPhenomXPrecVersion( lalparams, pflag )
        #     alert('pflag is effectively True and set to %i'%pflag)
        # else:
        #     alert('pflag is to the effect of False, so pflag will simply not be set in the lasagna (aka laldict)')
        
        '''Only set the final spin flag is the user desires a non-trivial value'''
        if fsflag:
            lalsim.SimInspiralWaveformParamsInsertPhenomXPFinalSpinMod( lalparams, fsflag )

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
def template_amp_phase(m1, m2, chi1_vec, chi2_vec, lm=(2,2),include_nu0=False,floor_dphi=True,**kwargs):
    
    # NOTE that mu4 is no longer to be used as it is completely degenerate with nu5 in PhenomX
    
    #
    import xcp
    from numpy import unwrap,angle,mean
    from positive import spline_diff,warning
    from numpy import ndarray
    
    #
    ell,emm = lm
    lmlist = [ (ell,emm) ]
    
    #
    def template_together_helper( f, mu1=0, mu2=0, mu3=0, mu4=0, nu4=0, nu5=0, nu6=0, zeta1=0, zeta2=0, nu0=0 ):
        
        # Calculate PhenomXPHM with the input deviations
        # NOTE that pflag=0 means that we use the default setting of PhenomXPHM as a reference model. NOTE that we try and except here becuase sometimes the optimization routines can stray outside of the accepted model domain thus causing LAL to throw an error
        
        # alert(nu0)
        if isinstance(nu0,(list,tuple,ndarray)):
            nu0 = nu0[0]
            
        # multipole_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, lmlist, m1, m2, chi1_vec, chi2_vec, pflag=0, mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4, nu0=nu0, nu4=nu4, nu5=nu5, nu6=nu6, zeta1=zeta1, zeta2=zeta2,**kwargs )
        
        try:
            # print('2>> ',*(mu1, mu2, mu3, mu4, nu0, nu4, nu5, nu6, zeta1, zeta2))
            multipole_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, lmlist, m1, m2, chi1_vec, chi2_vec, pflag=0, mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4, nu0=nu0, nu4=nu4, nu5=nu5, nu6=nu6, zeta1=zeta1, zeta2=zeta2,**kwargs )
        except:
            print('n0 = ',nu0)
            warning('Something went wrong with the standard evaluation:')
            multipole_dict = xcp.get_phenomxphm_coprecessing_multipoles( f, lmlist, m1, m2, chi1_vec, chi2_vec, pflag=0,**kwargs )
        
        # 
        complex_strain = multipole_dict[ell,ell]
            
        # Given the complex FD waveform, compute its amplitude
        amplitude = abs(complex_strain)
        # Given the complex FD waveform, compute its phase derivative
        complex_strain = multipole_dict[ell,ell]
        phase = unwrap( angle(complex_strain) )
        phase_derivative = spline_diff(f,phase,k=5)
        
        # Find min phase derivative to subtract off
        
        # determine the desired value  
        mask = (f>0.03*ell/2)&(f<0.12*ell/2)
        min_phase_derivative = min( phase_derivative[ mask ] )
        
        # Adjust phase derivative 
        if floor_dphi:
            phase_derivative -= min_phase_derivative 
        
        #
        return amplitude,phase_derivative
        
    # We need to acknowledge here that the (2,2) moment is not sensitive to mu4 (and maybe other inputs). We do this by prototyping the waveform function to explicitely ignore mu4 as an input 
    
    if not include_nu0:
        if ell == 2:
            # no mu4 
            template_together = lambda f, mu1=0, mu2=0, mu3=0, nu4=0, nu5=0, nu6=0, zeta1=0, zeta2=0: template_together_helper( f, mu1=mu1, mu2=mu2, mu3=mu3, mu4=0, nu0=0, nu4=nu4, nu5=nu5, nu6=nu6, zeta1=zeta1, zeta2=zeta2 )
        else:
            # 
            template_together = lambda f, mu1=0, mu2=0, mu3=0, mu4=0, nu4=0, nu5=0, nu6=0, zeta1=0, zeta2=0: template_together_helper( f, mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4, nu0=0, nu4=nu4, nu5=nu5, nu6=nu6, zeta1=zeta1, zeta2=zeta2 )
    else:
        if ell == 2:
            # no mu4 
            template_together = lambda f, mu1=0, mu2=0, mu3=0, nu4=0, nu5=0, nu6=0, zeta1=0, zeta2=0, nu0=0: template_together_helper( f, mu1=mu1, mu2=mu2, mu3=mu3, mu4=0, nu0=nu0, nu4=nu4, nu5=nu5, nu6=nu6, zeta1=zeta1, zeta2=zeta2 )
        else:
            # 
            template_together = lambda f, mu1=0, mu2=0, mu3=0, mu4=0, nu4=0, nu5=0, nu6=0, zeta1=0, zeta2=0, nu0=0: template_together_helper( f, mu1=mu1, mu2=mu2, mu3=mu3, mu4=mu4, nu0=nu0, nu4=nu4, nu5=nu5, nu6=nu6, zeta1=zeta1, zeta2=zeta2 )
        
        
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
def determine_data_fitting_region_legacy_2(raw_data_array, threshold=0.015, lm=(2,2), plot=False, simname=None, floor_dphi=True, f_lim = None, smooth_dphi = False, fring=None ):
    '''
    Given version2 data array, determine fitting region dynamically.
    This function assumes clean data within the default or input values of fmin and fmax, and then uses the phase derivate to determine new fmin and fmax values that are ultimately used to define a fitting region.
    '''

    # Import usefuls
    from numpy import argmin, log, arange
    from positive import smooth, find, lim, smoothest_part_by_threshold,findpeaks
    if plot:
        from matplotlib.pyplot import figure, plot, show, axhline, xlim, ylim

    # DETERMINE DATA FITTING REGION
    # ---
    
    #
    l,m = lm
    
    # 0. Select and unpack
    f,amp_fd,dphi_fd,alpha,beta,gamma = data
    
    #
    if smooth_dphi:
        error('dont do this; not fully tested')
        dphi_fd = smooth(dphi_fd,width=20).answer
    
    #
    if f_lim is None:
        
        # 1. Find smoothest part
        pre_mask = (f>0) #& (f<0.12)
        mask_1     = smoothest_part_by_threshold( dphi_fd[pre_mask], threshold=threshold, smooth_width=20, plot=False )
        dphi_fd_1  = smooth( dphi_fd[pre_mask][mask_1], width=30 ).answer
        f_1 = f[pre_mask][mask_1]
        
        
        #
        shift_1 = 0
        if l==3:
            if simname in 'q4a08t30dPm5p5dRm47_T_96_360':
                shift_1 = 0.015
            if simname in 'q8a08t30dPm9':
                shift_1 = 0.015
            if simname in 'q8a08t60Ditm45dr075_96_360':
                shift_1 = -0.008
            if simname in 'q2_a10_a28_ph0_th30':
                shift_1 = 0.008
        
        # 2. Handle unstable end behavior and remask
        peaks,peak_locations = findpeaks( dphi_fd_1 )
        using_mask_1 = False
        if len(peak_locations):
            mask_2    = range(0,peak_locations[-1]+1)
            dphi_fd_2 = dphi_fd_1[mask_2]
            f_2 = f_1[mask_2]
            mask_3 = f_2 < (min(0.12*l*0.5,0.152) + shift_1)
            dphi_fd_3 = dphi_fd_2[ mask_3 ]
        else:
            # alert('Using mask_1')
            using_mask_1 = True
            mask_2 = mask_1
            dphi_fd_3 = dphi_fd_1
        
        # 3. Determine location of the lorentzian min
        try:
            if using_mask_1:
                lorentzian_mindex = mask_1[argmin( dphi_fd_3 )]
            else:
                lorentzian_mindex = mask_1[mask_2[argmin( dphi_fd_3 )]]
        except:
            from matplotlib.pyplot import figure, plot, show, axhline, xlim, ylim
            from positive import error
            figure()
            #print('>> ',len(mask_2),mask_2[argmin( dphi_fd_3 )])
            plot( dphi_fd_3 )
            plot( argmin( dphi_fd_3 ), dphi_fd_3[argmin( dphi_fd_3 )], marker = 'o', ms=8 )
            show()
            error('something went wrong')
        dphi_lorentzian_min = min( dphi_fd_3 )
        
        # 4. Use the lorentzian_mindex to define the start and end of the fitting region
        f_lorentzian_min = f[pre_mask][ lorentzian_mindex ]
        # f_lorentzian_min = f[mask_1][mask_2][ lorentzian_mindex ]
        f_min = f_lorentzian_min * 0.2 
        # f_min = max(f_lorentzian_min * 0.22, 0.018) 
        f_max = f_lorentzian_min + 0.012# 0.018 # 0.012
        calibration_mask = arange(len(f))[(f>=f_min) & (f<=f_max)]
    
    else:
        
        #
        f_min,f_max = f_lim    
        calibration_mask = arange(len(f))[(f>=f_min) & (f<=f_max)]
        smoothed_dphi = smooth(dphi_fd,width=30).answer
        dphi_lorentzian_min = min(smoothed_dphi)
        f_lorentzian_min = f[ argmin( smoothed_dphi ) ]
    
    # 5. Select region for output 
    calibration_data = data.T[ calibration_mask ]
    # 6. Smooth phase derivative data and subtract min
    calibration_data.T[2] = dphi_fd[calibration_mask] - ( dphi_lorentzian_min if floor_dphi else 0 )
    # calibration_data.T[2] = smooth( dphi_fd, width=10 ).answer[calibration_mask] - dphi_lorentzian_min
    # calibration_data.T[2] = smooth( calibration_data.T[2], width=30 ).answer
    

    #
    return calibration_data, dphi_lorentzian_min, f_min, f_max, f_lorentzian_min


# Function to determine version2 data fitting region
def determine_data_fitting_region(raw_data_array, fring, lm=(2,2), floor_dphi=True, plot=False, simname=None, ax=None ):
    '''
    Given version4 data array, determine fitting region dynamically.
    This function assumes clean data within the default or input values of fmin and fmax, and then uses the phase derivate to determine new fmin and fmax values that are ultimately used to define a fitting region.
    '''

    # Import usefuls
    from numpy import argmin, log, arange, unwrap, mean
    from positive import smooth, find, lim, smoothest_part_by_threshold,findpeaks,limy
    if plot:
        from matplotlib.pyplot import plot,figure,figaspect,axvline,axhline,xlim,ylim,xscale,show,sca

    #
    ll,mm = lm
    
    #
    if simname is None: simname = ''
    
    #
    f,amp,raw_dphi,alpha,beta,gamma = raw_data_array
    
    #
    dphi_shifted = unwrap(raw_dphi,discont=4000)
    
    shift_mask = (f>0.014*ll*0.5) & (f<fring)
    dphi_reshifted = dphi_shifted - mean(dphi_shifted[shift_mask]-raw_dphi[shift_mask])
    
    # --- Case specific toggles --- #
    #
    smooth_width = 40
    upscale = 1.2
    downscale = 0.014
    out_smoothed_data = False
    #
    if ll==2:
        #
        if 'q1a08t60' in simname:
            smooth_width = 60
        if 'q1a06t30dPm35' in simname:
            smooth_width = 10
            downscale = 0.019
        if 'q2_a10_a28_ph0_th60' in simname:
            upscale = 1.1
        if 'q1a08t60' in simname:
            smooth_width = 90
        if 'q1a06t120' in simname:
            smooth_width = 20
            downscale = 0.019
    if ll==3:
        if 'q4a08t60dPm3dRm250' in simname:
            upscale = 1.05
    #
    # ----------------------------- #
    
    #
    dphi = dphi_reshifted
    s_dphi = smooth(dphi,width=smooth_width).answer

    #
    f0,f1 = downscale*ll*0.5, fring*upscale
    mask_1 = (f>f0) & (f<f1)
    
    #
    dphi_lorentzian_min = min(s_dphi[mask_1])
    f_lorentzian_min    = f[mask_1][ argmin(s_dphi[mask_1]) ]
    
    #
    f_min = f_lorentzian_min*0.2
    f_max = f_lorentzian_min+0.012
    
    #
    calibration_mask =  (f>f_min) & (f<f_max)
    
    #
    if plot:
        
        #
        if ax is None:
            figure( figsize=2*figaspect(0.618) )
        else:
            sca(ax)
        plot( f, dphi )
        plot( f, s_dphi, color='k', ls='-', lw=2, alpha=0.2 )
        plot( f_lorentzian_min, dphi_lorentzian_min, marker='o', ms=8, mfc='none', color='k', zorder=100, mew=2 )
        
        axvline( f0, c='g',ls='-')
        axvline( f1, c='g',ls='-')
        
        axvline( f_min, c='r',ls='--')
        axvline( f_max, c='r',ls='--')
        
        axvline( fring, c='k',ls='--')
        
        xlim( f_min*0.9, f_max/0.9)
        ylim( limy(f[calibration_mask],dphi[calibration_mask],dilate=0.1) )
        
        xscale('log')
        
        if ax is None: show()
    
    #
    if out_smoothed_data:
        calibration_dphi = s_dphi[ calibration_mask ]
    else:
        calibration_dphi = dphi[ calibration_mask ]
    
    #
    offset =  dphi_lorentzian_min if floor_dphi else 0 
    calibration_dphi -= offset
    # if ll==3: alert(offset)
    
    # Select region for output 
    calibration_data = raw_data_array.T[ calibration_mask ]
    # Store final dphi
    calibration_data.T[2] = calibration_dphi
    

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
def advanced_gmvx_plot( fit_object, ll, mm ):
    
    '''
    Advanced gloss atop mvpolyfit.plot and mvrfit.plot
    '''
    
    from matplotlib.pyplot import subplots, plot, xlabel, ylabel, title, sca, gca, figaspect, tight_layout
    from numpy import cos,sin,array,around,ones_like,sort,pi,linspace,mod
    from positive import eta2q,q2eta,eta2delta,rgb
    from glob import glob
    from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict
    
    # Load and unpack physical parameter space
    first_raw_domain = loadtxt(data_dir+'version2/fit_initial_binary_parameters.txt')
    theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_ = first_raw_domain.T
    
    #
    if mod(mm,2) == 0:
        raw_domain = first_raw_domain
    else:
        mask = around( m1/m2 ) != 1
        raw_domain = first_raw_domain[ mask, : ]
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
    set_fig_ax = set_fig_ax.flatten()
    tight_layout(w_pad=4,h_pad=4)
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
            case_eta   = linspace( min(_eta)*0.6,max(_eta),1000 ) 
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
    tight_layout(w_pad=4,h_pad=4)
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
            # case_theta   = linspace( min(_theta),max(_theta),1000 ) # 
            case_theta   = linspace( 0,pi,1500 ) # 
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
    num_figs = len(theta_set)*len(eta_set)
    a1_set_figs,set_fig_ax = subplots( len(theta_set), len(eta_set), figsize=5*array([ len(eta_set),len(theta_set) ]) )
    set_fig_ax = set_fig_ax.flatten()
    tight_layout(w_pad=4,h_pad=4)
    ax_counter = 0

    #
    for k,_theta in enumerate(theta_set):
        for _eta in eta_set:

            #
            eta_mask = (_eta==eta_point)
            theta_mask = (_theta==theta_point)
            mask = theta_mask & eta_mask

            #
            _a1 = a1_point[mask]

            #
            case_a1   = linspace( 0,1,100 ) # 
            case_u     = cos(_theta) * ones_like(case_a1) # cos(case_theta)
            case_eta   = _eta * ones_like(case_a1)
            case_delta = eta2delta( case_eta )

            #
            case_domain = array([case_u,case_eta,case_delta,case_a1]).T
            case_range = fit_object.eval(case_domain)
            opt_range  = fit_object.eval(fit_object.domain[mask,:])

            #
            # sca(ax[0])
            # ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color='k', ls='--' ) # colors[k]
            # ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color = 'tab:blue' if _a1==a1_set[0] else 'red' )
            
            #
            sca( set_fig_ax[ax_counter] ); ax_counter += 1
            plot( a1[mask], fit_object.range[mask] if hasattr(fit_object,'range') else fit_object.scalar_range[mask], marker='o',ls='none',color='k'  )
            plot( a1[mask], opt_range, marker='o',ms=10,mfc='none', color='forestgreen',ls='none'  )
            plot( case_a1, case_range, ls='-', color='forestgreen' )
            title( r'$\theta_{\mathrm{LS}}=%1.0f$, $q=%1.2f$'%(_theta*180.0/pi,eta2q(_eta)) )
            xlabel(r'$a_1$')
            ylabel(r'$\%s$'%fit_object.labels['python'][0])
            
    #
    return summary_fig, eta_set_figs, theta_set_figs, a1_set_figs
 

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