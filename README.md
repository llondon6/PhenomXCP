# PhenomXCP
High-level development of, PhenomXCP, PhenomXAS and PhenomXHM with some feratures tuned to Numerical Relativity simulations. 

-- Lionel London, londonl@mit.edu, l.t.london@uva.nl, pilondon2@gmail.com

## Notes and Action Items (12/21)

* We proceed tuning the low-level ringdown frequency, thus allowing changes to propagate to various aspects of the PhenomXAS layer

* Both amplitude and phase ansatz depend on the ringdown frequency -- if this frequency is to be tuned, it must either be tuned for only one, or both simultaneously

- [ ] Decide on approach for tuning `fring` and `fdamp`

* 

## Notes and Action Items (11/21)

**General idea for tuning**: Use `XLALSimIMRPhenomXASGenerateFD` as primary python interface (see previous notes). It takes as input a `LALDict`, `lalparams`, which can be setup to store parameters for amplitude and phase deviations. 
 * The key file is `LALSimIMRPhenomX.c` (for the l=m=2 multipole)
 * `int IMRPhenomXASGenerateFD` is the key low level function which constructs the model
 * Lines `523-539` generate structs for amp and phase model coefficients. Of particular interest are the functions `IMRPhenomXGetAmplitudeCoefficients` and `IMRPhenomXGetPhaseCoefficients`
 * As was done with previous model iterations, changes to the **ringdown frequency and damping time** of the l=m=2 multipole can be implemented via `lib/LALSimIMRPhenomX_precession.c` around line `393`. In that setting, deviations must be funneled into `pWF`. **Note however**, that we may ultimately want to input the ringdown frequency deviations downstream such that e.g. the collocation points are not affected. 
 * There are **two routes** for modifying the ringdown frequency:
    - **Route A:** Modify the core `fRing` and `fDamp` parameters as defined below. This will change all downstream parameters related to RD frequencies and dampings for `l=m=2`
    ```python
    LALSimIMRPhenomX_precession.c:302:    printf("fring  (non-prec)  : %e\n",pWF->fRING);
    LALSimIMRPhenomX_precession.c:370:  pWF->fRING     = evaluate_QNMfit_fring22(pWF->afinal) / (pWF->Mfinal);
    LALSimIMRPhenomX_precession.c:376:    printf("fring  (prec)  : %e\n",pWF->fRING);
    LALSimIMRPhenomX_precession.c:391:  pWF->fRINGCP = pWF->fRING;  // Default value of effecting RD frequency
    LALSimIMRPhenomX_precession.c:395:    const REAL8 fRING22 = pWF->fRING;
    LALSimIMRPhenomX_precession.c:404:    pWF->fRINGEffShiftDividedByEmm = (1.0-cos(pWF->betaRD)) * ( fRING22  -  fRING21 );
    LALSimIMRPhenomX_precession.c:407:    // pWF->fRINGCP = fRING22 - mm * (1.0-cos(1.5+pWF->betaRD)) * ( fRING22  -  fRING21 );
    LALSimIMRPhenomX_precession.c:409:    pWF->fRINGCP = fRING22 - emm * pWF->fRINGEffShiftDividedByEmm;
    LALSimIMRPhenomX_precession.c:414:      // pWF->fRING = pWF->fRINGCP;  
    LALSimIMRPhenomX_precession.c:419:    pWF->fRING = pWF->fRINGCP; 
    LALSimIMRPhenomX_precession.c:424:      printf("fring22 (prec)       : %e\n",pWF->fRING);
    LALSimIMRPhenomX_precession.c:425:      printf("fring22 (coprec)     : %e\n",pWF->fRINGCP);
    ```
    - **Route B:** Modify select uses of the ringdown frequency. Collocation points will be determined by EZH's effective frequency, be lines below should be considered for modification. **Note** that the lines below are only approximate recommendations based on the output of `grep`. In particular, we may want to change fpeak but not change collocation points.
    ```python
    LALSimIMRPhenomX_intermediate.c:1338:  double frd = pWF->fRING; // The default behavior
    LALSimIMRPhenomX_intermediate.c:1342:   frd = pWF->fRINGCP;
    LALSimIMRPhenomX_intermediate.c:1402:  double frd = pWF->fRING; // The default behavior
    LALSimIMRPhenomX_intermediate.c:1406:   frd = pWF->fRINGCP;
    LALSimIMRPhenomX_internals.c:611:	pAmp->fAmpRDMin     = IMRPhenomX_Ringdown_Amp_22_PeakFrequency(pAmp->gamma2,pAmp->gamma3,pWF->fRING,pWF->fDAMP,pWF->IMRPhenomXRingdownAmpVersion);
    LALSimIMRPhenomX_internals.c:618:	pAmp->gamma1 = ( pAmp->v1RD / (pWF->fDAMP * pAmp->gamma3) ) * (F1*F1 - 2.0*F1*pWF->fRING + pWF->fRING*pWF->fRING + pWF->fDAMP*pWF->fDAMP*pAmp->gamma3*pAmp->gamma3)
    LALSimIMRPhenomX_internals.c:619:													* exp( ( (F1 - pWF->fRING) * pAmp->gamma2 ) / (pWF->fDAMP * pAmp->gamma3) );
    LALSimIMRPhenomX_ringdown.c:155:  // double dfr = ff - pWF->fRING;
    LALSimIMRPhenomX_ringdown.c:160:  double dfr = ff - pWF->fRING; // The default behavior
    LALSimIMRPhenomX_ringdown.c:164:    dfr = ff - pWF->fRINGCP;
    LALSimIMRPhenomX_ringdown.c:198:  // double frd   = pWF->fRING;
    LALSimIMRPhenomX_ringdown.c:203:  double frd = pWF->fRING; // `the default behavior
    LALSimIMRPhenomX_ringdown.c:207:    frd = pWF->fRINGCP;
    LALSimIMRPhenomX_ringdown.c:472:  // double frd     = pWF->fRING;
    LALSimIMRPhenomX_ringdown.c:477:  double frd = pWF->fRING; // The default behavior
    LALSimIMRPhenomX_ringdown.c:481:    frd = pWF->fRINGCP;
    LALSimIMRPhenomX_ringdown.c:521:  // double frd      = pWF->fRING;
    LALSimIMRPhenomX_ringdown.c:526:  double frd = pWF->fRING; // The default behavior
    LALSimIMRPhenomX_ringdown.c:530:    frd = pWF->fRINGCP;
    ```
* **Strategy**: Route A is simpler, so we will try this first. Under this route, here are the parameters to be modified:
  * Ringdown frequency and damping (LALSimIMRPhenomX_precession.c ~ Line 400). The EZH effective frequency will be directly modified. 
  *  

- [x] Identify key model parameters to change, and where to change them.
- [x] Write tuning appropriate python wrapper for `IMRPhenomXASGenerateFD`: given system masses and spins, output function parameterized only be deviations from base model
- [x] Initiate basic tuning tests with standard optimization routines
 
## Notes and Action Items (8/21) 

* Move waveform flags for CP to IMRPhenomXHM_SetHMWaveformVariables

- [x] Choose c file in which to add functions for modified ringdown frequencies 
- [x] Add functions for modified (i.e. effective a la EZH) ringdown frequencies 
- [x] Make PhenomX option that incorporates effective ringdown frequencies 
- [x] Make python wrapper for PhenomX+options that results in PhenomXCPv0 = PhenomXAS+EffectiveRingdownFreqs. PhenomXCPv0 may be defined by separate amplitude and phase functions which require beta as input.
- [x] Compare PhenomXCPv0 to NR simulations

## Notes and Questions

* Plan for version 0 -- IMRPhenomXHM_SetHMWaveformVariables

* NOTE that `XPHM` calls `XAS` via `IMRPhenomXASGenerateFD` in `LALSimIMRPhenomXPHM.c`

* **NOTE:** Versions of PhenomXCP will only be available through XLAL functions, as it is expected that they will not be of use for e.g. parameter estimation through SimInspiralChooseFDWaveform. As a result, XCP model version will, in some sense, not be proper LAL approximants. Let's say that they will be "meta approximants". 

## Related code resources 

* Fork of LALSuite -- https://git.ligo.org/lionel.london/lalsuite/-/tree/hack-phenomx
* Jonathan's LALSuite fork -- https://git.ligo.org/jonathan.thompson/lalsuite/-/blob/o773b/lalsimulation/lib/LALSimIMRPhenomX_PNR_internals.c

## Related pages

 * LALSuite fork where development C code for PhenomXCP lives: https://git.ligo.org/lionel.london/lalsuite
 * Review of PhenomXHM: https://git.ligo.org/waveforms/reviews/imrphenomx/-/wikis/home
 * PhenomXAS paper: https://arxiv.org/pdf/2001.11412.pdf
 * PhenomXHM paper: https://arxiv.org/pdf/2001.10914.pdf
 
 
 ## Outline of / Plan for notebooks and scripts
 
 Most items will have an ipython notebook and python script. The notebooks will be used to develop the python scripts. The scripts will be executed in series or alone. Together, all scripts (really 2-8) in order generate the tuned model from scratch. Script 1 is a testbed for a version of PhenomX that has an effective ringdown frequency given by Eleanor's model.
 
* *0_preliminaries* -- basics package setup and testing 
 
* *1_PhenomXCPv0* -- LAL version of PhenomXCP that only has EZH effective ringdown frequencies 

* *2_collect_nr_data*

* *3_collect_metadata*

* *4_determine_fitting_regions*

* *5_fit_waveforms*

* *6_fit_hyperparameters*

* *7_document_fits*

* *8_matches*
