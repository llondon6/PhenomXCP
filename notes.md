
# Update Feb. 2023

**1.** Cecilio has confirmed that the low frequency relative phases of `XPHM` differ from those of `XHM`. He suggested that a time-shift alignment parameter be lowered --  **I tried this, but it did not work.** 

**2.** ***Questions***:
- [x] Does my comparison of low frequency behavior depend incorrectly on in-plane spin (somehow)? **Would be nice to clarify, but ultimately a moot point. See below.**
- [x] Is it viable to **force** a phase and time-shift of PNR-XCP such that the PhenomXHM low frequency phase behavior is recovered?  **Others agree. Try it.**
- [ ] Is it necessary and viable to restructure XCP so that the problematic dependence on the ringdown frequency is avoided. I suspect that this would still require a post-fact alignment step since the phase's morphology changes amid heuristic alignment frequencies.  **If forcing alignment with XHM works, then no.**

**3.** ***Thoughts on forcing alignment with the `XHM` phase & phase derivative***:
* This would ostensibly require computing the `XHM` phase and phase derivative at some reference frequency. 
```C
\* The XCP phase derivatives will be set to the respective XHM values at
this frequency. Then, the values of the XCP phases will be aligned with
respective XHM values at this frequency.*\
  
f_inpiral_align = 0.001;
```
* For that, one would need to generate the `XHM` parameter structure and then store it, or pre-compute the related quantities, and them apply them when needed. **Let's use the strcuts in-place**, and try to pre-compute the parameters needed for the desired alignment. 
* Since computation of the `XHM` structure is computationally cheap, this need not be a "no-go" step? Most likely not. The code already does this many times.

**4.** ***Implementation notes for forced phase alignment***

* Development will happen on branch `pnr-review-candidate-ll-dev-2`
* Items needed:
  - Routines to compute single frequency phase and phase derivative values for `XHM` without PNR deviations
  - The `XHM` waveforms, amplitude and phase structs need to be generated and stored upon model initialization. They should be stored in the `XCP22` waveforms struct or within the `XCPHM` waveforms struct 
  
**5.** ***Implementation Action Items***

- [ ] Review and make short notes on how final `XHM` phase quantities are computed: ***Which structs are needed? Which functions?***
  * Can `SimIMRPhenomXHMPhase` be used to get the absolute (*i.e.* not relative) value of the phase?  Perhaps not, and instead something like `XLALSimIMRPhenomXHMFrequencySequenceOneMode` is needed?
  * Is there an equivalent function for the phase derivative?
- [ ] Add fields to `XHM` and `XAS` structs the are inputs / data-holders for forced phase alignment.
  * Define a phase alignment frequency: `f_inpiral_align`
- [ ] Generate `XHM` $(2,2)$ and HM structs upon initialization of `XCP` structs. Do this separately for $(2,2)$ and HM.

**6.** ***Useful Notes / Other Action Items***

* For generating waveforms at a **single frequency**: `XLALSimIMRPhenomXHMFrequencySequenceOneMode` uses `XLALSimIMRPhenomXASFrequencySequence` and `IMRPhenomXHMGenerateFDOneMode`. **But these are likely not the functions / perspectives that one wants to use.** 
* In `lib/LALSimIMRPhenomXHM_internals.c`, the `XAS` phase is constructed at a single frequency point
```
pWFHM->phiref22 = -1./pWF22->eta*IMRPhenomX_Phase_22(pWF22->MfRef, 
&powers_of_MfRef, pPhase22, pWF22) - pWFHM->timeshift*pWF22->MfRef
 - pWFHM->phaseshift + 2.0*pWF22->phi0 + LAL_PI_4;
```
Here, note that `pWFHM->phiref22` is the $(2,2)$ moment's phase at the reference frequency `pWF22->MfRef`. I should be able to use the same logic to set an absolute value of the `XCP` phase at some low frequency.
- [ ] For $(2,2)$, reconstruct the absolute value of the $(2,2)$ phase in the same way that `pWFHM->phiref22` is computed. Zero the `XCP` phase at some low reference frequency, and then add the $(2,2)$ phase from `XAS`.
* Assuming that works, what should one do for the higher moments?
  * **Note** that `IMRPhenomXHM_Phase_noModeMixing` takes float input, so by construction it may be called at a single frequency point. So it is fit for purpose.
  * When `IMRPhenomXHM_Phase_noModeMixing` is called in `IMRPhenomXHMGenerateFDOneMode`, a factor of $(-1)^\ell$ is applied to the amplitude. This has the effect of adding $\ell \,\pi$ to the phase. The same is effecgtively done in `IMRPhenomXHM_MultiMode2`. So the factor should be added when `XCP` is aligned with `XHM`.
- [ ] Don't forget to add $\ell \pi$ to the overall phase of the `XHM` moment when forcing `XCP` to have its value.
* Are there any other modifications made to the `XHM` phases that must be taken into account? **The answer appears to be NO up to how the moments are added in `IMRPhenomXHMFDAddMode`**.
* Ok. Now what about the phase derivative??
* According to `IMRPhenomXHM_Phase_noModeMixing` the phase derivative during inspiral should be
```C
DPhiIns + pPhase->C1INSP 
``` 
where `DPhiIns` is an evaluation of `IMRPhenomXHM_Inspiral_Phase_Ansatz`
- [ ] Given the points above, show that you can set the `XCP` phase and phase derivative to zero at a chosen alignment frequency. Plot for $(2,2)$ and $(3,3)$ **(Fig. 1)**. This and following figures have 4 panels.
- [ ] Add `XHM` to the plot **(Fig. 2)**.
- [ ] Add `XHM` to the plot **(Fig. 3)**.
- [ ] Align `XCP` with `XHM` **(Fig. 4)**.
- [ ] Convince yourself that you are done.

# Summary of my annoyances with the PhenomXHM phase

**1.** Their clever way of imposing relative phases (`4.13` of the [XHM paper](https://arxiv.org/pdf/2001.10914.pdf)) is superfluous: it's entirely tetrad dependent, and is ultimately indistinguishable from changes in extrinsic variables -- ***polarization angle and orbital phase***

**2.** The XHM paper says that they impose the above condition in the limit of zero frequency, but really they impose this at `0.6` times the `MECO` frequency, which is high enough to *not* be in the regime where the relative phases are constant -- this means that their attempt to impose a tetrad convention on the relative phase is sensitive to the frequency dependence (e.g shape) of the phase in late inspiral. So even if `4.13` is correct, XHM is not actually guaranteed to set the desired condition on the relative phase. This is particularly true for `XPHM` (and so `PNR`), but not so true for `XHM`

  * **Possible action:**  When the tuned coprecessing model is active, use a different fraction of MECO. Perhaps `falign = 0.4 * fMECO`?

**3.** I'm not sure that `4.13` as stated is correct. If one re-derives this equation according to the description in their paper, then it appears that one cannot ignore the overall relative phases from PN (i.e. at zero frequency, the `PN` moments have non-zero relative phase). It's true that my statement in parenthesis is also tetrad dependent, and it's not obvious to me that (as e.g. Mishra et al have argued in private correspondence) that there exists a tetrad for which the PN amplitudes have no relative phase at leading order.

  * **Possible action:** Since the tetrad convention appears to be degenerate with changes in the 

**4.** Okay, the zero frequency limit is perhaps not the correct limit in which to think about this to begin within -- the zero frequency limit is `0 PN`, and we know that most of the interesting PN effect turn on later. In particular, as evident in the XHM paper's appendix, PN amplitudes pick up complex terms at `2PN`, meaning that `2PN` effects impact both relative phase and time-shift.

**5.** Regarding the coprecessing model, I've concluded that PhenomXHM is actually behaving *exactly as it is coded to*, and that what it's coded to do is closer to incorrect than correct. I don't feel very comfortable coming out (to the XHM devs) and saying the XHM is incorrect, so I am ruminating on the evidence, until it is plainly undeniable (this in itself if adjacent to "fixing" `XHM`). Doing a proper comparison between PN, NR and the models is at least as hard as making hybrids -- likely harder since hybrid construction can be simplified by assuming that NR is the most correct of all possibilities. 

# Getting phase directly from the XPHM waveform generator

* This will be done with `SimIMRPhenomXPHMFrequencySequence` since it is used at the python level for e.g. tuning. At that level, individual modes are extracted using the `SimInspiralCreateModeArray()`

* Note that ***we will assume the following***:

    1. Individual multipole moments are referenced externally via `SimInspiralCreateModeArray()`
    2. When the user wants the phase of an individual multipole moment with $(\ell,m)$, then this will be concurrent with the user outputting only the non-twised-up moments via the option below:
    ``` python
        # Tell the model to return the coprecessing mode -- only works on our development branches
        lalsim.SimInspiralWaveformParamsInsertPhenomXReturnCoPrec(lalparams, 1)
    ```
    3. We will also assume that no Multibanding is used, i.e.
    ```Python
        # Turn off multibanding
        lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(lalparams, 0)
        lalsim.SimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(lalparams, 0)
    ```
    4. We will assume that the user is only interested in moments that have approximately no spheroidal mixing.
    5. We will assume that the user is only interested in moments with $m>0$. This is important because the code changes noted below do not effect the `PhenomXHM` approach for mapping $m>0$ moments to $m<0$ ones.

* Within `SimIMRPhenomXPHMFrequencySequence`, the coprecessing moment is called `htildelm`
* `htildelm` is defined as follows
    1. For $(\ell,m)=(2,2)$
    ```C
        status = IMRPhenomXASGenerateFD(&htildelm, freqs, pWF, lalParams);
    ```
    2. For all other moments 
    ```C
        status = IMRPhenomXHMGenerateFDOneMode(&htildelm, freqs, pWF, ell, emmprime, lalParams);
    ```

## Strategy for getting phase

**(a)** The two functions, `IMRPhenomXASGenerateFD` and `IMRPhenomXHMGenerateFDOneMode`, must be modified such that, if an appropriate option is passed via the `LALDict`, then `htildelm` holds the multipole moment phase. 
**(b)** These modification will be backwards-compatible and in-place, meaning that if-else conditions will be used for implementation, and no new functions will be created. Since the phase option of interest will only be used for testing and diagnostics, there is no need for the evaluation to be speed-optimized.
**(c)** The new `LALDict` option will be `PhenomXOnlyReturnPhase`. This is implemented in 'PhenomX_internals.c' via
```C
wf->PhenomXOnlyReturnPhase = XLALSimInspiralWaveformParamsLookupPhenomXOnlyReturnPhase(LALParams);
``` 
where `XLALSimInspiralWaveformParamsLookupPhenomXOnlyReturnPhase` and related functions are defined in the usual way in `lib/LALSimInspiralWaveformParams.h` and `lib/LALSimInspiralWaveformParams.c`. By default, `wf->PhenomXOnlyReturnPhase` will be $0$ (i.e. Off).

## Key Code Changes

* In `IMRPhenomXASGenerateFD`, the new code is

    ```C
      if ( pWF->PhenomXOnlyReturnPhase ) {
        //
        ((*htilde22)->data->data)[jdx] = phi;
      } else {
        /* Reconstruct waveform: h(f) = A(f) * Exp[I phi(f)] */
        ((*htilde22)->data->data)[jdx] = Amp0 * powers_of_Mf.m_seven_sixths * amp * cexp(I * phi);
      }
    ```
    
* In `IMRPhenomXHMGenerateFDOneMode`, the new code is 

    ```C
      if ( pWF->PhenomXOnlyReturnPhase ) {
        //
        ((*htildelm)->data->data)[idx+offset] = phi;
      } else {
        /* Reconstruct waveform: h_l-m(f) = A(f) * Exp[I phi(f)] */
        ((*htildelm)->data->data)[idx+offset] = Amp0 * amp * cexp(I * phi);
      }
    ```
    
* In `SimIMRPhenomXPHMFrequencySequence` where `pWF->IMRPhenomXReturnCoPrec == 1` is checked, we now add e.g.
    ```C
      if( pWF->IMRPhenomXReturnCoPrec == 1 )
      {
          // Do not twist up
          hplus  =  0.5 * hlmcoprec;
          hcross = -0.5 * I * hlmcoprec;
          
          // 
          if( pWF->PhenomXOnlyReturnPhase ){
          // Set hplus to phase (as will be stored in hlmcoprec) and hcross to zero
          hplus  = hlmcoprec;
          hcross = 0;
          }
          
      }
    ```

## Final Comments

* Once implemented, this option appears to reveal a possible inconsistency in either the interaction of various setting, or between XHM and XPHM, whereby within `get_phenomxphm_coprecessing_multipoles`, options `'5-xhm-tuning'` and `'4-xhm'` do not yield consistent results, with the latter yielding phases that are clearly incorrect. 