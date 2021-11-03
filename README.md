# PhenomXCP
High-level development of, PhenomXCP, PhenomXAS and PhenomXHM with some feratures tuned to Numerical Relativity simulations. 

-- Lionel London, londonl@mit.edu, l.t.london@uva.nl, pilondon2@gmail.com

## Action Items

* Move waveform flags for CP to IMRPhenomXHM_SetHMWaveformVariables

1. Choose c file in which to add functions for modified ringdown frequencies 
2. Add functions for modified (i.e. effective a la EZH) ringdown frequencies 
3. Make PhenomX option that incorporates effective ringdown frequencies 
4. Make python wrapper for PhenomX+options that results in PhenomXCPv0 = PhenomXAS+EffectiveRingdownFreqs. PhenomXCPv0 may be defined by separate amplitude and phase functions which require beta as input.
5. Compare PhenomXCPv0 to NR simulations

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
