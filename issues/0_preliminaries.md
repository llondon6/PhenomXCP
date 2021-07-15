# Notes for issue 0_preliminaries

*Here we will keep track of lessons learned from this issue as well as some related ideas.*

## General notes

* A simple way to modify PhenomX, and perhaps a required way based on the model's structure, is to directly modify the ringdown frequency. This may be a worthy test bed for other modifications.

* We wish to first modify the PhenomX infrastructure such that **PhenomXAS** can be modified in the following operative sense: Tuning to NR requires a way of passing modified parameters from python to C. For PhenomDCP this was done entirely in python, but used two "layers": a package containing the python version of PhenomD, and another package containing code for the modeling graft ie PhenomDCP. The later created an abstraction of PhenomD, given some set of modified parameters (deviations away from original parameters). One can think of this as "PhenomDCP-OnTheFly", and a global optimum of these on the fly models was ultimately labeled PhenomDCP. We need to do something similar here, but we will not be able to initiate an abstraction (i.e. class object) in the way that was done for PhenomDCP. Instead, a less efficient but practical route will likely be to create a python wrapper which simply takes in *all* model inputs, including new parameters.

* It seems that the C code does not really treat **PhenomXAS** on par with other PhenomX models; instead, **PhenomXAS is somewhat confusingly often simply referred to as PhenomX** 

## Action items:
Action items marked as complete have associated results in the issue notebook.

- [x] Install lalsuite (see notes [here](https://github.com/llondon6/positive/blob/master/docs/notes/install_lalsuite_locally.md))  
- [x] Call default PhenomX model from python
- [x] Compare PhenomXHM to NR cases (write preliminary code for this purpose)
- [x] Create test modification to PhenomXAS
- [ ] Verify that test modification has the desired effect on the waveform
 * Which files should be modified?
 * At what level should modifications be passes in C, and is this level currently visible in Python?

## Lessons learned

* To get PhenomX to print in debugging mode, one must compile with `PHENOMXHMDEBUG` passed at the config stage via `CFLAGS` : 
```bash
./configure --prefix=${CONDA_PREFIX} --enable-swig-python\
            --disable-lalstochastic --disable-lalxml \
            --disable-lalinference --disable-laldetchar \
            --disable-lalapps --disable-lalframe \
            --disable-lalmetaio --disable-lalburst \
            --disable-lalinspiral CFLAGS="-g -D PHENOMXHMDEBUG"
```

* **Adding options to LALDict is a faff** as it requires the following convoluted steps which must be followed for every new option one wishes to add:
    1. Open `lalsimulation/lib/LALSimInspiralWaveformParams.c`, and then open `lalsimulation/lib/LALSimInspiralWaveformParams.h`
    2. In `lalsimulation/lib/LALSimInspiralWaveformParams.c` prototype three functions. **For example** if we wish to add a field named `CPFlag`, then the following functions should be added in their respective section within the file (i.e. your new functions should be surrounded by like functions)
        a. `int XLALSimInspiralWaveformParamsInsertPhenomXCPFlag(LALDict *params, INT4 value);`
        b. `INT4 XLALSimInspiralWaveformParamsLookupPhenomXCPFlag(LALDict *params);`
        c. `int XLALSimInspiralWaveformParamsPhenomXCPFlagIsDefault(LALDict *params);`
    Here, **note** that the first declaration encodes to type of the dictionary entry.
    3. Similarly, in `lalsimulation/lib/LALSimInspiralWaveformParams.c` prototype three functions:
        a. `DEFINE_INSERT_FUNC(PhenomXCPFlag, INT4, "CPFlag", 0);`
        b. `DEFINE_LOOKUP_FUNC(PhenomXCPFlag, INT4, "CPFlag", 0);`
        c. `DEFINE_ISDEFAULT_FUNC(PhenomXCPFlag, INT4, "CPFlag", 0);`
    **Note** here that the type of the entry is always the second input. 

## PhenomX files and why they are useful 

* *LALSimIMRPhenomX_internals.c* :  