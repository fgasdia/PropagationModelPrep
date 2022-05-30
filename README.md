# PropagationModelPrep.jl

Parse `.json` files which describe a VLF Earth-ionosphere waveguide propagation scenario (from the [LongwaveModePropagator.jl](https://github.com/fgasdia/LongwaveModePropagator.jl)) and generate and run files for other models, such as the EMP2D FDTD model and LWPC.

This package is specific to the models and compute resources available to me, but can serve as a template to others. It is not meant for "general use".

## ComputeJobs

The `ComputeJob` is an `abstract type` to represent compute jobs with specific file paths and resources. Concrete subtypes include specific computers (such as `Summit` at CU Boulder) or generic `Local` computers or `LocalParallel` computers with resources for running propagation models in parallel.

---

This package provides submodules `.EMP2D` and `.LWPC`. They both produce `BasicOutput`-formatted .json files from LongwaveModePropagator.jl.

## EMP2D

EMP2D is an FDTD propagation model ((Marshall, 2017)[https://doi.org/10.1109/TAP.2017.2758392]) which is usually run with OpenMP. A typical ~2000 km VLF propagation path may take 2 hours on a modern 4 core machine.

The `EMP2D` submodule provides `struct`s to generate input files and read output files from LWPC. These are similar to the MATLAB files written by Bob Marshall for setup and analysis of EMP2D results. There may be minor differences in default parameters to ensure consistency between the inputs to EMP2D and LWPC.

We provide for both `linear_exponential` and `ricker` source waveforms, but the `ricker` source has not yet been successfully verified.

See `EMP2D.Inputs` to define model inputs.

## LWPC

`src/preplwpc.jl` is a script which automatically sets up LWPC to be run in parallel on a single computer. It copies the original LWPC folder and changes hard-coded paths so that each copy of LWPC only writes and reads within its own folder.

This submodule can automatically run LWPC jobs in parallel. Depending on procesor and drive read speed, it may be necessary to adjust the `sleeptime` argument, but I recommend trying the default first.

Model inputs can be specified as `ExponentialInput` or `BatchInput{ExponentialInput}`. For example,

```julia
batch = BatchInput{ExponentialInput}()
# fill in `batch`...
computejob = LocalParallel("jobname", ".", "C:\\LWPCv21\\lwpm.exe", numexe, 90)
output = LWPC.run(batch, computejob)  # returns a `BatchOutput{BasicOutput}`
```
