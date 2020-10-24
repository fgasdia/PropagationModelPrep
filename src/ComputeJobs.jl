"""
    ComputeJob

An abstract type for different computer clusters and job managers.

`runname` does not need to be set because it will be automatically overwritten by the JSON
file name.

If `rundir` path does not end in a directory with the same name as `runname`, such a
directory will be created.

All `ComputeJob`s should have `runname`, `rundir`, and `exefile` fields at a minimum.
"""
abstract type ComputeJob end

"""
    ParallelComputeJob

An abstract type for `ComputeJob`s that execute in parallel (by OpenMP or Julia's
Distributed).

`ParallelComputeJob` types should have a `numnodes` field to represent the number of
"processors" (or threads, etc.) over which the job will run. The exact meaning of `numnodes`
can vary between concrete types.
"""
abstract type ParallelComputeJob <: ComputeJob end

"""
    Summit <: ParallelComputeJob

`walltime` should be specified as a `String` of the form `"00:57:00"` to represent 57
minutes of expected total runtime.
"""
mutable struct Summit <: ParallelComputeJob
    runname::String
    rundir::String
    exefile::String
    numnodes::Int
    walltime::String
end

"""
    Local <: ComputeJob

Mutable struct holding `runname`, `rundir`, and `exefile` strings refering to paths on the
local machine.
"""
mutable struct Local <: ComputeJob
    runname::String
    rundir::String
    exefile::String
end

mutable struct LocalParallel <: ParallelComputeJob
    runname::String
    rundir::String
    exefile::String
    numnodes::Int
end
