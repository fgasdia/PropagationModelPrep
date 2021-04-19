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

function Base.copy(cj::ComputeJob)
    T = typeof(cj)
    newcj = T()
    for fn in fieldnames(T)
        setfield!(newcj, fn, getfield(cj, fn))
    end
    return newcj
end

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

    Summit() = new()
end
function Summit(runname, rundir, exefile, numnodes, walltime)
    cj = Summit()
    cj.runname = runname
    cj.rundir = rundir
    cj.exefile = exefile
    cj.numnodes = numnodes
    cj.walltime = walltime
    return cj
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

    Local() = new()
end
function Local(runname, rundir, exefile)
    cj = Local()
    cj.runname = runname
    cj.rundir = rundir
    cj.exefile = exefile
    return cj
end

mutable struct LocalParallel <: ParallelComputeJob
    runname::String
    rundir::String
    exefile::String
    numnodes::Int
    walltime::Int  # seconds, time limit per run

    LocalParallel() = new()
end
function LocalParallel(runname, rundir, exefile, numnodes, walltime)
    cj = LocalParallel()
    cj.runname = runname
    cj.rundir = rundir
    cj.exefile = exefile
    cj.numnodes = numnodes
    cj.walltime = walltime

    return cj
end
