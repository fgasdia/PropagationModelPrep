"""
    ComputeJob

An abstract type for different computer clusters and job managers.

`runname` does not need to be set because it will be automatically overwritten
by the JSON file name.

If `rundir` path does not end in a directory with the same name as `runname`,
such a directory will be created.

All `ComputeJob`s should have `runname`, `rundir`, and `exefile` fields at a
minimum.
"""
abstract type ComputeJob end

"""
walltime "00:57:00" is 57 minutes
"""
mutable struct Summit <: ComputeJob
    runname::String
    rundir::String
    numnodes::Int
    walltime::String
    exefile::String
end

mutable struct LocalOMP <: ComputeJob
    runname::String
    rundir::String
    numnodes::Int
    exefile::String
end

mutable struct LairLwpc <: ComputeJob
    runname::String
    rundir::String
    exefile::String
end
