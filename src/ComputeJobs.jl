"""
    ComputeJob

An abstract type for different computer clusters and job managers.

All `ComputeJob`s should have `runname`, `rundir`, and `exefile` fields at a
minimum.
"""
abstract type ComputeJob end

"""
walltime "00:57:00" is 57 minutes
"""
struct Summit <: ComputeJob
    runname::String
    rundir::String
    numnodes::Int
    walltime::String
    exefile::String
end

struct LocalOMP <: ComputeJob
    runname::String
    rundir::String
    numnodes::Int
    exefile::String
end


function writeshfile(s::LocalOMP)

    runname = s.runname
    rundir = s.rundir
    numnodes = s.numnodes
    exefile = s.exefile
    exefile = basename(exefile)

    shfile = joinpath(rundir, runname*".sh")
    endline = "\n"

    open(shfile, "w") do f
        write(f, "#!/bin/sh\n", endline)
        write(f, endline)
        write(f, "rm -f $rundir/output_K.dat", endline)
        write(f, "rm -f $rundir/output_E.dat", endline)
        write(f, "rm -f $rundir/output_D.dat", endline)
        write(f, "rm -f $rundir/output_H.dat", endline)
        write(f, "rm -f $rundir/output_J.dat", endline)
        write(f, "rm -f $rundir/output_T.dat", endline)
        write(f, "rm -f $rundir/output_O.dat", endline)
        write(f, "rm -f $rundir/output_S.dat", endline)
        write(f, "rm -f $rundir/Probe.dat", endline)
        write(f, "rm -f $rundir/elve.dat", endline)
        write(f, "rm -f $rundir/sferic.dat", endline)
        write(f, endline)
        write(f, "export OMP_NUM_THREADS=$numnodes", endline)
        write(f, endline)
        write(f, joinpath(rundir,exefile), endline)
    end

    return shfile
end

function writeshfile(s::Summit)

    runname = s.runname
    rundir = s.rundir
    numnodes = s.numnodes
    walltime = s.walltime
    exefile = s.exefile
    exefile = basename(exefile)

    shfile = joinpath(rundir, runname*".sh")
    endline = "\n"

    open(shfile, "w") do f
        write(f, "#!/bin/sh\n", endline)
        write(f, "#SBATCH --job-name=$runname", endline)
        write(f, "#SBATCH --partition=shas", endline)
        write(f, "#SBATCH --nodes=1", endline)
        write(f, "#SBATCH --ntasks=$numnodes", endline)
        write(f, "#SBATCH --time=$walltime", endline)
        write(f, endline)
        write(f, "rm -f $rundir/output_K.dat", endline)
        write(f, "rm -f $rundir/output_E.dat", endline)
        write(f, "rm -f $rundir/output_D.dat", endline)
        write(f, "rm -f $rundir/output_H.dat", endline)
        write(f, "rm -f $rundir/output_J.dat", endline)
        write(f, "rm -f $rundir/output_T.dat", endline)
        write(f, "rm -f $rundir/output_O.dat", endline)
        write(f, "rm -f $rundir/output_S.dat", endline)
        write(f, "rm -f $rundir/Probe.dat", endline)
        write(f, "rm -f $rundir/elve.dat", endline)
        write(f, "rm -f $rundir/sferic.dat", endline)
        write(f, endline)
        write(f, "module purge", endline)
        write(f, "module load intel", endline)
        write(f, endline)
        write(f, "export OMP_NUM_THREADS=$numnodes", endline)
        write(f, endline)
        write(f, joinpath(rundir,exefile), endline)
    end

    return shfile
end

function runjob(s::LocalOMP, shfile)
    jobname = read(`$shfile`, String)

    println(jobname)

    return nothing
end

function runjob(s::Summit, shfile)
    jobname = read(`sbatch $shfile`, String)
    jobid = strip(jobname)

    println("Job $jobid submitted!\n")

    return nothing
end
