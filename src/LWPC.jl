module LWPC

using JSON3, CSV, DataFrames

using ..PropagationModelPrep
using ..PropagationModelPrep: rounduprange, unwrap!
import ..LWMS


"""
    run(file, computejob::ComputeJob, submitjob=true)

Generate the input files and attempt to run the emp2d code for the scenario
described by `file` as `computejob`.

Optionally provide an inputs::Inputs() struct. Otherwise default values are used.
"""
function run(file, computejob::ComputeJob; submitjob=true)
    isfile(file) || error("$file is not a valid file name")

    s = LWMS.parse(file)

    if computejob.runname != s.name
        @info "Updating computejob runname to $(s.name)"
        computejob.runname = s.name
    end

    # Get abspath in case a relpath is provided
    origrundir = splitpath(abspath(computejob.rundir))[end]
    if origrundir != computejob.runname
        # Check if computejob rundir path ends with a directory called runname
        rundir = joinpath(computejob.rundir, computejob.runname)
        @info "Running in $rundir/"
    else
        rundir = computejob.rundir
    end

    if !isdir(rundir)
        # Create rundir if it doesn't exist
        @info "Creating $rundir/"
        mkpath(rundir)
    end

    build(s, computejob)

    if submitjob
        runjob(computejob)
    end

    return nothing
end

"""
    build(s::LWMS.BasicInput, computejob::ComputeJob, inputs::Inputs)

This is essentially a "private" function that sets default parameters for LWPC
and generates the necessary input files.
"""
function build(s::LWMS.BasicInput, computejob::ComputeJob)
    rundir = computejob.rundir

    writeinp(s, computejob)
    writendx(s, computejob)

    return nothing
end

function writeinp(s::LWMS.BasicInput, computejob::ComputeJob)
    exepath = computejob.exefile
    lwpcpath, exename = splitdir(exepath)
    runname = computejob.runname

    freq = s.frequency/1e3  # in kHz
    max_range = trunc(Int, rounduprange(maximum(s.output_ranges))/1e3)  # convert to km
    diffrange = diff(s.output_ranges)
    drange = trunc(Int, diffrange[1]/1e3)  # in km
    length(unique(round.(diffrange))) == 1 || @info "Using drange = $drange km"

    if length(s.segment_ranges) == 1
        homogeneous = true
        beta = only(s.betas)
        hprime = only(s.hprimes)
        ionosphere = "homogeneous exponential $beta $hprime"
    else
        homogeneous = false
        ionosphere = "range exponential $runname"
    end

    endline = "\n"
    open(joinpath(lwpcpath, "cases", runname*".inp"), "w") do f
        write(f, "file-mds    Output\\", endline)
        write(f, "file-lwf    Output\\", endline)
        write(f, "file-prf    cases\\", endline)
        write(f, "file-ndx    cases\\", endline)
        write(f, "case-id     $runname", endline)
        write(f, "tx          $runname", endline)
        write(f, "tx-data     vlftx  $freq  0.0  0.0  100.0  0.0  0.000  0.0", endline)
        write(f, "receivers   0.0000  0.0000", endline)
        write(f, "range-max   $max_range", endline)
        write(f, "ionosphere  $ionosphere", endline)
        write(f, "rx-data     vertical 0", endline)
        write(f, "lwf-vs-distance $max_range $drange", endline)
        write(f, "mc-options  full-wave 50 true", endline)
        write(f, "print-lwf   1", endline)
        write(f, "lwflds", endline)
        write(f, "preseg", endline)
        for i in eachindex(s.segment_ranges)
            r = trunc(Int, s.segment_ranges[i]/1e3)  # dist in km
            b_az = s.b_az[i]  # deg east of north
            b_dip = s.b_dip[i]  # deg from horizontal
            b_mag = s.b_mag[i]*1e4  # XXX: supposedly Tesla, but probably Gauss
            gsigma = s.ground_sigmas[i]
            gepsr = s.ground_epsr[i]
            beta = s.betas[i]
            hprime = s.hprimes[i]

            write(f, " $r,$b_az,$b_dip,$b_mag,-1,$gsigma,$gepsr,-1,$beta,$hprime", endline)
        end
        write(f, " 40000", endline)
        write(f, "start", endline)
        write(f, "quit")
    end

    return nothing
end

function writendx(s::LWMS.BasicInput, computejob::ComputeJob)
    exepath = computejob.exefile
    lwpcpath, exename = splitdir(exepath)
    runname = computejob.runname

    endline = "\n"
    open(joinpath(lwpcpath, "cases", runname*".ndx"), "w") do f
        for i in eachindex(s.segment_ranges)
            r = s.segment_ranges[i]/1e3  # dist in km
            beta = s.betas[i]
            hprime = s.hprimes[i]

            write(f, "$r $beta $hprime", endline)
        end
    end

    return nothing
end

function runjob(computejob::Local)
    exepath = computejob.exefile
    lwpcpath, exename = splitdir(exepath)
    runname = computejob.runname

    # Delete any files that might already exist in LWPC Output folder for same runname
    output_files = readdir(joinpath(lwpcpath, "Output"))
    for file in output_files
        if splitext(basename(file))[1] == runname
            rm(joinpath(lwpcpath, "Output", file), force=true)
        end
    end

    rm(joinpath(lwpcpath, "cases", runname*".log"), force=true)

    input = joinpath("cases", runname)

    # cd to lwpcpath, run lwpc, and return here
    origdir = pwd()
    cd(lwpcpath)
    output = read(`$exename $input`, String)
    cd(origdir)

    println("Job submitted!\n")
    println(output)

    return nothing
end

function readlog(file)
    # Search for first and last line of data
    lines = readlines(file)
    firstdataline = findfirst(startswith.(lines, "  dist   amplitude  phase")) + 1
    lastdataline = findfirst(startswith.(lines, "nc nrpt bearng")) - 1
    skiplastlines = length(lines) - lastdataline

    # Read log file
    raw = CSV.File(file; skipto=firstdataline, footerskip=skiplastlines,
                   delim=' ', ignorerepeated=true, header=false)

    # Stack the columns together
    dist = vcat(raw.Column1, raw.Column4, raw.Column7)
    amp = vcat(raw.Column2, raw.Column5, raw.Column8)
    phase = vcat(raw.Column3, raw.Column6, raw.Column9)

    # Clean up end of the last column
    delidxs = ismissing.(dist)  # assuming these are at the end
    delidxs[findfirst(delidxs)-1] = true  # last valid entry is at receiver distance
    deleteat!(dist, delidxs)
    deleteat!(amp, delidxs)
    deleteat!(phase, delidxs)

    return dist, amp, phase
end

function process(jsonfile, computejob::ComputeJob)
    exepath = computejob.exefile
    lwpcpath, exename = splitdir(exepath)
    runname = computejob.runname
    rundir = computejob.rundir

    if !isfile(jsonfile)
        @info "$jsonfile not found. Assuming run name is $pathname."
        name = pathname
        description = "LWPC results"
        datetime = Dates.now()
    else
        s = LWMS.parse(jsonfile)
        name = s.name
        description = s.description
        datetime = s.datetime
    end

    dist, amp, phase = readlog(joinpath(lwpcpath, "cases", runname*".log"))
    dist *= 1e3  # convert to m

    output = LWMS.BasicOutput()
    output.name = name
    output.description = description
    output.datetime = datetime

    # Strip last index of each field because they're 0 (and then inf)
    output.output_ranges = round.(dist, digits=-3)  # fix floating point, rounded to nearest km
    output.amplitude = amp
    output.phase = phase

    # Save output
    json_str = JSON3.write(output)

    open(joinpath(rundir,name*"_lwpc.json"), "w") do f
        write(f, json_str)
    end

    return output
end

end  # module
