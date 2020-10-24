module LWPC

using Distributed
using JSON3, CSV, DataFrames, Printf

using ..PropagationModelPrep
using ..PropagationModelPrep: rounduprange, unwrap!, LWMS
using ..LWMS

"""
    run(file, computejob::ComputeJob, submitjob=true)

Generate the input files and attempt to run the emp2d code for the scenario
described by `file` as `computejob`.

Optionally provide an inputs::Inputs() struct. Otherwise default values are used.
"""
function run(file, computejob::ComputeJob; submitjob=true)
    isfile(file) || error("$file is not a valid file name")

    s = LWMS.parse(file)

    rundir = computejob.rundir

    if !isdir(rundir)
        # Create rundir if it doesn't exist
        newpath = abspath(rundir)
        @info "Creating $newpath"
        mkpath(newpath)
    end

    if s isa BatchInput
        buildrunjob(s, computejob)
    else
        build(s, computejob)
        submitjob && runjob(computejob)
    end

    return nothing
end


function build(s::BasicInput, computejob::ComputeJob)
    if computejob.runname != s.name
        @info "Updating computejob runname to $(s.name)"
        computejob.runname = s.name
    end

    writeinp(s, computejob)
    writendx(s, computejob)

    return nothing
end

function writeinp(s::BasicInput, computejob::ComputeJob)
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
            # NOTE: According to LWPC manual, b_mag is Tesla, but it appears to actually
            # be Gauss. Also, rounded up because exactly 0 is not supported

            r = round(Int, s.segment_ranges[i]/1e3)  # dist in km
            b_az = rad2deg(s.b_azs[i])  # deg east of north
            b_dip = rad2deg(s.b_dips[i])  # deg from horizontal
            b_mag = round(s.b_mags[i]*1e4, digits=6, RoundUp)
            gsigma = round(s.ground_sigmas[i], digits=6)
            gepsr = s.ground_epsrs[i]
            beta = round(s.betas[i], digits=3)
            hprime = round(s.hprimes[i], digits=3)

            preseg_str = @sprintf(" %d,%.1f,%.1f,%.6f,-1,%.6f,%d,,%.3f,%.3f",
                r, b_az, b_dip, b_mag, gsigma, gepsr, beta, hprime)

            write(f, preseg_str, endline)
        end
        write(f, " 40000", endline)
        write(f, "start", endline)
        write(f, "quit")
    end

    return nothing
end

function writendx(s::BasicInput, computejob::ComputeJob)
    exepath = computejob.exefile
    lwpcpath, exename = splitdir(exepath)
    runname = computejob.runname

    endline = "\n"
    open(joinpath(lwpcpath, "cases", runname*".ndx"), "w") do f
        for i in eachindex(s.segment_ranges)
            r = round(Int, s.segment_ranges[i]/1e3)  # dist in km
            beta = round(s.betas[i], digits=3)
            hprime = round(s.hprimes[i], digits=3)

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
            outfilepath = joinpath(lwpcpath, "Output", file)
            isfile(outfilepath) && rm(outfilepath)
        end
    end

    logfilepath = joinpath(lwpcpath, "cases", runname*".log")
    isfile(logfilepath) && rm(logfilepath)

    input = joinpath("cases", runname)

    # cd to lwpcpath, run lwpc, and return here
    origdir = pwd()
    cd(lwpcpath)
    output = read(`$exename $input`, String)
    cd(origdir)

    return output
end

function buildrunjob(s::BatchInput{BasicInput}, computejob::LocalParallel)
    nprocs() == computejob.numnodes+1 || @warn "computejob.numnodes does not match nprocs(). Using $(nprocs()) processes."

    _buildrunjob = function (s)
        pid = myid()

        exefile, exeext = splitext(computejob.exefile)
        lwpcpath, exefilename = splitdir(exefile)

        newlwpcpath = lwpcpath*"_"*string(pid)
        newexefile = joinpath(newlwpcpath, exefilename*string(pid)*exeext)

        cj = Local(s.name, newlwpcpath, newexefile)

        build(s, cj)
        runjob(cj)

        tstart = time()  # seconds
        goodlog = false
        while !goodlog && (time() - tstart < 60)
            try
                o = readlog(joinpath(newlwpcpath, "cases", s.name*".log"))
                goodlog = true
            catch
                continue
            end
        end

        sleep(0.2)
        return nothing
    end

    pmap(_buildrunjob, s.inputs)
end

function readlog(file)
    # Search for first and last line of data
    lines = readlines(file)
    firstdataline = findfirst(startswith.(lines, "  dist   amplitude  phase")) + 1
    lastdataline = findfirst(startswith.(lines, "nc nrpt bearng")) - 1
    skiplastlines = length(lines) - lastdataline

    # Read log file
    raw = CSV.File(file; skipto=firstdataline, footerskip=skiplastlines,
                   delim=' ', ignorerepeated=true, header=false,
                   silencewarnings=true)

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

    # If phase gets above 9999 deg in log file, there is no space between amp and phase
    if count(ismissing.(amp)) != count(ismissing.(phase))
        for i in eachindex(phase)
            if ismissing(phase[i])
                phase[i] = parse(Float64, amp[i][end-9:end])  # phase is 10000.0000
                amp[i] = parse(Float64, amp[i][1:end-10])
            end
        end
        # Other elements in the same column will also be string type
        for i in eachindex(amp)
            if amp[i] isa String
                amp[i] = parse(Float64, amp[i])
            end
        end
    end

    return dist, amp, phase
end

function process(jsonfile, computejob::ComputeJob)
    s = LWMS.parse(jsonfile)

    return process(s, computejob)
end

function process(s::BasicInput, computejob::ComputeJob)
    exepath = computejob.exefile
    lwpcpath, exename = splitdir(exepath)
    runname = computejob.runname
    rundir = computejob.rundir

    name = s.name
    description = s.description
    datetime = s.datetime

    dist, amp, phase = readlog(joinpath(lwpcpath, "cases", runname*".log"))
    dist *= 1e3  # convert to m

    output = BasicOutput()
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

function process(jsonfile, computejob::LocalParallel)
    s = LWMS.parse(jsonfile)

    return process(s, computejob)
end

function process(s::BatchInput{BasicInput}, computejob::LocalParallel)
    lwpcfile, ext = splitext(computejob.exefile)
    lwpcpath, lwpcfilename = splitdir(lwpcfile)

    batch = BatchOutput{BasicOutput}()
    batch.name = s.name
    batch.description = s.description
    batch.datetime = Dates.now()

    for n in 1:nprocs()
        newlwpcpath = lwpcpath*"_"*string(n)

        for i in eachindex(s.inputs)
            runname = s.inputs[i].name
            logfile = joinpath(newlwpcpath, "cases", runname*".json")

            if isfile(logfile)
                dist, amp, phase = readlog(logfile)
                dist *= 1e3  # convert to m

                output = BasicOutput()
                output.name = name
                output.description = description
                output.datetime = datetime

                # Strip last index of each field because they're 0 (and then inf)
                output.output_ranges = round.(dist, digits=-3)  # fix floating point, rounded to nearest km
                output.amplitude = amp
                output.phase = phase

                push!(batch.outputs, output)
            end
        end
    end

    # Save output
    json_str = JSON3.write(batch)

    open(joinpath(computejob.rundir, computejob.name*"_lwpc.json"), "w") do f
        write(f, json_str)
    end

    return batch
end

end  # module
