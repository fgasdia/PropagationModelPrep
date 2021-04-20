module LWPC

using Printf, Random
using JSON3, CSV, DataFrames

using ..PropagationModelPrep
using ..PropagationModelPrep: LMP
using ..LMP

mutable struct ProcessInfo
    inputidx::Int
    process::Union{Base.Process,Nothing}
    starttime::Float64

    ProcessInfo() = new(0,nothing,0.0)
end
start!(p::ProcessInfo) = (p.starttime = time())
elapsed(p::ProcessInfo) = time() - p.starttime


"""
    run(file, computejob; submitjob=true, savefile=true)
    run(input, computejob; submitjob=true, savefile=true)

Generate the input files and run the LWPC code for the scenario described by `file` with the
compute parameters defined by `computejob`. If `file` does not contain a `BatchInput` and
`submitjob` is `false`, then the input files are created but LWPC is not run.

The output is written to a `.json` file if `savefile` is `true`.

If `computejob` is a `LocalParallel`, it is assumed that there exist directories named
"C:\\LWPCv21_0" to "C:\\LWPCv21_N" where "N" is `computejob.numnodes`. In each directory
it is assumed there is an associated executable "lwpm0.exe" to "lwpmN.exe". The `computejob`
should have `exefile = C:\\LWPCv21\\lwpm.exe`.
"""
function run(file::String, computejob; submitjob=true, savefile=true)
    isfile(file) || error("$file is not a valid file name")

    s = LMP.parse(file)

    run(s, computejob; submitjob=submitjob, savefile=savefile)
end

function run(input, computejob; submitjob=true, savefile=true)
    rundir = computejob.rundir

    if !isdir(rundir)
        # Create rundir if it doesn't exist
        newpath = abspath(rundir)
        @info "Creating $newpath"
        mkpath(newpath)
    end

    output = build_runjob(input, computejob; submitjob=submitjob)

    if savefile
        json_str = JSON3.write(output)

        open(joinpath(rundir, input.name*"_lwpc.json"), "w") do f
            write(f, json_str)
        end
    end

    return output
end

"""
    build(input, computejob)

Write LWPC `.inp` and `.ndx` files.
"""
function build(input::BasicInput, computejob)
    if computejob.runname != input.name
        @info "Updating computejob runname to $(input.name)"
        computejob.runname = input.name
    end

    writeinp(input, computejob)
    writendx(input, computejob)

    return nothing
end

"""
    randtransmittername()

Return a random ten character `'a':'z'` string to be used as a transmitter name.

If the scenario transmitter frequency doesn't match up with an existing transmitter name
in LWPC's xmtr.lis, LWPC will break.
"""
randtransmittername() = randstring('a':'z', 10)

"""
    writeinp(input, computejob)
"""
function writeinp(input::BasicInput, computejob)
    exepath = computejob.exefile
    lwpcpath, exename = splitdir(exepath)
    runname = computejob.runname

    freq = input.frequency/1e3  # in kHz
    max_range = ceil(Int, maximum(input.output_ranges)/1e3)  # convert to km
    diffrange = diff(input.output_ranges)
    drange = trunc(Int, diffrange[1]/1e3)  # in km
    length(unique(round.(diffrange))) == 1 || @info "Using drange = $drange km"

    if length(input.segment_ranges) == 1
        homogeneous = true
        beta = only(input.betas)
        hprime = only(input.hprimes)
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
        write(f, "tx-data     $(randtransmittername())  $freq  0.0  0.0  100.0  0.0  0.000  0.0", endline)
        write(f, "receivers   0.0000  0.0000", endline)
        write(f, "range-max   $max_range", endline)
        write(f, "ionosphere  $ionosphere", endline)
        write(f, "rx-data     vertical 0", endline)
        write(f, "lwf-vs-distance $max_range $drange", endline)
        write(f, "mc-options  full-wave 50 true", endline)
        write(f, "print-lwf   1", endline)
        write(f, "lwflds", endline)
        write(f, "preseg", endline)
        for i in eachindex(input.segment_ranges)
            # NOTE: According to LWPC manual, b_mag is Tesla, but it appears to actually
            # be Gauss. Also, rounded up because exactly 0 is not supported

            r = round(Int, input.segment_ranges[i]/1e3)  # dist in km
            b_az = rad2deg(input.b_azs[i])  # deg east of north
            b_dip = rad2deg(input.b_dips[i])  # deg from horizontal
            b_mag = round(input.b_mags[i]*1e4, digits=6, RoundUp)
            gsigma = round(input.ground_sigmas[i], digits=6)
            gepsr = input.ground_epsrs[i]
            beta = round(input.betas[i], digits=3)
            hprime = round(input.hprimes[i], digits=3)

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

"""
    writendx(input, computejob)
"""
function writendx(input::BasicInput, computejob)
    exepath = computejob.exefile
    lwpcpath, exename = splitdir(exepath)
    runname = computejob.runname

    endline = "\n"
    open(joinpath(lwpcpath, "cases", runname*".ndx"), "w") do f
        for i in eachindex(input.segment_ranges)
            r = round(Int, input.segment_ranges[i]/1e3)  # dist in km
            beta = round(input.betas[i], digits=3)
            hprime = round(input.hprimes[i], digits=3)

            write(f, "$r $beta $hprime", endline)
        end
    end

    return nothing
end

"""
    deletefiles(lwpcpath, runname)
    deletefiles(computejob)

Delete any files that might already exist in LWPC "Output" folder and the `.log` file for
`runname`.
"""
function deletefiles(lwpcpath, runname)
    output_files = readdir(joinpath(lwpcpath, "Output"))
    for file in output_files
        if splitext(basename(file))[1] == runname
            outfilepath = joinpath(lwpcpath, "Output", file)
            isfile(outfilepath) && rm(outfilepath)
        end
    end

    logfilepath = joinpath(lwpcpath, "cases", runname*".log")
    isfile(logfilepath) && rm(logfilepath)
end
deletefiles(computejob) = deletefiles(splitdir(computejob.exefile)[1], computejob.runname)

"""
    runjob(computejob)

Run the `computejob` and return the associated `Process`.
"""
function runjob(computejob::Local)
    exepath = computejob.exefile
    lwpcpath, exename = splitdir(exepath)
    runname = computejob.runname

    deletefiles(lwpcpath, runname)

    inputname = joinpath("cases", runname)

    cmd = Cmd(`$exename $inputname`; dir=lwpcpath)
    process = Base.run(cmd; wait=false)  # run asynchronously

    return process
end

"""
    build_runjob(input::BasicInput, computejob; submitjob=true)

Construct the LWPC files for `input`, and if `submitjob` is true, run LWPC and return results
as a `BasicOutput`.

No matter if `computejob` is parallel or not, `input::BasicInput` will be run with a single
process.
"""
function build_runjob(input::BasicInput, computejob; submitjob=true)
    exefile, exeext = splitext(computejob.exefile)
    lwpcpath, exefilename = splitdir(exefile)

    build(input, computejob)

    if submitjob
        output = BasicOutput()
        output.name = input.name
        output.description = input.description
        output.datetime = input.datetime

        t0 = time()
        process = runjob(computejob)
        completed = false
        while !completed && (time() - t0 < computejob.walltime)
            if process_exited(process)
                dist, amp, phase = readlog(joinpath(lwpcpath, "cases", input.name*".log"))
                dist *= 1e3  # convert to m
                phase .= deg2rad.(phase)  # convert from deg to rad

                # Strip last index of each field because they're 0 (and then inf)
                output.output_ranges = round.(dist, digits=-3)  # fix floating point, rounded to nearest km
                output.amplitude = amp
                output.phase = phase

                completed = true
            else
                sleep(0.1)
            end
        end
        if !completed
            @warn "LWPC time limit exceeded"

            output.output_ranges = [NaN]
            output.amplitude = [NaN]
            output.phase = [NaN]
        end
    else
        output = nothing
    end

    return output
end

"""
    build_runjob(inputs::BatchInput, computejob; submitjob=true)

For each of `inputs.inputs`, construct the input files and run LWPC, returning results as a
`BatchOutput{BasicOutput}`.

If `computejob` is a `LocalParallel`, it is assumed that there exist directories named
"C:\\LWPCv21_0" to "C:\\LWPCv21_N" where "N" is `computejob.numnodes`. In each directory
it is assumed there is an associated executable "lwpm0.exe" to "lwpmN.exe". The `computejob`
should have `exefile = C:\\LWPCv21\\lwpm.exe`.

!!! note
    The `submitjob` argument is ignored if `computejob` is `LocalParallel`.
"""
function build_runjob(inputs::BatchInput{BasicInput}, computejob::LocalParallel; submitjob=true)
    exefile, exeext = splitext(computejob.exefile)
    lwpcpath, exefilename = splitdir(exefile)

    numinputs = length(inputs.inputs)

    batch = BatchOutput{BasicOutput}()
    batch.name = inputs.name
    batch.description = inputs.description
    batch.datetime = inputs.datetime
    batch.outputs = Vector{BasicOutput}(undef, numinputs)
    
    completed = falses(numinputs)
    processes = Tuple(ProcessInfo() for i in 1:computejob.numnodes)
    for (inputidx, input) in enumerate(inputs.inputs)
        submitted = false
        while !submitted
            for (pid, proc) in enumerate(processes)
                newlwpcpath = lwpcpath*"_"*string(pid-1)
                newexefile = joinpath(newlwpcpath, exefilename*string(pid-1)*exeext)
                ii = proc.inputidx

                if isnothing(proc.process)
                    # Process is available, can assign immediately
                    @debug "Process $pid is `nothing`"

                    cj = Local(input.name, newlwpcpath, newexefile, computejob.walltime)

                    build(input, cj)
                    proc.inputidx = inputidx
                    start!(proc)
                    proc.process = runjob(cj)
                    submitted = true
                    @debug "Input $inputidx submitted"
                    break
                elseif process_exited(proc.process)
                    # Process is completed
                    @debug "Process $pid has status exited"

                    dist, amp, phase = readlog(joinpath(newlwpcpath, "cases", inputs.inputs[ii].name*".log"))
                    dist *= 1e3  # convert to m
                    phase .= deg2rad.(phase)  # convert from deg to rad

                    output = BasicOutput()
                    output.name = inputs.inputs[ii].name
                    output.description = inputs.inputs[ii].description
                    output.datetime = inputs.inputs[ii].datetime

                    # Strip last index of each field because they're 0 (and then inf)
                    output.output_ranges = round.(dist, digits=-3)  # fix floating point, rounded to nearest km
                    output.amplitude = amp
                    output.phase = phase

                    batch.outputs[ii] = output
                    completed[ii] = true
                    @debug "Input $ii is completed"

                    # Start next process
                    cj = Local(input.name, newlwpcpath, newexefile, computejob.walltime)
                    build(input, cj)
                    proc.inputidx = inputidx
                    start!(proc)
                    proc.process = runjob(cj)
                    submitted = true
                    @debug "Input $inputidx submitted"
                    break
                elseif elapsed(proc) > computejob.walltime
                    # Process timed out
                    @warn "LWPC time limit exceeded"
                    @debug "Process $pid has exceeded walltime"
                    
                    kill(proc.process)

                    output = BasicOutput()
                    output.name = inputs.inputs[ii].name
                    output.description = inputs.inputs[ii].description
                    output.datetime = inputs.inputs[ii].datetime

                    output.output_ranges = [NaN]
                    output.amplitude = [NaN]
                    output.phase = [NaN]

                    batch.outputs[ii] = output
                    completed[ii] = true
                    @debug "Input $ii is completed"

                    # Start the next process
                    cj = Local(input.name, newlwpcpath, newexefile, computejob.walltime)
                    build(input, cj)
                    proc.inputidx = inputidx
                    start!(proc)
                    proc.process = runjob(cj)
                    submitted = true
                    @debug "Input $inputidx submitted"
                    break
                end
            end
            sleep(0.1)
        end
    end

    @debug "All inputs submitted"
    @debug "$(count(completed)) inputs completed"

    # Remaining results
    while any(!, completed)
        for (pid, proc) in enumerate(processes)
            isnothing(proc.process) && continue

            newlwpcpath = lwpcpath*"_"*string(pid-1)
            ii = proc.inputidx
            if !completed[ii] && process_exited(proc.process)
                dist, amp, phase = readlog(joinpath(newlwpcpath, "cases", inputs.inputs[ii].name*".log"))
                dist *= 1e3  # convert to m
                phase .= deg2rad.(phase)  # convert from deg to rad

                output = BasicOutput()
                output.name = inputs.inputs[ii].name
                output.description = inputs.inputs[ii].description
                output.datetime = inputs.inputs[ii].datetime

                # Strip last index of each field because they're 0 (and then inf)
                output.output_ranges = round.(dist, digits=-3)  # fix floating point, rounded to nearest km
                output.amplitude = amp
                output.phase = phase

                batch.outputs[ii] = output
                completed[ii] = true
            elseif elapsed(proc) > computejob.walltime
                # Process timed out
                @warn "LWPC time limit exceeded"
                    
                kill(proc.process)

                output = BasicOutput()
                output.name = inputs.inputs[ii].name
                output.description = inputs.inputs[ii].description
                output.datetime = inputs.inputs[ii].datetime

                output.output_ranges = [NaN]
                output.amplitude = [NaN]
                output.phase = [NaN]

                batch.outputs[ii] = output
                completed[ii] = true
            end
        end
        sleep(0.1)
    end

    @debug "Final: All $(count(completed)) inputs completed"

    return batch
end

function build_runjob(inputs::BatchInput{BasicInput}, computejob::Local; submitjob=true)
    exefile, exeext = splitext(computejob.exefile)
    lwpcpath, exefilename = splitdir(exefile)

    batch = BatchOutput{BasicOutput}()
    batch.name = inputs.name
    batch.description = inputs.description
    batch.datetime = inputs.datetime
    batch.outputs = Vector{BasicOutput}(undef, length(inputs.inputs))
    
    for i in eachindex(inputs.inputs)
        cj = Local(inputs.inputs[i].name, computejob.rundir, computejob.exefile, computejob.walltime)
        o = build_runjob(inputs.inputs[i], cj; submitjob=submitjob)
        batch.outputs[i] = o
    end

    return batch
end

"""
    readlog(file)

Return vectors `(dist, amp, phase)` read from the log named `file`.

These are the raw values (no unit conversions) from LWPC.
"""
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

end  # module
