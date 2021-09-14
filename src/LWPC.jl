module LWPC

using Printf, Random
using JSON3, CSV
using ProgressLogging

using ..PropagationModelPrep
using ..PropagationModelPrep: LMP
using ..LMP

mutable struct ProcessInfo
    inputidx::Int
    process::Union{Base.Process,Nothing}
    starttime::Float64

    ProcessInfo() = new(0, nothing, 0.0)
end
start!(p::ProcessInfo) = (p.starttime = time())
elapsed(p::ProcessInfo) = time() - p.starttime


"""
    run(file, computejob; submitjob=true, savefile=true, sleeptime=0.1)
    run(input, computejob; submitjob=true, savefile=true, sleeptime=0.1)

Generate the input files and run the LWPC code for the scenario described by `file` with the
compute parameters defined by `computejob`. If `file` does not contain a `BatchInput` and
`submitjob` is `false`, then the input files are created but LWPC is not run.

The output is written to a `.json` file if `savefile` is `true`.

`sleeptime` is the amount of time to `sleep` in seconds in between submitting jobs to
`batch_runjob`.

If `computejob` is a `LocalParallel`, it is assumed that there exist directories named
"C:\\LWPCv21_0" to "C:\\LWPCv21_N" where "N" is `computejob.numnodes`. In each directory
it is assumed there is an associated executable "lwpm0.exe" to "lwpmN.exe". The `computejob`
should have `exefile = C:\\LWPCv21\\lwpm.exe`.
"""
function run(file::String, computejob; submitjob=true, savefile=true, sleeptime=0.1)
    isfile(file) || error("$file is not a valid file name")

    s = LMP.parse(file)

    run(s, computejob; submitjob, savefile, sleeptime)
end

function run(input, computejob; submitjob=true, savefile=true, sleeptime=0.1)
    rundir = computejob.rundir

    if !isdir(rundir)
        # Create rundir if it doesn't exist
        newpath = abspath(rundir)
        @info "Creating $newpath"
        mkpath(newpath)
    end

    output = build_runjob(input, computejob; submitjob, sleeptime)

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
function build(input::ExponentialInput, computejob)
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
    writeinp(input::ExponentialInput, computejob)

Write `.inp` file for LWPC using an `ExponentialInput` from LongwaveModePropagator.jl.

Hardcoded to use a transmitter power of 100 kW to match LongwaveModePropagator when run
using an `ExponentialInput`.
"""
function writeinp(input::ExponentialInput, computejob)
    length(input.output_ranges) == 1 &&
        throw(ArgumentError("`output_ranges` with length 1 is not currently supported."))

    exepath = computejob.exefile
    lwpcpath, _ = splitdir(exepath)
    runname = computejob.runname

    freq = input.frequency/1e3  # in kHz
    max_range = ceil(Int, maximum(input.output_ranges)/1e3)  # convert to km
    diffrange = diff(input.output_ranges)
    drange = trunc(Int, diffrange[1]/1e3)  # in km
    length(unique(round.(diffrange))) == 1 || @info "Using drange = $drange km"
    length(0:drange:max_range) > 1000 && throw(ArgumentError("output_ranges is too long"))

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
        write(f, "tx-data     $(randtransmittername())  $freq  0.0  0.0  1.0  0.0  0.000  0.0", endline)
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
            b_az = rad2deg(mod2pi(input.b_azs[i]))  # deg east of north
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
function writendx(input::ExponentialInput, computejob)
    exepath = computejob.exefile
    lwpcpath, _ = splitdir(exepath)
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
    build_runjob(input::ExponentialInput, computejob; submitjob=true, sleeptime=0.1)

Construct the LWPC files for `input`, and if `submitjob` is true, run LWPC and return results
as a `BasicOutput`.

No matter if `computejob` is parallel or not, `input::ExponentialInput` will be run with a single
process.
"""
function build_runjob(input::ExponentialInput, computejob; submitjob=true, sleeptime=0.1)
    exefile, _ = splitext(computejob.exefile)
    lwpcpath, _ = splitdir(exefile)

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
                sleep(sleeptime)
                dist, amp, phase = readlog(joinpath(lwpcpath, "cases", input.name*".log"))
                dist *= 1e3  # convert to m
                phase .= deg2rad.(phase)  # convert from deg to rad

                # Strip last index of each field because they're 0 (and then inf)
                output.output_ranges = round.(dist, digits=-3)  # fix floating point, rounded to nearest km
                output.amplitude = amp
                output.phase = phase

                completed = true
            else
                sleep(sleeptime)
            end
        end
        if !completed
            @warn "LWPC time limit exceeded. $process"

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
    build_runjob(batchinput::BatchInput, computejob; submitjob=true, sleeptime=0.1)

For each of `batchinput.inputs`, construct the input files and run LWPC, returning results as a
`BatchOutput{BasicOutput}`.

If `computejob` is a `LocalParallel`, it is assumed that there exist directories named
"C:\\LWPCv21_0" to "C:\\LWPCv21_N" where "N" is `computejob.numnodes`. In each directory
it is assumed there is an associated executable "lwpm0.exe" to "lwpmN.exe". The `computejob`
should have `exefile = C:\\LWPCv21\\lwpm.exe`. It is also recommended to update the contents
of the file "C:\\LWPCv21_N\\lwpcDAT.loc" in each of "N" to "C:\\LWPCv21_N\\Data\\".
".

`sleeptime` is the amount of time to `sleep` in seconds between submitting jobs.

!!! note
    The `submitjob` argument is ignored if `computejob` is `LocalParallel`.
"""
function build_runjob(batchinput::BatchInput{ExponentialInput}, computejob::LocalParallel;
    submitjob=true, sleeptime=0.1)

    exefile, exeext = splitext(computejob.exefile)
    lwpcpath, exefilename = splitdir(exefile)

    inputs = batchinput.inputs
    numinputs = length(inputs)

    batch = BatchOutput{BasicOutput}()
    batch.name = batchinput.name
    batch.description = batchinput.description
    batch.datetime = batchinput.datetime
    batch.outputs = Vector{BasicOutput}(undef, numinputs)
    
    completed = falses(numinputs)
    processes = Tuple(ProcessInfo() for i in 1:computejob.numnodes)

    @withprogress name="Batch LWPC" begin
    i = 1  # i: inputidx, ii: proc.inputidx
    while i <= numinputs
        for (pid, proc) in enumerate(processes)
            i > numinputs && break
            
            newlwpcpath = lwpcpath*"_"*string(pid-1)
            newexefile = joinpath(newlwpcpath, exefilename*string(pid-1)*exeext)
            ii = proc.inputidx

            if isnothing(proc.process)
                # Process is available, can assign immediately
                @debug "Process $(pid-1) is `nothing`"

                cj = Local(inputs[i].name, newlwpcpath, newexefile, computejob.walltime)
                build(inputs[i], cj)
                proc.inputidx = i
                start!(proc)
                proc.process = runjob(cj)
                @debug "Input $i submitted"
                i += 1
            elseif process_exited(proc.process)
                # Process is completed
                @debug "Process $(pid-1) has status exited with proc.inputidx: $ii"
                @debug "Trying to read: $(joinpath(newlwpcpath, "cases", inputs[ii].name*".log"))"

                sleep(0.05)
                dist, amp, phase = readlog(joinpath(newlwpcpath, "cases", inputs[ii].name*".log"))
                dist *= 1e3  # convert to m
                phase .= deg2rad.(phase)  # convert from deg to rad

                output = BasicOutput()
                output.name = inputs[ii].name
                output.description = inputs[ii].description
                output.datetime = inputs[ii].datetime

                # Strip last index of each field because they're 0 (and then inf)
                output.output_ranges = round.(dist, digits=-3)  # fix floating point, rounded to nearest km
                output.amplitude = amp
                output.phase = phase

                batch.outputs[ii] = output
                completed[ii] = true
                @logprogress count(completed)/numinputs
                @debug "Input $ii is completed"

                # Start next process
                cj = Local(inputs[i].name, newlwpcpath, newexefile, computejob.walltime)
                build(inputs[i], cj)
                proc.inputidx = i
                start!(proc)
                proc.process = runjob(cj)
                @debug "Input $i submitted"
                i += 1
            elseif elapsed(proc) > computejob.walltime
                # Process timed out
                @warn "LWPC time limit exceeded. $(proc.process)"
                @debug "Process $(pid-1) has exceeded walltime with proc.inputidx: $ii"
                
                kill(proc.process)

                output = BasicOutput()
                output.name = inputs[ii].name
                output.description = inputs[ii].description
                output.datetime = inputs[ii].datetime

                output.output_ranges = [NaN]
                output.amplitude = [NaN]
                output.phase = [NaN]

                batch.outputs[ii] = output
                completed[ii] = true
                @logprogress count(completed)/numinputs
                @debug "Input $ii is completed"

                # Start the next process
                cj = Local(inputs[i].name, newlwpcpath, newexefile, computejob.walltime)
                build(inputs[i], cj)
                proc.inputidx = i
                start!(proc)
                proc.process = runjob(cj)
                @debug "Input $i submitted"
                i += 1
            end
        end
        sleep(sleeptime)  # somehow without this the first process "appears" to time out
    end

    @debug "$(count(completed)) inputs completed"
    @debug "All inputs submitted"

    # Remaining results
    while any(!, completed)
        for (pid, proc) in enumerate(processes)
            isnothing(proc.process) && continue

            newlwpcpath = lwpcpath*"_"*string(pid-1)
            ii = proc.inputidx
            if !completed[ii] && process_exited(proc.process)
                sleep(sleeptime)
                dist, amp, phase = readlog(joinpath(newlwpcpath, "cases", inputs[ii].name*".log"))
                dist *= 1e3  # convert to m
                phase .= deg2rad.(phase)  # convert from deg to rad

                output = BasicOutput()
                output.name = inputs[ii].name
                output.description = inputs[ii].description
                output.datetime = inputs[ii].datetime

                # Strip last index of each field because they're 0 (and then inf)
                output.output_ranges = round.(dist, digits=-3)  # fix floating point, rounded to nearest km
                output.amplitude = amp
                output.phase = phase

                batch.outputs[ii] = output
                completed[ii] = true
                @logprogress count(completed)/numinputs
            elseif elapsed(proc) > computejob.walltime
                # Process timed out
                @warn "LWPC time limit exceeded. $(proc.process)"
                    
                kill(proc.process)

                output = BasicOutput()
                output.name = inputs[ii].name
                output.description = inputs[ii].description
                output.datetime = inputs[ii].datetime

                output.output_ranges = [NaN]
                output.amplitude = [NaN]
                output.phase = [NaN]

                batch.outputs[ii] = output
                completed[ii] = true
                @logprogress count(completed)/numinputs
            end
        end
        sleep(sleeptime)
    end
    end  # withprogress 

    @debug "Final: All $(count(completed)) inputs completed"

    return batch
end

function build_runjob(batchinput::BatchInput{ExponentialInput}, computejob::Local; submitjob=true)
    batch = BatchOutput{BasicOutput}()
    batch.name = batchinput.name
    batch.description = batchinput.description
    batch.datetime = batchinput.datetime
    batch.outputs = Vector{BasicOutput}(undef, length(batchinput.inputs))
    
    for i in eachindex(batchinput.inputs)
        cj = Local(batchinput.inputs[i].name, computejob.rundir, computejob.exefile, computejob.walltime)
        o = build_runjob(batchinput.inputs[i], cj; submitjob=submitjob)
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
