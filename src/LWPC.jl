module LWPC

using JSON3, CSV, DataFrames

using ..PropagationModelPrep
import ..LWMS


"""
    run(file, computejob::ComputeJob, inputs=false; submitjob=true)

Generate the input files and attempt to run the emp2d code for the scenario
described by `file` as `computejob`.

Optionally provide an inputs::Inputs() struct. Otherwise default values are used.
"""
function run(file, computejob::ComputeJob; inputs=nothing, submitjob=true)
    isfile(file) || error("$file is not a valid file name")

    s = LWMS.parse(file)

    if computejob.runname != s.name
        @info "Updating computejob runname to $(s.name)"
        computejob.runname = s.name
    end

    if splitdir(computejob.rundir)[2] != computejob.runname
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

    build(s, computejob, inputs)

    if submitjob
        # rm("Output\\")
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

    writeinp(s, computejob; path=rundir)
    writendx(s, computejob; path=rundir)

    return nothing
end

function writeinp(s::LWMS.BasicInput, computejob::ComputeJob; path="")
    runname = computejob.runname
    freq = s.frequency
    max_range = rounduprange(maximum(s.output_ranges))/1e3  # convert to km
    diffrange = diff(s.output_ranges)
    drange = diffrange[1]/1e3  # in km
    length(unique(round(diff, digits=3))) == 1 || @info "Using drange = $drange km"

    endline = "\n"
    open(joinpath(path, runname*".inp")) do f
        write(f, "file-mds    Output\\", endline)
        write(f, "file-lwf    Output\\", endline)
        write(f, "file-prf    cases\\", endline)
        write(f, "file-ndx    cases\\", endline)
        write(f, "case-id     $runname", endline)
        write(f, "tx          $runname", endline)
        write(f, "tx-data     vlftx  $freq  0.0  0.0  100.0  0.0  0.000  0.0", endline)
        write(f, "receivers   0.0000  0.0000", endline)
        write(f, "range-max   $max_range", endline)
        write(f, "ionosphere  range exponential $runname", endline)
        write(f, "rx-data     vertical 0", endline)
        write(f, "lwf-vs-distance $max_range $drange", endline)
        write(f, "mc-options  full-wave 50 true", endline)
        write(f, "print-lwf   1", endline)
        write(f, "lwflds", endline)
        write(f, "preseg", endline)
        for i in eachindex(s.segment_ranges)
            r = s.segment_ranges[i]/1e3  # dist in km
            b_az = s.b_az[i]  # deg east of north
            b_dip = s.b_dip[i]  # deg from horizontal
            b_mag = s.b_mag[i]*1e4  # XXX: supposedly Tesla, but maybe Gauss?
            gsigma = s.ground_sigmas[i]
            gepsr = s.ground_epsr[i]
            beta = s.betas[i]
            hprime = s.hprimes[i]

            write(f, " $r,$b_az,$b_dip,$b_mag,,$gsigma,$gepsr,,$beta,$hprime", endline)
        end
        write(" 40000", endline)
        write("start", endline)
        write("quit")
    end

    return nothing
end

function writendx(s::LWMS.BasicInput, computejob::ComputeJob; path="")
    runname = computejob.runname

    endline = "\n"
    open(joinpath(path, runname*".ndx")) do f
        for i in eachindex(s.segment_ranges)
            r = s.segment_ranges[i]/1e3  # dist in km
            beta = s.betas[i]
            hprime = s.hprimes[i]

            write(f, "$r $beta $hprime", endline)
        end
    end

    return nothing
end

function runjob(computejob::LairLwpc)
    out = read(`lwms`, String)

    println("Job submitted!\n")
    println(out)

    return nothing
end

function readlog(file)
    raw = CSV.File(file;
                   skipto=40, delim=' ', ignorerepeated=true, header=false)

    dist = vcat(raw.Column1, raw.Column4, raw.Column7)
    amp = vcat(raw.Column2, raw.Column5, raw.Column8)
    phase = vcat(raw.Column3, raw.Column6, raw.Column9)

    return dist, amp, phase
end

function process(path)
    fullpath = abspath(path)
    pathname = splitpath(fullpath)[end]  # proxy for filename

    jsonfile = joinpath(fullpath, pathname*".json")

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

    dist, amp, phase = readlog(joinpath(fullpath, name*".log"))
    dist *= 1e3  # convert to m

    delidxs = ismissing.(dist)  # assuming these are at the end
    delidxs[findfirst(delidxs)-1] = true  # last valid entry is at receiver distance
    deleteat!(dist, delidxs)
    deletat!(amp, delidxs)
    deleteat!(phase, delidxs)

    output = LWMS.BasicOutput()
    output.name = name
    output.description = description
    output.datetime = datetime

    # Strip last index of each field because they're 0 (and then inf)
    output.output_ranges = round.(dist, digits=3)  # fix floating point
    output.amplitude = amp
    output.phase = phase

    # Save output
    json_str = JSON3.write(output)

    open(joinpath(path,name*"_lwpc.json"), "w") do f
        write(f, json_str)
    end

    return output
end

end  # module
