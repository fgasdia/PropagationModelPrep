using Test, Dates, Statistics
using JSON3
using LongwaveModePropagator
const LMP = LongwaveModePropagator

using PropagationModelPrep

function generatejson()
    # Waveguide definition
    segment_ranges = [0]
    hprimes = [75]
    betas = [0.32]
    b_mags = fill(50e-6, length(segment_ranges))
    b_dips = fill(π/2.0, length(segment_ranges))
    b_azs = fill(0.0, length(segment_ranges))
    ground_sigmas = [0.001]
    ground_epsrs = [15]

    # Transmitter
    frequency = 24e3

    # Outputs
    output_ranges = collect(0:5e3:2500e3)

    input = LMP.ExponentialInput()
    input.name = "homogeneous1"
    input.description = "homogeneous ionosphere"
    input.datetime = Dates.now()

    input.segment_ranges = segment_ranges
    input.hprimes = hprimes
    input.betas = betas
    input.b_mags = b_mags
    input.b_dips = b_dips
    input.b_azs = b_azs
    input.ground_sigmas = ground_sigmas
    input.ground_epsrs = ground_epsrs
    input.frequency = frequency
    input.output_ranges = output_ranges

    json_str = JSON3.write(input)

    open("$(input.name).json","w") do f
        write(f, json_str)
    end

    return nothing
end

function generate_batchbasic()
    N = 2
    rep(v) = repeat(v, 1, N)

    # Waveguide definition
    nsegments = 4
    segment_ranges = [0, 500e3, 1000e3, 1500e3]
    b_mags = fill(50e-6, nsegments)
    b_dips = fill(π/2, nsegments)
    b_azs = fill(0.0, nsegments)
    ground_sigmas = [0.001, 0.001, 0.0005, 0.0001]
    ground_epsrs = [4, 4, 10, 10]

    hprimes = rep([72.0, 74, 75, 76])
    betas = rep([0.28, 0.30, 0.32, 0.35])

    # Transmitter
    frequency = 24e3

    # Outputs
    output_ranges = collect(0:20e3:2000e3)

    binput = BatchInput{ExponentialInput}()
    binput.name = "batchbasic"
    binput.description = "Test BatchInput with ExponentialInput"
    binput.datetime = Dates.now()

    inputs = Vector{ExponentialInput}(undef, N)
    for i in eachindex(inputs)
        input = ExponentialInput()

        input.name = "batchbasic"
        input.description = "ExponentialInput $i"
        input.datetime = binput.datetime
        input.segment_ranges = segment_ranges
        input.hprimes = hprimes[:,i]
        input.betas = betas[:,i]
        input.b_mags = b_mags
        input.b_dips = b_dips
        input.b_azs = b_azs
        input.ground_sigmas = ground_sigmas
        input.ground_epsrs = ground_epsrs
        input.frequency = frequency
        input.output_ranges = output_ranges

        inputs[i] = input
    end

    binput.inputs = inputs

    json_str = JSON3.write(binput)

    open("batchbasic.json","w") do f
        write(f, json_str)
    end

    return nothing
end

function largebatchinput()
    N = 24

    # Waveguide definition
    segment_ranges = [0.0]
    b_mags = [50e-6]
    b_dips = [π/2]
    b_azs = [0.0]
    ground_sigmas = [0.001]
    ground_epsrs = [10]

    hprimes = [72.0]
    betas = [0.30]

    # Transmitter
    frequency = 24e3

    # Outputs
    output_ranges = collect(0:20e3:2000e3)

    binput = BatchInput{ExponentialInput}()
    binput.name = "largebatchinput"
    binput.description = "Test BatchInput with ExponentialInput"
    binput.datetime = Dates.now()

    inputs = Vector{ExponentialInput}(undef, N)
    for i in eachindex(inputs)
        input = ExponentialInput()

        input.name = "batch"
        input.description = "ExponentialInput $i"
        input.datetime = binput.datetime
        input.segment_ranges = segment_ranges
        input.hprimes = hprimes
        input.betas = betas
        input.b_mags = b_mags
        input.b_dips = b_dips
        input.b_azs = b_azs
        input.ground_sigmas = ground_sigmas
        input.ground_epsrs = ground_epsrs
        input.frequency = frequency
        input.output_ranges = output_ranges

        inputs[i] = input
    end
    binput.inputs = inputs

    return binput
end

function test_unwrap()
    # Generate parametric spiral
    t = range(0, 6π, length=201)
    x = t/π .* cos.(t)
    y = t/π .* sin.(t)

    # Extract phase angle of spiral
    p = atan.(y,x)

    @test !issorted(p)

    PropagationModelPrep.unwrap!(p)
    @test issorted(p)
end

function mismatchedrunnames()
    computejob = Summit("homogeneous2", "homogeneous1", "dummy_exe", 12, "01:00:00")
    EMP2D.run("homogeneous1.json", computejob; submitjob=false)
end

function newrundir()
    computejob = Summit("homogeneous1", "homogeneous2", "dummy_exe", 12, "01:00:00")
    EMP2D.run("homogeneous1.json", computejob; submitjob=false)
end

function filename_error()
    computejob = Summit("homogeneous1", "homogeneous1", "dummy_exe", 12, "01:00:00")
    EMP2D.run("homogeneous999.json", computejob; submitjob=false)
end

function test_emp2dinputs()
    computejob = Summit("homogeneous1", "homogeneous1", "dummy_exe", 12, "01:00:00")
    inputs = EMP2D.Inputs(6366e3, 110e3, 50e3, 200, 100, 4000e3, 100, [24e3])
    EMP2D.run("homogeneous1.json", computejob; inputs=inputs, submitjob=false)

    testinputs = EMP2D.readinputs("homogeneous1")
    for field in fieldnames(EMP2D.Inputs)
        @test getfield(inputs, field) == getfield(testinputs, field)
    end
end

function test_lwpclocal()
    scenarioname = "homogeneous1"
    computejob = Local(scenarioname, ".", "C:\\LWPCv21\\lwpm.exe", 90)

    exepath = computejob.exefile
    lwpcpath, _ = splitdir(exepath)

    logfile = joinpath(lwpcpath, "cases", scenarioname*".log")
    inpfile = joinpath(lwpcpath, "cases", scenarioname*".inp")
    ndxfile = joinpath(lwpcpath, "cases", scenarioname*".ndx")
    lwffile = joinpath(lwpcpath, "Output", scenarioname*".lwf")
    mdsfile = joinpath(lwpcpath, "Output", scenarioname*".mds")

    output = LWPC.run(scenarioname*".json", computejob)
    @test output isa BasicOutput

    # Compare to LMP
    loutput = propagate(scenarioname*".json")
    ares = mean(abs, output.amplitude .- loutput.amplitude) < 0.2
    pres = rad2deg(mean(abs, output.phase .- loutput.phase)) < 1

    if !ares && !pres
        @warn "Comparison to LongwaveModePropagator failed. This is expected for LongwaveModePropagator `v0.2`."
    end
    @test ares
    @test pres

    # Test `deletefiles`
    @test isfile(logfile)
    @test isfile(lwffile)
    @test isfile(mdsfile)
    @test isfile(inpfile)
    @test isfile(ndxfile)
    LWPC.deletefiles(computejob)
    @test !isfile(logfile)
    @test !isfile(lwffile)
    @test !isfile(mdsfile)
    @test isfile(inpfile)  # shouldn't be deleted
    @test isfile(ndxfile)  # shouldn't be deleted

    # Test `build` (and indirectly, `writendx` and `writeinp`)
    rm(logfile; force=true)
    s = LMP.parse(scenarioname*".json")
    LWPC.build(s, computejob)
    @test isfile(ndxfile)
    @test isfile(inpfile)
    rm(inpfile; force=true)
    rm(ndxfile; force=true)

    # Test `runjob`
    LWPC.deletefiles(computejob)
    LWPC.build(s, computejob)
    process = LWPC.runjob(computejob)
    sleep(5)
    d, a, p = LWPC.readlog(logfile)
    @test length(d) == length(a) == length(p) == length(s.output_ranges)

    # Test `build_runjob`
    rm(inpfile; force=true)
    rm(ndxfile; force=true)
    output = LWPC.build_runjob(s, computejob; submitjob=true)
    @test isfile(logfile)
    d, a, p = LWPC.readlog(logfile)
    @test output.output_ranges/1000 == d
    @test output.amplitude == a
    @test output.phase == deg2rad.(p)

    rm(inpfile; force=true)
    rm(ndxfile; force=true)
    LWPC.deletefiles(computejob)
    output = LWPC.build_runjob(s, computejob; submitjob=false)
    @test isnothing(output)
    @test isfile(inpfile)
    @test isfile(ndxfile)
end

function test_lwpclocal_batch()
    scenarioname = "batchbasic"
    computejob = Local(scenarioname, ".", "C:\\LWPCv21\\lwpm.exe", 90)

    exepath = computejob.exefile
    lwpcpath, _ = splitdir(exepath)

    inpfile = joinpath(lwpcpath, "cases", scenarioname*".inp")
    ndxfile = joinpath(lwpcpath, "cases", scenarioname*".ndx")

    s = LMP.parse(scenarioname*".json")

    # Test `build_runjob`
    rm(inpfile; force=true)
    rm(ndxfile; force=true)
    o = LWPC.build_runjob(s, computejob; submitjob=true)
    @test o isa BatchOutput
    @test length(o.outputs) == 2

    for i in eachindex(s.inputs)
        d, a, p = o.outputs[i].output_ranges, o.outputs[i].amplitude, o.outputs[i].phase
        @test length(d) == length(a) == length(p) == length(s.inputs[i].output_ranges)
    end
end

function test_lwpclocalparallel_batch(numnodes)
    s = largebatchinput()
    scenarioname = s.name
    computejob = LocalParallel(scenarioname, ".", "C:\\LWPCv21\\lwpm.exe", numnodes, 90)

    # Test `build_runjob`
    o = LWPC.build_runjob(s, computejob; submitjob=true)
    @test o isa BatchOutput
    @test length(o.outputs) == 24

    for i in eachindex(s.inputs)
        d, a, p = o.outputs[i].output_ranges, o.outputs[i].amplitude, o.outputs[i].phase
        @test length(d) == length(a) == length(p) == length(s.inputs[i].output_ranges)
    end

    # Test `run` against `build_runjob`
    o2 = LWPC.run(s, computejob)

    for field in (:name, :description, :datetime)
        @test getfield(o, field) == getfield(o2, field)
    end

    for i in eachindex(o.outputs)
        for field in fieldnames(BasicOutput)
            @test getfield(o.outputs[i], field) == getfield(o2.outputs[i], field)
        end
    end
end

@testset "PropagationModelPrep" begin
    generatejson()
    generate_batchbasic()

    @testset "Utils" begin
        @test PropagationModelPrep.rounduprange(2314.2e3) == 4000e3
        test_unwrap()
    end

    @testset "EMP2D" begin
        @info "Testing EMP2D"
        test_emp2dinputs()

        @test_logs (:info,
            "Updating computejob runname to homogeneous1") mismatchedrunnames()

        # @info prints with "\" but joinpath has "\\" we effectively "strip" it with printf
        rundir = joinpath("homogeneous2", "homogeneous1", "")
        @test_logs (:info, "Running in $rundir") (:info, "Creating $rundir") newrundir()
        @test_throws ErrorException filename_error()
    end


    if Sys.iswindows() && isfile("C:\\LWPCv21\\lwpm.exe")
        @testset "LWPC" begin
            @info "Testing LWPC"
            test_lwpclocal()
            test_lwpclocal_batch()
        end
    end

    # Cleanup
    rm("homogeneous1", force=true, recursive=true)
    rm("homogeneous2", force=true, recursive=true)
    isfile("batchbasic.json") && rm("batchbasic.json")
    isfile("homogeneous1.json") && rm("homogeneous1.json")
    isfile("homogeneous1_output.json") && rm("homogeneous1_output.json")
    isfile("homogeneous1_lwpc.json") && rm("homogeneous1_lwpc.json")
    isfile("largebatchinput_lwpc.json") && rm("largebatchinput_lwpc.json")
end
