using Test, Dates, Distributed
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

    input = LMP.BasicInput()
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

    binput = BatchInput{BasicInput}()
    binput.name = "batchbasic"
    binput.description = "Test BatchInput with BasicInput"
    binput.datetime = Dates.now()

    inputs = Vector{BasicInput}(undef, N)
    for i in eachindex(inputs)
        input = BasicInput()

        input.name = "$i"
        input.description = "BasicInput $i"
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
    computejob = Local(scenarioname, ".", "C:\\LWPCv21\\lwpm.exe")
    LWPC.run(scenarioname*".json", computejob)
    LWPC.process(scenarioname*".json", computejob)
end

@testset "PropagationModelPrep" begin
    generatejson()
    generate_batchbasic()

    @testset "Utils" begin
        @test PropagationModelPrep.rounduprange(2314.2e3) == 4000e3
        test_unwrap()
    end

    # TODO: Toggle EMP2D and LWPC tests
    @testset "EMP2D" begin
        test_emp2dinputs()

        @test_logs (:info,
            "Updating computejob runname to homogeneous1") mismatchedrunnames()

        # @info prints with "\" but joinpath has "\\" we effectively "strip" it with printf
        rundir = joinpath("homogeneous2", "homogeneous1", "")
        @test_logs (:info, "Running in $rundir") (:info, "Creating $rundir") newrundir()
        @test_throws ErrorException filename_error()
    end

    @testset "LWPC" begin
        # @test test_lwpclocal()
        # TODO
    end

    # Cleanup
    rm("homogeneous1", force=true, recursive=true)
    rm("homogeneous2", force=true, recursive=true)
    isfile("batchbasic.json") && rm("batchbasic.json")
    isfile("homogeneous1.json") && rm("homogeneous1.json")
    isfile("homogeneous1_lwpc.json") && rm("homogeneous1_lwpc.json")
end
