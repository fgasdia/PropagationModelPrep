using Test, Dates
using JSON3
using LongwaveModeSolver
const LWMS = LongwaveModeSolver

using PropagationModelPrep

function generatejson()
    # Waveguide definition]
    segment_ranges = [0]
    hprimes = [75]
    betas = [0.32]
    b_mag = fill(50e-6, length(segment_ranges))  # TODO: rename to be consistent (plural?)
    b_dip = fill(90.0, length(segment_ranges))
    b_az = fill(0.0, length(segment_ranges))
    ground_sigmas = [0.001]
    ground_epsr = [15]

    # Transmitter
    frequency = 24e3

    # Outputs
    output_ranges = collect(0:5e3:2500e3)

    input = LWMS.BasicInput()
    input.name = "homogeneous1"
    input.description = "homogeneous ionosphere"
    input.datetime = Dates.now()

    input.segment_ranges = segment_ranges
    input.hprimes = hprimes
    input.betas = betas
    input.b_mag = b_mag
    input.b_dip = b_dip
    input.b_az = b_az
    input.ground_sigmas = ground_sigmas
    input.ground_epsr = ground_epsr
    input.frequency = frequency
    input.output_ranges = output_ranges

    json_str = JSON3.write(input)

    open("$(input.name).json","w") do f
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

    issorted(p) == false

    PropagationModelPrep.unwrap!(p)
    issorted(p) == true || return false

    return true
end

function test_mismatchedrunnames()
    computejob = Summit("homogeneous2", "homogeneous1", 12, "01:00:00", "dummy_exe")
    EMP2D.run("homogeneous1.json", computejob; submitjob=false)
end

function test_newrundir()
    computejob = Summit("homogeneous1", "homogeneous2", 12, "01:00:00", "dummy_exe")
    EMP2D.run("homogeneous1.json", computejob; submitjob=false)
end

function test_filename_error()
    computejob = Summit("homogeneous1", "homogeneous1", 12, "01:00:00", "dummy_exe")
    EMP2D.run("homogeneous999.json", computejob; submitjob=false)
end

function test_emp2dinputs()
    computejob = Summit("homogeneous1", "homogeneous1", 12, "01:00:00", "dummy_exe")
    inputs = EMP2D.Inputs(6366e3, 110e3, 50e3, 200, 100, 4000e3, 100, [24e3])
    EMP2D.run("homogeneous1.json", computejob; inputs=inputs, submitjob=false)

    testinputs = EMP2D.readinputs("homogeneous1")
    for field in fieldnames(EMP2D.Inputs)
        getfield(inputs, field) == getfield(testinputs, field) || return false
    end
    return true
end

function delete_tmpdirs()
    rm("homogeneous1", force=true, recursive=true)
    rm("homogeneous2", force=true, recursive=true)
end

@testset "PropagationModelPrep" begin
    @info "Generating json input file"
    generatejson()

    @testset "Utils" begin
        @test PropagationModelPrep.rounduprange(2314.2e3) == 4000e3
        @test test_unwrap()
    end

    @testset "EMP2D" begin
        @test test_emp2dinputs()

        @test_logs (:info,
            "Updating computejob runname to homogeneous1") test_mismatchedrunnames()
        @test_logs (:info,
            "Running in homogeneous2/homogeneous1/") (:info,
            "Creating homogeneous2/homogeneous1/") test_newrundir()
        @test_throws ErrorException test_filename_error()
    end

    # Cleanup
    delete_tmpdirs()
end
