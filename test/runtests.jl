using Test, Dates
using JSON3
using LongwaveModeSolver
const LWMS = LongwaveModeSolver

using PropagationModelPrep

function generatehomogeneous()
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

# Generate files
generatehomogeneous()

# Read files
computejob = Summit("homogeneous1", ".", 12, "01:00:00", "/projects/emp/emp2/emp2d")
emp2d("homogeneous1.json", computejob; submitjob=false)

# Run file
computejob = LocalOMP("homogeneous1", ".", 2, "../emp2d")
emp2d("homogeneous1.json", computejob)

@testset "PropagationModelPrep" begin
    #
end
