using Test, Dates, Distributed
using JSON3
using LongwaveModeSolver
const LWMS = LongwaveModeSolver

using PropagationModelPrep

function generatejson()
    # Waveguide definition
    segment_ranges = [0]
    hprimes = [75]
    betas = [0.32]
    b_mag = fill(50e-6, length(segment_ranges))  # TODO: rename to be consistent (plural?)
    b_dip = fill(90.0, length(segment_ranges))
    b_az = fill(0.0, length(segment_ranges))
    ground_sigmas = [0.001]
    ground_epsrs = [15]

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
    input.b_mags = b_mag
    input.b_dips = b_dip
    input.b_azs = b_az
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

function test_lwpclocal()
    scenarioname = "homogeneous1"
    computejob = Local(scenarioname, ".", "C:\\LWPCv21\\lwpm.exe")
    LWPC.run(scenarioname*".json", computejob)
    LWPC.process(scenarioname*".json", computejob)
end

function test_lwpclocalparallel()
    scenarioname = "batchbasic"
    computejob = LocalParallel(scenarioname, ".", "C:\\LWPCv21\\lwpm.exe", 4)

    LWPC.run(scenarioname*".json", computejob)
end

function test_lwpcparallellocal()
    nprocs() < 4 && addprocs(4 - nprocs())

    scenarioname = ""
    computejob = LocalParallel(scenarioname, ".", "C:\\LWPCv21\\lwpm.exe")

    jobs = RemoteChannel(()->Channel{Int}(32));
    results = RemoteChannel(()->Channel{Tuple}(32));

    @everywhere function do_work(jobs, results) # define work function everywhere
        while true
            job_id = take!(jobs)
            exec_time = 5rand()
            put!(results, (job_id, exec_time, myid()))
        end
    end

    function make_jobs(n)
        for i in 1:n
            put!(jobs, i)
        end
    end

    n = 12

    @async make_jobs(n) # feed the jobs channel with "n" jobs

    for p in workers() # start tasks on the workers to process requests in parallel
        remote_do(do_work, p, jobs, results)
    end

    @elapsed while n > 0 # print out results
        job_id, exec_time, where = take!(results)
        println("$job_id finished in $(round(exec_time; digits=2)) seconds on worker $where")
        global n = n - 1
    end

end

@testset "PropagationModelPrep" begin
    generatejson()
    generate_batchbasic()

    @testset "Utils" begin
        @test PropagationModelPrep.rounduprange(2314.2e3) == 4000e3
        @test test_unwrap()
    end

    # TODO: Toggle EMP2D and LWPC tests
    @testset "EMP2D" begin
        @test test_emp2dinputs()

        @test_logs (:info,
            "Updating computejob runname to homogeneous1") test_mismatchedrunnames()
        @test_logs (:info,
            "Running in homogeneous2/homogeneous1/") (:info,
            "Creating homogeneous2/homogeneous1/") test_newrundir()
        @test_throws ErrorException test_filename_error()
    end

    @testset "LWPC" begin
        @test test_lwpclocal()
    end

    # Cleanup
    rm("homogeneous1", force=true, recursive=true)
    rm("homogeneous2", force=true, recursive=true)
    isfile("batchbasic.json") && rm("batchbasic.json")
end
