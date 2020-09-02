module EMP2D

using DSP

using ..PropagationModelPrep
import ..LWMS

export emp2d

#==
EMP2D STRUCTS
==#

struct Defaults
    taur::Float64
    tauf::Float64
    I0::Float64
    Ic::Float64
    fcut::Float64
    sourcealt::Float64
    rsspeed::Float64
end
Defaults() = Defaults(20e-6, 50e-6, 10e3, 0, 60000, 5000, -0.75*LWMS.C0)

"""
    Inputs

Fields of `Inputs.dat` for `emp2d-slimfork.cpp`
"""
mutable struct Inputs
    Re::Float64
    dopml_top::Int
    dopml_wall::Int
    doionosphere::Int
    savefields::Vector{Int}  # always size 6
    groundmethod::Int
    maxalt::Float64
    stepalt::Float64
    dr0::Float64
    dr1::Float64
    dr2::Float64
    nground::Int
    range::Float64
    drange::Float64
    dt::Float64
    tsteps::Int
    sig::Float64
    sigm::Float64
    numfiles::Int
    planet::Int
    decfactor::Int
    doDFT::Int
    numDFTfreqs::Int
    DFTfreqs::Vector{Float64}
    read2Dionosphere::Int

    Inputs() = new()
end

"""

"""
function Inputs(Re, maxalt, stepalt, dr0, dr1, dr2, range, drange, DFTfreqs,
    dt=1e-7, decfactor=2, savefields=[1, 0, 0, 0, 0, 0])

    # decfactor: decimate outputs before writing (for 100m, probably want this to be 4 or 5)
    # savefields = [E J H N/A N/A N/A]
    # dt = 1e-7 for D-region up to 150 km, grids >= 100 m

    s = Inputs()

    # Defaults
    setfield!(s, :dopml_top, 1)
    setfield!(s, :dopml_wall, 1)
    setfield!(s, :doionosphere, 1)
    setfield!(s, :groundmethod, 1)  # SIBC
    setfield!(s, :nground, 0)
    setfield!(s, :sig, 0.0)
    setfield!(s, :sigm, 0.0)
    setfield!(s, :numfiles, 50) # in part, used to determine log file steps
    setfield!(s, :planet, 0) # Earth
    setfield!(s, :doDFT, 1)
    setfield!(s, :read2Dionosphere, 1)

    # Arguments
    setfield!(s, :Re, Re)
    setfield!(s, :maxalt, maxalt)  # usually 110e3 m
    setfield!(s, :stepalt, stepalt)  # 50e3 m
    setfield!(s, :dr0, convert(Float64,dr0))  # 100 m
    setfield!(s, :dr1, convert(Float64, dr1))  # 500 m
    setfield!(s, :dr2, convert(Float64, dr2))  # 250 m
    setfield!(s, :range, convert(Float64, range))  # m
    setfield!(s, :drange, convert(Float64, drange))  # 500 m
    setfield!(s, :dt, dt)

    tstepcoeff = 1.1
    maxdist = sqrt(range^2 + maxalt^2)
    tsteps = floor(Int, tstepcoeff*maxdist/LWMS.C0/dt)
    setfield!(s, :tsteps, tsteps)

    length(savefields) == 6 || error("`savefields` must be length 6")
    setfield!(s, :savefields, savefields)
    setfield!(s, :decfactor, decfactor)
    setfield!(s, :numDFTfreqs, length(DFTfreqs))
    setfield!(s, :DFTfreqs, DFTfreqs)

    return s
end

mutable struct Source
    nalt_source::Int
    nt_source::Int
    source::Matrix{Float64}

    Source() = new()
end

function Source(inputs::Inputs)
    s = Source()

    source = create_emp_source(Defaults(), inputs)
    nalt_source, nt_source = size(source)

    setfield!(s, :nalt_source, nalt_source)
    setfield!(s, :nt_source, nt_source)
    setfield!(s, :source, source)

    return s
end

mutable struct Ground
    gsigma::Vector{Float64}
    gepsilon::Vector{Int}

    Ground() = new()
end

#==
FUNCTIONS TO WRITE FILES
==#

function writeinputs(s::Inputs; path="")
    open(joinpath(path,"inputs.dat"), "w") do f
        for field in fieldnames(Inputs)
            write(f, getfield(s, field))
        end
    end
end

function writesource(s::Source; path="")
    open(joinpath(path,"source.dat"), "w") do f
        for field in fieldnames(Source)
            write(f, getfield(s, field))
        end
    end
end

function writeground(s::Ground; path="")
    open(joinpath(path,"ground.dat"), "w") do f
        for field in fieldnames(Ground)
            write(f, getfield(s, field))
        end
    end
end

function writebfield(Bmag, in::Inputs; path="")
    thmax = in.range/in.Re
    dth = in.drange/in.Re
    hh = round(Int, thmax/dth) + 1

    Br = fill(convert(Float64, Bmag), hh)
    Bt = zeros(hh)
    Bp = zeros(hh)

    open(joinpath(path,"B0.dat"), "w") do f
        write(f, Br)
        write(f, Bt)
        write(f, Bp)
    end
end

function writene(ne; path="")
    open(joinpath(path,"ne.dat"), "w") do f
        write(f, ne)
    end
end

function writeni(ne; path="")
    ni = copy(ne)
    ni[ni .< 100e6] .= 100e6

    open(joinpath(path,"ni.dat"), "w") do f
        write(f, ni)
    end
end

# function writend(in::Inputs)
#     r, dr = generate_rdr(in)
#     nd = MSISatmosphere((r-in.Re)/1000)
#     ndt = nd.total*1e6
#
#     open("nd.dat", "w") do f
#         write(f, ndt)
#     end
# end

function writenu(nu; path="")
    open(joinpath(path,"nu.dat"), "w") do f
        write(f, nu)
    end
end

#==
MAIN FUNCTIONS
==#

function emp2d(file::AbstractString, computejob::ComputeJob; runjob=true)
    isfile(file) || ispath(file) || error("$file is not a valid file name")

    s = LWMS.parse(file)
    buildandrun(s, computejob; runjob=runjob)

    return nothing
end

"""
    buildandrun()

With default (coarse) Inputs.
"""
function buildandrun(s::LWMS.BasicInput, computejob::ComputeJob; runjob=true)

    all(s.b_dip .≈ 90) || @warn "Segment magnetic field is not vertical"
    length(unique(s.b_mag)) == 1 || @warn "Magnetic field is not homogeneous"

    ne_threshold = 3e9

    # Useful values
    max_range = last(s.output_ranges)
    # round up to nearest thousand km and go 1000 km beyond that
    max_range = round(max_range+1000e3, digits=-3, RoundUp)

    drange = 500  # m
    Nrange = round(Int, max_range/drange) + 1
    rangevec = range(0, max_range, length=Nrange)

    inputs = Inputs(LWMS.EARTH_RADIUS, 110e3, 50e3, 100, 500, 250, max_range, drange, [s.frequency])
    r, dr = generate_rdr(inputs)
    altitudes = r .- LWMS.EARTH_RADIUS

    # Fill in values for each waveguide segment
    ne = Matrix{Float64}(undef, length(r), Nrange)
    nu = similar(ne)
    gsigma = Vector{Float64}(undef, Nrange)
    gepsilon = Vector{Int}(undef, Nrange)

    for i in eachindex(s.segment_ranges)
        segment_begin_idx = findfirst(x->x==s.segment_ranges[i], rangevec)
        if i == lastindex(s.segment_ranges)
            segment_end_idx = Nrange
        else
            segment_end_range = s.segment_ranges[i+1]
            segment_end_idx = findfirst(x->x==segment_end_range, rangevec)
        end

        # Electron density profile
        neprofile = LWMS.waitprofile.(altitudes, s.hprimes[i], s.betas[i], cutoff_low=50e3, threshold=3e9)
        neprofile[neprofile .> ne_threshold] .= ne_threshold
        ne[:,segment_begin_idx:segment_end_idx] .= neprofile

        # Electron collision profile
        nuprofile = LWMS.electroncollisionfrequency.(altitudes)
        nu[:,segment_begin_idx:segment_end_idx] .= nuprofile

        # Ground profile
        gsigma[segment_begin_idx:segment_end_idx] .= s.ground_sigmas[i]
        gepsilon[segment_begin_idx:segment_end_idx] .= s.ground_epsr[i]
    end

    source = Source(inputs)
    ground = Ground()
    ground.gsigma = gsigma
    ground.gepsilon = gepsilon

    rundir = computejob.rundir
    exefile = computejob.exefile

    writeinputs(inputs, path=rundir)
    writesource(source, path=rundir)
    writeground(ground, path=rundir)
    writebfield(s.b_mag[1], inputs, path=rundir)
    writene(ne, path=rundir)
    writeni(ne, path=rundir)
    writenu(nu, path=rundir)

    shfile = writeshfile(computejob)

    if runjob
        run(`cp $exefile $rundir`)

        jobname = read(`sbatch $shfile`, String)
        jobid = strip(jobname)

        println("Job $jobid submitted!\n")
    end

    return nothing
end

function create_emp_source(def::Defaults, in::Inputs)
    Iin = zeros(in.tsteps)

    # Current waveform: linear rise, exp decay
    for t  = 1:in.tsteps
        if t*in.dt < def.taur
            Iin[t] = def.I0*t*in.dt/def.taur
        else
            Iin[t] = def.I0*exp(-(t*in.dt-def.taur)^2/def.tauf^2) + def.Ic*(1 - exp(-(t*in.dt-def.taur)^2/def.tauf^2))
        end
    end

    # Filter waveform
    fc = def.fcut/(1/in.dt/2);  # cutoff frequency
    order = 41; # filter order

    window = hamming(order)
    filter = digitalfilter(Lowpass(fc), FIRWindow(window));

    Iin2 = filt(filter, Iin);
    Iin2 .= Iin2./maximum(Iin2).*def.I0  # rescaled to peak current

    # Fix delays
    Iin = [Iin2; zeros(2*round(Int,def.sourcealt/abs(def.rsspeed)/in.dt))]

    Iin[Iin .< 1e-20] .= 1e-20

    # Spatial variation
    satop = in.nground + floor(Int, def.sourcealt/in.dr1) + 1

    source = permutedims(repeat(Iin, 1, satop))

    return source
end

function generate_rdr(in::Inputs)
    rr = convert(Int, in.stepalt/in.dr1 + (in.maxalt - in.stepalt)/in.dr2 + 1)

    r = Vector{Float64}(undef, rr)
    dr = Vector{Float64}(undef, rr-1)

    r[1] = in.Re
    for i = 2:rr
        if r[i-1] < (in.Re + in.stepalt)
            r[i] = r[i-1] + in.dr1
        else
            r[i] = r[i-1] + in.dr2
        end
        dr[i-1] = r[i] - r[i-1]
    end
    return r, dr
end

end  # module
