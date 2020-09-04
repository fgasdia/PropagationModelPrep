module EMP2D

using DSP

using ..PropagationModelPrep
import ..LWMS

export emp2d

"""
    Defaults

Struct to hold some default values for EMP2D. Instantiated as `Defaults()`.

Currently this only holds values related to the emp source.
"""
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

Mutable struct that holds fields of `Inputs.dat` for `emp2d-slimfork.cpp`.
"""
mutable struct Inputs
    Re::Float64
    dopml_top::Int32
    dopml_wall::Int32
    doionosphere::Int32
    savefields::Vector{Int32}  # always size 6
    groundmethod::Int32
    maxalt::Float64
    stepalt::Float64
    dr0::Float64
    dr1::Float64
    dr2::Float64
    nground::Int32
    range::Float64
    drange::Float64
    dt::Float64
    tsteps::Int32
    sig::Float64
    sigm::Float64
    numfiles::Int32
    planet::Int32
    decfactor::Int32
    doDFT::Int32
    numDFTfreqs::Int32
    DFTfreqs::Vector{Float64}
    read2Dionosphere::Int32

    Inputs() = new()
end

"""
    Inputs(Re, maxalt, stepalt, dr1, dr2, range, drange, DFTfreqs,
        dt=1e-7, decfactor=2, savefields=[0, 0, 0, 0, 0, 0])

Outer constructor for the `Inputs` mutable struct.

This constructor fills in default values for the additional fields that are not
arguments of this function. All arguments will be `convert`ed to the appropriate
type as necessary.

Arguments:

    - `Re` is Earth radius
    - `maxalt` is altitude of top of ionosphere (typically 110e3)
    - `stepalt` is altitude of change of grid resolution (typically 50e3)
    - `dr1` is ionosphere grid cell resolution below `stepalt` (typically 500)
    - `dr2` is ionosphere grid cell resolution above `stepalt` (typically 250)
    - `range` is the maximum horizontal ground range to run the simulation to
    - `drange` is the grid cell resolution along the range direction (typically 500)
    - `DFTfreqs` is a vector of frequencies to examine
    - `dt=1e-7` is the timestep
    - `decfactor=2` is the decimation factor for output fields. Larger values
        are greater decimation
    - `savefields=[0, 0, 0, 0, 0, 0]` turns on output for E, J, and H when any
        entry is greater than 0. The last 3 fields are not supported.

!!! note

    All units are SI (mks).
"""
function Inputs(Re, maxalt, stepalt, dr1, dr2, range, drange, DFTfreqs,
    dt=Float64(1e-7), decfactor=Int32(2), savefields=Int32[1, 0, 0, 0, 0, 0])

    # decfactor: decimate outputs before writing (for 100m, probably want this to be 4 or 5)
    # savefields = [E J H N/A N/A N/A]
    # dt = 1e-7 for D-region up to 150 km, grids >= 100 m

    s = Inputs()

    # Defaults
    setfield!(s, :dopml_top, Int32(1))
    setfield!(s, :dopml_wall, Int32(1))
    setfield!(s, :doionosphere, Int32(1))
    setfield!(s, :groundmethod, Int32(1))  # SIBC
    setfield!(s, :dr0, Float64(100))  # not actually used b/c nground = 0
    setfield!(s, :nground, Int32(0))
    setfield!(s, :sig, Float64(0))
    setfield!(s, :sigm, Float64(0))
    setfield!(s, :numfiles, Int32(50)) # in part, used to determine log file steps
    setfield!(s, :planet, Int32(0)) # Earth
    setfield!(s, :doDFT, Int32(1))
    setfield!(s, :read2Dionosphere, Int32(1))

    # Arguments
    setfield!(s, :Re, convert(Float64, Re))
    setfield!(s, :maxalt, convert(Float64, maxalt))  # usually 110e3 m
    setfield!(s, :stepalt, convert(Float64, stepalt))  # 50e3 m
    setfield!(s, :dr1, convert(Float64, dr1))  # 500 m
    setfield!(s, :dr2, convert(Float64, dr2))  # 250 m
    setfield!(s, :range, convert(Float64, range))  # m
    setfield!(s, :drange, convert(Float64, drange))  # 500 m
    setfield!(s, :dt, convert(Float64, dt))

    tstepcoeff = 1.1
    maxdist = sqrt(range^2 + maxalt^2)
    tsteps = floor(Int32, tstepcoeff*maxdist/LWMS.C0/dt)
    setfield!(s, :tsteps, Int32(tsteps))

    length(savefields) == 6 || error("`savefields` must be length 6")
    if any(x->x>0, savefields[4:6])
        @info "Setting `savefields[4:6]` to 0"
        savefields[4:6] .= 0
    end
    setfield!(s, :savefields, convert(Vector{Int32}, savefields))
    setfield!(s, :decfactor, Int32(decfactor))
    setfield!(s, :numDFTfreqs, Int32(length(DFTfreqs)))
    setfield!(s, :DFTfreqs, convert(Vector{Float64}, DFTfreqs))

    return s
end

"""
    Source

Mutable struct for the source.dat input to emp2d-slimfork.cpp.
"""
mutable struct Source
    nalt_source::Int32
    nt_source::Int32
    source::Matrix{Float64}

    Source() = new()
end

"""
    Source(inputs::Inputs)

Outer constructor for the `Source` struct. This will automatically generate an
emp source suitable for emp2d given an `Inputs` struct.
"""
function Source(inputs::Inputs)
    s = Source()

    source = create_emp_source(Defaults(), inputs)
    nalt_source, nt_source = size(source)

    setfield!(s, :nalt_source, Int32(nalt_source))
    setfield!(s, :nt_source, Int32(nt_source))
    setfield!(s, :source, convert(Matrix{Float64}, source))

    return s
end

"""
    Ground

Mutable struct to hold ground parameters for emp2d.
"""
mutable struct Ground
    gsigma::Vector{Float64}
    gepsilon::Vector{Float64}

    Ground() = new()
end

########

"""
    writeinputs(s::Inputs; path="")

Write inputs.dat at `path` given `s`.
"""
function writeinputs(s::Inputs; path="")
    open(joinpath(path,"inputs.dat"), "w") do f
        for field in fieldnames(Inputs)
            write(f, getfield(s, field))
        end
    end
end

"""
    writesource(s::Source; path="")

Write source.dat at `path` given `s`.
"""
function writesource(s::Source; path="")
    open(joinpath(path,"source.dat"), "w") do f
        for field in fieldnames(Source)
            write(f, getfield(s, field))
        end
    end
end

"""
    writeground(s::Ground; path="")

Write ground.dat at `path` given `s`.
"""
function writeground(s::Ground; path="")
    open(joinpath(path,"ground.dat"), "w") do f
        for field in fieldnames(Ground)
            write(f, getfield(s, field))
        end
    end
end

"""
    writebfield(Bmag, in::Inputs; path="")

Write B0.dat at `path` given magnetic field magnitude `Bmag` and `in`.

`emp2d-slimfork.cpp` only supports vertically oriented magnetic fields.
"""
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

"""
    writene(ne::Array{Float64}; path="")

Write ne.dat at `path` given electron density array `ne`.

`ne` will typically be a two-dimensional array of size (altitude grid cells,
range grid cells). It can also be flattened vector of the same length.
"""
function writene(ne::Array{Float64}; path="")
    open(joinpath(path,"ne.dat"), "w") do f
        write(f, ne)
    end
end

"""
    writeni(ni::Array{Float64}; path="")

Write ni.dat at `path` given ion density array `ni`. `ni` will typically be the
same as `ne`.

`ni` will typically be a two-dimensional array of size (altitude grid cells,
range grid cells). It can also be flattened vector of the same length.
"""
function writeni(ni::Array{Float64}; path="")
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

"""
    writenu(nu::Array{Float64}; path="")

Write nu.dat at `path` given electron collision frequency array `nu`.

Ion collision frequency is calculated as `nu/100` in `emp2d-slimfork.cpp`.

`nu` will typically be a two-dimensional array of size (altitude grid cells,
range grid cells). It can also be flattened vector of the same length.
"""
function writenu(nu::Array{Float64}; path="")
    open(joinpath(path,"nu.dat"), "w") do f
        write(f, nu)
    end
end

########

"""
    emp2d(file::AbstractString, computejob::ComputeJob, inputs=false; submitjob=true)

Generate the input files and attempt to run the emp2d code for the scenario
described by `file` as `computejob`.

Optionally provide an inputs::Inputs() struct. Otherwise default values are used.
"""
function emp2d(file::AbstractString, computejob::ComputeJob; inputs=nothing, submitjob=true)
    isfile(file) || error("$file is not a valid file name")

    s = LWMS.parse(file)

    if computejob.runname != s.name
        @info "Updating computejob runname to $(s.name)"
        computejob.runname = s.name
    end

    if splitdir(computejob.rundir)[2] != computejob.runname
        # Check if computejob rundir path ends with a directory called runname
        rundir = joinpath(computejob.rundir, computejob.runname)
        @info "Running in $rundir"
    else
        rundir = computejob.rundir
    end

    if !isdir(rundir)
        # Create rundir if it doesn't exist
        @info "Creating $rundir"
        mkpath(rundir)
    end

    if isnothing(inputs)
        max_range = last(s.output_ranges)
        # round up to nearest thousand km and go 1000 km beyond that
        max_range = round(max_range+1000e3, digits=-3, RoundUp)

        inputs = Inputs(LWMS.EARTH_RADIUS, 110e3, 50e3, 500, 250, max_range, 500, [s.frequency])
    end

    shfile = build(s, computejob, inputs)

    if submitjob
        runjob(computejob, shfile)
    end

    return nothing
end

"""
    buildandrun(s::LWMS.BasicInput, computejob::ComputeJob, inputs::Inputs)

This is essentially a "private" function that sets default parameters for emp2d
and generates the necessary input files.
"""
function build(s::LWMS.BasicInput, computejob::ComputeJob, inputs::Inputs)

    all(s.b_dip .≈ 90) || @warn "Segment magnetic field is not vertical"
    length(unique(s.b_mag)) == 1 || @warn "Magnetic field is not homogeneous"

    ne_threshold = 3e9

    # Useful values
    max_range = inputs.range
    drange = inputs.drange
    Nrange = round(Int, max_range/drange) + 1
    rangevec = range(0, max_range, length=Nrange)

    r, dr = generate_rdr(inputs)  # somewhat confusing: r ≂̸ range, it is height
    altitudes = r .- LWMS.EARTH_RADIUS

    # Fill in values for each waveguide segment
    ne = Matrix{Float64}(undef, length(r), Nrange)
    nu = similar(ne)
    gsigma = Vector{Float64}(undef, Nrange)
    gepsilon = similar(gsigma)

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

    run(`cp $exefile $rundir`)

    return shfile
end

"""
    create_emp_source(def::Defaults, in::Inputs)

Create a standard emp source as a current waveform with linear rise and exponential
decay, which is then filtered with a lowpass filter according to the parameters in
`def`.
"""
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

"""
    generate_rdr(in::Inputs)

Create the `r` and `dr` vectors for the radial grids which extend along the
altitude direction beginning at `in.Re` and extending up until `in.maxalt`.
"""
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
