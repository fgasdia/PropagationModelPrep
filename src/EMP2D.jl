module EMP2D

using Dates
using DSP, JSON3
using Interpolations

using ..PropagationModelPrep
using ..PropagationModelPrep: rounduprange, unwrap!, LWMS
using ..LWMS

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
    doioniz::Int32
    doelve::Int32
    dodetach::Int32
    dotransmitter::Int32
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
    camdist::Float64
    camalt::Float64
    elvesteps::Int32
    numfiles::Int32
    planet::Int32
    decfactor::Int32
    nprobes::Int32
    prober::Vector{Int32}
    probet::Vector{Int32}
    dogwave::Int32
    gwavemag::Float64
    gwavemaxalt::Float64
    gwavekh::Float64
    nonlinearstart::Int32
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
    - `savefields=[0, 0, 0, 0, 0, 0]` turns on output for E, J, H, K, and D when
        an entry is greater than 0. The last field (6) is not supported.

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
    # Re
    setfield!(s, :dopml_top, Int32(1))
    setfield!(s, :dopml_wall, Int32(1))
    setfield!(s, :doionosphere, Int32(1))
    setfield!(s, :doioniz, Int32(0))
    setfield!(s, :doelve, Int32(0))
    setfield!(s, :dodetach, Int32(0))
    setfield!(s, :dotransmitter, Int32(0))
    # savefields
    setfield!(s, :groundmethod, Int32(1))  # SIBC
    # maxalt
    # stepalt
    setfield!(s, :dr0, Float64(0))  # not actually used b/c nground = 0
    # dr1
    # dr2
    setfield!(s, :nground, Int32(0))
    # range
    # drange
    # dt
    # tsteps
    setfield!(s, :sig, Float64(0))
    setfield!(s, :sigm, Float64(0))
    setfield!(s, :camdist, Float64(0))
    setfield!(s, :camalt, Float64(0))
    setfield!(s, :elvesteps, Int32(0))
    setfield!(s, :numfiles, Int32(50)) # used to determine log file steps
    setfield!(s, :planet, Int32(0)) # Earth
    # decfactor
    setfield!(s, :nprobes, Int32(1))
    setfield!(s, :prober, Int32[0])
    setfield!(s, :probet, Int32[0])
    setfield!(s, :dogwave, Int32(0))
    setfield!(s, :gwavemag, Float64(0))
    setfield!(s, :gwavemaxalt, Float64(0))
    setfield!(s, :gwavekh, Float64(0))
    setfield!(s, :nonlinearstart, Int32(0))
    setfield!(s, :doDFT, Int32(1))
    # numDFTfreqs
    # DFTfreqs
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
    halfchlength::Int32
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
    nt_source, nalt_source = size(source)

    setfield!(s, :nalt_source, Int32(nalt_source))
    setfield!(s, :nt_source, Int32(nt_source))
    setfield!(s, :halfchlength, Int32(0))  # not used by emp2d
    setfield!(s, :source, source)

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

"""
    Phasor{T}

Store `amp` and `phase` of a complex value.
"""
struct Phasor{T}  # Float64 or Vector{Float64}
    amp::T
    phase::T
end

"""
    DFTFields{T}

Mutable struct of DFT electromagnetic fields of type `Phasor{T}`, as well as
distance vector `dist` and frequencies `DFTfreqs`.

Declare a new `DFTFields` of type `T` as `DFTFields(T)`.
"""
mutable struct DFTFields{T}
    dist::Vector{Float64}
    DFTfreqs::Vector{Float64}
    Er::Phasor{T}
    Et::Phasor{T}
    Ep::Phasor{T}
    Hr::Phasor{T}
    Ht::Phasor{T}
    Hp::Phasor{T}

    DFTFields(T::Type) = new{T}()
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
    readinputs(path)

Read the `inputs.dat` file located in directory `path`.
"""
function readinputs(path)
    inputs = Inputs()

    open(joinpath(path, "inputs.dat")) do f
        for field in fieldnames(Inputs)
            if field == :savefields
                tmp = fieldtype(Inputs, field)(undef, 6)
                for i = 1:6
                    v = read(f, eltype(fieldtype(Inputs, field)))
                    tmp[i] = v
                end
                setfield!(inputs, field, tmp)
            elseif field == :prober
                nprobes = getfield(inputs, :nprobes)
                tmp = fieldtype(Inputs, field)(undef, nprobes)
                for i = 1:nprobes
                    v = read(f, eltype(fieldtype(Inputs, field)))
                    tmp[i] = v
                end
                setfield!(inputs, field, tmp)
            elseif field == :probet
                nprobes = getfield(inputs, :nprobes)
                tmp = fieldtype(Inputs, field)(undef, nprobes)
                for i = 1:nprobes
                    v = read(f, eltype(fieldtype(Inputs, field)))
                    tmp[i] = v
                end
                setfield!(inputs, field, tmp)
            elseif field == :DFTfreqs
                numDFTfreqs = getfield(inputs, :numDFTfreqs)
                tmp = fieldtype(Inputs, field)(undef, numDFTfreqs)
                for i = 1:numDFTfreqs
                    v = read(f, eltype(fieldtype(Inputs, field)))
                    tmp[i] = v
                end
                setfield!(inputs, field, tmp)
            else
                v = read(f, fieldtype(Inputs, field))
                setfield!(inputs, field, v)
            end
        end
    end

    return inputs
end

"""
    writesource(s::Source; path="")

Write source.dat at `path` given `s`.
"""
function writesource(s::Source; path="")
    open(joinpath(path,"source.dat"), "w") do f
        for field in fieldnames(Source)
            if field == :source
                write(f, permutedims(getfield(s, :source)))  # for c++
            else
                write(f, getfield(s, field))
            end
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
        write(f, permutedims(ne))  # for c++
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
        write(f, permutedims(ni))  # for c++
    end
end

"""
    writenu(nu::Array{Float64}; path="")

Write nu.dat at `path` given electron collision frequency array `nu`.

Ion collision frequency is calculated as `nu/100` in `emp2d-slimfork.cpp`.

`nu` will typically be a two-dimensional array of size (altitude grid cells,
range grid cells). It can also be flattened vector of the same length.
"""
function writenu(nu::Array{Float64}; path="")
    open(joinpath(path,"nu.dat"), "w") do f
        write(f, permutedims(nu))  # for c++
    end
end

########

"""
    run(file, computejob::ComputeJob, inputs=false; submitjob=true)

Generate the input files and attempt to run the emp2d code for the scenario
described by `file` as `computejob`.

Optionally provide an inputs::Inputs() struct. Otherwise default values are used.
"""
function run(file, computejob::ComputeJob; inputs=nothing, submitjob=true)
    isfile(file) || error("$file is not a valid file name")

    s = LWMS.parse(file)

    if isnothing(inputs)
        max_range = maximum(s.output_ranges)
        max_range = rounduprange(max_range)

        inputs = Inputs(LWMS.EARTH_RADIUS, 110e3, 50e3, 500, 250, max_range, 500, [s.frequency])
    end

    shfiles = build(s, computejob, inputs)

    if submitjob
        submitandrun(shfiles)
    end

    return nothing
end

function submitandrun(shfile, computejob::ComputeJob)
    rundir, shfile = splitdir(shfile)
    exefile = computejob.exefile
    exename = splitdir(exefile)[2]
    newexefile = joinpath(rundir, exename)
    isfile(newexefile) || cp(exefile, newexefile)
    runjob(computejob, shfile)
end

function submitandrun(shfiles::AbstractVector, computejob::ComputeJob)
    exefile = computejob.exefile
    exename = splitdir(exefile)[2]

    for shfile in shfiles
        rundir, shfile = splitdir(shfile)
        runname, shext = splitext(shfile)

        newexefile = joinpath(rundir, exename)
        isfile(newexefile) || cp(exefile, newexefile)
        computejob.runname = runname
        runjob(computejob, shfile)
    end
end

function prepbuild(s, computejob::ComputeJob, inputs::Inputs)
    if computejob.runname != s.name
        @info "Updating computejob runname to $(s.name)"
        computejob.runname = s.name
    end

    # Get abspath in case a relpath is provided
    origrundir = splitpath(abspath(computejob.rundir))[end]
    if origrundir != computejob.runname
        # Check if computejob rundir path ends with a directory called runname
        rundir = joinpath(computejob.rundir, computejob.runname)*"/"
        @info "Running in $rundir"
    else
        rundir = computejob.rundir*"/"
    end

    if !isdir(rundir)
        # Create rundir if it doesn't exist
        @info "Creating $rundir"
        mkpath(rundir)
    end

    return nothing
end

function build(s::BatchInput, computejob::ComputeJob, inputs::Inputs)
    shfiles = Vector{String}(undef, length(s.inputs))
    for i in eachindex(s.inputs)
        shfiles[i] = build(s.inputs[i], computejob, inputs)
    end
    return shfiles
end

"""
    build(s::BasicInput, computejob::ComputeJob, inputs::Inputs)

This is essentially a "private" function that sets default parameters for emp2d
and generates the necessary input files.
"""
function build(s::BasicInput, computejob::ComputeJob, inputs::Inputs)
    prepbuild(s, computejob, inputs)

    all(s.b_dips .≈ π/2) || @warn "Segment magnetic field is not vertical"
    length(unique(s.b_mags)) == 1 || @warn "Magnetic field is not homogeneous"

    # Useful values
    r, dr, th = generategrid(inputs)  # r is along radial direction (height)
    altitudes = r .- inputs.Re

    rr = length(r)
    hh = length(th)

    max_range = inputs.range
    drange = inputs.drange
    rangevec = range(0, max_range, length=hh)

    # Fill in values for each waveguide segment
    ne = Matrix{Float64}(undef, rr, hh)
    nu = similar(ne)
    gsigma = Vector{Float64}(undef, hh)
    gepsilon = similar(gsigma)

    for i in eachindex(s.segment_ranges)
        # Find closest match of segment range to rangevec
        segment_begin_idx = argmin(abs.(rangevec .- s.segment_ranges[i]))
        if i == lastindex(s.segment_ranges)
            segment_end_idx = hh
        else
            segment_end_range = s.segment_ranges[i+1]
            segment_end_idx = argmin(abs.(rangevec .- segment_end_range))
        end

        # Electron density profile
        neprofile = LWMS.waitprofile.(altitudes, s.hprimes[i], s.betas[i], cutoff_low=50e3, threshold=3e9)
        ne[:,segment_begin_idx:segment_end_idx] .= neprofile

        # Electron collision profile
        nuprofile = LWMS.electroncollisionfrequency.(altitudes)
        nu[:,segment_begin_idx:segment_end_idx] .= nuprofile

        # Ground profile
        gsigma[segment_begin_idx:segment_end_idx] .= s.ground_sigmas[i]
        gepsilon[segment_begin_idx:segment_end_idx] .= s.ground_epsrs[i]
    end

    source = Source(inputs)
    ground = Ground()
    ground.gsigma = gsigma
    ground.gepsilon = gepsilon

    rundir = computejob.rundir

    writeinputs(inputs, path=rundir)
    writesource(source, path=rundir)
    writeground(ground, path=rundir)
    writebfield(s.b_mags[1], inputs, path=rundir)
    writene(ne, path=rundir)
    writeni(ne, path=rundir)
    writenu(nu, path=rundir)

    shfile = writeshfile(computejob)

    return shfile
end

"""
    build(s::TableInput, computejob::ComputeJob, inputs::Inputs)

This is essentially a "private" function that sets default parameters for emp2d
and generates the necessary input files.
"""
function build(s::TableInput, computejob::ComputeJob, inputs::Inputs)
    prepbuild(s, computejob, inputs)

    all(s.b_dips .≈ π/2) || @warn "Segment magnetic field is not vertical"
    length(unique(s.b_mags)) == 1 || @warn "Magnetic field is not homogeneous"

    # Useful values
    r, dr, th = generategrid(inputs)  # r is along radial direction (height)
    altitudes = r .- inputs.Re

    rr = length(r)
    hh = length(th)

    max_range = inputs.range
    drange = inputs.drange
    rangevec = range(0, max_range, length=hh)

    # Fill in values for each waveguide segment
    ne = Matrix{Float64}(undef, rr, hh)
    nu = similar(ne)
    gsigma = Vector{Float64}(undef, hh)
    gepsilon = similar(gsigma)

    for i in eachindex(s.segment_ranges)
        # Find closest match of segment range to rangevec
        segment_begin_idx = argmin(abs.(rangevec .- s.segment_ranges[i]))
        if i == lastindex(s.segment_ranges)
            segment_end_idx = hh
        else
            segment_end_range = s.segment_ranges[i+1]
            segment_end_idx = argmin(abs.(rangevec .- segment_end_range))
        end

        # Electron density profile
        density_itp = LinearInterpolation(s.altitude, s.density[i], extrapolation_bc=Line())
        ne[:,segment_begin_idx:segment_end_idx] .= density_itp(altitudes)

        # Electron collision profile
        collision_itp = LinearInterpolation(s.altitude, s.collision_frequency[i], extrapolation_bc=Line())
        nu[:,segment_begin_idx:segment_end_idx] .= collision_itp(altitudes)

        # Ground profile
        gsigma[segment_begin_idx:segment_end_idx] .= s.ground_sigmas[i]
        gepsilon[segment_begin_idx:segment_end_idx] .= s.ground_epsrs[i]
    end

    source = Source(inputs)
    ground = Ground()
    ground.gsigma = gsigma
    ground.gepsilon = gepsilon

    rundir = computejob.rundir

    writeinputs(inputs, path=rundir)
    writesource(source, path=rundir)
    writeground(ground, path=rundir)
    writebfield(s.b_mags[1], inputs, path=rundir)
    writene(ne, path=rundir)
    writeni(ne, path=rundir)
    writenu(nu, path=rundir)

    shfile = writeshfile(computejob)

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

    source = repeat(Iin, 1, satop)

    return source
end

"""
    generategrid(in)

Calculate `r`, `dr`, and `th` vectors for the radial and theta grids which
effectively extend along the altitude direction from `in.Re` to `in.maxalt` and
the propagation direction, respectively.

!!! note
    `rr` and `hh` can be obtained from `length(r)` and `length(th)`, respectively.
"""
function generategrid(in::Inputs)
    # Set up radial array
    rr = trunc(Int32, in.stepalt/in.dr1 + (in.maxalt - in.stepalt)/in.dr2 + 1 + in.nground)

    r = Vector{Float64}(undef, rr)
    dr = Vector{Float64}(undef, rr-1)

    r[1] = in.Re - in.nground*in.dr0
    if in.nground > 0
        for i = 2:nground+1
            r[i] = r[i-1] + in.dr0
            dr[i-1] = r[i] - r[i-1]
        end
    end
    for i = in.nground+2:rr
        if r[i-1] < (in.Re + in.stepalt)
            r[i] = r[i-1] + in.dr1
        else
            r[i] = r[i-1] + in.dr2
        end
        dr[i-1] = r[i] - r[i-1]
    end

    # Set up theta array
    thmax = in.range/in.Re
    thmax > π && (thmax = π)
    dth = in.drange / in.Re
    hh = round(Int32, thmax/dth) + 1

    th = Vector{Float64}(undef, hh)
    th[1] = 0
    for i = 2:hh
        th[i] = th[i-1] + dth
    end

    return r, dr, th
end

"""
    readdfts(path)

Read `inputs.dat` and `dfts.dat` in directory `path` and return a `DFTFields`.
"""
function readdfts(path)
    inputs = readinputs(path)
    r, dr, th = generategrid(inputs)
    rr = length(r)
    hh = length(th)

    f = open(joinpath(path, "dft.dat"),"r")

    numDFTfreqs = read(f, Int32)
    DFTfreqs = Vector{Float64}(undef, numDFTfreqs)
    read!(f, DFTfreqs)

    Er = Matrix{Float64}(undef, 2*numDFTfreqs, hh)
    Et = similar(Er)
    Ep = similar(Er)
    Hr = similar(Er)
    Ht = similar(Er)
    Hp = similar(Er)

    read!(f, Er)
    read!(f, Et)
    read!(f, Ep)
    read!(f, Hr)
    read!(f, Ht)
    read!(f, Hp)

    close(f)

    out = DFTFields(Vector{Float64})
    out.dist = th*inputs.Re
    out.DFTfreqs = DFTfreqs

    for m = 1:numDFTfreqs
        for (field, vals) in Dict(:Er=>Er, :Et=>Et, :Ep=>Ep, :Hr=>Hr, :Ht=>Ht, :Hp=>Hp)
            tmp = complex.(vals[2*m-1,:], vals[2*m,:])
            setfield!(out, field, Phasor(abs.(tmp), unwrap!(angle.(tmp))))
        end
    end

    return out
end

"""
    process(path)

Read `inputs.dat` and `dfts.dat` in directory `path` and return a `BasicOutput`
with amplitude in dB μV/m and phase in degrees corrected for dispersion.
"""
function process(path)
    dfts = readdfts(path)
    inputs = readinputs(path)

    # TODO: support multiple frequencies
    length(dfts.DFTfreqs) > 1 && @info "Only using DFTfreq $(out.DFTfreqs[1]) Hz"
    freq = first(dfts.DFTfreqs)
    freq_kHz = freq/1e3

    dist = dfts.dist
    dist_km = dist./1e3

    drange_km = inputs.drange/1e3

    amp = 20*log10.(dfts.Er.amp/1e-6)  # dB μV/m
    amp .-= maximum(amp)

    phase = rad2deg.(dfts.Er.phase) + (360*freq/LWMS.C0*dist)

    # 2D FDTD numerical dispersion correction
    p = [0.00136588, -0.026108, 0.154035]
    k = p[1]*freq_kHz^2 + p[2]*freq_kHz + p[3]  # deg/km²
    phase += k*dist_km*drange_km^2

    # Create output
    fullpath = abspath(path)
    pathname = splitpath(fullpath)[end]  # proxy for filename

    jsonfile = joinpath(fullpath, pathname*".json")

    if !isfile(jsonfile)
        @info "$jsonfile not found. Assuming run name is $pathname."
        name = pathname
        description = "emp2d results"
        datetime = Dates.now()
    else
        s = LWMS.parse(jsonfile)
        name = s.name
        description = s.description
        datetime = s.datetime
    end

    output = BasicOutput()
    output.name = name
    output.description = description
    output.datetime = datetime

    # Strip last index of each field because they're 0 (and then inf)
    output.output_ranges = round.(dist[1:end-1], digits=3)  # fix floating point
    output.amplitude = amp[1:end-1]
    output.phase = phase[1:end-1]

    # Save output
    json_str = JSON3.write(output)

    open(joinpath(path,name*"_emp2d.json"), "w") do f
        write(f, json_str)
    end

    return output
end

function writeshfile(s::LocalParallel)
    runname = s.runname
    rundir = s.rundir
    numnodes = s.numnodes
    exefile = s.exefile
    exefile = basename(exefile)

    shfile = joinpath(rundir, runname*".sh")
    endline = "\n"

    open(shfile, "w") do f
        write(f, "#!/bin/sh", endline)
        write(f, endline)
        write(f, "rm -f $rundir/output_K.dat", endline)
        write(f, "rm -f $rundir/output_E.dat", endline)
        write(f, "rm -f $rundir/output_D.dat", endline)
        write(f, "rm -f $rundir/output_H.dat", endline)
        write(f, "rm -f $rundir/output_J.dat", endline)
        write(f, "rm -f $rundir/output_T.dat", endline)
        write(f, "rm -f $rundir/output_O.dat", endline)
        write(f, "rm -f $rundir/output_S.dat", endline)
        write(f, "rm -f $rundir/Probe.dat", endline)
        write(f, "rm -f $rundir/elve.dat", endline)
        write(f, "rm -f $rundir/sferic.dat", endline)
        write(f, endline)
        write(f, "export OMP_NUM_THREADS=$numnodes", endline)
        write(f, endline)
        write(f, joinpath(rundir,exefile), endline)
    end

    return shfile
end

function writeshfile(s::Summit)
    runname = s.runname
    rundir = s.rundir
    numnodes = s.numnodes
    walltime = s.walltime
    exefile = s.exefile
    exefile = basename(exefile)

    shfile = joinpath(rundir, runname*".sh")
    endline = "\n"

    open(shfile, "w") do f
        write(f, "#!/bin/sh", endline)
        write(f, "#SBATCH --job-name=$runname", endline)
        write(f, "#SBATCH --partition=shas", endline)
        write(f, "#SBATCH --nodes=1", endline)
        write(f, "#SBATCH --ntasks=$numnodes", endline)
        write(f, "#SBATCH --time=$walltime", endline)
        write(f, endline)
        write(f, "rm -f $rundir/output_K.dat", endline)
        write(f, "rm -f $rundir/output_E.dat", endline)
        write(f, "rm -f $rundir/output_D.dat", endline)
        write(f, "rm -f $rundir/output_H.dat", endline)
        write(f, "rm -f $rundir/output_J.dat", endline)
        write(f, "rm -f $rundir/output_T.dat", endline)
        write(f, "rm -f $rundir/output_O.dat", endline)
        write(f, "rm -f $rundir/output_S.dat", endline)
        write(f, "rm -f $rundir/Probe.dat", endline)
        write(f, "rm -f $rundir/elve.dat", endline)
        write(f, "rm -f $rundir/sferic.dat", endline)
        write(f, endline)
        write(f, "module purge", endline)
        write(f, "module load intel", endline)
        write(f, endline)
        write(f, "export OMP_NUM_THREADS=$numnodes", endline)
        write(f, endline)
        write(f, joinpath(rundir,exefile), endline)
    end

    return shfile
end

function runjob(s::LocalParallel, shfile)
    jobname = read(`$shfile`, String)

    println(jobname)

    return nothing
end

function runjob(s::Summit, shfile)
    # jobname = read(`sbatch $shfile`, String)
    # jobid = strip(jobname)

    # println("Job $jobid submitted!\n")

    # TEMP: Running Julia from inside Singularity isn't allowing this command
    println("Please run the command:")
    println("sbatch $shfile")

    return nothing
end

end  # module
