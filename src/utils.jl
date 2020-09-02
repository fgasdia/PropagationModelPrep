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
        write(f, "#!/bin/sh\n", endline)
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
