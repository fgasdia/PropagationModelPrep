nums = 0:15

for n in nums
    rootdir = "C:\\LWPCv21_$n"
    @info rootdir

    # Set lwpcDir before building
    lines = readlines(joinpath(rootdir, "setLWPC.cmd"); keep=true)
    for i in eachindex(lines)
        if startswith(lines[i], "set lwpcDir")
            lines[i] = "set lwpcDir=C:\\LWPCv21_$n\\"
        end
    end
    open(joinpath(rootdir, "setLWPC.cmd"), "w") do writer
        for line in lines
            write(writer, line)
        end
    end

    # Change hard-coded path
    lines = readlines(joinpath(rootdir, "library", "lwpc_dat_loc.for"); keep=true)
    for i in eachindex(lines)
        if occursin("C:\\LWPCv21\\a.loc", lines[i])
            lines[i] = replace(lines[i], "C:\\LWPCv21\\a.loc" => joinpath(rootdir, "a.loc"))
        end
    end
    open(joinpath(rootdir, "library", "lwpc_dat_loc.for"), "w") do writer
        for line in lines
            write(writer, line)
        end
    end
    open(joinpath(rootdir, "a.loc"), "w") do writer
        write(writer, joinpath(rootdir, "data", ""))
    end

    # Build LWPC
    run(`cmd /c cd $rootdir '&&' buildLWPC.cmd`)

    # Update executable name
    oldexefile = joinpath(rootdir, "lwpm.exe")
    newexefile = joinpath(rootdir, "lwpm$n.exe")
    if isfile(oldexefile)
        mv(oldexefile, newexefile)
    end
end
