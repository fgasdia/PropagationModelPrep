nums = 0:15

for n in nums
    rootdir = "C:\\LWPCv21_$n"

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

    # Build LWPC
    run(`cmd /c cd $rootdir '&&' buildLWPC.cmd`)

    # Update executable name
    oldexefile = joinpath(rootdir, "lwpm.exe")
    newexefile = joinpath(rootdir, "lwpm$n.exe")
    if isfile(oldexefile)
        mv(oldexefile, newexefile)
    end
end

