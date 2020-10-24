module LocalParallelModule

using ..LWMS
# export batchrunjobs, make_jobs

function batchrunjobs(jobs, results)
    while true
        s = take!(jobs)
        pid = myid()

        exefile, exeext = splitext(computejob.exefile)
        exepath, exefilename = splitdir(exefile)
        newexepath = exepath*"_"*string(pid)
        newexefile = joinpath(newexepath, exefilename*string(pid)*exeext)

        cj = Local(s.name, newexepath, newexefile)
        build(s, cj)
        runjob(cj)
        put!(results, s.name)
    end
end

end
