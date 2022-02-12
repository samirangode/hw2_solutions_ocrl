

const NUMQUESTIONS = 3

# Create a module for each question
for i in 1:NUMQUESTIONS 
@eval module $(Symbol("Q" * string(i)))
include("autograder.jl")
getname() = string(split(string(@__MODULE__), ".")[end])
grade() = Autograder.gradequestion(getname()) 
checktestsets(solutiondir=joinpath(@__DIR__, "..")) = Autograder.checktestsets(getname(), solutiondir)
end
end

solutiondir = get(ARGS, 1, joinpath(@__DIR__, ".."))

# Grade all of the questions
modules = [@eval $(Symbol("Q" * string(i))) for i = 1:NUMQUESTIONS]
points,ts = modules[1].grade()
points,ts = modules[2].grade()
points,ts = modules[3].grade()


results = map(modules) do mod
    mod.checktestsets(solutiondir)
    mod.grade()[1]
end