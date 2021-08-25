println("##############################################################################")
println("START COMPILE ################################################################")
println("##############################################################################")
using Pkg
#Pkg.add(["HTTP", "LibSSH2_jll", "IJulia"])

#projectenv = ARGS[1] #"/home/clara/Data/GAU/Work/Projects/Patchwork"
projectdirectory = ARGS[1] #"/home/clara/Data/GAU/Work/Projects/Patchwork"
precompiled = projectdirectory * "/precompiled.jl"
outdirectory = ARGS[2] #"compiled"

println("##############################################################################")
println("MAKE OUTDIR ##################################################################")
println("##############################################################################")
#this block is new: 
try
    mkpath(outdirectory)
catch e
    println("WARNING for" * outdirectory * ": " * e)
    println("Using the present working directory " * pwd() * " as outdir.")
    global outdirectory = pwd()
end

println("##############################################################################")
println("ACTIVATE ENV #################################################################")
println("##############################################################################")
Pkg.activate(projectdirectory)
println("##############################################################################")
println("ACTIVATED ####################################################################")
println("##############################################################################")


Pkg.add("PackageCompiler")
#Pkg.add(["HTTP", "LibSSH2_jll", "IJulia"])
using PackageCompiler

#Pkg.activate(projectdirectory)
#println("environment activated")
#Pkg.instantiate()
#println("instantiated")

println("##############################################################################")
println("CREATE APP ###################################################################")
println("##############################################################################")
create_app(projectdirectory, outdirectory; precompile_statements_file = precompiled, 
           force = true)
println("##############################################################################")
println("DONE ######## ################################################################")
println("##############################################################################")


