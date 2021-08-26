println("##############################################################################")
println("START COMPILE ################################################################")
println("##############################################################################")
using Pkg
Pkg.add("PackageCompiler")
using PackageCompiler

#projectenv = ARGS[1] #"/home/clara/Data/GAU/Work/Projects/Patchwork"
projectdirectory = ARGS[1]                      # the conda build source directory SRC_DIR
precompiled = projectdirectory * "/precompiled.jl"
outdirectory = ARGS[2]                          # the conda build build/patchwork directory

println("##############################################################################")
println("MAKE OUTDIR ##################################################################")
println("##############################################################################")

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


