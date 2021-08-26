using Pkg
using PackageCompiler

projectdirectory = ARGS[1]                      # the conda build source directory SRC_DIR
precompiled = projectdirectory * "/precompiled.jl"
outdirectory = ARGS[2]                          # the conda build build/patchwork directory

try
    mkpath(outdirectory)
catch e
    println("WARNING for" * outdirectory * ": " * e)
    println("Using the present working directory " * pwd() * " as outdir.")
    global outdirectory = pwd()
end

Pkg.activate(projectdirectory)

create_app(projectdirectory, outdirectory; precompile_statements_file = precompiled, 
           force = true, app_name = "patchwork")