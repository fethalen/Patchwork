# this script can be run on its own, e.g. from inside the Subsample directory: 
# julia --trace-compile="src/precompiled.jl" src/Subsample.jl --contigs CONTIGSFILE --reference REFERENCEFILE
# julia src/compile.jl . src/subsample_precompiled.jl build


import Pkg 
Pkg.add(name="PackageCompiler", version="1.7.1")
Pkg.instantiate()
using PackageCompiler

println("Initializing...")
projectdirectory = ARGS[1]
precompiled = ARGS[2]
outdirectory = ARGS[3]

println("Starting compilation process...")
create_app(projectdirectory, outdirectory; precompile_statements_file = precompiled, force = true, app_name="subsample")
println("Compilation finished.")
println("Creating symlink to executable...")
try
	run(`cd $outdirectory`)
	run(`ln -s bin/subsample subsample`)
catch e
	println("Failed to create symlink to executable. See below for error message. ")
	println(e)
end
println("Done.")