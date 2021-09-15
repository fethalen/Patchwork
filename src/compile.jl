# this script can be run on its own, e.g. from inside the Patchwork directory: 
# julia --trace-compile="src/precompiled.jl" src/Patchwork.jl --contigs CONTIGSFILE --reference REFERENCEFILE
# julia src/compile.jl . src/precompiled.jl ../patchwork
# or run it as part of conda build, from inside the Patchwork recipe directory
# that contains the meta.yaml, build.sh, Project.toml and Manifest.toml files:
# (no precompilation necessary) 
# conda build .

import Pkg
# ERROR occurred in bioconda build test; proposed solution was:
#Pkg.add("PackageCompiler") # HERE
##############################################################
using PackageCompiler

println("START") # HERE
projectdirectory = ARGS[1]
precompiled = ARGS[2]
outdirectory = ARGS[3]

println("ACTIVATE") # HERE
Pkg.activate(projectdirectory)

println("CREATE") # HERE
create_app(projectdirectory, outdirectory; precompile_statements_file = precompiled, 
           force = true, app_name = "patchwork")