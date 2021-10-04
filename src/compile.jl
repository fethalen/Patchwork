# this script can be run on its own, e.g. from inside the Patchwork directory: 
# julia --trace-compile="src/precompiled.jl" src/Patchwork.jl --contigs CONTIGSFILE --reference REFERENCEFILE
# julia src/compile.jl . src/precompiled.jl ../patchwork
# or run it as part of conda build, from inside the Patchwork recipe directory
# that contains the meta.yaml, build.sh, Project.toml and Manifest.toml files:
# (no precompilation necessary) 
# conda build .

println("START")
projectdirectory = ARGS[1]
precompiled = ARGS[2]
outdirectory = ARGS[3]

import Pkg 
println("ACTIVATE")
Pkg.activate(projectdirectory)
Pkg.add(["ArgParse", "PackageCompiler"])
using PackageCompiler

println("CREATE")
create_app(projectdirectory, outdirectory; precompile_statements_file = precompiled, 
           force = true, app_name = "patchwork")