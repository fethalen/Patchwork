#!/usr/bin/env julia

const BASEDIR    = "/" * joinpath(split(@__FILE__, "/")[1:end-1]...)
const PROJECTDIR = abspath(joinpath(BASEDIR, ".."))
const COMPILEDIR = abspath(joinpath(BASEDIR, "compiled"))

# Remove compiled directory if it already exists
if !isdir(COMPILEDIR)
    mkdir(COMPILEDIR)
    @info "Created directory $COMPILEDIR."
else
    rm(COMPILEDIR, recursive=true, force=true)
    @info "Removed pre-existing compiled directory $COMPILEDIR."
end

# Resolve dependencies and activate project
using Pkg
try
    Pkg.activate(PROJECTDIR)
    Pkg.instantiate()
catch e
    @warn "Could not read project's dependencies! ($e)"
end

# Build executable
Pkg.activate(PROJECTDIR)
using PackageCompiler
create_app(
    PROJECTDIR,
    COMPILEDIR,
    force=true,
    incremental=true,
    include_transitive_dependencies=false,
    filter_stdlibs = false,
    cpu_target="native"
)
