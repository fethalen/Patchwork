#!/usr/bin/env julia

const BUILDDIR = @__DIR__
const PROJECTDIR = dirname(BUILDDIR)
const COMPILEDIR = joinpath(BUILDDIR, "compiled")

using Pkg

# Instantiate the build tooling and application in their separate environments.
Pkg.activate(BUILDDIR)
Pkg.instantiate()
Pkg.activate(PROJECTDIR)
Pkg.instantiate()
Pkg.activate(BUILDDIR)

using PackageCompiler

# Remove a previous build only after both environments are ready.
if isdir(COMPILEDIR)
    rm(COMPILEDIR, recursive = true, force = true)
    @info "Removed pre-existing compiled directory $COMPILEDIR."
end

create_app(
    PROJECTDIR,
    COMPILEDIR,
    executables = ["patchwork" => "julia_main"],
    force = true,
    incremental = true,
    include_transitive_dependencies = false,
    filter_stdlibs = false,
    cpu_target = "native"
)
