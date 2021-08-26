#!/usr/bin/env bash

compile_file="${RECIPE_DIR}/compile.jl"
build_dir="${SRC_DIR}/build"

mkdir -p "${build_dir}/patchwork"
cd "${build_dir}" || return

#julia -e 'using Pkg; Pkg.add(["PackageCompiler"])'
julia "${compile_file}" "${RECIPE_DIR}" "${build_dir}/patchwork"

mkdir -p "${PREFIX}/bin"
cp -r ${build_dir}/* "${PREFIX}/bin"
#install -d ${PREFIX}/bin
#install gff ${PREFIX}/bin/