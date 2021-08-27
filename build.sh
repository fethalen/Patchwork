#!/usr/bin/env bash

# script for compiling:
compile_file="${SRC_DIR}/src/compile.jl"
# file with precompilation statements: 
precompiled_file="${SRC_DIR}/src/precompiled.jl"
# test files for precompilation run, if necessary:
contigs="${RECIPE_DIR}/test/07673_lcal.fa"
reference="${RECIPE_DIR}/test/07673_Alitta_succinea.fa"
# directory to contain the finished build:
build_dir="${SRC_DIR}/build"
# inside build_dir, contains directories artifacts, bin with executable patchwork, and lib:
patchwork_dir="patchwork-${PKG_VERSION}"

# precompile, if necessary: 
if [ ! -f "${precompiled_file}" ]
then 
    julia --trace-compile="${precompiled_file}" "${SRC_DIR}/src/Patchwork.jl" --contigs "${contigs}" --reference "${reference}"
fi

# this is not very nice, but you need the .toml files for compilation and conda only copies 
# the stuff inside the "src"-directory of the repo to its own working directory SRC_DIR
# so make sure you have the .toml files in the RECIPE_DIR together with meta.yaml and build.sh
cp "${RECIPE_DIR}/Project.toml" "${SRC_DIR}"
cp "${RECIPE_DIR}/Manifest.toml" "${SRC_DIR}"

mkdir -p "${build_dir}/${patchwork_dir}"
# compile Patchwork: 
julia "${compile_file}" "${SRC_DIR}" "${precompiled_file}" "${build_dir}/${patchwork_dir}"

# bundle together for conda packaging: 
mkdir -p "${PREFIX}/bin"
cp -r ${build_dir}/* "${PREFIX}/bin"
ln -s "${PREFIX}/bin/${patchwork_dir}/bin/patchwork" "${PREFIX}/bin/patchwork"