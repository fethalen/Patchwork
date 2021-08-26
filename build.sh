#!/usr/bin/env bash

compile_file="${RECIPE_DIR}/compile.jl"
build_dir="${SRC_DIR}/build"
patchwork_dir="patchwork-${PKG_VERSION}"

mkdir -p "${build_dir}/${patchwork_dir}"
cd "${build_dir}" || return

#julia -e 'using Pkg; Pkg.add(["PackageCompiler"])'
julia "${compile_file}" "${RECIPE_DIR}" "${build_dir}/${patchwork_dir}"

mkdir -p "${PREFIX}/bin"
cp -r ${build_dir}/* "${PREFIX}/bin"
ln -s "${PREFIX}/bin/${patchwork_dir}/bin/patchwork" "${PREFIX}/bin/patchwork"
#ln -s "${PREFIX}/bin/patchwork/bin/patchwork" "${PREFIX}/bin/patchwork"
#chmod u+x ${PREFIX}/bin/patchwork