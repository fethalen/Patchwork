#!/usr/bin/env bash
#
# Compiles the Patchwork software. Requires that Julia is installed.

set -o errtrace
set -o errexit
set -o nounset
set -o pipefail

# Show line numbers and function names when debugging
export PS4='+ ${LINENO}:${FUNCNAME[0]:-}() '

readonly VERSION="0.1.0"
PATCHWORK_SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
readonly PATCHWORK_BASEDIR="${PATCHWORK_SCRIPTDIR}/.."
readonly PATCHWORK_SOURCEDIR="${PATCHWORK_BASEDIR}/src"
readonly PATCHWORK_TESTDIR="${PATCHWORK_BASEDIR}/test"
readonly CONTIGS="${PATCHWORK_TESTDIR}/c_australis_k55_node_82332.fa"
readonly REFERENCE="${PATCHWORK_TESTDIR}/a_succinea_m_30.fa"

usage() { printf "%s" "\
usage:
  compile_patchwork.sh [--help] [--version] --outdir PATH

description:
  Compile Patchwork and store the program in the provided output directory.

options:
  -h, --help           display this help message and exit
  -v, --version        display the version number and exit
  -o, --outdir PATH    path to the directory where the program will be stored
                       (required)
"
  exit 1
}

# Display the provided error message before exiting with the provided status
# (default: 1).
error() {
  local message="$1"
  local status="${2-1}" # default exit status: 1
  echo "$message"
  exit "$status"
}

# Display the version number and exit.
version_info() {
  echo "$VERSION"
  exit 1
}

# Parse user-provided arguments
[[ "$#" -lt 1 ]] && usage
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
   -h | --help)
      usage
      ;;
    -v | --version)
      version_info
      ;;
    -o | --outdir)
      outdir="$2"
      shift # past argument
      shift # past value
      ;;
    --) # end of all options
      break
      ;;
    -*) # unknown option
      error "unknown option: $key"
      ;;
    *) # end of options
      [[ -n "$basedir" ]] && error "unrecognized positional argument: $key"
      basedir="$key"
      shift # past argument
      ;;
  esac
done

compile_patchwork() {
  julia --trace-compile="${PATCHWORK_SOURCEDIR}/precompiled.jl"\
    "${PATCHWORK_SOURCEDIR}/Patchwork.jl"\
    --contigs "$CONTIGS"\
    --reference "$REFERENCE"\
    --output-dir "$test_run_out"
}

test_run() {
  julia "${PATCHWORK_SOURCEDIR}/compile.jl"\
    "$PATCHWORK_BASEDIR"\
    "${PATCHWORK_SOURCEDIR}/precompiled.jl"\
    "$outdir"
}

main() {
  test_run_out=$(mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir')
  julia "${PATCHWORK_SCRIPTDIR}/install_packages.jl"
  compile_patchwork
  test_run
  rm -rf test_run_out
}

main "$@"
