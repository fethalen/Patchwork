#!/usr/bin/env bash
#
# Run Patchwork on the provided set of directories.

set -o errtrace
set -o errexit
set -o nounset
set -o pipefail

# Show line numbers and function names when debugging
export PS4='+ ${LINENO}:${FUNCNAME[0]:-}() '

readonly VERSION="0.1.0"
readonly COMBINED_DIR="markers_all_species"
basedir=''
reference=''
outdir='multi_patchwork_out'
files_only=0
combine_only=0
skip_combine=0
id_offset=3
kmer_offset=1

usage() { printf "%s" "\
usage:
  multi_patchwork.sh [--help] [--version] DIRECTORY

description:
  Run Patchwork on the provided set of files and combine the outputs into one.

recommended usage:
  - Display input files using '--files-only' to verify paths and species IDs
  - Once your happy with the input, do a normal run. Remember that
    '--reference' and 'DIRECTORY' are always required
  - _If_ you need to perform multiple runs, set '--skip-combine'
  - After performing _multiple runs_, set '--combine-only' to combine outputs
    separately

options:
  miscellaneous:
    -h, --help          display this help message and exit
    -v, --version       display the version number and exit

  input/output control:
    -r, --reference     path to a FASTA file containing protein sequences to
                        match with (REQUIRED)
    -o, --outdir        save all output to this directory
                        (default: ${outdir})

  offsets:
    -k, --kmer-offset     get K-mer size from the nth parent directory
                          (default: "$kmer_offset")
    -i, --id-offset       get ID from the nth parent directory (default: "$id_offset")

  flow control:
    -f, --files-only    display Patchwork input files and exit
    -c, --combine-only  do not run Patchwork but combine files found in output
    -s, --skip-combine  do not combine output files (do it later using
                        \`--combine-only\`)
"
  exit 1
}

# Display the version number and exit.
version_info() {
  echo "$VERSION"
  exit 1
}

# Display the provided error message before exiting with the provided status
# (default: 1).
error() {
  local message="$1"
  local status="${2-1}" # default exit status: 1
  echo "ega_wrapper: error: $message"
  exit "$status"
}

# Takes the path to a directory and an integer as an input. Go to the provided
# directory's parent directory X times.
up() {
  local dir="$1"
  local levels="$2"
  [[ ! -d "$dir" && ! -f "$dir" ]] && error "provided directory not found: $dir"
  while (( levels > 0 ))
  do
    dir=$( dirname "$dir" )
    (( levels = levels - 1 ))
  done
  echo "$dir"
}

# This function takes an array of assembly paths and returns the subset of
# those array with the smallest contigs.
get_smallest_contigs() {
  local assemblies=("${@}")
  local smallest_contigs=()
  local last_kmer_size=""
  local last_id=""
  local best_choice=""

  for assembly in "${assemblies[@]}"
  do
    kmer_size=$( basename "$( up "$assembly" "$kmer_offset" )" | tr -d 'K' )
    id=$( basename "$( up "$assembly" "$id_offset" )" )

    if [[ "$last_id" && "$last_kmer_size" ]]
    then
      if (( "$kmer_size" < "$last_kmer_size" ))
      then
        best_choice="$assembly"
      fi
      if [[ "$best_choice" && "$last_id" != "$id" ]]
      then
        smallest_contigs=("${smallest_contigs[@]}" "${best_choice}")
        best_choice="$assembly"
      fi
    else
      best_choice="$assembly"
    fi

    last_kmer_size="$kmer_size"
    last_id="$id"
  done
  if [[ "$best_choice" ]]
  then
    smallest_contigs=("${smallest_contigs[@]}" "${best_choice}")
  fi

  echo "${smallest_contigs[@]}"
}

# Takes a filename and a species name as an input. Inserts the species ID given
# by the 3rd parent directory of the file. If an at symbol ('@') already exist
# then don't do anything.
insert_species_prefix() {
  local filename="$1"
  local id=''
  id=$( basename "$( up "$filename" "$id_offset" )" )
  first_line=$(head -n 1 "${filename}" )
  if [[ ! "$first_line" =~ "@" ]]; then
      sed -i "s/>/>${id}@/g" "$filename"
  fi
}

get_sequence_ids() {
  local filename="$1"
  local seq_ids=()
  readarray -d '' seq_ids < <( grep '^>' "$filename" | sed 's/>[^>]*@//' )
  echo "${seq_ids[@]}"
}

combine_output() {
  local multi_patchwork_out="$1"
  marker_count=0
  for seq_id in $( get_sequence_ids "$reference" )
  do
    rm -f "${outdir}/${COMBINED_DIR}/${seq_id}.fa"
    readarray -d '' matches <\
      <(find "$multi_patchwork_out" -type f -name "${seq_id}.fa" -print0)
    for match in "${matches[@]}"
    do
      cat "$match" >> "${outdir}/${COMBINED_DIR}/${seq_id}.fa"
      (( marker_count++ ))
    done
  done
  echo "$marker_count"
}

run_patchwork() {
  local assemblies=(${@}) # doesn't split correctly when quoted
  local id=""
  for assembly in "${assemblies[@]}"
  do
    id=$( basename "$( up "$assembly" "$id_offset" )" )
    patchwork_output="${outdir}/${id}_patchwork_out"
    patchwork\
      --contigs "$assembly"\
      --reference "$reference"\
      --output-dir "$patchwork_output"\
      --overwrite
  done
}

# Returns all FASTQ files found in the provided `path`
fastq_files() {
  local path="$1"
  find "$path" -type f \( --iname \*.fastq.gz -o -iname \*.fq.gz \) |\
    sort
}

[[ $# -lt 1 ]] && usage

# Parse user-provided arguments
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
    -r | --reference)
      reference="$2"
      shift # past argument
      shift # past value
      ;;
    -f | --files-only)
      files_only=1
      shift # past argument
      ;;
    -c | --combine-only)
      combine_only=1
      shift # past argument
      ;;
    -s | --skip-combine)
      skip_combine=1
      shift # past argument
      ;;
    -k | --kmer-offset)
      kmer_offset="$2"
      shift # past argument
      shift # past value
      ;;
    -i | --id-offset)
      id_offset="$2"
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

main() {
  # Display an error message if no arguments were provided
  [[ -z "$basedir" ]] && error "missing mandatory argument DIRECTORY"
  [[ -z "$reference" ]] && error "missing mandatory argument --reference"
  [[ ! -d "$basedir" ]] && error "provided directory not found: $basedir"
  [[ ! -f "$reference" ]] && error "provided file not found: $reference"
  if (( ! files_only ))
  then
    mkdir -p "$outdir"
    mkdir -p "${outdir}/${COMBINED_DIR}"
  fi
  readarray -d '' assemblies <\
    <(find "$basedir" -type f -name "final_contigs.fasta" -not -wholename "*/work/*" -print0 )
  local smallest_contigs=( "$( get_smallest_contigs "${assemblies[@]}" )" )

  if (( files_only ))
  then
    echo "These files would be used for input when \`--files-only\` is unset:"
    for assembly in ${smallest_contigs[@]}
    do
      id=$( basename "$( up "${assembly}" "$id_offset" )" )
      first_line=$(head -n 1 "${assembly}" )
      echo -e "\nFilepath:     ${assembly}"
      echo "Species ID:   ${id}"
      echo "First header: ${first_line}"
    done
    exit 0
  elif (( combine_only ))
  then
    echo "Combining output files..."
    marker_count=$( combine_output "$outdir" )
    echo "Concatenated ${marker_count} markers into ${outdir}"
    exit 0
  fi

  # Insert species ID if there is none (follows format '>SPECIES@SEQUENCE_ID')
  echo "Inserting species names to query sequences..."
  for filename in ${smallest_contigs[@]}; do
    insert_species_prefix "$filename"
  done

  echo "Running Patchwork..."
  run_patchwork "${smallest_contigs[@]}"

  if (( skip_combine ))
  then
    echo "\`--skip-combine\` is set, Patchwork output combination is skipped"
    exit 0
  fi

  echo "Combining output files..."
  marker_count=$( combine_output "$outdir" )
  echo "Concatenated ${marker_count} markers into ${outdir}"
}

main
