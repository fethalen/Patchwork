#!/usr/bin/env bash

readonly VERSION='0.0.1'
basedir=''
separator='@'
fasta_extension='fas'

usage() { printf "%s" "\
usage:
  insert_gene_name.sh [--help] [--version] DIRECTORY

description:
  For each file in the provided directory, insert the filename into each into
  each sequence headers in all FASTA files within that directory.

example:
  16s.fas:
    >Drosophila_melanogaster -> >Drosophila_melanogaster@16s
  18s.fas:
    >Drosophila_melanogaster -> >Drosophila_melanogaster@18s

options:
  -h, --help          display this help message and exit
  -v, --version       display the version number and exit
  -s, --separator     use this string to separate the species name from the
                      sequence ID (default: ${separator})
  -e, --extension     search for FASTA files with this extension (default: ${fasta_extension})
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
  echo "insert_gene_name: error: $message"
  exit "$status"
}

insert_into_headers() {
  local string="$1"
  local filename="$2"
  local separator="$3"
  sed -i "s/\(^>.*\)/\1${separator}${string}/g" "$filename"
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
    -s | --separator)
      separator="$2"
      shift # past argument
      shift # past value
      ;;
    -e | --extension)
      fasta_extension="$2"
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
  readarray -d '' sequence_files <\
    <(find "$basedir" -type f -name "*.${fasta_extension}" -print0 )
  echo "Found ${#sequence_files[@]} files with extension '${fasta_extension}'"
  # insert_into_headers "test" "${sequence_files[1]}" "$separator"
  for filename in "${sequence_files[@]}"
  do
    gene_name=$( basename "${filename%\."${fasta_extension}"}" )
    insert_into_headers "$gene_name" "$filename" "$separator"
  done
}

main
