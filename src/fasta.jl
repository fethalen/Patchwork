# Provides utilities and types for working with FASTA files

using FASTX
using BioSequences

"""
    fastafiles(directory::String, extensions::Array{String, 1})::Array{String, 1}

Takes the path to a directory and a list of acceptable FASTA filetype
extensions as an input and returns all of the FASTA files with that filetype
extension within the provided directory in form of a names.
"""
function fastafiles(directory::String, extensions::Array{String, 1})::Array{String, 1}
    absdir = abspath(directory)
    isdir(absdir) || error(*("the directory ", directory, " does not exist"))
    files = []

    for path in readdir(directory, join=true)
        if isfile(path)
            filetype_ext = strip(last(splitext(path)), '.')
            if filetype_ext in extensions
                push!(files, path)
            end
        end
    end

    return files
end

"""
    readmsa(fastafile::String)::MultipleSequenceAlignment

Takes the path to a multiple sequence alignment (MSA) as an input and returns a
MultipleSequenceAlignment object
"""
function readmsa(fastafile::String, delimiter::Char)::MultipleSequenceAlignment
    abs_fastafile = abspath(fastafile)
    isfile(abs_fastafile) || error(*("cannot locate file ", fastafile))
    record = FASTA.Record()
    msa = MultipleSequenceAlignment(abs_fastafile)

    open(FASTA.Reader, fastafile) do reader
        while !eof(reader)
            read!(reader, record)
            alignment = SequenceRecord(FASTA.identifier(record), FASTA.sequence(record),
                                        delimiter)
            addalignment!(msa, alignment)
        end
    end
    return msa
end

"""
    selectsequences(fastafile, identifiers; delimiter)

Searches within a `fastafile` for sequence records with identifiers that matches the
provided `identifiers`. The provided `delimiter` separates the OTU from the sequence
name (set to '@' by default).
"""
function selectsequences(fastafile::AbstractString, 
                         identifiers::Array{String};
                         delimiter='@'::Char)::MultipleSequenceAlignment
    abs_fastafile = abspath(fastafile)
    isfile(abs_fastafile) || error(*("cannot locate file ", fastafile))
    record = FASTA.Record()
    msa = MultipleSequenceAlignment(abs_fastafile)

    open(FASTA.Reader, fastafile) do reader
        while !eof(reader)
            read!(reader, record)
            identifier = FASTA.identifier(record)
            if identifier in identifiers
                alignment = SequenceRecord(identifier, FASTA.sequence(record),
                                           delimiter)
                addalignment!(msa, alignment)
            end
        end
    end
    return msa
end