# Provides utilities and types for working with FASTA files

using FASTX
using BioSequences
using BioGenerics

"""
    fastafiles(directory::String, extensions::Array{String, 1})::Array{String, 1}

Takes the path to a directory and a list of acceptable FASTA filetype
extensions as an input and returns all of the FASTA files with that filetype
extension within the provided directory in form of a names.
"""
function fastafiles(
    directory::String,
    extensions::Array{String, 1}
)::Array{String, 1}
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
function readmsa(fastafile::String)::MultipleSequenceAlignment
    abs_fastafile = abspath(fastafile)
    isfile(abs_fastafile) || error(*("cannot locate file ", fastafile))
    record = FASTA.Record()
    msa = MultipleSequenceAlignment(abs_fastafile)

    open(FASTA.Reader, fastafile) do reader
        while !eof(reader)
            read!(reader, record)
            alignment = SequenceRecord(FASTA.identifier(record), FASTA.sequence(record))
            addalignment!(msa, alignment)
        end
    end
    return msa
end

"""
    get_fullseq(fastafile::String)::SequenceRecord

Returns the first and only the first sequence within the provided FASTA file as a
`SequenceRecord`.
"""
function get_fullseq(fastafile::String)::SequenceRecord
    abs_fastafile = abspath(fastafile)
    isfile(abs_fastafile) || error(*("cannot locate file ", fastafile))
    record = FASTA.Record()
    msa = MultipleSequenceAlignment(abs_fastafile)

    open(FASTA.Reader, fastafile) do reader
        while !eof(reader)
            read!(reader, record)
            firstrecord = SequenceRecord(FASTA.identifier(record), FASTA.sequence(record))
            return firstrecord
        end
    end
end

"""
    get_fullseq(fastafile::String, id::SequenceIdentifier)::SequenceRecord

Returns the first and only the first sequence within the provided FASTA file whose
identifier matches `id` as a `SequenceRecord`.
"""
function get_fullseq(fastafile::String, id::String)::SequenceRecord
    abs_fastafile = abspath(fastafile)
    isfile(abs_fastafile) || error(*("cannot locate file ", fastafile))
    record = FASTA.Record()
    msa = MultipleSequenceAlignment(abs_fastafile)

    open(FASTA.Reader, fastafile) do reader
        while !eof(reader)
            read!(reader, record)
            if isfilled(record)
                isequal(FASTA.identifier(record), id) && return SequenceRecord(
                        FASTA.identifier(record), FASTA.sequence(record))
            end
        end
    end

    return SequenceRecord()
end

"""
    selectsequence(fastafile, identifier)

Searches within a `fastafile` for a sequence record that matches the
provided `identifier`. `nothing` is return if no matching records were found.
"""
function selectsequence(
    fastafile::AbstractString,
    identifier::String
)
    abs_fastafile = abspath(fastafile)
    isfile(abs_fastafile) || error(*("cannot locate file ", fastafile))
    record = FASTA.Record()
    open(FASTA.Reader, fastafile) do reader
        while !eof(reader)
            read!(reader, record)
            currentid = FASTA.identifier(record)
            if currentid == identifier
                return SequenceRecord(currentid, FASTA.sequence(record))
            end
        end
    end
end

"""
    selectsequences(fastafile, identifiers)

Searches within a `fastafile` for sequence records with identifiers that matches the
provided `identifiers`.
"""
function selectsequences(
    fastafile::AbstractString,
    identifiers::Array{String}
)::MultipleSequenceAlignment
    abs_fastafile = abspath(fastafile)
    isfile(abs_fastafile) || error(*("cannot locate file ", fastafile))
    record = FASTA.Record()
    msa = MultipleSequenceAlignment(abs_fastafile)

    open(FASTA.Reader, fastafile) do reader
        while !eof(reader)
            read!(reader, record)
            identifier = FASTA.identifier(record)
            if identifier in identifiers
                alignment = SequenceRecord(identifier, FASTA.sequence(record))
                addalignment!(msa, alignment)
            end
        end
    end
    return msa
end
