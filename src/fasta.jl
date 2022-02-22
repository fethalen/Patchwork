# Provides utilities and types for working with FASTA files

using FASTX
using CodecZlib
using BioSequences
using BioGenerics

const FASTQEXTENSIONS = ["fq", "fastq"]

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
    readmsa(file::String, removeduplicates::Bool=true)::MultipleSequenceAlignment

Takes the path to a multiple sequence alignment (MSA) as an input and returns a
MultipleSequenceAlignment object. If `removeduplicates` is set, duplicate sequences
will be removed before returning the object.
"""
function readmsa(
    file::String;
    removeduplicates::Bool=true,
    bysequence::Bool=true,
    byid::Bool=true
)::MultipleSequenceAlignment
    absfile = abspath(file)
    isfile(absfile) || error(*("cannot locate file ", file))
    msa = MultipleSequenceAlignment(absfile)

    removeduplicates && !bysequence && !byid && error("Please choose whether duplicates " *
        "should be removed based on sequence or ID.")

    if isfastafile(file)
        record = FASTA.Record()
        reader = isgzipcompressed(file) ?
            FASTA.Reader(GzipDecompressorStream(open(file))) : FASTA.Reader(open(file, "r"))
    elseif isfastqfile(file)
        record = FASTQ.Record()
        reader = isgzipcompressed(file) ?
            FASTQ.Reader(GzipDecompressorStream(open(file))) : FASTQ.Reader(open(file, "r"))
    else
        error("incorrect file type.")
    end

    while !eof(reader)
        read!(reader, record)
        alignment = SequenceRecord(FASTX.identifier(record), FASTX.sequence(record))
        # remove duplicated sequences (including reverse complement) and or IDs while reading
        if length(msa) >= 1 && removeduplicates
            (bysequence && (in(alignment.sequencedata, map(aln -> aln.sequencedata, msa.sequences)) ||
            in(reverse_complement(alignment.sequencedata), map(aln -> aln.sequencedata, msa.sequences))) ||
            byid && in(alignment.id.id, map(aln -> aln.id.id, msa.sequences))) && continue
            push!(msa, alignment)
        else
            push!(msa, alignment)
        end
    end
    close(reader)
    #removeduplicates && remove_duplicates!(msa, bysequence, byid) # not necessary when removing duplicates while reading
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

function selectsequence(
    fastafile::AbstractString,
    identifier::SequenceIdentifier
)
    abs_fastafile = abspath(fastafile)
    isfile(abs_fastafile) || error(*("cannot locate file ", fastafile))
    record = FASTA.Record()
    open(FASTA.Reader, fastafile) do reader
        while !eof(reader)
            read!(reader, record)
            currentid = FASTA.identifier(record)
            if isequal(currentid, identifier.id)
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
                push!(msa, alignment)
            end
        end
    end
    return msa
end

function isgzipcompressed(file::AbstractString)
    occursin(".", file) || return false
    extension = rsplit(file, ".", limit=2)[2]
    isequal(extension, "gz") && return true
    return false
end

function isfastafile(
    path::AbstractString,
    ext::AbstractVector{String}=FASTAEXTENSIONS
)::Bool
    # splits = split(path, ".")
    # length(splits) > 1 && last(splits) in ext && return true
    # length(splits) > 2 && splits[lastindex(splits)-1] in ext && isgzipcompressed(path) && return true
    # return false
	reader = isgzipcompressed(path) ? FASTA.Reader(GzipDecompressorStream(open(path))) : FASTA.Reader(open(path))
    record = FASTA.Record()
    try # in case the file is a tmp file without extension
        read!(reader, record)
        close(reader)
    catch e
        close(reader)
        return false
    end
    return true
end

# function clean_tmpfasta(file::AbstractString; bysequence::Bool=true, byid::Bool=true)
#     isfasta = isfastafile(file)
#     isfasta || isfastqfile(file) || error("Wrong file type: $file")

#     tmpfile, tmpio = mktemp()
#     type = isfasta ? FASTA : FASTQ
#     reader = type.Reader(open(file))
#     writer = type.Writer(tmpio)
#     record = type.Record()
#     lastrecord = type.Record()

#     if eof(reader)
#         close(reader)
#         close(writer)
#         return tmpfile
#     end
#     read!(reader, lastrecord)
#     write(writer, lastrecord)
#     while !eof(reader)
#         read!(reader, record)
#         (bysequence && isequal(FASTX.sequence(record), FASTX.sequence(lastrecord)) ||
#             byid && isequal(FASTX.identifier(record), FASTX.identifier(lastrecord))) && continue
#         write(writer, record)
#         lastrecord = record
#     end
#     close(reader)
#     close(writer)
#     return tmpfile
# end

"""
    countsequences(path)

Returns the number of sequence records (i.e., the number of `>`s) found in the provided
`path`.
"""
function countsequences(path::AbstractString)::Int
    !isfile(path) && error("path not found or not a file: $path")
    count = 0
    for line in readlines(path)
        if first(line) == '>'
            count += 1
        end
    end
    return count
end

