function isfastqfile(
    file::AbstractString,
    ext::AbstractVector{String}=FASTQEXTENSIONS
)::Bool
    # splits = split(file, ".")
    # length(splits) > 1 && last(splits) in ext && return true
    # length(splits) > 2 && splits[lastindex(splits)-1] in ext && isgzipcompressed(file) && return true
    # return false
    reader = isgzipcompressed(file) ? FASTQ.Reader(GzipDecompressorStream(open(file))) : FASTQ.Reader(open(file))
    record = FASTQ.Record()
    try # in case the file is a tmp file without extension
        read!(reader, record)
    catch e
        close(reader)
        return false
    end
    close(reader)
    return true
end

function fastq2fasta(
    infile::AbstractString;
    ext::Vector{String}=FASTQEXTENSIONS,
    removeduplicates::Bool=true,
    bysequence::Bool=true,
    byid::Bool=true
)::AbstractString
    isfastqfile(infile, ext) || error("Wrong file format; input must be a fastq file.")
    msa = readmsa(infile; removeduplicates=removeduplicates, bysequence=bysequence, byid=byid)
    return mktemp_fasta(msa)
end

function fastq2fasta(
    infiles::Vector{String};
    ext::Vector{String}=FASTQEXTENSIONS,
    removeduplicates::Bool=true,
    bysequence::Bool=true,
    byid::Bool=true
)::Vector{AbstractString}
    output = String[]
    for file in infiles
        outfile = fastq2fasta(file; ext=ext, removeduplicates=removeduplicates, bysequence=bysequence, byid=byid)
        push!(output, outfile)
    end
    return output
end

function splitfile(path::AbstractString; recordsperfile::Int64=1000000)
    isfastqfile(path) || error("Incorrect file type.")

    files = Tuple{String, Int64}[]
    record = FASTQ.Record()
    reader = isgzipcompressed(path) ? FASTQ.Reader(GzipDecompressorStream(open(path))) : FASTQ.Reader(open(path, "r"))

    count = 1
    msa = MultipleSequenceAlignment()
    while !eof(reader)
        if count > recordsperfile # done reading `recordsperfile` records into msa
            #msa = sort(msa)
            remove_duplicates!(msa) # includes sorting by sequence and by ID
            tmpfile = mktemp_fasta(msa) # write
            # push!(files, tmpfile)
            push!(files, (tmpfile, length(msa))) # keep the file length
            msa = MultipleSequenceAlignment() # init
            count = 1
        end
        read!(reader, record)
        push!(msa, SequenceRecord(FASTX.identifier(record), FASTX.sequence(record)))
        count += 1
    end

    if !isempty(msa) # eof(reader)
        tmpfile = mktemp_fasta(msa)
        push!(files,(tmpfile, length(msa)))
    end

    return files # empty if input file empty
end

function combinefiles(files::AbstractVector{String})
    file, io = mktemp() # deduplicated (output) fasta file
    writer = FASTA.Writer(io)
    # stores tuples of record and index of the reader it came from:
    seqrecords = repeat([(SequenceRecord(), 0)], length(files))
    record = FASTA.Record()
    # one reader for each infile:
    readers = Vector{FASTA.Reader}(undef, length(files))
    # index of the reader from which the next record should be read:
    nextreader_index = 0
    for (i, f) in enumerate(files) # init readers (open each file) and seqrecords
        readers[i] = FASTA.Reader(open(f))
        eof(readers[i]) && continue # empty input file; should not happen but who knows...
        read!(readers[i], record)
        seqrecord = SequenceRecord(FASTA.identifier(record), FASTA.sequence(record))
        seqrecords[i] = (seqrecord, i)
    end
    filter!(tuple -> tuple[2] != 0, seqrecords) # keep only initialised (non-empty infiles)
    if isempty(seqrecords)
        close.(readers)
        close(writer)
        return file
    end

    sort!(seqrecords, by = tuple -> tuple[1]) # sort by sequence
    mappedsequences = map(tuple -> tuple[1].sequencedata, seqrecords)
    while !all(eof.(readers)) # > 1 stream still open
        # first record should be written to outfile (--> keep outfile sorted)
        record = FASTA.Record(seqrecords[1][1].id.id, seqrecords[1][1].sequencedata)
        # the reader it came from should provide the next record...
        nextreader_index = seqrecords[1][2]
        write(writer, record)
        deleteat!(seqrecords, 1)
        # ...but if eof, continue: write first record, get nextreader_index, delete record.
        # At some point, all mentions of this particular index will have been removed as the
        # records get written to the output and deleted from the array.
        eof(readers[nextreader_index]) && continue
        # if not eof, read next record...
        read!(readers[nextreader_index], record)
        seqrecord = SequenceRecord(FASTA.identifier(record), FASTA.sequence(record))
        # ...and add to array only if it doesn't yet contain any of sequence, revcomp or ID
        if !(in(seqrecord.sequencedata, mappedsequences) ||
            in(reverse_complement(seqrecord.sequencedata), mappedsequences) ||
            in(seqrecord.id.id, map(tuple -> tuple[1].id.id, seqrecords)))
            push!(seqrecords, (seqrecord, nextreader_index))
            sort!(seqrecords, by = tuple -> tuple[1]) # keep array sorted at all times
            mappedsequences = map(tuple -> tuple[1].sequencedata, seqrecords)
        end # else, skip this record
    end
    # all readers eof --> sort and write stored records
    sort!(seqrecords, by = tuple -> tuple[1])
    for seqrecord in seqrecords
        record = FASTA.Record(seqrecord[1].id.id, seqrecord[1].sequencedata)
        write(writer, record)
    end
    close.(readers) # close all readers
    close(writer)
    return file
end
