using FASTX

function splitfile(path::AbstractString; recordsperfile::Int64=1000000)
    isfastqfile(path) || error("Incorrect file type.")

    files = []
    record = FASTQ.Record()
    reader = isgzipcompressed(path) ? FASTQ.Reader(GzipDecompressorStream(open(file))) : FASTQ.Reader(open(file, "r"))

    count = 1
    msa = MultipleSequenceAlignment()
    # tmpfile, tmpio = mktemp()
    # writer = FASTQ.Writer(tmpio)
    while !eof(reader)
        if count > recordsperfile
            msa = sort(msa)
            remove_duplicates!(msa)
            tmpfile = mktemp_fasta(msa)
            push!(files, tmpfile)
            # close(writer)
            # tmpfile, tmpio = mktemp()
            # writer = FASTQ.Writer(tmpio)
            msa = MultipleSequenceAlignment()
            count = 1
        end
        read!(reader, record)
        push!(msa, SequenceRecord(FASTX.identifier(record), FASTX.sequence(record)))
        # write(writer, record)
        # flush(writer)
        count += 1
    end
    
    count <= recordsperfile && push!(files, tmpfile)
    # close(writer)
    # close(reader)
    return files
end

function combinefiles(files::AbstractVector{String})
    file, io = mktemp()
    writer = FASTA.Writer(io)
    seqrecords = repeat([(SequenceRecord(), 0)], length(files))
    record = FASTA.Record()
    readers = Vector{FASTA.Reader}(undef, length(files))
    nextreader_index = 0
    for (i, f) in enumerate(files)
        readers[i] = FASTA.Reader(open(f))
        eof(readers[i]) && continue
        read!(readers[i], record)
        seqrecord = SequenceRecord(FASTA.identifier(record), FASTA.sequence(record))
        seqrecords[i] = (seqrecord, i)
    end
    filter!(tuple -> tuple[2] != 0, seqrecords)
    while length(readers) >= 1 && any(!eof.(readers))
        sort!(seqrecords, by = tuple -> tuple[1])
        record = FASTA.Record(seqrecords[1][1].id.id, seqrecords[1][1].sequencedata)
        nextreader_index = seqrecords[1][2]
        write(writer, record)
        deleteat!(seqrecords, 1)
        if eof(readers[nextreader_index]) 
            close(readers, nextreader_index)
            deleteat!(readers, nextreader_index)
            continue
        end
        read!(readers[nextreader_index], record)
        seqrecord = SequenceRecord(FASTA.identifier(record), FASTA.sequence(record))
        push!(seqrecords, seqrecord)
    end
    sort!(seqrecords, by = tuple -> tuple[1])
    for seqrecord in seqrecords
        record = FASTA.Record(seqrecord.id.id, seqrecord.sequencedata)
        write(writer, record)
    end
    length(readers) == 0 || close.(readers) # REMOVE if @assert doesn't error
    @assert length(readers) != 0
    close(writer)
    return file
end