using BioAlignments
using FASTX

function cleanfiles(files::String...)
    for file in files
        isfile(file) && run(pipeline(`rm $file`))
    end
end

function write_alignmentfile(file::AbstractString, id::SequenceIdentifier, contigs::Int64,
                             alignment::BioAlignments.PairwiseAlignment)
    subjectlength = length(alignment.b)
    querylength = length(alignment.a.seq)
    exactmatches = count_matches(alignment)
    mismatches = count_mismatches(alignment)
    #unspecifed = count(alignment, OP_MATCH)
    deletions = count_deletions(alignment)
    aln_occupancy = Patchwork.occupancy(alignment)
    open(file, "a") do io
        #if filesize(file) == 0
        #    print(io, "reference_ID\treference_length\tquery_length\tcontigs\tmatches\t" * 
        #              "mismatches\tdeletions\toccupancy\n")
        #end
        print(io, "reference_ID: " * id.id * "\t")
        print(io, "reference_length" * string(subjectlength) * "\t")
        print(io, "query_length" * string(querylength) * "\t")
        print(io, "contigs: " * string(contigs) * "\t")
        print(io, "matches: " * string(exactmatches) * "\t")
        print(io, "mismatches: " * string(mismatches) * "\t")
        print(io, "deletions: " * string(deletions) * "\t")
        print(io, "occupancy: " * string(aln_occupancy) * "\n")
        print(io, "\n")
        print(io, alignment)
        print(io, "\n")
    end
end

function write_fasta(file::AbstractString, id::SequenceIdentifier, 
                     alignment::PairwiseAlignment)
    queryname = otupart(id)
    fastawriter = FASTA.Writer(open(file, "a"))
    write(fastawriter, FASTA.Record(queryname, alignment.a.seq))
    close(fastawriter)
end
