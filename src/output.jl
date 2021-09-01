using BioAlignments
using FASTX

const WIDTH = 80

function cleanfiles(files::String...)
    for file in files
        isfile(file) && run(pipeline(`rm $file`))
    end
end

function warn_overwrite()::String
    println("WARNING: Found output files from a previous run.")
    println("Do you want to overwrite them? (y/n): ")
    answer = readline()
    while !isequal(answer, "y") && !isequal(answer, "n")
        println("Please type either 'y' for yes or 'n' for no.")
        answer = readline()
    end
    return answer
end

function write_alignmentfile(file::AbstractString, id::SequenceIdentifier, contigs::Int64,
                             alignment::BioAlignments.PairwiseAlignment, index::Int64)
    count = string(index) * ". "    
    subjectlength = length(alignment.b)
    querylength = length(alignment.a.seq)
    exactmatches = count_matches(alignment)
    mismatches = count_mismatches(alignment)
    deletions = count_deletions(alignment)
    aln_occupancy = Patchwork.occupancy(alignment)
    open(file, "a") do io
        print(io, count * repeat('-', WIDTH - length(count)) * "\n")
        print(io, "\n")
        print(io, "Reference ID:        " * id.id * "\n")
        print(io, "Reference Length:    " * string(subjectlength) * "\n")
        print(io, "Query Length:        " * string(querylength) * "\n")
        print(io, "Contigs:             " * string(contigs) * "\n")
        print(io, "Matches:             " * string(exactmatches) * "\n")
        print(io, "Mismatches:          " * string(mismatches) * "\n")
        print(io, "Deletions:           " * string(deletions) * "\n")
        print(io, "Occupancy:           " * string(aln_occupancy) * "\n")
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
