using BioAlignments
using FASTX

const WIDTH = 80

function cleanfiles(paths::String...)
    for path in paths
        isfile(path) && rm(path, force = true)
        isdir(path) && rm(path, force = true, recursive = true)
    end
end

function warn_overwrite()::String
    print("WARNING: found output from a previous run, overwrite old files? (y/n):")
    answer = readline()
    while !isequal(answer, "y") && !isequal(answer, "n")
        println("Please answer 'y' for yes or 'n' for no")
        answer = readline()
    end
    return answer
end

function write_alignmentfile(file::AbstractString, id::AbstractString, contigs::Int64,
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
        print(io, "Reference ID:        " * id * "\n")
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
                     alignment::PairwiseAlignment, delimiter::Char)
    queryname = ""
    if delimiter in id.id
        queryname = otupart(id, delimiter)
    end
    fastawriter = FASTA.Writer(open(file, "a"))
    write(fastawriter, FASTA.Record(queryname, alignment.a.seq))
    close(fastawriter)
end
