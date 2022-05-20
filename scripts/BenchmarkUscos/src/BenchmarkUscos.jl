# Example usage: julia BenchmarkUscos.jl --query '/media/feli/16tb_storage/Dimorphilus_gyrociliatus/spades_assembly_x_h_robusta_buscos_patchwork_out/d_gyrociliatus_x_h_robusta_no_stops.fa'\
# --reference /media/feli/16tb_storage/Dimorphilus_gyrociliatus/d_gyrociliatus_busco/d_gyrociliatus_single_copy_busco_seqs.faa\
# --output-dir test\
# --threads 28

"""
    Consensus

This package provides a quick and easy way of comparing BUSCOs obtained
from different datasets.
"""
module BenchmarkUscos

export
    otupart, sequencepart, splitdescription,otupart, sequencepart

using ArgParse
using CSV
using DataFrames

include("alignment.jl")
include("diamond.jl")
include("sequenceidentifier.jl")
include("simplestats.jl")

const WIDTH = 80

function warn_overwrite()::String
    print("WARNING: found output from a previous run, overwrite old files? (y/n):")
    answer = readline()
    while !isequal(answer, "y") && !isequal(answer, "n")
        println("Please answer 'y' for yes or 'n' for no")
        answer = readline()
    end
    return answer
end

function write_alignmentfile(
    index::Int,
    outfile::AbstractString,
    seqid::AbstractString,
    refid::AbstractString,
    row::Vector{Float64},
    alignment::BioAlignments.PairwiseAlignment
)
    count = string(index) * ". "
    open(outfile, "a") do io
        print(io, count * repeat('-', WIDTH - length(count)) * "\n")
        print(io, "\n")
        print(io, "Query ID:         " * seqid * "\n")
        print(io, "Reference ID:     " * refid * "\n")
        print(io, "Query length:     " * string(Int(row[1])) * "\n")
        print(io, "Reference length: " * string(Int(row[2])) * "\n")
        print(io, "Distance:         " * string(Int(row[7])) * "\n")
        print(io, "# matches:        " * string(Int(row[3])) * "\n")
        print(io, "# mismatches:     " * string(Int(row[4])) * "\n")
        print(io, "# insertions:     " * string(Int(row[5])) * "\n")
        print(io, "# deletions:      " * string(Int(row[6])) * "\n")
        print(io, "% query coverage: " * string(row[8]) * "\n")
        print(io, "% identity:       " * string(row[9]) * "\n")
        print(io, "\n")
        print(io, alignment)
        print(io, "\n")
    end
end

function ArgParse.parse_item(::Type{T}, argument::AbstractString) where T <: AbstractVector
    convert(T, split(argument, " "))
end

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

function uniquehits(blasthits::Array{DiamondSearchResult,1})
    hitdata = DataFrame(blasthits)
    transform!(hitdata, :subjectid => ByRow(id -> id.id) => :subjectid)
    transform!(hitdata, :queryid => ByRow(id -> id.id) => :queryid)
    sort!(hitdata, [:subjectid, :bitscore])
    return combine(groupby(hitdata, :subjectid), last)
end

function parse_parameters()
    overview = """
    Compare universal single-copy orthologs (USCOs) across different datasets.
    """
    settings = ArgParseSettings(description=overview,
                                version = "0.1.0",
                                add_version = true)
    @add_arg_table! settings begin
        "--query"
            help = "Path to a sequence file in FASTA format"
            required = true
            arg_type = String
            metavar = "PATH"
        "--reference"
            help = "Path to the file containing the single-copy orthologs in FASTA format"
            default = "../data/helobdella_robusta_uscos.faa"
            arg_type = String
            metavar = "PATH"
        "--output-dir"
            help = "Write output files to this directory"
            required = true
            arg_type = String
            default = "patchwork_output"
            metavar = "PATH"
        "--threads"
            help = "Number of threads to utilize (default: all available)"
            default = Sys.CPU_THREADS
            arg_type = Int64
            metavar = "NUMBER"
        "--overwrite"
            help = "Overwrite output from previous runs without warning"
            action = :store_true
    end

    return ArgParse.parse_args(settings)
end

function main()
    args = parse_parameters()
    queryseqs = args["query"]
    referenceseqs = args["reference"]
    outdir = args["output-dir"]

    if isdir(outdir)
        if !args["overwrite"]
            answer = warn_overwrite()
            isequal(answer, "n") && return
        end
        rm(outdir, recursive = true)
    end
    mkdir(outdir)

    !isfile(queryseqs) && error("query file not found: $queryseqs")
    !isfile(referenceseqs) && println("reference file not found: $referenceseqs")

    blastout = diamond_blastp(
        queryseqs, diamond_makeblastdb(referenceseqs, outdir), outdir)
    refseqs_count = countsequences(referenceseqs)

    alignmentstats = DataFrame(
        query_len = Int64[],
        reference_len = Int64[],
        no_matches = Int64[],
        no_mismatches = Int64[],
        no_insertions = Int64[],
        no_deletions = Int64[],
        distance = Int64[],
        query_cover = Float64[],
        percent_identity = Float64[]
    )

    sort!(alignmentstats, :percent_identity, rev = true)
    alignmentfile = outdir * "/" * "alignments.txt"

    for (count, hitrow) in enumerate(eachrow(uniquehits(blastout)))
        (alignment, alignmentscore) = scorealignment(hitrow.full_querysequence, hitrow.full_subjectsequence)
        push!(alignmentstats, alignmentscore)
        write_alignmentfile(count, alignmentfile, hitrow.queryid, hitrow.subjectid,
            alignmentscore, alignment.aln)
    end

    plot_percentident(alignmentstats.percent_identity, refseqs_count, outdir * "/percent_identity.png")
    plot_querycover(alignmentstats.query_cover, outdir * "/query_coverage.png")

    summary = select(describe(alignmentstats), Not([:nmissing, :eltype]))
    statsfile = outdir * "/" * "marker_stats.csv"
    summaryfile = outdir * "/" * "summary_stats.csv"
    println(summary)

    CSV.write(statsfile, alignmentstats)
    CSV.write(summaryfile, summary)
end

main()

end # module
