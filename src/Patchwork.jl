# julia --trace-compile=precompiled.jl Patchwork.jl --contigs "../test/07673_lcal.fa" --reference "../test/07673_Alitta_succinea.fa" --diamond-flags "--frameshift 15 --ultra-sensitive" --output-dir "../test/patchwork-output"
# diamond blastx --query 07673_dna.fa --db 07673_Alitta_succinea.fa --outfmt 6 qseqid qseq full_qseq qstart qend qframe sseqid sseq sstart send cigar pident bitscore --out diamond_results.tsv --frameshift 15

module Patchwork

export
    # types
    SequenceIdentifier,
    SequenceRecord,
    MultipleSequenceAlignment,
    DiamondSearchResult,
    AlignedRegion,
    AlignedRegionCollection,

    # alignedregion
    cleancigar, query_isreverse, query2subject, subject2query, subject2fullsubject,
    subject2fullquery, fullsubject2subject, fullsubject2query, fullsubject2fullquery,
    query2fullquery, fullsubject_queryboundaries, subject_queryboundaries, slicealignment,
    query_leftposition, query_rightposition, subject_leftposition, subject_rightposition,
    subjectinterval, queryinterval, queryframe, isordered, precedes, #isoverlapping,
    samerange, identifier, sameid, samesequence, merge, leftintersection, rightintersection,
    overlap, beforeoverlap, afteroverlap, shadows, equalrange, totalrange, longest,
    isnucleotide, isaminoacid,

    # alignedregioncollection
    sameids, uniquesequences, has_uniquesequences, hasoverlaps, show_subjectintervals,
    mergeoverlaps, isnucleotide, isaminoacid, queryids,

    # alignment
    DEFAULT_SCOREMODEL, pairalign_global, pairalign_local, order, realign,

    # alignmentconcatenation
    createbridgealignment, concatenate, countmatches, occupancy, maskgaps, countgaps,

    # checkinput
    MATRICES, GAPS, GAPDEFAULTS, #checkmatrixtype, getmatrixtype, setmatrixname,
    #read_custommatrix, getmatrix, setgapopen, setgapextend, checkgappenalty,
    #checkdiamondflags, checkmakedbflags, setpatchworkflags, setdiamondflags,
    #collectdiamondflags,

    # diamond
    FIELDS, OUTPUT_FORMAT, readblastTSV, writeblastTSV, diamond_blastx, diamond_makeblastdb,
    queryids, subjectids, isdiamonddatabase,

    # fasta
    FASTQEXTENSIONS, fastafiles, readmsa, get_fullseq, selectsequence, isfastafile, #isfastqfile,
    isgzipcompressed, #fastq2fasta,

    # fastq
    isfastqfile, fastq2fasta, splitfile, combinefiles,

    # filtering
    #remove_duplicates,

    # multiplesequencealignment
    addalignment, removealignment, hasgaps, otus, otufrequencies, countotus, coverage,
    equal_length, gapmatrix, gapfrequencies, mktemp_fasta, remove_duplicates,
    remove_duplicates!, pool,

    # output
    #WIDTH, cleanfiles, warn_overwrite, write_alignmentfile, write_fasta,

    # sequenceidentifier
    splitdescription, otupart, sequencepart,

    # sequencerecord
    missingdata, gappositions, nongap_range, has_compoundregions, fillmissing, compoundregions,
    compoundranges, nongap_ranges

using Base: Bool, Int64, func_for_method_checked, DEFAULT_COMPILER_OPTS, Cint
using ArgParse
using BioAlignments
using CSV
using FASTX
using DataFrames
using Statistics

include("sequenceidentifier.jl")
include("sequencerecord.jl")
include("multiplesequencealignment.jl")
include("filtering.jl")
include("diamond.jl")
include("alignedregion.jl")
include("alignedregioncollection.jl")
include("alignment.jl")
include("alignmentconcatenation.jl")
include("checkinput.jl")
include("fasta.jl")
#include("fastq.jl")
include("output.jl")

const EMPTY = String[]
const DIAMONDFLAGS = ["--ultra-sensitive"]
const MIN_DIAMONDVERSION = "2.0.3"
const MATRIX = "BLOSUM62"
const ALIGNMENTOUTPUT = "alignments.txt"
const FASTAOUTPUT = "query_sequences"
const DEFAULT_FASTA_EXT = ".fas"
const DIAMONDOUTPUT = "diamond_out"
const STATSOUTPUT = "sequence_stats"

"""
    printinfo()

Print program title and basic information.
"""
function printinfo(
    diamondversion::AbstractString,
    threads::Int64,
    outputdir::AbstractString
)
    about = """
       ◣
      ◢◤  P A T C H W O R K
    ◢◤
    ◥  ◣  Developers: Felix Thalén & Clara G. Köhne
      ◢◤  Contact   : <felix.thalen@uni-goettingen.de>
    ◢◤    Wiki      : github.com/fethalen/patchwork/wiki
    ◥  ◣  Cite      : In prep.
      ◢◤  DIAMOND v.: $diamondversion
    ◢◤    Threads   : $threads
    ◥     Output    : $outputdir
    """
    println(about)
    return
end

"""
    min_diamondversion(version)

Returns `true` if DIAMOND is installed with a version number equal to or higher than the
provided `minversion`. Returns false if the version number is lower or DIAMOND is not found at
all.
"""
function min_diamondversion(minversion::AbstractString)::Bool
    version = get_diamondversion()
    isempty(version) && return false
    diamondversion_vector = split(version, ".")
    minversion_vector = split(minversion, ".")
    for (v, r) in zip(diamondversion_vector, minversion_vector)
        version = parse(Int64, v)
        required = parse(Int64, r)
        version < required && return false
        version > required && return true
    end
    return true
end

function get_diamondversion()::AbstractString
    try
        run(`diamond --version`, wait=false)
    catch
        return ""
    end
    versioncmd = read(`diamond --version`, String)
    return last(split(versioncmd))
end

function ArgParse.parse_item(::Type{T}, argument::AbstractString) where T <: AbstractVector
    convert(T, split(argument, " "))
end

function parse_parameters()
    overview = """
    Alignment-based retrieval and concatenation of phylogenetic markers from
    whole-genome sequencing data
    """
    settings = ArgParseSettings(description=overview,
                                version = "0.1.0",
                                add_version = true)
    @add_arg_table! settings begin
        "--contigs"
            help = "Path to one or more sequences in FASTA format"
            #required = true #Not required if reference is a .tsv
            arg_type = String
            nargs = '*'
            metavar = "PATH"
        "--reference"
            #help = "Either (1) a path to one or more sequences in FASTA format, (2) a
            #        subject database (DIAMOND or BLAST database), or (3) a DIAMOND output
            #        file in tabular format. For (3), set the `--tabular` flag."
            help = "A path to one or more sequences in FASTA format. Additionally, you can
                    also provide a DIAMOND output file in tabular format (use --search-results)
                    or a DIAMOND or BLAST database (use --database)."
            required = true
            arg_type = String
            nargs = '+'
            metavar = "PATH"
        "--search-results"
            help = "Provide a DIAMOND output file in tabular format. The first line of
                    such files is considered to be the header. Please adhere to Patchwork's
                    DIAMOND output format:
                    6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send
                    evalue bitscore qframe sseq seq."
            arg_type = String
            nargs = '?'
            metavar = "PATH"
        "--database"
            help = "Provide a subject DIAMOND or BLAST database to search against."
            arg_type = String
            nargs = '?'
            metavar = "PATH"
        "--output-dir"
            help = "Write output files to this directory"
            arg_type = String
            default = "patchwork_output"
            metavar = "PATH"
        "--diamond-flags"
            help = "Flags sent to DIAMOND"
            arg_type = Vector{String}
            default = DIAMONDFLAGS
            metavar = "LIST"
        "--fasta-extension"
            help = "Filetype extension used for output FASTA files"
            arg_type = String
            default = DEFAULT_FASTA_EXT
            metavar = "STRING"
        "--makedb-flags"
            help = "Flags sent to DIAMOND makedb"
            arg_type = Vector{String}
            default = EMPTY
            metavar = "LIST"
        "--matrix"
            help = "Set scoring matrix"
            arg_type = String
            metavar = "NAME"
        "--custom-matrix"
            help = "Use a custom scoring matrix"
            arg_type = String
            metavar = "PATH"
        "--gapopen"
            help = "Set gap open penalty (positive integer)"
            arg_type = Int64
            metavar = "NUMBER"
        "--gapextend"
            help = "Set gap extension penalty (positive integer)"
            arg_type = Int64
            metavar = "NUMBER"
        "--retain-stops"
            help = "Do not remove stop codons (`*`) in the output sequences"
            action = :store_true
        "--threads"
            help = "Number of threads to utilize (default: all available)"
            default = Sys.CPU_THREADS
            arg_type = Int64
            metavar = "NUMBER"
        "--species-delimiter"
            help = "Used to distinguish the OTU from the rest in sequence IDs"
            default = '@'
            arg_type = Char
            metavar = "CHARACTER"
        "--wrap-column"
            help = "Wrap output sequences at this column number (default: no wrap)"
            default = 0
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
    if !min_diamondversion(MIN_DIAMONDVERSION)
        error("Patchwork requires \'diamond\' with a version number above
               $MIN_DIAMONDVERSION to run")
    end
    printinfo(get_diamondversion(), args["threads"], args["output-dir"])

    setpatchworkflags!(args)
    setdiamondflags!(args)

    if length(args["reference"]) == 1
        references_file = args["reference"][1]              # 1 .fa
    else                                                    # multiple .fa files
        references_file = mktemp_fasta(pool(args["reference"]; removeduplicates=false))
    end
    # TODO: queries doesn't need to be stored as a MSA all the time!
    # You just need 1 fasta file for DIAMOND.
    queries = pool(args["contigs"])                          # MultipleSequenceAlignment
    outdir = args["output-dir"]
    alignmentoutput = outdir * "/" * ALIGNMENTOUTPUT
    fastaoutput = outdir * "/" * FASTAOUTPUT
    diamondoutput = outdir * "/" * DIAMONDOUTPUT
    statsoutput = outdir * "/" * STATSOUTPUT
    statistics = DataFrame(id = String[],
                           length_reference = Int[],
                           length_query = Int[],
                           regions = Int[],
                           contigs = Int[],
                           matches = Int[],
                           mismatches = Int[],
                           deletions = Int[],
                           occupancy = Float64[])

    if (isfile(alignmentoutput) || isdir(statsoutput) && !isempty(readdir(statsoutput))
        || (isdir(fastaoutput) && !isempty(readdir(fastaoutput))))
        if !args["overwrite"]
            answer = warn_overwrite()
            isequal(answer, "n") && return
        end
        cleanfiles(alignmentoutput, statsoutput, fastaoutput)
    end

    map(mkpath, [diamondoutput, fastaoutput, statsoutput])

    if isnothing(args["search-results"])
        if isempty(queries)
            println("Please provide one or more query files if not running with
                `--search-results` mode.")
            return
        end
        if isnothing(args["database"])
            println("Building DIAMOND database")
            reference_db = diamond_makeblastdb(references_file, outdir, args["makedb-flags"])
        else
            reference_db = args["database"]
        end
        diamondparams = collectdiamondflags(args)
        println("Aligning query sequences against reference database")
        diamondsearch = diamond_blastx(queries, reference_db, outdir, diamondparams)
    else
        diamondsearch = args["search-results"]
    end

    println("Merging overlapping hits")
    allhits = readblastTSV(diamondsearch)
    referenceids = unique(subjectids(allhits))

    for (index, referenceid) in enumerate(referenceids)
        diamondhits = filter(hit -> isequal(referenceid, hit.subjectid), allhits)
        @assert length(diamondhits) != 0
        writeblastTSV(*(diamondoutput, "/", referenceid.id, ".tsv"), diamondhits; header = true)
        #regions = AlignedRegionCollection(get_fullseq(args["reference"][index]), diamondhits)
        regions = AlignedRegionCollection(selectsequence(references_file, referenceid), diamondhits)
        # assuming all queries belong to same species:
        mergedregions = mergeoverlaps(regions)
        concatenation = concatenate(mergedregions, args["species-delimiter"])
        finalalignment = maskgaps(concatenation).aln

        write_alignmentfile(alignmentoutput, referenceid, length(regions), finalalignment, index)
        write_fasta(
            *(fastaoutput, "/", sequencepart(referenceid), args["fasta-extension"]),
            regions.records[1].queryid,
            finalalignment
        )
        stats_row = [
            mergedregions.referencesequence.id.id,
            length(mergedregions.referencesequence),
            length(finalalignment.a.seq),
            length(mergedregions), # only the merged regions
            length(regions), # all contigs retrieved by DIAMOND search
            BioAlignments.count_matches(finalalignment),
            BioAlignments.count_mismatches(finalalignment),
            BioAlignments.count_deletions(finalalignment),
            round(occupancy(finalalignment), digits=2)
        ]
        push!(statistics, stats_row)
    end

    CSV.write(*(statsoutput, "/statistics.csv"), statistics, delim="\t")
    averages = DataFrame(mean_length_query = Float64[],
                         mean_no_regions = Float64[],
                         mean_no_contigs = Float64[],
                         mean_no_matches = Float64[],
                         mean_no_mismatches = Float64[],
                         mean_no_deletions = Float64[],
                         mean_occupancy = Float64[])
    push!(averages, map(col -> mean(col), eachcol(select(statistics,
          Not([:id, :length_reference])))))
    CSV.write(*(statsoutput, "/average.csv"), averages, delim="\t")
    println("\nStatistics:")
    println("  Description                Value")
    output = ["Mean query length", "Mean no. of regions", "Mean no. of contigs",
              "Mean no. of matches", "Mean no. of mismatches", "Mean no. of deletions",
              "Mean occupancy"]
    for (name, value) in zip(output, averages[1, :])
        println("  ", name, repeat(" ", 25 - length(name)), "  ", round(value, digits=2))
    end
end

function julia_main()::Cint
    try
        main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end

    return 0
end

if length(ARGS) >= 2
    julia_main()
end

end # module
