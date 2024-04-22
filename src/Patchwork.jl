# julia --trace-compile=precompiled.jl Patchwork.jl --contigs "../test/07673_lcal.fa" --reference "../test/07673_Alitta_succinea.fa" --diamond-flags "--frameshift 15 --ultra-sensitive" --output-dir "../test/patchwork-output"
# julia --project=. src/Patchwork.jl --contigs "test/07673_lcal.fa" --reference "test/07673_Alitta_succinea.fa" --frameshift 15 --ultra-sensitive --output-dir "test/patchwork-output" --overwrite
# diamond blastx --query 07673_dna.fa --db 07673_Alitta_succinea.fa --outfmt 6 qseqid qseq full_qseq qstart qend qframe sseqid sseq sstart send cigar pident bitscore --out diamond_results.tsv --frameshift 15

module Patchwork

using ArgParse
using Base: Bool, Int64, func_for_method_checked, DEFAULT_COMPILER_OPTS, Cint
using BioAlignments
using BioSequences
using BioSymbols
using CSV
using CodecZlib
using DataFrames
using FASTX
using PackageCompiler
using PrettyTables
using Statistics

import Base
import BioGenerics
import Pkg
import Plots
import Random
import UnicodePlots

include("sequenceidentifier.jl")
include("sequencerecord.jl")
include("fasta.jl")
include("fastq.jl")
include("multiplesequencealignment.jl")
include("diamond.jl")
include("alignedregion.jl")
include("alignedregioncollection.jl")
include("alignment.jl")
include("filtering.jl")
include("alignment.jl")
include("alignmentconcatenation.jl")
include("checkinput.jl")
include("output.jl")
include("plotting.jl")
include("sliding.jl")
include("clustering.jl")

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
    subjectinterval, queryinterval, queryframe, isordered, precedes, isoverlapping,
    samerange, identifier, sameid, samesequence, merge, leftintersection, rightintersection,
    overlap, beforeoverlap, afteroverlap, shadows, equalrange, totalrange, longest,
    isnucleotide, isaminoacid,

    # alignedregioncollection
    sameids, uniquesequences, has_uniquesequences, hasoverlaps, show_subjectintervals,
    mergeoverlaps, isnucleotide, isaminoacid, queryids,

    # alignment
    DEFAULT_SCOREMODEL, pairalign_global, pairalign_local, order, realign,
    gapexcluded_identity,

    # alignmentconcatenation
    createbridgealignment, concatenate, countmatches, occupancy, maskgaps, countgaps,

    # checkinput
    MATRICES, GAPS, GAPDEFAULTS, DMND_ITERATE, #checkmatrixtype, getmatrixtype, setmatrixname,
    #read_custommatrix, getmatrix, setgapopen, setgapextend, checkgappenalty,
    #checkdiamondflags, checkmakedbflags, setpatchworkflags, setdiamondflags,
    #collectdiamondflags,

    # diamond
    FIELDS, OUTPUT_FORMAT, readblastTSV, frameshift, writeblastTSV, diamond_blastx,
    diamond_makeblastdb, queryids, subjectids, isdiamonddatabase,

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
    remove_duplicates!, pool, cat,

    # output
    #WIDTH, cleanfiles, warn_overwrite, write_alignmentfile, write_fasta,

    # sequenceidentifier
    splitdescription, otupart, sequencepart,

    # sequencerecord
    missingdata, gappositions, nongap_range, has_compoundregions, fillmissing, compoundregions,
    compoundranges, nongap_ranges

const EMPTY = String[]
const DIAMONDFLAGS = ["--iterate", "--evalue", "0.001"]
const MIN_DIAMONDVERSION = "2.0.10" # for --iterate option
const MATRIX = "BLOSUM62"
const ALIGNMENTOUTPUT = "untrimmed_alignments.txt"
const TRIMMEDALIGNMENT_OUTPUT = "trimmed_alignments.txt"
const FASTAOUTPUT = "protein_query_sequences"
const DNAFASTAOUTPUT = "nucleotide_query_sequences"
const DEFAULT_FASTA_EXT = ".fas"
const DIAMONDOUTPUT = "diamond_out"
const STATSOUTPUT = "sequence_stats"
const PLOTSOUTPUT = "plots"
const RULER = repeat('─', 74)
const VERSION = "0.6.4"
# Default scoremodel taken from DIAMOND's defaults for BLOSUM62
const DEFAULT_SCOREMODEL = BioAlignments.AffineGapScoreModel(BioAlignments.BLOSUM62,
    gap_open = -11, gap_extend = -1)
const GZIP_MAGIC_NUMBER = UInt8[0x1f, 0x8b]

"""
    printinfo()

Print program title and basic information.
"""
function printinfo(
    diamondversion::AbstractString,
    threads::Int64
)
    about = """
    P A T C H W O R K v$VERSION
    - Developers: Felix Thalén & Clara G. Köhne       - Cite (DOI): 10.1093/gbe/evad227
    - Contact   : <felix.thalen@cardio-care.ch>       - DIAMOND v.: $diamondversion
    - Wiki      : github.com/fethalen/patchwork/wiki  - Threads   : $threads
    $RULER"""
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
        run(`diamond --version`, wait = false)
    catch
        return ""
    end
    versioncmd = read(`diamond --version`, String)
    return last(split(versioncmd))
end

function getpercentage(
    seqsin::Int,
    seqsout::Int
)::Float64
    percent_markersout = 0.0
    if seqsout > 0
        percent_markersout = round((seqsin / seqsout) * 100, digits = 1)
    end
    return percent_markersout
end

function ArgParse.parse_item(::Type{T}, argument::AbstractString) where {T<:AbstractVector}
    convert(T, split(argument, " "))
end

function parse_parameters()
    overview = """
    Alignment-based retrieval and concatenation of phylogenetic markers from
    whole-genome sequencing data
    """
    settings = ArgParseSettings(description = overview,
        version = VERSION,
        add_version = true)

    add_arg_group!(settings, "input/output")
    @add_arg_table! settings begin
        "--contigs"
        help = """PATH to 1+ nucleotide sequence files in FASTA or FASTQ format. Can be GZip
            compressed."""
        arg_type = String
        nargs = '+'
        metavar = "PATH"
        "--reference"
        help = "PATH to 1+ amino acid sequence files in the FASTA format."
        required = true
        arg_type = String
        nargs = '+'
        metavar = "PATH"
        "--search-results"
        help = """PATH to a tabular DIAMOND output file, with one header line in format:
            6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send
            evalue bitscore qframe sseq seq."""
        arg_type = String
        metavar = "PATH"
        "--database"
        help = "Path to a subject DIAMOND or BLAST database to search against."
        arg_type = String
        metavar = "PATH"
        "--matrix"
        help = "Specifies the NAME of the scoring matrix"
        arg_type = String
        metavar = "NAME"
        "--custom-matrix"
        help = "PATH to a custom scoring matrix"
        arg_type = String
        metavar = "PATH"
        "--species-delimiter"
        help = "Set the CHARACTER used to separate the OTU from the rest in sequence IDs"
        default = '@'
        arg_type = Char
        metavar = "CHARACTER"
        "--fasta-extension"
        help = "Filetype extension used for output FASTA files"
        arg_type = String
        default = DEFAULT_FASTA_EXT
        metavar = "STRING"
        "--wrap-column"
        help = "Wrap output sequences at column NUMBER. 0 = no wrap"
        default = 0
        arg_type = Int64
        metavar = "NUMBER"
        "--no-plots"
        help = "Do not include plots"
        action = :store_true
        "--output-dir"
        help = "Write output files to this directory PATH"
        arg_type = String
        default = "patchwork_output"
        metavar = "PATH"
        "--overwrite"
        help = "Overwrite old content in the output directory"
        action = :store_true
    end # Input/output

    add_arg_group!(settings, "DIAMOND BLASTX")
    @add_arg_table! settings begin
        "--query-gencode"
        help = """Genetic code used for translation of query sequences. A list
            of possible values can be found on the NCBI website. Standard Code is
            used by default"""
        arg_type = Int64
        metavar = "NUMBER"
        "--strand"
        help = """Specifies the strand of the query. Possible values are:
            'both', 'plus', and 'minus'. Both strands are searched by default"""
        arg_type = String
        metavar = "STRING"
        "--min-orf"
        help = """DIAMOND ignores translated sequences with smaller open reading frames.
            Default is: disabled for sequences smaller than 30, 20 fro sequences smaller
            than 100, and 40 otherwise. Set to 1 to disable"""
        arg_type = Int64
        metavar = "NUMBER"
        "--fast"
        help = "Set DIAMOND sensitivity mode to 'fast'."
        action = :store_true
        "--mid-sensitive"
        help = "Set DIAMOND sensitivity mode to 'mid-sensitive'."
        action = :store_true
        "--sensitive"
        help = "Set DIAMOND sensitivity mode to 'sensitive'."
        action = :store_true
        "--more-sensitive"
        help = "Set DIAMOND sensitivity mode to 'more-sensitive'."
        action = :store_true
        "--very-sensitive"
        help = "Set DIAMOND sensitivity mode to 'very-sensitive'."
        action = :store_true
        "--ultra-sensitive"
        help = "Set DIAMOND sensitivity mode to 'ultra-sensitive'."
        action = :store_true
        # ATTENTION: error when combining iterate and frameshift (GitHub issue #593)
        "--iterate" # iterate with list requires diamond v. 2.0.12
        help = """Set DIAMOND option --iterate. In version 2.0.12 or higher, you can
            optionally specify a space-separated list of sensitivity modes to iterate over.
            Allowed values are 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive',
            'very-sensitive', 'ultra-sensitive', 'default' and none"""
        arg_type = String
        default = DMND_ITERATE
        metavar = "MODE"
        nargs = '*'
        "--frameshift"
        help = """Allow frameshift in DIAMOND and set frameshift penalty. Without this
            option, frameshift is disabled entirely"""
        arg_type = Int64
        metavar = "NUMBER"
        "--evalue"
        help = "Only report DIAMOND hits with lower e-values than the given value"
        arg_type = Float64
        metavar = "NUMBER"
        "--min-score"
        help = "Only report DIAMOND hits with bitscores >= the given value. Overrides
            the --evalue option"
        arg_type = Float64
        metavar = "NUMBER"
        "--max-target-seqs"
        help = """The maximum number of subject sequences that DIAMOND may report per
            query. Default is 25; setting it to 0 will report all hits"""
        arg_type = Int64
        metavar = "NUMBER"
        "--top"
        help = """Discard DIAMOND hits outside the given percentage range of the top
            alignment score. This option overrides --max-target-seqs"""
        arg_type = Int64
        metavar = "NUMBER"
        "--max-hsps"
        help = """Maximum number of HSPs DIAMOND may report per target sequence for each
            query. Default is reporting only the highest-scoring HSP. Setting this option
            to 0 will report all alternative HSPs"""
        arg_type = Int64
        metavar = "NUMBER"
        "--id"
        help = "Discard DIAMOND hits with less sequence identity than the given percentage"
        arg_type = Float64
        metavar = "PERCENTAGE"
        "--query-cover"
        help = "Discard DIAMOND hits with less query cover than the given percentage"
        arg_type = Float64
        metavar = "PERCENTAGE"
        "--subject-cover"
        help = "Discard DIAMOND hits with less subject cover than the given percentage"
        arg_type = Float64
        metavar = "PERCENTAGE"
        "--masking"
        help = """Set the DIAMOND mode for repeat masking. Note that, contrary to DIAMOND
            default (tantan masking enabled), Patchwork disables masking by default
            (default: 0)! Set to 1 to enable tantan masking, or to 2 to enable default
            BLASTP SEG masking. Note that the latter requires a DIAMOND version >= 2.0.12."""
        arg_type = Int64
        metavar = "MODE"
        default = 0 # TODO: default 0 or default 1? (which one works better for Patchwork?)
        "--block-size"
        help = "Billions of sequence letters to be processed at a time. A larger block size
            leads to increased performance at the expense of disk and memory usage. Values
            >20 are not recommended."
        arg_type = Float64
        metavar = "NUMBER"
        default = 2.0
    end # DIAMOND BLASTX

    add_arg_group!(settings, "alignment")
    @add_arg_table! settings begin
        "--len"
        help = "Discard DIAMOND hits shorter than the provided NUMBER"
        arg_type = Int64
        metavar = "NUMBER"
        "--gapopen"
        help = "Set the gap open penalty to this positive NUMBER"
        arg_type = Int64
        metavar = "NUMBER"
        "--gapextend"
        help = "Set the gap extension penalty to this positive NUMBER"
        arg_type = Int64
        metavar = "NUMBER"
        "--retain-stops"
        help = "Do not remove stop codons (`*`) in the output sequences"
        action = :store_true
        "--retain-ambiguous"
        help = "Do not remove ambiguous characters from the output sequences"
        action = :store_true
    end

    add_arg_group!(settings, "sliding window")
    @add_arg_table! settings begin
        "--no-trimming"
        help = "Skip sliding window-based trimming of alignments"
        action = :store_true
        "--window-size"
        help = "Specifices the NUMBER of positions to average across"
        arg_type = Int64
        default = 4
        metavar = "NUMBER"
        "--required-distance"
        help = "Specifies the average distance required"
        arg_type = Float64
        default = -7.0
        metavar = "NUMBER"
    end

    add_arg_group!(settings, "resources")
    @add_arg_table! settings begin
        "--threads"
        help = "Number of threads to utilize"
        default = Sys.CPU_THREADS
        arg_type = Int64
        metavar = "NUMBER"
    end

    return ArgParse.parse_args(settings)
end

function main()
    args = parse_parameters()

    if !min_diamondversion(MIN_DIAMONDVERSION)
        error("Patchwork requires DIAMOND with a version number above $MIN_DIAMONDVERSION")
    end
    printinfo(get_diamondversion(), args["threads"])

    setpatchworkflags!(args)

    if length(args["reference"]) == 1
        references_file = args["reference"][1]              # 1 .fa
    else                                                    # multiple .fa files
        references_file = mktemp_fasta(pool(args["reference"]; removeduplicates = false))
    end
    refseqs_count = countsequences(references_file)

    contigpaths = args["contigs"]
    contigsformat = sequencefiles_format(contigpaths)
    println("Sequence data file format detected: ", contigsformat)

    # Concatenate multiple sequence files if necessary.
    if length(contigpaths) > 1
        println("Concatenating query sequence files...")
        queries = cat(contigpaths)
    elseif length(contigpaths) == 1
        queries = only(contigpaths)
    end

    contigscount = countrecords(queries, contigsformat)

    outdir = args["output-dir"]
    alignmentoutput = outdir * "/" * ALIGNMENTOUTPUT
    fastaoutput = outdir * "/" * FASTAOUTPUT
    dnafastaoutput = outdir * "/" * DNAFASTAOUTPUT
    diamondoutput = outdir * "/" * DIAMONDOUTPUT
    statsoutput = outdir * "/" * STATSOUTPUT
    plotsoutput = outdir * "/" * PLOTSOUTPUT
    statistics = DataFrame(id = String[],
        reference_len = Int[],
        query_len = Int[],
        regions = Int[],
        contigs = Int[],
        matches = Int[],
        mismatches = Int[],
        deletions = Int[],
        query_coverage = Float64[],
        identity = Float64[])
    trimmedalignment_output = outdir * "/" * TRIMMEDALIGNMENT_OUTPUT

    if (isfile(alignmentoutput) || isdir(statsoutput) && !isempty(readdir(statsoutput))
        || (isdir(fastaoutput) && !isempty(readdir(fastaoutput)))
        || (isdir(dnafastaoutput) && !isempty(readdir(dnafastaoutput))))
        if !args["overwrite"]
            answer = warn_overwrite()
            isequal(answer, "n") && return
        end
        cleanfiles(alignmentoutput, statsoutput, fastaoutput, dnafastaoutput)
    end

    map(mkpath, [diamondoutput, fastaoutput, dnafastaoutput, statsoutput, plotsoutput])

    if isnothing(args["search-results"])
        if isempty(queries)
            println("Please provide one or more query files if not running with
                `--search-results` mode.")
            return
        end
        if isnothing(args["database"])
            println("Building DIAMOND database...")
            reference_db = diamond_makeblastdb(references_file, outdir)
        else
            reference_db = args["database"]
        end
        diamondparams = collectdiamondflags(args)
        println("Aligning query sequences against reference database...")
        diamondsearch = diamond_blastx(queries, reference_db, outdir, diamondparams)
    else
        diamondsearch = args["search-results"]
    end

    println("Merging overlapping hits...")
    allhits = readblastTSV(diamondsearch)
    referenceids = unique(subjectids(allhits))

    if !args["no-trimming"]
        println("Trimming alignments...")
    end

    hit = allhits[1]
    queryseqs_count = 0
    for (index, referenceid) in enumerate(referenceids)
        queryseqs_count += 1
        diamondhits = filter(hit -> isequal(referenceid, hit.subjectid), allhits)
        @assert length(diamondhits) != 0
        writeblastTSV(*(diamondoutput, "/", referenceid.id, ".tsv"), diamondhits; header = true)
        regions = AlignedRegionCollection(selectsequence(references_file, referenceid), diamondhits)
        # assuming all queries belong to same species:

        mergedregions = mergeoverlaps(regions)
        concatenation, dna_concatenation = concatenate(mergedregions, args["species-delimiter"])

        # Mask inserts
        maskedalignment, maskeddna = maskgaps(concatenation, dna_concatenation)

        # Mask stop codons and ambiguous characters
        maskedalignment, finaldna = maskalignment(maskedalignment.aln, maskeddna, DEFAULT_SCOREMODEL,
            args["retain-stops"], args["retain-ambiguous"])

        finalalignment = maskedalignment.aln
        write_alignmentfile(alignmentoutput, referenceid, length(mergedregions), finalalignment, index)

        # Alignment trimming
        if !args["no-trimming"]
            finalalignment, finaldna = slidingwindow(finalalignment, finaldna, args["window-size"],
                args["required-distance"], DEFAULT_SCOREMODEL)
            write_alignmentfile(trimmedalignment_output, referenceid, length(mergedregions), finalalignment, index)
        end

        write_fasta(
            *(fastaoutput, "/", sequencepart(referenceid), args["fasta-extension"]),
            regions.records[1].queryid,
            finalalignment.a.seq
        )
        write_fasta(
            *(dnafastaoutput, "/", sequencepart(referenceid), args["fasta-extension"]),
            regions.records[1].queryid,
            finaldna
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
            round(occupancy(finalalignment) * 100., digits = 2),
            gapexcluded_identity(BioAlignments.count_matches(finalalignment),
                                 BioAlignments.count_mismatches(finalalignment))
        ]
        push!(statistics, stats_row)
    end

    println(RULER)
    percentmarkers = getpercentage(queryseqs_count, refseqs_count)
    println("# query sequences in: ", contigscount)
    println("# reference sequences in: ", refseqs_count)
    println("# alignments out: ", queryseqs_count, " (", percentmarkers, "%)")

    CSV.write(*(statsoutput, "/statistics.csv"), statistics, delim = ",")
    statssummary = select(describe(select(statistics, Not(:id))), Not([:nmissing, :eltype]))
    pretty_table(statssummary, nosubheader = true)

    if !args["no-plots"]
        plot_querycover(statistics.query_coverage, *(plotsoutput, "/query_coverage.png"))
        plot_percentident(statistics.identity, refseqs_count, *(plotsoutput, "/percent_identity.png"))
    end

    CSV.write(*(statsoutput, "/average.csv"), statssummary, delim = ",")
end

function julia_main()::Cint
    try
        @time main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end

    return 0
end

if length(ARGS) >= 1
    julia_main()
end

end # module
