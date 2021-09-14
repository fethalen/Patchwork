# julia --trace-compile=precompiled.jl Patchwork.jl --contigs "../test/07673_lcal.fa" --reference "../test/07673_Alitta_succinea.fa" --diamond-flags "--frameshift 15 --ultra-sensitive" --output-dir "../test/patchwork-output"
# diamond blastx --query 07673_dna.fa --db 07673_Alitta_succinea.fa --outfmt 6 qseqid qseq full_qseq qstart qend qframe sseqid sseq sstart send cigar pident bitscore --out diamond_results.tsv --frameshift 15

module Patchwork

# ERROR occurred in bioconda build test; proposed solution was:
import Pkg
Pkg.add("ArgParse")
##############################################################

using Base: Bool, Int64, func_for_method_checked, DEFAULT_COMPILER_OPTS, Cint
using ArgParse
using BioAlignments
using CSV
using FASTX
using DataFrames

include("alignment.jl")
include("alignedregion.jl")
include("alignedregioncollection.jl")
include("alignmentconcatenation.jl")
include("checkinput.jl")
include("diamond.jl")
include("fasta.jl")
include("multiplesequencealignment.jl")
include("output.jl")
include("sequencerecord.jl")

const FASTAEXTENSIONS = ["aln", "fa", "fn", "fna", "faa", "fasta", "FASTA"]
const DIAMONDDB = "dmnd"
const EMPTY = String[]
# --evalue defaults to 0.001 in DIAMOND
# --threads defaults to autodetect in DIAMOND
const DIAMONDFLAGS = ["--ultra-sensitive"]
# const FRAMESHIFT = 15
# const DIAMONDMODE = "--ultra-sensitive"
const MIN_DIAMONDVERSION = "2.0.3"
const MATRIX = "BLOSUM62"

const ALIGNMENTOUTPUT = "alignments.txt"
const FASTAOUTPUT = "queries_out.fa"
const DIAMONDOUTPUT = "diamond_results.tsv"
const DATABASE = "database.dmnd"
const STATSOUTPUT = "stats.csv"

"""
    printinfo()

Print basic information about this program.
"""
function printinfo()
    about = """
    P A T C H W O R K
    Developed by: Felix Thalén and Clara G. Köhne
    © Dept. for Animal Evolution and Biodiversity, University of Göttingen

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
function min_diamondversion(minversion::AbstractString)
    try
        run(`diamond --version`)
    catch
        return false
    end
    versioncmd = read(`diamond --version`, String)
    diamondversion_vector = split(last(split(versioncmd)), ".")
    minversion_vector = split(minversion, ".")
    for (v, r) in zip(diamondversion_vector, minversion_vector)
        version = parse(Int64, v)
        required = parse(Int64, r)
        version < required && return false
        version > required && return true
    end
    return true
end

function ArgParse.parse_item(::Type{T}, argument::AbstractString) where T <: AbstractVector
    convert(T, split(argument, " "))
end

function parse_parameters()
    overview = """
    Alignment-based Exon Retrieval and Concatenation with Phylogenomic
    Applications
    """
    settings = ArgParseSettings(description=overview,
                                version = "0.1.0",
                                add_version = true)
    @add_arg_table! settings begin
        "--contigs"
            help = "Path to one or more sequences in FASTA format"
            required = true
            arg_type = String
            metavar = "PATH"
        "--reference"
            help = "Either (1) a path to one or more sequences in FASTA format, (2) a 
                    subject database (DIAMOND or BLAST database), or (3) a DIAMOND output 
                    file in tabular format."
            required = true
            arg_type = String
            metavar = "PATH"
        #"--database"
        #    help = "When specified, \"--reference\" points to a DIAMOND/BLAST database"
        #    arg_type = Bool
        #    action = :store_true
        #"--tabular"
        #    help = "When specified, \"--reference\" points to a tabular DIAMOND output 
        #            file generated in a previous Patchwork run."
        #    arg_type = Bool
        #    action = :store_true
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
        "--threads"
            help = "Number of threads to utilize (default: all available)"
            default = Sys.CPU_THREADS
            arg_type = Int64
            metavar = "NUMBER"
        "--wrap-column"
            help = "Wrap output sequences at this column number (default: no wrap)"
            default = 0
            arg_type = Int64
            metavar = "NUMBER"
    end

    return ArgParse.parse_args(settings)
end

function main()
    args = parse_parameters()
    printinfo()
    if !min_diamondversion(MIN_DIAMONDVERSION)
        error("Patchwork requires \'diamond\' with a version number above
               $MIN_DIAMONDVERSION to run")
    end

    setpatchworkflags!(args)
    setdiamondflags!(args)
    reference = args["reference"]
    query = args["contigs"]
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

    mkpath(outdir)
    if isfile(alignmentoutput) || isfile(fastaoutput)   # !isempty(readdir(outdir))
        answer = warn_overwrite()
        isequal(answer, "n") && return
        cleanfiles(alignmentoutput, fastaoutput)        # if answer == 'y'
    end
    #if !args["tabular"]
    println("Creating DIAMOND database...")
    reference_db = diamond_makeblastdb(reference, outdir, args["makedb-flags"])
    diamondparams = collectdiamondflags(args)
    index = 1 # dummy count for working with only one query file
    # in case of multiple query files: pool first? else process each file separately: 
    #for (index, query) in enumerate(queries)
    println("Performing DIAMOND BLASTX search...")
    diamondsearch = diamond_blastx(query, reference_db, outdir, diamondparams)
    println("DIAMOND BLASTX search done.")
    println("Doing Patchwork Magic...")
    diamondhits = readblastTSV(diamondsearch)
    writeblastTSV(diamondoutput, diamondhits; header = true)
    regions = AlignedRegionCollection(get_fullseq(reference), diamondhits)
    referencename = regions.referencesequence.id
    mergedregions = mergeoverlaps(regions)
    concatenation = concatenate(mergedregions)
    finalalignment = maskgaps(concatenation).aln
    println("Patchwork Magic done.")
    println("Saving output...")
    write_alignmentfile(alignmentoutput, referencename, length(regions), finalalignment, index)
    # only one query species allowed in regions!: 
    write_fasta(fastaoutput, regions.records[1].queryid, finalalignment)
    stats_row = [mergedregions.referencesequence.id.id,
                 length(mergedregions.referencesequence),
                 length(finalalignment.a.seq),
                 length(mergedregions),
                 length(unique(map(region -> region.queryid.id, mergedregions))),
                 BioAlignments.count_matches(finalalignment),
                 BioAlignments.count_mismatches(finalalignment),
                 BioAlignments.count_deletions(finalalignment),
                 round(occupancy(finalalignment), digits=2)]
    push!(statistics, stats_row)
    CSV.write(statsoutput, statistics, delim = "\t")
    #end
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
