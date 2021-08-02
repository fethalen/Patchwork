module Patchwork

using Base: Bool, Int64, func_for_method_checked, DEFAULT_COMPILER_OPTS
using ArgParse
using DataFrames

include("alignment.jl")
include("alignedregion.jl")
include("alignedregioncollection.jl")
include("alignmentconcatenation.jl")
include("diamond.jl")
include("fasta.jl")
include("sequencerecord.jl")
include("multiplesequencealignment.jl")

const FASTAEXTENSIONS = ["aln", "fa", "fn", "fna", "faa", "fasta", "FASTA"]
const MAKEBLASTDB_FLAGS = ["--threads", Sys.CPU_THREADS]
const DIAMONDFLAGS = ["--evalue", 0.001, "--frameshift", 15, "--threads",
                      "--ultra-sensitive", Sys.CPU_THREADS]
const MIN_DIAMONDVERSION = "2.0.3"
const MATRIX = "BLOSUM62"
const GAPOPEN = 11
const GAPEXTEND = 1

"""
    printinfo()

Print basic information about this program.
"""
function printinfo()
    about = """
    P A T C H W O R K
    Developed by: Felix Thalen and Clara Köhne
    Dept. for Animal Evolution and Biodiversity, University of Göttingen

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
    diamondversion = last(split(versioncmd))
    return diamondversion >= minversion
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
            help = "Either (1) a path to one or more sequences in FASTA format or (2) a
                    subject database (set --database)"
            required = true
            arg_type = String
            metavar = "PATH"
        # TODO: Doesn't have to be a flag, check filetype extension instead
        "--database"
            help = "When specified, \"--reference\" points to a DIAMOND/BLAST database"
            arg_type = Bool
            action = :store_true
        "--output-dir"
            help = "Write output files to this directory"
            arg_type = String
            default = "patchwork_output"
            metavar = "PATH"
        "--diamond-flags"
            help = "Flags sent to DIAMOND"
            arg_type = Vector
            default = DIAMONDFLAGS
            metavar = "LIST"
        "--matrix"
            help = "Set scoring matrix"
            arg_type = String
            default = "BLOSUM62"
            metavar = "NAME"
        "--custom-matrix"
            help = "Use a custom scoring matrix"
            arg_type = String
            default = "BLOSUM62"
            metavar = "PATH"
        "--gapopen"
            help = "Set gap open penalty (positive integer)"
            arg_type = Int64
            default = 11
            metavar = "NUMBER"
        "--gapextend"
            help = "Set gap extension penalty (positive integer)"
            arg_type = Int64
            default = 1
            metavar = "NUMBER"
        "--seq-type"
            help = "Type of input alignments (nucleotide/aminoacid; default: autodetect)"
            default = "autodetect"
            arg_type = String
            metavar = "NAME"
        "--species-delimiter"
            help = "Used to separate the species name from the rest of the identifier in
                    FASTA records (default: @)"
            default = '@'
            arg_type = Char
            metavar = "CHARACTER"
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

    subject = "test/07673_Alitta_succinea.fa"
    query = "/media/feli/Storage/nereidid_data/2020-01-02_allgenetics/ceratonereis_australis/spades_assembly/K125/Ceratonereis_australis_k125_spades_assembly/final_contigs.fasta"
    subject_db = Patchwork.diamond_makeblastdb(subject, MAKEBLASTDB_FLAGS)
    # diamondresults = Patchwork.diamond_blastx(query, subject_db, DIAMONDFLAGS)
    diamondresults = "test/c_australis_x_07673.tsv"
    hits = Patchwork.readblastTSV(diamondresults)
    # querymsa = Patchwork.selectsequences(query, Patchwork.queryids(hits, speciesdelimiter))
    # referenceseq = Patchwork.readmsa(subject)
    # The full sequences are required for evaluating the subject <-> query methods
    full_subjectseq = Patchwork.get_fullseq(subject)
    regions = Patchwork.AlignedRegionCollection(full_subjectseq, hits)
    while Patchwork.hasoverlaps(regions)
        regions = Patchwork.mergeoverlapping(regions)
    end
    # a = regions[12]
    # b = regions[13]
    # aqueryseq = Patchwork.selectsequence(query, a.queryid.id)
    # bqueryseq = Patchwork.selectsequence(query, b.queryid.id)
end

end # module
