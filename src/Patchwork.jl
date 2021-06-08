module Patchwork

using Base: Bool
using ArgParse
using DataFrames

include("alignedregion.jl")
include("alignedregioncollection.jl")
include("diamond.jl")
include("fasta.jl")
include("sequencerecord.jl")
include("multiplesequencealignment.jl")

"""
    commandexists(command)

Returns `true` if the provided command exists within the current path and throw an error
otherwise.
"""
function commandexists(command::AbstractString)
    cmdexists = `command -v $command`
    try
        run(cmdexists)
    catch
        return false
    end
    return true
end

function parse_parameters()
    overview = """
    Concatenates two or more sequence alignments into a single supermatrix
    and provides detailed information about that supermatrix.
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
        "--database"
            help = "When specified, \"--reference\" points to a DIAMOND/BLAST database"
            arg_type = Bool
            action = :store_true
        "--output-dir"
            help = "Write output files to this directory"
            arg_type = String
            default = "patchwork_output"
            metavar = "PATH"
        "--fasta-extensions"
            help = "Filetype extensions used to detect alignment FASTA files (default:
                    [\"aln\", \"fa\", \"fn\", \"fna\", \"faa\", \"fasta\", \"FASTA\"])"
            arg_type = Array{String, 1}
            default = ["aln", "fa", "fn", "fna", "faa", "fasta", "FASTA"]
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
            help = "Set gap open penalty (non-negative integer)"
            arg_type = UInt
            default = 11
            metavar = "NUMBER"
        "--gapextend"
            help = "Set gap extension penalty (non-negative integer)"
            arg_type = UInt
            default = 1
            metavar = "NUMBER"
        "--mid-sensitive"
            help = "Run DIAMOND in mid-sensitive mode"
            arg_type = Bool
            action = :store_true
        "--sensitive"
            help = "Run DIAMOND in sensitive mode"
            arg_type = Bool
            action = :store_true
        "--more-sensitive"
            help = "Run DIAMOND in more sensitive mode"
            arg_type = Bool
            action = :store_true
        "--very-sensitive"
            help = "Run DIAMOND in very sensitive mode"
            arg_type = Bool
            action = :store_true
        "--ultra-sensitive"
            help = "Run DIAMOND in ultra-sensitive mode"
            arg_type = Bool
            action = :store_true
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
            arg_type = UInt
            metavar = "NUMBER"
        "--wrap-column"
            help = "Wrap output sequences at this column number (default: no wrap)"
            arg_type = UInt
            metavar = "NUMBER"
    end

    return ArgParse.parse_args(settings)
end

function main()
    # args = parse_parameters()
    println("Patchwork")
    if commandexists("blastx")
        blastengine = "blastx"
    else
        error("\'diamond\' not found in current path and is required")
    end

    matrix = "BLOSUM62"
    gapopen = 11
    gapextend = 1
    speciesdelimiter = '@'
    subject = "test/07673_Alitta_succinea.fa"
    query = "/media/feli/Storage/phylogenomics/1st_wo_nextera/ceratonereis_australis/spades_assembly/K125/Ceratonereis_australis_k125_spades_assembly/final_contigs.fasta"
    # subject_db = Patchwork.diamond_makeblastdb(subject, ["--threads", Sys.CPU_THREADS])
    diamondresults = Patchwork.diamond_blastx(query, subject_db, ["--threads", Sys.CPU_THREADS])
    # blastresults = Patchwork.readblastTSV(diamondresults)
    hits = Patchwork.readblastTSV("test/c_australis_x_07673.tsv")
    # querymsa = Patchwork.selectsequences(query, Patchwork.queryids(hits, speciesdelimiter))
    referenceseq = Patchwork.readmsa(subject, speciesdelimiter)
    regions = Patchwork.AlignedRegionCollection(hits)
    uniqueregions = Patchwork.uniquesequences(regions)
end

end # module
