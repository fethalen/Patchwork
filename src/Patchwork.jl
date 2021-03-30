module Patchwork

using ArgParse
using Crayons
using DataFrames

include("alignedregion.jl")
include("alignedregioncollection.jl")
include("blast.jl")
include("fasta.jl")
include("sequencerecord.jl")
include("multiplesequencealignment.jl")

"""
    print_logo()

Prints the Patchwork logo to stdout.
"""
function print_logo()
    r = string(Crayons.Crayon(foreground = :red))
    g = string(Crayons.Crayon(foreground = :green))
    y = string(Crayons.Crayon(foreground = :yellow))
    b = string(Crayons.Crayon(foreground = :blue))
    c = string(Crayons.Crayon(foreground = :cyan))
    e = string(Crayon(foreground = :dark_gray))
    res = string(Crayons.Crayon(reset = true))
    und = string(Crayons.Crayon(underline = true))
    bol = string(Crayons.Crayon(bold = true))

    str = """
    $(r)██$(e)██$(y)██$(g)██$(res)
    $(e)██$(c)██$(g)██$(r)██$(res) $(bol)Patchwork$(res)
    $(y)██$(g)██$(r)██$(c)██$(res)
    """

    print(str)
end

"""
    commandexists(command)

Returns `true` if the provided command exists within the current path and throw an error
otherwise.
"""
function commandexists(command::AbstractString)
    cmdexists = `which $command`
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
        "--query-alignments"
            help = "Path to one or more query alignments in FASTA format"
            required = true
            metavar = "PATH"
        "--subject"
            help = "Path to one or more subject sequences or a subject database"
            required = true
            metavar = "PATH"
        "--output"
            help = "Store outputs in this folder"
            default = "patchworks_output"
            metavar = "PATH"
        "--blast-engine"
            help = "Which program to use for performing the BLAST search (diamond/blast; 
                    default: try diamond, fall back to blast)"
            default = "diamond"
            metavar = "PATH"
        "--database"
            help = "When true, subject holds a path to a BLAST database"
            action = :store_true
        "--extensions"
            help = "Filetype extensions used to detect alignment FASTA files (default:
                    [\"aln\", \"fa\", \"fn\", \"fna\", \"faa\", \"fasta\", \"FASTA\"])"
            arg_type = Array{String, 1}
            default = ["aln", "fa", "fn", "fna", "faa", "fasta", "FASTA"]
            metavar = "LIST"
        "--seq-type"
            help = "Type of input alignments (default: auto-detect; valid types:
                    nucleotide, aminoacid, and auto)"
            default = "auto"
            arg_type = String
            metavar = "NAME"
        "--translate"
            help = "Translate nucleotide sequence into amino acid sequences"
            action = :store_true
        "--reverse-translate"
            help = "Reverse translate amino acid sequences into nucleotide sequences"
            action = :store_true
        "--species-delimiter"
            help = "Used to separate species names from sequence identifiers in FASTA records
                    (default: species@identifier)"
            default = "@"
            arg_type = Char
            metavar = "CHARACTER"
        "--threads"
            help = "Number of threads that will be utilized (default: all available)"
            default = Sys.CPU_THREADS
            arg_type = Int
        "--wrap-column"
            help = "Wrap output sequences at the provided column number"
            arg_type = UInt
            metavar = "NUMBER"
    end

    return ArgParse.parse_args(settings)
end

function main()
    # args = parse_parameters()
    print_logo()
    if ! commandexists("mafft")
        error("Cannot find command \'mafft\' in current path")
    end
    if commandexists("diamond")
        blastengine = "diamond"
    elseif commandexists("blastx")
        blastengine = "blastx"
    else
        error("\'diamond\' or \'blastx\' must be in the current path to run Patchwork")
    end
    println("BLAST-engine: $blastengine")

    query = "/media/feli/Storage/phylogenomics/1st_wo_nextera/ceratonereis_australis/spades_assembly/K125/Ceratonereis_australis_k125_spades_assembly/final_contigs.fasta"
    subject = "/home/feli/ownCloud/projects/Patchwork/test/07673_Alitta_succinea.fa"
    subject_db = diamond_makeblastdb(subject, ["--threads", Sys.CPU_THREADS])
    hits = diamond_blastx(query, subject_db, ["--threads", Sys.CPU_THREADS])
    queryids = convert(Array, select(DataFrames.DataFrame(hits), :queryid))
    querymsa = selectsequences(query, queryids)

    # merge!(queryalignment, results)
    # querysubject_aln = mafft_linsi(queryalignment, ["--thread", Sys.CPU_THREADS])
    # regions = AlignedRegionCollection(querymsa,hits)
    # # uniqueregions = uniquesequences(regions)

    # # alittahits = filter(result -> result.queryotu == "ASUC", results)
    # merge!(queryalignment, alittahits)
    # querysubject_aln = mafft_linsi(queryalignment, ["--thread", Sys.CPU_THREADS])
    # regions = AlignedRegionCollection(queryalignment, hits)
end

end # module
