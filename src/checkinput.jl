using BioAlignments
using BioSymbols

#include("diamond.jl")

const MATRICES = Dict("BLOSUM45"=>BLOSUM45, "BLOSUM50"=>BLOSUM50,"BLOSUM62"=>BLOSUM62,
                      "BLOSUM80"=>BLOSUM80, "BLOSUM90"=>BLOSUM90, "PAM30"=>PAM30,
                      "PAM70"=>PAM70, "PAM250"=>PAM250)
const GAPS = Dict("BLOSUM45"=>[(10, 3), (11, 3), (12, 3), (13, 3), (12, 2), (13, 2),
                               (14, 2), (15, 2), (16, 2), (16, 1), (17, 1), (18, 1),
                               (19, 1)],
                  "BLOSUM50"=>[(9, 3),(10, 3), (11, 3), (12, 3), (13, 3), (12, 2), (13, 2),
                              (14, 2), (15, 2), (16, 2), (15, 1), (16, 1), (17, 1),
                              (18, 1), (19, 1)],
                  "BLOSUM62"=>[(6, 2), (7, 2), (8, 2), (9, 2), (10, 2), (11, 1), (9, 1),
                               (10, 1), (11, 1), (12, 1), (13, 1)],
                  "BLOSUM80"=>[(6, 2), (7, 2), (8, 2), (9, 2), (9, 1), (10, 1), (11, 1),
                               (13, 2), (25, 2)],
                  "BLOSUM90"=>[(6, 2), (7, 2), (8, 2), (9, 2), (9, 1), (10, 1), (11, 1)],
                  "PAM30"=>[(5, 2), (6, 2), (7, 2), (8, 1), (9, 1), (10, 1)],
                  "PAM70"=>[(6, 2), (7, 2), (8, 2), (9, 1), (10, 1), (11, 1)],
                  "PAM250"=>[(11, 3), (12, 3), (13, 3), (14, 3), (15, 3), (13, 2), (14, 2),
                             (15, 2), (16, 2), (17, 2), (17, 1), (18, 1), (19, 1), (20, 1),
                             (21, 1)])
const GAPDEFAULTS = Dict("BLOSUM45"=>(14, 2), "BLOSUM50"=>(13, 2), "BLOSUM62"=>(11, 1),
                         "BLOSUM80"=>(10, 1), "BLOSUM90"=>(10, 1), "PAM30"=>(9, 1),
                         "PAM70"=>(10, 1), "PAM250"=>(14, 2))
const AVAIL_DIAMOND = ["iterate", "frameshift", "evalue", "min-score", #"sensitivity",
    "max-target-seqs", "top", "max-hsps", "id", "query-cover", "subject-cover", 
    "query-gencode", "strand", "min-orf", "masking"]
const DIAMONDSTRANDS = ["both", "plus", "minus"]
const DIAMONDMODES = ["fast", "mid-sensitive", "sensitive", "more-sensitive", 
    "very-sensitive", "ultra-sensitive"]
const DMND_ITERATE = ["PATCHWORK_OFF"]
const DIAMONDMASK = ["0", "1", "seg"]

# FUNCTIONS ###############################################################################

# Matrix ##################################################################################

"""
    checkmatrixtype(args::Dict{String, Any})

Ensure that at most one of `--matrix` and `--custom-matrix` was provided as command line
argument. The name supplied to the argument `--matrix` must correspond to one of
`BioAlignements`'s pre-defined matrices (BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90,
PAM30, PAM70, or PAM250). In case a custom matrix is used, check that `--lambda` and `--K`
were set in `--diamond-flags`.
"""
function checkmatrixtype(args::Dict{String, Any})
    if !isnothing(args["matrix"]) && !isnothing(args["custom-matrix"])
        error("You can only supply one of --matrix and --custom-matrix")
    elseif !isnothing(args["custom-matrix"]) && (!("--lambda" in args["diamond-flags"])
                                                 || !("--K" in args["diamond-flags"]))
        error("DIAMOND custom scoring matrices require setting the --lambda and",
              " --K options.")
    elseif !isnothing(args["matrix"]) && !(args["matrix"] in keys(MATRICES))
        error("Invalid scoring matrix: $matrix.")
    end
end

"""
    getmatrixtype(args::Dict{String, Any})::String

Get the type of matrix the user supplied, i.e. `matrix` or `custom-matrix`.
"""
function getmatrixtype(args::Dict{String, Any})::String
    !isnothing(args["custom-matrix"]) && return "custom-matrix"
    return "matrix"
end

"""
    setmatrixname!(args::Dict{String, Any})

If neither `--matrix` nor `--custom-matrix` were supplied, set the substitution matrix to
the default `MATRIX`. For matrices provided with `--matrix`, convert the name to capitals
for internal compatibility with other functions using the same `args`.
"""
function setmatrixname!(args::Dict{String, Any})
    matrixtype = getmatrixtype(args)
    if !isequal("custom-matrix", matrixtype)
        if isnothing(args[matrixtype])
            args[matrixtype] = MATRIX
        else
            args[matrixtype] = uppercase(args[matrixtype])
        end
    end
end

"""
    read_custommatrix(path::AbstractString)::SubstitutionMatrix

Parse the file provided in `--custom-matrix` as an amino acid substitution matrix,
creating a new `BioAlignments.SubstitutionMatrix` object from the data.
"""
function read_custommatrix(path::AbstractString)::SubstitutionMatrix
    return BioAlignments.parse_ncbi_submat(BioSymbols.AminoAcid, path)
end

# matrixname = args[getmatrixtype(args)] and matrixtype = getmatrixtype(args)
"""
    getmatrix(matrixname::AbstractString, matrixtype::AbstractString)::SubstitutionMatrix

Get the `BioAlignments.SubstitutionMatrix` object corresponding to the substitution matrix.
If `--custom-matrix` was set, this function will parse the provided file as an amino acid
substitution matrix and create a new `BioAlignments.SubstitutionMatrix` object. Otherwise,
the matrix must correspond to one of `BioAlignments`'s pre-defined matrices (BLOSUM45,
BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, PAM30, PAM70, or PAM250).
"""
function getmatrix(
    matrixname::AbstractString,
    matrixtype::AbstractString
)::SubstitutionMatrix
    isequal("custom-matrix", matrixtype) && return read_custommatrix(matrixname)
    return MATRICES[matrixname]
end

# Gap Penalties ###########################################################################

"""
    setgapopen!(args::Dict{String, Any})

Set the gap-opening penalty to the default value, depending on the matrix provided in
`--matrix`, if no gap-opening penalty was explicitely provided via `--gapopen`.
If necessary, convert the penalty to a positive integer. This function throws an error if
`--custom-matrix` was set without providing a gap-opening penalty.
"""
function setgapopen!(args::Dict{String, Any})
    if !isnothing(args["gapopen"])
        args["gapopen"] = abs(args["gapopen"])
        return
    end
    matrixtype = getmatrixtype(args)
    if isequal("matrix", matrixtype)
        args["gapopen"] = GAPDEFAULTS[args[matrixtype]][1]
    else
        error("Custom scoring matrices require setting the --gapopen option.")
    end
end

"""
    setgapextend!(args::Dict{String, Any})

Set the gap-extension penalty to the default value, depending on the matrix provided in
`--matrix`, if no gap-extension penalty was explicitely provided via `--gapextend`.
If necessary, convert the penalty to a positive integer. This function throws an error if
`--custom-matrix` was set without providing a gap-extension penalty.
"""
function setgapextend!(args::Dict{String, Any})
    if !isnothing(args["gapextend"])
        args["gapextend"] = abs(args["gapextend"])
        return
    end
    matrixtype = getmatrixtype(args)
    if isequal("matrix", matrixtype)
        args["gapextend"] = GAPDEFAULTS[args[matrixtype]][2]
    else
        error("Custom scoring matrices require setting the --gapextend option.")
    end
end

"""
    checkgappenalty(matrix::AbstractString, gapopen::Int64, gapextend::Int64)

Ensure that the gap penalties are compatible with the matrix.
"""
function checkgappenalty(
    matrix::AbstractString,
    gapopen::Int64,
    gapextend::Int64
)
    if matrix in keys(MATRICES)
        allowedgaps = GAPS[matrix]
        if isnothing(findfirst(isequal((gapopen, gapextend)), allowedgaps))
            error("Invalid gap open and/or gap extend scores.")
        end
    end
end

# Diamond Flags ###########################################################################
#
#"""
#    set_diamondframeshift!(args::Dict{String, Any})
#
#If no `--frameshift` argument was provided in `--diamond-flags`, set the `--frameshift`
#to the default value.
#"""
#function set_diamondframeshift!(args::Dict{String, Any})
#    if !("--frameshift" in args["diamond-flags"])
#        push!(args["diamond-flags"], "--frameshift", FRAMESHIFT)
#    end
#end

#"""
#    set_diamondmode!(args::Dict{String, Any})
#
#If none of `--fast`, `--mid-sensitive`, `--sensitive`, `--more-sensitive`, `very-sensitive`
#or `--ultra-sensitive` was set in `--diamond-flags`, set the mode to the default value.
#"""
#function set_diamondmode!(args::Dict{String, Any})
#    if !("--fast" in args["diamond-flags"]
#         || "--mid-sensitive" in args["diamond-flags"]
#         || "--sensitive" in args["diamond-flags"]
#         || "--more-sensitive" in args["diamond-flags"]
#         || "--very-sensitive" in args["diamond-flags"]
#         || "--ultra-sensitive" in args["diamond-flags"])
#        push!(args["diamond-flags"], DIAMONDMODE)
#    end
#end

# """
#     checkdiamondflags(args::Dict{String, Any})

# Ensure that no conflicts are caused by command line arguments provided in `--diamond-flags`,
# i.e. that parameters used and set by `Patchwork` are not also provided to `DIAMOND` in a
# separate place.
# """
# function checkdiamondflags(args::Dict{String, Any})
#     flags = args["diamond-flags"]
#     #isempty(flags) && return
#     patchworkflags = ["matrix", "custom-matrix", "gapopen", "gapextend", "threads"]
#     for flag in patchworkflags
#         if flag in flags
#             error("Please use Patchwork's own --$flag option instead.")
#         end
#     end
#     if "--query" in flags || "-q" in flags
#         error("You should not explicitly provide a DIAMOND query.",
#               "The file provided to Patchwork as --contigs will be used as DIAMOND query.")
#     end
#     if "--db" in flags || "-d" in flags
#         error("You should not explicitly provide a DIAMOND database.",
#               "The file or database provided with --reference will be used by DIAMOND.",
#               "If you are providing a database, please make sure to also set --database.")
#     end
#     if "--outfmt" in flags || "-f" in flags
#         error("The DIAMOND output format is built into Patchwork and cannot be changed.")
#     end
#     if "--out" in flags || "-o" in flags
#         error("The DIAMOND results file is built into Patchwork and cannot be changed.")
#     end
# end

# function checkmakedbflags(args::Dict{String, Any})
#     flags = args["makedb-flags"]
#     isempty(flags) && return
#     if "--in" in flags
#         error("Please use Patchwork's --reference option instead. ",
#               "It willl be used to construct your DIAMOND database.")
#     end
#     if "--db" in flags || "-d" in flags
#         error("The creation of your DIAMOND database will be internally handled. ",
#               "You should not explicitly provide an output database filename.")
#     end
#     if "--threads"  in flags
#         error("Please use Patchwork's own --threads option instead.")
#     end
# end

# Check and Collect all DIAMOND Flags in an Array #########################################

"""
    setpatchworkflags!(args::Dict{String, Any})

Ensure that no conflicts are caused by combinations of, or missing, command line arguments.
Set gap-opening and gap-extension penalties correctly, modifying `args`. This function
throws an error if any conflicts are detected.
"""
function setpatchworkflags!(args::Dict{String, Any})     # Run this fct. before the next.
    checkmatrixtype(args)
    setmatrixname!(args)
    setgapopen!(args)
    setgapextend!(args)
    checkgappenalty(args[getmatrixtype(args)], args["gapopen"], args["gapextend"])
end

# """
#     setdiamondflags!(args::Dict{String, Any})

# Ensure that no conflicts are caused by combinations of, or missing, command line arguments.
# Set `DIAMOND` frameshift and mode correctly, modifying `args`. This function throws an
# error if any conflicts are detected.
# """
# function setdiamondflags!(args::Dict{String, Any})      # Run this fct. before the next.
#     checkdiamondflags(args)
#     #set_diamondframeshift!(args)
#     #set_diamondmode!(args)
#     checkmakedbflags(args)
#     if length(args["reference"]) > 1 || !isdiamonddatabase(args["reference"]...)
#         push!(args["makedb-flags"], "--threads", string(args["threads"]))
#         #isfastafile(args["reference"]) && push!(args["makedb-flags"], "-d",
#         #            args["output-dir"] * DATABASE)
#     end
# end

# """
#     collectdiamondflags(args::Dict{String, Any})::Vector{String}

# Build the complete vector of parameters for running `DIAMOND`.
# """
# function collectdiamondflags(args::Dict{String, Any})::Vector{String}
#     @assert !isnothing(args["gapopen"]) && !isnothing(args["gapextend"]) """set gap
#     penalties with setpatchworkflags!(args) before calling this function"""
#     return [args["diamond-flags"]..., "--" * getmatrixtype(args), args[getmatrixtype(args)],
#             "--gapopen", string(args["gapopen"]), "--gapextend", string(args["gapextend"]),
#             "--threads", string(args["threads"])]
# end

function collectdiamondflags(args::Dict{String, Any})::Vector{String}
    #setpatchworkflags!(args)
    @assert !isnothing(args["gapopen"]) && !isnothing(args["gapextend"]) """set gap
        penalties with setpatchworkflags!(args) before calling this function"""
    diamondflags = ["--" * getmatrixtype(args), args[getmatrixtype(args)],
        "--gapopen", string(args["gapopen"]), "--gapextend", string(args["gapextend"]),
        "--threads", string(args["threads"])]

    if !isnothing(args["strand"])
        !in(args["strand"], DIAMONDSTRANDS) && error("Only values 'both', 'plus', and 'minus' 
            are allowed for the --strand option.")
    end
    # if !isnothing(args["sensitivity"])
    #     !in(args["sensitivity"], DIAMONDMODES) && error("""Only values 'fast', 'mid-sensitive', 
    #         'sensitive', 'more-sensitive', 'very-sensitive', and 'ultra-sensitive' are 
    #         allowed for the --sensitivity option.""")
    #     push!(diamondflags, "--" * args["sensitivity"])
    # end
    if !isequal(args["iterate"], DMND_ITERATE) 
        if !isempty(args["iterate"])
            !min_diamondversion("2.0.12") && error("""DIAMOND version must be at least 2.0.12
            to support providing a list of sensitivities with the --iterate option.""")
            # if in("default", args["iterate"]) && length(args["iterate"]) > 1
            #     error("'default' may not be used in combination with other values when setting 
            #         --iterate.")
            # end
            # if !isequal(args["iterate"][1], "default")
            for val in args["iterate"]
                (!in(val, DIAMONDMODES) && !isequal(val, "default")) && error("""The only 
                    values allowed for the --iterate option are 'default', a combination of 
                    the values 'fast', 'mid-sensitive', 'sensitive', 'more-sensitive', 
                    'very-sensitive', and 'ultra-sensitive', provided as a space-separated 
                    list, or nothing at all.""")
            end
            # end
        end
        push!(diamondflags, "--iterate", args["iterate"]...)
    end
    if args["masking"] < 0 || args["masking"] > 2
        error("Only values '0', '1', and '2' are allowed for the --masking option.")
    end
    push!(diamondflags, "--masking", DIAMONDMASK[args["masking"] + 1]) # 1-based indexing

    sens = false
    for opt in DIAMONDMODES
        if args[opt] 
            if !sens
                push!(diamondflags, "--" * opt)
                sens = true
            else
                error("Please specify <= 1 sensitivity mode for DIAMOND.")
            end
        end
    end
    for opt in AVAIL_DIAMOND
        # in(opt, ["sensitivity", "iterate", "masking"]) && continue
        # if !isnothing(args[opt])
        #     push!(diamondflags, "--" * opt, string(args[opt]))
        # end
        in(opt, ["iterate", "masking"]) && continue
        if !isnothing(args[opt])
            push!(diamondflags, "--" * opt, string(args[opt]))
        end
    end
    # println(diamondflags)
    return diamondflags
end
