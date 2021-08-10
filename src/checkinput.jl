using BioAlignments
using BioSymbols

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

# FUNCTIONS ###############################################################################

# if isfastafile(args["reference"]), build DIAMOND database before running DIAMOND
# NOW ALSO FOUND IN DIAMOND.JL
function isfastafile(path::AbstractString)::Bool
    splits = split(path, ".")
    if length(splits) > 1
        extension = last(splits)
        if extension in FASTAEXTENSIONS
            return true
        end
    end
    return false
end

function isdiamonddatabase(path::AbstractString)::Bool
    splits = split(path, ".")
    if length(splits) > 1
        extension = last(splits)
        if isequal(extension, DIAMONDDB)
            return true
        end
    end
    return false
end

# Matrix ##################################################################################

"""
    checkmatrixtype(args::Dict{String, Any})

Ensure that only one of `--matrix` and `--custom-matrix` was supplied as command line 
argument. Also check that `--lambda` and `--K` were set in `--diamond-flags` in case a 
custom matrix is used.
"""
function checkmatrixtype(args::Dict{String, Any})
    if "matrix" in keys(args) && "custom-matrix" in keys(args)
        error("You can only supply one of --matrix and --custom-matrix")
    elseif "custom-matrix" in keys(args)
        if !("--lambda" in args["diamond-flags"]) || !("--K" in args["diamond-flags"])
            error("DIAMOND custom scoring matrices require setting the --lambda and", 
                  " --K options.")
        end
    end
end

"""
    getmatrixtype(args::Dict{String, Any})::String

Get the type of matrix the user supplied, i.e. `matrix` or `custom-matrix`.
"""
function getmatrixtype(args::Dict{String, Any})::String
    if "custom-matrix" in keys(args)
        return "custom-matrix"
    else
        return "matrix"
    end
end

"""
    getmatrixname(args::Dict{String, Any})::String

Get the name (in case the parameter `--matrix` was set) or the path (in case 
`--custom-matrix` was set) of the matrix the user provided on the command line. 
The name supplied with the argument `--matrix` must correspond to one of `BioAlignements`'s 
pre-defined matrices (BLOSUM45, BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, PAM30, PAM70, or 
PAM250). If neither `--matrix` nor `--custom-matrix` were supplied, the default 
`Patchwork.MATRIX` is used.
"""
function getmatrixname(args::Dict{String, Any})::String
    matrixtype = getmatrixtype(args)
    if isequal("custom-matrix", matrixtype)
        matrix = args[matrixtype]
    else                                    # --matrix default BLOSUM62 or user-specified
        matrix = uppercase(args[matrixtype])
        if !(matrix in keys(MATRICES))
            error("Invalid scoring matrix: $matrix.")
        end
    end
    return matrix
end

"""
    read_custommatrix(path::AbstractString)::SubstitutionMatrix

Parse the file provided in `--custom-matrix` as an amino acid substitution matrix, 
creating a new `BioAlignments.SubstitutionMatrix` object from the data. 
"""
function read_custommatrix(path::AbstractString)::SubstitutionMatrix
    return BioAlignments.parse_ncbi_submat(BioSymbols.AminoAcid, path)
end

# matrixname = getmatrixname(args) and matrixtype = getmatrixtype(args)
"""
    getmatrix(matrixname::AbstractString, matrixtype::AbstractString)::SubstitutionMatrix

Get the `BioAlignments.SubstitutionMatrix` object corresponding to the substitution matrix. 
If `--custom-matrix` was set, this function will parse the provided file as an amino acid 
substitution matrix and create a new `BioAlignments.SubstitutionMatrix` object. Otherwise, 
the matrix must correspond to one of `BioAlignments`'s pre-defined matrices (BLOSUM45, 
BLOSUM50, BLOSUM62, BLOSUM80, BLOSUM90, PAM30, PAM70, or PAM250).
"""
function getmatrix(matrixname::AbstractString, matrixtype::AbstractString)::SubstitutionMatrix
    if isequal("custom-matrix", matrixtype)
        return read_custommatrix(matrixname)            # matrixname is a path
    else
        return MATRICES[matrixname]
    end
end

# Gap Penalties ###########################################################################

"""
    setgapopen!(args::Dict{String, Any})

Set the gap-opening penalty to the default value, depending on the matrix provided in 
`--matrix`, if no gap-opening penalty was explicitely provided via `--gapopen`. 
If necessary, convert the penalty to a positive integer. This function throws an error if 
`--custom-matrix` was set without providing a gap-opening penalty.
"""
function setgapopen!(args::Dict{String, any})
    if "gapopen" in keys(args)
        args["gapopen"] = abs(args["gapopen"])
        return
    end
    matrixtype = getmatrixtype(args)
    if isequal("matrix", matrixtype)
        args["gapopen"] = GAPS[args[matrixtype]][1]
    elseif isequal("custom-matrix", matrixtype)
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
function setgapextend!(args::Dict{String, any})
    if "gapextend" in keys(args)
        args["gapextend"] = abs(args["gapextend"])
        return
    end
    matrixtype = getmatrixtype(args)
    if isequal("matrix", matrixtype)
        args["gapextend"] = GAPS[args[matrixtype]][2]
    elseif isequal("custom-matrix", matrixtype)
        error("Custom scoring matrices require setting the --gapextend option.")
    end
end

"""
    checkgappenalty(matrix::AbstractString, gapopen::Int64, gapextend::Int64)

Ensure that the gap penalties are compatible with the matrix.
"""
function checkgappenalty(matrix::AbstractString, gapopen::Int64, gapextend::Int64)
    if matrix in keys(MATRICES)
        allowedgaps = GAPS[matrix]
        if isnothing(findfirst(isequal((gapopen, gapextend)), allowedgaps))
            error("Invalid gap open and/or gap extend scores.")
        end
    end
end

# Diamond Flags ###########################################################################

"""
    set_diamondframeshift!(args::Dict{String, Any})

If no `--frameshift` argument was provided in `--diamond-flags`, set the `--frameshift` 
to the default value.
"""
function set_diamondframeshift!(args::Dict{String, Any})
    if !("--frameshift" in args["diamond-flags"])
        push!(args["diamond-flags"], "--frameshift", 15)
    else
        args["diamond-flags"] = ["--frameshift", 15]
    end
end

"""
    set_diamondmode!(args::Dict{String, Any})

If none of `--fast`, `--sensitive`, `--more-sensitive`, or `--ultra-sensitive` was set in 
`--diamond-flags`, set the mode to the default value.
"""
function set_diamondmode!(args::Dict{String, Any})
    if !("--fast" in args["diamond-flags"] 
         || "--sensitive" in args["diamond-flags"] 
         || "--more-sensitive" in args["diamond-flags"] 
         || "--ultra-sensitive" in args["diamond-flags"])   
        push!(args["diamond-flags"], "--ultra-sensitive")
    else
        args["diamond-flags"] = ["--ultra-sensitive"]
    end
end

"""
    checkdiamondflags(args::Dict{String, Any})

Ensure that no conflicts are caused by command line arguments provided in `--diamond-flags`, 
i.e. that parameters used and set by `Patchwork` are not also provided to `DIAMOND` in a 
separate place. 
"""
function checkdiamondflags(args::Dict{String, Any})
    flags = args["diamond-flags"]
    patchworkflags = ["--matrix", "--custom-matrix", "--gapopen", "--gapextend", 
                      "--threads"]
    for flag in patchworkflags
        if flag in flags
            error("Please use Patchwork's own $flag option instead.")
        end
    end
    if "--query" in flags || "-q" in flags
        error("You should not explicitly provide a DIAMOND query.",
              "The file provided to Patchwork as --contigs will be used as DIAMOND query.")
    end
    if "--db" in flags || "-d" in flags
        error("You should not explicitly provide a DIAMOND database.",
              "The file or database provided with --reference will be used by DIAMOND.",
              "If you are providing a database, please make sure to also set --database.")
    end
    if "--outfmt" in flags || "-f" in flags
        error("The DIAMOND output format is built into Patchwork and cannot be changed.")
    end
    if "--out" in flags || "-o" in flags
        error("The DIAMOND results file is built into Patchwork and cannot be changed.")
    end
end

# makedb flags can't include threads bc handled by patchwork
function checkmakedbflags(args::Dict{String, Any})
    flags = args["makedb-flags"]
    if "--db" in flags || "-d" in flags
        error("The creation of your DIAMOND database will be internally handled.", 
              "You should not explicitly provide a DIAMOND database filename.")
    if "--threads"  in flags
        error("Please use Patchwork's own --threads option instead.")
    end
end

# Check and Collect all DIAMOND Flags in an Array #########################################

"""
    setpatchworkflags!(args::Dict{String, Any})

Ensure that no conflicts are caused by combinations of, or missing, command line arguments. 
Set gap-opening and gap-extension penalties correctly, modifying `args`. This function 
throws an error if any conflicts are detected. 
"""
function setpatchworkflags!(args::Dict{String, Any})     # Run this fct. before the next.
    checkmatrixtype(args)
    setgapopen!(args)
    setgapextend!(args)
    checkgappenalty(getmatrixname(args), args["gapopen"], args["gapextend"])
end

"""
    setdiamondflags!(args::Dict{String, Any})

Ensure that no conflicts are caused by combinations of, or missing, command line arguments. 
Set `DIAMOND` frameshift and mode correctly, modifying `args`. This function throws an 
error if any conflicts are detected. 
"""
function setdiamondflags!(args::Dict{String, Any})      # Run this fct. before the next.
    checkdiamondflags(args)
    set_diamondframeshift!(args)
    set_diamondmode!(args)
    if isequal(getmatrixtype(args), "custom-matrix")
        checkmakedbflags(args)
        push!(args["makedb-flags"], "--threads", args["threads"])
    end
end

"""
    collectdiamondflags(args::Dict{String, Any})::Vector{String}

Build the complete vector of parameters for running `DIAMOND`. 
"""
function collectdiamondflags(args::Dict{String, Any})::Vector{String}
    @assert "gapopen" in keys(args) && "gapextend" in keys(args) """set gap penalties with 
    setpatchworkflags!(args) before calling this function"""
    return [args["diamond-flags"]..., "--" * getmatrixtype(args), getmatrixname(args), 
            "--gapopen", args["gapopen"], "--gapextend", args["gapextend"], "--threads", 
            args["--threads"]]
end