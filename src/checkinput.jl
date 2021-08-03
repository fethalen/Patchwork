using BioAlignments
using BioSymbols

cd("src")
include("diamond.jl")

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

# Matrix ##################################################################################

function getmatrixtype(args::Dict{String, Any})::AbstractString
    if "matrix" in keys(args) && "custom-matrix" in keys(args)
        error("You can only supply one of --matrix and --custom-matrix")
    elseif "custom-matrix" in keys(args)
        if !("--lambda" in args["diamond-flags"]) || !("--K" in args["diamond-flags"])
            error("Custom scoring matrices require setting the --lambda and --K options.")
        end
        return "--custom-matrix"
    else
        return "--matrix"
    end
end

function getmatrixname(args::Dict{String, Any})::AbstractString
    if "custom-matrix" in keys(args)
        matrix = args["custom-matrix"]
    else                                    # --matrix default BLOSUM62 or user-specified
        matrix = uppercase(args["matrix"])
        if !(matrix in keys(MATRICES))
            error("Invalid scoring matrix: $matrix.")
        end
    end
    return matrix
end

# supporting only amino acid alignmnents!
function read_custommatrix(path::AbstractString)::SubstitutionMatrix
    return parse_ncbi_submat(AminoAcid, path)
end

# matrixname = getmatrixname(args) and matrixtype = getmatrixtype(args)
function getmatrix(matrixname::AbstractString, matrixtype::AbstractString)::SubstitutionMatrix
    if matrixtype == "--matrix"
        return MATRICES[matrixname]
    else
        return read_custommatrix(matrixname)
    end
end

# Gap Penalties ###########################################################################

function setgapopen!(args::Dict{String, any})
    if "gapopen" in keys(args)
        args["gapopen"] = abs(args["gapopen"])
    end
    if "matrix" in keys(args)
        args["gapopen"] = GAPS[args["matrix"]][1]
    elseif "custom-matrix" in keys(args)
        error("Custom scoring matrices require setting the --gapopen option.")
    end
end

function setgapextend!(args::Dict{String, any})
    if "gapextend" in keys(args)
        args["gapextend"] = abs(args["gapextend"])
    end
    if "matrix" in keys(args)
        args["gapextend"] = GAPS[args["matrix"]][2]
    elseif "custom-matrix" in keys(args)
        error("Custom scoring matrices require setting the --gapextend option.")
    end
end

function checkgappenalty(matrix::AbstractString, gapopen::Int64, gapextend::Int64)
    if matrix in keys(MATRICES)
        allowedgaps = GAPS[matrix]
        if isnothing(findfirst(isequal((gapopen, gapextend)), allowedgaps))
            error("Invalid gap open and/or gap extend scores.")
        end
    end
end

# Diamond Flags ###########################################################################

function set_diamondframeshift!(args::Dict{String, Any})
    if !("--frameshift" in args["diamond-flags"])
        push!(args["diamond-flags"], "--frameshift", 15)
    end
end

function set_diamondmode!(args::Dict{String, Any})
    if !("--fast" in args["diamond-flags"] 
         || "--sensitive" in args["diamond-flags"] 
         || "--more-sensitive" in args["diamond-flags"] 
         || "--ultra-sensitive" in args["diamond-flags"])   
        push!(args["diamond-flags"], "--ultra-sensitive")
    end
end

function checkdiamondflags(args::Dict{String, Any})
    flags = args["diamond-flags"]
    patchworkflags = ["--matrix", "--custom-matrix", "--gapopen", "--gapextend"]
    for flag in patchworkflags
        if flag in flags
            error("Please use Patchwork's own $flag option instead.")
        end
    end
    if "--query" in flags || "-q" in flags
        error("""You should not explicitly provide a DIAMOND query. 
        The file provided to Patchwork as --contigs will be used as DIAMOND query.""")
    end
    if "--db" in flags || "-d" in flags
        error("""You should not explicitly provide a DIAMOND database. 
        The file or database provided to Patchwork as --reference will be used by DIAMOND.
        If you are providing a database, please make sure to also set --database.""")
    end
    if "--outfmt" in flags || "-f" in flags
        error("The DIAMOND output format is built into Patchwork and cannot be changed.")
    end
    if "--out" in flags || "-o" in flags
        error("The DIAMOND results file is built into Patchwork and cannot be changed.")
    end
end

# Collect all DIAMOND Flags in an Array ###################################################

function collectdiamondflags(args::Dict{String, Any})::Vector{Any}
    checkdiamondflags(args)                 # Could also be done outside of this function
    matrixtype = getmatrixtype(args)
    matrix = getmatrixname(args)
    setgapopen!(args)                       # Do this outside of this function ?
    setgapextend!(args)
    checkgappenalty(matrix, args["gapopen"], args["gapextend"])
    set_diamondframeshift!(args)
    set_diamondmode!(args)
    return [args["diamond-flags"]..., matrixtype, matrix, "--gapopen", args["gapopen"], 
            "--gapextend", args["gapextend"]]
end