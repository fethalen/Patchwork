# Provides I/O utilities and basic types for working with the general feature
# format (GFF)

import Base
import GFF3
import DataFrames
#using GenomicFeatures
using BioCore

TESTDIR = "/home/clara/Desktop/SHK-Job_Bleidorn/Projects/BioFmtSpecimens/GFF3/"
TESTFILES = [TESTDIR * "au9_scaffold_subset.gff3", TESTDIR * "directives.gff3", TESTDIR * "TAIR10.part.gff.bgz"]
TESTRECORD = GFF3.Record("CCDS1.1\tCCDS\tgene\t801943\t802434\t.\t-\t.\tID=LINC00115;Name=LINC00115")

columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

# change the variable names to something shorter
# throw errors instead of filtering

###############################################################################################################
################################### converting GFF3.Records to other types ####################################
###############################################################################################################

# change that to individual data fields
mutable struct GFF3Record
    record::Tuple

    function GFF3Record()
        return new(Tuple([]))
    end

    function GFF3Record(gff3record::GFF3.Record)
        result = AbstractVector{Any}(split(GFF3.content(gff3record), '\t'))
        for i in eachindex(result)
            if result[i] == "."
                result[i] = missing # or nothing ? or don't change at all ? 
            elseif i in 4:5
                result[i] = parse(Int64, result[i])
            elseif i == 6
                result[i] = parse(Float64, result[i])
            # else
            #   result[i] = String(result[i])
            elseif i == 9
                result[i] = GFF3.attributes(gff3record)
            end
        end
        return new(Tuple(result))
    end
end

function Base.getindex(gff3record::GFF3Record, index::Integer)
    return getindex(gff3record.record, index)
end
function Base.getindex(gff3record::GFF3Record, indexes::UnitRange{Int})
    result = []
    for i in indexes
        push!(result, getindex(gff3record, i))
    end
    return Tuple(result)
end

Base.length(gff3record::GFF3Record) = length(gff3record.record)
Base.firstindex(gff3record::GFF3Record) = 1
Base.lastindex(gff3record::GFF3Record) = length(gff3record)
Base.eachindex(gff3record::GFF3Record) = Base.OneTo(lastindex(gff3record))

function Base.isempty(gff3record::GFF3Record)
    return isempty(gff3record.record)
end

function Base.iterate(gff3record::GFF3Record)
    return iterate(gff3record.record)
end

function Base.iterate(gff3record::GFF3Record, i::Int)
    return iterate(gff3record.record, i)
end

# convert vector or tuple back to gff3-record --> for gff3-writer
function Base.convert(::Type{GFF3.Record}, gff3record::GFF3Record) 
    return GFF3.Record(String(gff3record))
end
function GFF3.Record(gff3record::GFF3Record)
    return convert(GFF3.Record, gff3record)
end

function Base.convert(::Type{String}, gff3record::GFF3Record)
    tmp = replace(join(gff3record[1:8], "\t"), "missing" => ".")
    attributes_tmp = []
    for (key, value) in gff3record[9]
        if !ismissing(value)
            if isa(value, AbstractVector{String})
                push!(attributes_tmp, join([key, join(value, ",")], "="))
            else
                push!(attributes_tmp, join([key, value], "="))
            end
        end
    end
    return join([tmp, join(attributes_tmp, ";")], "\t")
end
function Base.String(gff3record::GFF3Record)
    return convert(String, gff3record)
end

function Base.convert(::Type{DataFrames.DataFrame}, gff3record::GFF3Record)
    result = DataFrames.DataFrame([gff3record.record])
    DataFrames.rename!(result, columns)
    return result
end
function DataFrames.DataFrame(gff3record::GFF3Record)
    return convert(DataFrames.DataFrame, gff3record)
end

###############################################################################################################
######################################### GFF3RecordCollection ################################################
###############################################################################################################

mutable struct GFF3RecordCollection
    records::Vector{GFF3Record}
end

function GFF3RecordCollection()
    return GFF3RecordCollection([])
end

function GFF3RecordCollection(reader::GFF3.Reader)
    result = GFF3RecordCollection(readfeatureGFF3(reader))
    filterGFF3input!(result)
    return result
end

function GFF3RecordCollection(file::AbstractString; idx=:auto) # Tabix index ??
    result = GFF3RecordCollection(readfeatureGFF3(file; idx))
    filterGFF3input!(result)
    return result
end

function GFF3RecordCollection(input::IO; idx=nothing) # Tabix index ??
    result = GFF3RecordCollection(readfeatureGFF3(input; idx))
    filterGFF3input!(result)
    return result
end

# read whole gff3 file: 
# what about directive and comment records ? 
# --> checkkind and checkfilled for each record before saving into result ? 
# Do you really nead to convert GFF3.Records to Tuples for creating a DF ? 
function readfeatureGFF3(reader::GFF3.Reader)
    #result = GFF3RecordCollection()
    result = []
    for record in reader
        if BioCore.isfilled(record) && GFF3.isfeature(record) # omit others; is checking type :feature really always necessary? 
            push!(result, GFF3Record(record))
        end
    end
    return result
end

# you lose the directives and comments here:
function readfeatureGFF3(file::AbstractString; idx=:auto) # Tabix index ??
    reader = GFF3.Reader(file; index=idx) # index = nothing ??
    result = readfeatureGFF3(reader)
    close(reader)
    return result
end

# ... and here, too: 
function readfeatureGFF3(input::IO; idx=nothing) # Tabix index ??
    reader = GFF3.Reader(input; index=idx) # index = nothing ??
    result = readfeatureGFF3(reader)
    close(reader)
    return result
end

function Base.length(gff3recordcollection::GFF3RecordCollection)
    return length(gff3recordcollection.records)
end

function Base.push!(gff3recordcollection::GFF3RecordCollection, gff3record::GFF3Record)
    return push!(gff3recordcollection.records, gff3record)
end

function Base.getindex(gff3recordcollection::GFF3RecordCollection, index::Integer)
    return getindex(gff3recordcollection.records, index)
end

Base.firstindex(gff3recordcollection::GFF3RecordCollection) = 1
Base.lastindex(gff3recordcollection::GFF3RecordCollection) = length(gff3recordcollection)
Base.eachindex(gff3recordcollection::GFF3RecordCollection) = Base.OneTo(lastindex(gff3recordcollection))

function Base.isempty(gff3recordcollection::GFF3RecordCollection)
    return isempty(gff3recordcollection.records)
end

# function Base.iterate(gff3recordcollection::GFF3RecordCollection) = iterate(gff3recordcollection.records)
function Base.iterate(gff3recordcollection::GFF3RecordCollection)
    if isempty(gff3recordcollection)
        return nothing
    else
        i = firstindex(gff3recordcollection)
        return getindex(gff3recordcollection, i), i + 1
    end
end

# function Base.iterate(gff3recordcollection::GFF3RecordCollection, i::Int) = iterate(gff3recordcollection.records, i)
function Base.iterate(gff3recordcollection::GFF3RecordCollection, i::Int)
    if i > lastindex(gff3recordcollection)
        return nothing
    else
        return getindex(gff3recordcollection, i), i + 1
    end
end

function Base.eltype(::Type{GFF3RecordCollection})
    return GFF3Record
end

function Base.deleteat!(gff3recordcollection::GFF3RecordCollection, i::Int)
    return deleteat!(gff3recordcollection.records, i)
end
function Base.popat!(gff3recordcollection::GFF3RecordCollection, i::Int)
    return popat!(gff3recordcollection.records, i)
end
function Base.pop!(gff3recordcollection::GFF3RecordCollection)
    return popat!(gff3recordcollection, lastindex(gff3recordcollection))
end

function Base.filter(f, gff3recordcollection::GFF3RecordCollection)
    return filter(f, gff3recordcollection.records)
end
function Base.filter!(f, gff3recordcollection::GFF3RecordCollection)
    return filter!(f, gff3recordcollection.records)
end

# better way ??
function Base.convert(::Type{DataFrames.DataFrame}, gff3recordcollection::GFF3RecordCollection)
    columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    result = []
    for gff3record in gff3recordcollection
        if !isempty(gff3record)
            push!(result, gff3record.record)
        end
    end
    result = DataFrames.DataFrame(result)
    DataFrames.rename!(result, columns)
    return result
end

function DataFrames.DataFrame(gff3recordcollection::GFF3RecordCollection)
    return convert(DataFrames.DataFrame, gff3recordcollection)
end

################################################################################################################
###################################### checking input data/ filtering ##########################################
################################################################################################################

# What happened before: You've loaded the whole file into a GFF3RecordColleciton.
testcollection = readfeatureGFF3(TESTFILES[1]; idx=nothing)
testcollection = GFF3RecordCollection(TESTFILES[1]; idx=nothing)
testdf = DataFrames.DataFrame(testcollection)

# type: must be a SO sequence_feature or is_a child - How do you check that? 

# start <= end
filter(["start", "end"] => (start, end_) -> start <= end_, testdf)
# for GFF3RecordCollection{GFF3.Record}: filter(rec -> (GFF3.seqstart(rec) <= GFF3.seqend(rec)), testcollection)
filter(rec -> (rec[4] <= rec[5]), testcollection)

# phase required in CDS features (and is_a children?)
filter(["type", "phase"] => (type, phase) -> !(type == "CDS" && phase === missing), testdf)
#filter(rec -> !(GFF3.featuretype(rec) == "CDS" && GFF3.phase(rec) === missing), testcollection)
filter(rec -> !(rec[3] == "CDS" && rec[8] === missing), testcollection)

# attributes ##################################################################################################
# alias, parent, dbxref and ontologyterm allow multiple values

# ID required for features with children; i.e. you can only describe a part_of relationship via the IDs
# --> if there's no ID, there's nothing to be done about it/ how are you supposed to know if they just forgot it?
function allids(gff3recordcollection::GFF3RecordCollection)
    allids = []
    for gff3record in gff3recordcollection
        recordattributes = Dict(gff3record[9])
        if "ID" in keys(recordattributes)
            @assert length(recordattributes["ID"]) == 1 "Can only assign up to one ID per record."
            push!(allids, recordattributes["ID"]...)
        end
    end
    return allids
end
function checkids(gff3recordcollection::GFF3RecordCollection)
    ids = allids(gff3recordcollection)
    for (i, id) in enumerate(ids)
        @assert length(id) == 1 "Can only assign up to one ID per record."
        ids[i] = id[1]
    end
    return ids
end

# parent: only part_of relationship (between sequence_features?), no cycles
# ID must be found within file (?)
function allparents(gff3recordcollection::GFF3RecordCollection)
    parents = []
    for gff3record in gff3recordcollection
        recordattributes = Dict(gff3record[9])
        if "Parent" in keys(recordattributes) && !(recordattributes["Parent"] in parents)
            push!(parents, recordattributes["Parent"]...)
        end
    end
    return parents
end
function checkparents(gff3recordcollection::GFF3RecordCollection)
    ids = allids(gff3recordcollection)
    parents = allparents(gff3recordcollection)
    for parentid in parents
        @assert parentid in ids "Parent $parentid must be found within scope of file"
    end
    return
end
checkparents(testcollection)

# right formatting of "target": "target_id start end [strand]" with strand = "+" or "-"
# target_id must be found in file (?)

# r"^.+ [0-9]+ [0-9]+( [+-])?$"
function alltargets(gff3recordcollection::GFF3RecordCollection)
    targets = []
    for gff3record in gff3recordcollection
        recordattributes = Dict(gff3record[9])
        if "Target" in keys(recordattributes) && !(recordattributes["Target"] in targets)
            push!(targets, recordattributes["Target"])
        end
    end
    return targets
end
function checktargets(gff3recordcollection::GFF3RecordCollection)
    ids = allids(gff3recordcollection)
    targets = alltargets(gff3recordcollection)
    for t in targets
        @assert length(t) == 1 "Can only assign up to one alignement target per record."
        @assert !isnothing(match(r"^.+ [0-9]+ [0-9]+( [+-])?$", t[1])) "Target $(t[1]) has wrong format."
        target = split(t[1], ' ')
        @assert target[1] in ids "Target $(target[1]) must be found within scope of file."
    end
    return
end

# gap: only when target was specified (?)
#filter(rec -> !(GFF3Attributes(rec).target === missing && !ismissing(GFF3Attributes(rec).gap)), testcollection.records)
function gap_hastarget(gff3record::GFF3Record)
    recordattributes = Dict(gff3record[9])
    attributekeys = keys(recordattributes)
    if "Gap" in attributekeys && !("Target" in attributekeys)
        return false
    else
        return true
    end
end
filter(rec -> gap_hastarget(rec), testcollection)

# check gap format: 
function checkgaps!(gff3recordcollection::GFF3RecordCollection)
    pattern = r"(D[0-9]+|F[0-9]+|I[0-9]+|M[0-9]+|R[0-9]+)( (D[0-9]+|F[0-9]+|I[0-9]+|M[0-9]+|R[0-9]+))*"
    filter!(gff3record -> gap_hastarget(gff3record), gff3recordcollection)
    for gff3record in gff3recordcollection
        recordattributes = Dict(gff3record[9])
        attributekeys = keys(recordattributes)
        if "Gap" in attributekeys
            gap = recordattributes["Gap"]
            @assert length(gap) == 1 "Can only specify up to one alignment per record/ target."
            @assert !isnothing(match(pattern, gap)) "Gap $gap has wrong format."
        end
    end
    return
end

# derives_from: must be found in file (?)
function allderivesfrom(gff3recordcollection::GFF3RecordCollection)
    derivesfrom = []
    for gff3record in gff3recordcollection
        recordattributes = Dict(gff3record[9])
        if "Derives_from" in keys(recordattributes)
            push!(derivesfrom, recordattributes["Derives_from"])
        end
    end
    return derivesfrom
end
function checkderivesfrom(gff3recordcollection::GFF3RecordCollection)
    ids = allids(gff3recordcollection)
    derivesfrom = allderivesfrom(gff3recordcollection)
    for originid in derivesfrom
        @assert length(originid) == 1 "Can only specify up to one Derives_from origin per record."
        @assert originid[1] in ids "Derives_from origin $originid must be found within scope of file"
    end
    return
end

# the whole check-up:
# think about order of checks! 
# filter or throw an error? 
function filterGFF3input!(gff3recordcollection::GFF3RecordCollection)
    #gff3dataframe = DataFrames.DataFrame(gff3recordcollection)

    # start <= end
    filter!(record -> (record[4] <= record[5]), gff3recordcollection)
    
    # phase
    filter!(record -> !(record[3] == "CDS" && record[8] === missing), gff3recordcollection)

    # attributes: 
    # unique IDs 
    allids(gff3recordcollection)
    # parents
    checkparents(gff3recordcollection)
    # targets
    checktargets(gff3recordcollection)
    # gaps
    checkgaps!(gff3recordcollection)
    # derives_from origins
    checkderivesfrom(gff3recordcollection)
end

filterGFF3input!(testcollection)
