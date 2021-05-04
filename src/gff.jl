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
TESTATTRIBUTES = GFF3Attributes(TESTRECORD)

columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
attributes = ["ID", "Name", "Alias", "Parent", "Target", "Gap", "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular"]

mutable struct GFF3Attributes
    id::Union{String, Missing}
    name::Union{String, Missing}
    alias::Union{AbstractVector{String}, Missing} # multiple values
    parent::Union{AbstractVector{String}, Missing} # multiple values
    target::Union{String, Missing}
    gap::Union{String, Missing} # RegEx?
    derivesfrom::Union{String, Missing}
    note::Union{String, Missing}
    dbcrossreference::Union{AbstractVector{String}, Missing} # multiple values
    ontologyterm::Union{AbstractVector{String}, Missing} # multiple values
    iscircular::Union{Bool, Missing}

    # attention: case sensitivity in attribute keys ! 
    function GFF3Attributes(gff3record::GFF3.Record)
        attributes = ["ID", "Name", "Alias", "Parent", "Target", "Gap", "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular"]
        allowsmultiplevalues = [3,4,9,10] # attribute data fields that allow multiple values
        result = Vector{Any}(missing, length(attributes))
        recordattributes = Dict(GFF3.attributes(gff3record))
        recordkeys = collect(keys(recordattributes))
        
        for (i, key) in enumerate(attributes)
            if key in recordkeys
                result[i] = recordattributes[key]
            end
            if !(i in allowsmultiplevalues) && !(result[i] === missing)
                result[i] = result[i][1]
            end
        end
        return new(result...)
    end
end 

#function Base.convert(::AbstractVector{Any}, gff3attributes::GFF3Attributes)
#    return [gff3attributes.id, gff3attributes.name, gff3attributes.alias, gff3attributes.parent, 
#            gff3attributes.target, gff3attributes.gap, gff3attributes.derivesfrom, gff3attributes.note, 
#            gff3attributes.dbcrossreference, gff3attributes.ontologyterm, gff3attributes.iscircular]
#end

# function to convert to list of tuples or pairs ? 
function Base.convert(::Dict, gff3attributes::GFF3Attributes)
    attributes = ["ID", "Name", "Alias", "Parent", "Target", "Gap", "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular"]
    return Dict(zip(attributes, [gff3attributes.id, gff3attributes.name, gff3attributes.alias, gff3attributes.parent, 
                                 gff3attributes.target, gff3attributes.gap, gff3attributes.derivesfrom, gff3attributes.note, 
                                 gff3attributes.dbcrossreference, gff3attributes.ontologyterm, gff3attributes.iscircular]))
end

function Base.Dict(gff3attributes::GFF3Attributes)
    return convert(Dict, gff3attributes)
end

# For :directive and :comment gff3-records ??
function Base.convert(::Type{AbstractVector{Any}}, gff3record::GFF3.Record) 
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
            result[i] = GFF3Attributes(gff3record)
        end
    end
    return result
end

function Base.AbstractVector(gff3record::GFF3.Record)
    return convert(AbstractVector{Any}, gff3record)
end

function Base.convert(::Type{Tuple}, gff3record::GFF3.Record)
    return Tuple(AbstractVector{Any}(gff3record))
end

function Base.Tuple(gff3record::GFF3.Record)
    return convert(Tuple{Any}, gff3record)
end

# convert vector or tuple back to gff3-record --> for gff3-writer
function Base.convert(::Type{GFF3.Record}, gff3record::Union{AbstractVector{Any}, Tuple{Any}}) 
    tmp = replace(join(gff3record[1:8], "\t"), "missing" => ".")
    attributes_tmp = []
    for (key, value) in Dict(gff3record[9]) 
        if !ismissing(value)
            if isa(AbstractVector{String}, value)
                push!(attributes_tmp, join([key, join(value, ",")], "="))
            else
                push!(attributes_tmp, join([key, value], "="))
            end
        end
    end
    return GFF3.Record(join(tmp, join(attributes_tmp, ";"), "\t"))
end

mutable struct GFF3RecordCollection
    records::Vector{GFF3.Record}
end

function GFF3RecordCollection()
    return GFF3RecordCollection([])
end

function Base.length(gff3recordcollection::GFF3RecordCollection)
    return length(gff3recordcollection.records)
end

function Base.push!(gff3recordcollection::GFF3RecordCollection, gff3record::GFF3.Record)
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

function Base.iterate(gff3recordcollection::GFF3RecordCollection)
    if isempty(gff3recordcollection)
        return nothing
    else
        i = firstindex(gff3recordcollection)
        return getindex(gff3recordcollection, i), i + 1
    end
end

function Base.iterate(gff3recordcollection::GFF3RecordCollection, i::Int)
    if i > lastindex(gff3recordcollection)
        return nothing
    else
        return getindex(gff3recordcollection, i), i + 1
    end
end

# better way ??
function Base.convert(::Type{DataFrames.DataFrame}, gff3recordcollection::GFF3RecordCollection)
    columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    result = []
    for record in gff3recordcollection
        if BioCore.isfilled(record)
            push!(result, Tuple(record))
        end
    end
    result = DataFrames.DataFrame(result)
    DataFrames.rename!(result, columns)
    return result
end

function DataFrames.DataFrame(gff3recordcollection::GFF3RecordCollection)
    return convert(DataFrames.DataFrame, gff3recordcollection)
end

# read whole gff3 file: 
# what about directive and comment records ? 
# --> checkkind and checkfilled for each record before saving into result ? 
# filter result DF for lines with invalid data 
# Do you really nead to convert GFF3.Records to Tuples for creating a DF ? 
function readfeatureGFF3(reader::GFF3.Reader)
    result = GFF3RecordCollection()
    for record in reader
        if(BioCore.isfilled(record)) # && GFF3.isfeature(record)) # omit others; is checking type :feature really always necessary? 
            push!(result, record)
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




# once you have a DataFrame from a gff3-file: check for errors in data (data cleaning)
# maybe for later:



function requires_phase(gff3record::GFF3.Record)
    GFF3.checkfilled(gff3record)
    GFF3.checkkind(gff3record, :feature)
    if GFF3.source(gff3record) == "CDS" && GFF3.phase(gff3record) != r"[012]"
        print("Phase is required for feature records.") # Throw an error!
    end
end

function requires_phase(gff3recordvector::AbstractVector{Any})
    gff3record = convert(GFF3.Record, gff3recordvector)
    return requires_phase(gff3record)
end

function hasvalidseqend(gff3record::GFF3.Record)
    GFF3.checkfilled(gff3record)
    GFF3.checkkind(gff3record, :feature)
    return (GFF3.seqend(gff3record) - GFF3.seqstart(gff3record)) >= 0 # missingerror if no start or end
end


for c in collection
    print(c)
end