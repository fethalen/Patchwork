# Provides basic types for working with the general feature format (GFF).

#import Base
import GFF3
import DataFrames
using GenomicFeatures
using BioCore

columns = [:seqid, :source, :type, :start, :end_, :score, :strand, :phase, :attributes]
attributelist = ["ID", "Name", "Alias", "Parent", "Target", "Gap", "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular"]

###############################################################################################################
################################################ GFF3Attributes ###############################################
###############################################################################################################

# only makes sense when used inside a GFF3Record...
"""
    struct GFF3Attributes

Datastructure which holds the attributes of a `GFF3Record`.
All data fields are optional and can contain the value `missing`.
- `id::String`: the ID assigned to the record; required for features with children.
- `name::String`: an arbitrary name.
- `alias::AbstractVector{String}`: an arbitrary number of aliases.
- `parent::AbstractVector{String}`: the parent feature ID(s); indicates Sequence Ontology part_of relationships.
- `target::String`: the ID of the alignment target.
- `gap::String`: the alignment of the target to this record; requires that `target` be specified.
- `derivesfrom::String`: the ID of the ancestral feature, when the relationship is not a part_of relationship;
                         e.g. for polycistronic genes.
- `note::String`: a free text.
- `dbcrossreference::AbstractVector{String}`: an arbitrary number of database cross references.
- `ontologyterm::AbstractVector{String}`: an arbitrary number of cross references to ontology terms.
- `iscircular::Bool`: indicates whether the feature is circular.
- `other`: any attributes that are not covered by the above; stored as a Vector of Pairs.
"""
mutable struct GFF3Attributes # mutable ? 
    id::Union{String, Missing}
    name::Union{String, Missing}
    alias::Union{AbstractVector{String}, Missing} # multiple values
    parent::Union{AbstractVector{String}, Missing} # multiple values
    target::Union{String, Missing}
    gap::Union{String, Missing}
    derivesfrom::Union{String, Missing}
    note::Union{String, Missing} 
    dbcrossreference::Union{AbstractVector{String}, Missing} # multiple values
    ontologyterm::Union{AbstractVector{String}, Missing} # multiple values
    iscircular::Union{Bool, Missing}
    other::Any # Union{AbstractVector{Pair{String, AbstractVector{String}}}, Missing}
end 

# attention: case sensitivity in attribute keys ! 
"""
    GFF3Attributes(record::GFF3.Record)

Get all the attributes from a GFF3.Record.
This enables easier access to the data fields as well as validity checks on the data.
"""
function GFF3Attributes(record::GFF3.Record)
    allowsmultiplevalues = [3,4,9,10,12]
    result = Vector{Any}(missing, length(attributelist) + 1)
    result[12] = [] # other
    recordattributes = Dict(GFF3.attributes(record))
    recordkeys = collect(keys(recordattributes))
    
    for key in recordkeys
        if key in attributelist
            index = findfirst(isequal(key), attributelist)
            if key == "Note"
                result[index] = join(recordattributes[key], ',')
            elseif !(index in allowsmultiplevalues)
                if length(recordattributes[key]) != 1 
                    error("Attribute $key only allows one value per record.")
                end
                result[index] = recordattributes[key][1]
            else
                result[index] = recordattributes[key]
            end
        else
            push!(result[12], key => recordattributes[key])
        end
    end

    if(isempty(result[12]))
        result[12] = missing
    end

    return GFF3Attributes(result...)
end

function Base.isempty(attributes::GFF3Attributes)
    return (ismissing(attributes.id) && ismissing(attributes.name) && ismissing(attributes.alias) &&
            ismissing(attributes.parent) && ismissing(attributes.target) && ismissing(attributes.gap) &&
            ismissing(attributes.derivesfrom) && ismissing(attributes.note) && ismissing(attributes.dbcrossreference) &&
            ismissing(attributes.ontologyterm) && ismissing(attributes.iscircular) && ismissing(attributes.other))
end

function Base.convert(::Type{String}, attributes::GFF3Attributes)
    if isempty(attributes)
        return ""
    end

    recordattributes = [attributes.id, attributes.name, attributes.alias,
                        attributes.parent, attributes.target, attributes.gap, 
                        attributes.derivesfrom, attributes.note, attributes.dbcrossreference, 
                        attributes.ontologyterm, attributes.iscircular]
    attributes_tmp = []

    for (key, value) in zip(attributelist, recordattributes)
        if !ismissing(value)
            if isa(value, AbstractVector{String})
                push!(attributes_tmp, join([key, join(value, ",")], "="))
            else
                push!(attributes_tmp, join([key, value], "="))
            end
        end
    end

    if !ismissing(attributes.other)
        for (key, value) in attributes.other 
            if isa(value, AbstractVector{String})
                push!(attributes_tmp, join([key, join(value, ",")], "="))
            else
                push!(attributes_tmp, join([key, value], "="))
            end
        end
    end

    return join(attributes_tmp, ";")
end

function Base.String(attributes::GFF3Attributes)
    return convert(String, attributes)
end

###############################################################################################################
################################################ GFF3Record ###################################################
###############################################################################################################

"""
    struct GFF3Record

A Datastructure which holds a single feature record in the GFF3 format.
All data fields are optional and can contain the value `missing`.
- `seqid::String`: the ID of the landmark used to establish the coordinate system for the current feature.
- `source::String`: the program or procedure used to generate the feature.
- `type::String`: the type of the feature; restricted to Sequence Ontology sequence_features and their is_a children.
- `start::Int64`: the start of the feature (given as a positive 1-based integer), relative to the landmark in `seqid`.
- `end_::Int64`: the end of the feature (given as a positive 1-based integer), relative to the landmark in `seqid`.
- `score::Float64`: the score of the feature; e-values recommended for sequence similarity features
                    and p-values for ab initio gene prediction features.
- `strand::GenomicFeatures.Strand`: the strand of the feature ('+', '-', or '?').
- `phase::Char`: the phase of the feature (0, 1, or 2); required for features of type "CDS".
- `attributes::GFF3Attributes`: optional additional information, see `GFF3Attributes`. 
"""
mutable struct GFF3Record
    seqid::Union{String, Missing}
    source::Union{String, Missing}
    type::Union{String, Missing}
    start::Union{Int64, Missing}
    end_::Union{Int64, Missing}
    score::Union{Float64, Missing}
    strand::Union{GenomicFeatures.Strand, Missing}
    phase::Union{Char, Missing}
    attributes::Union{GFF3Attributes, Missing}
end

"""
    GFF3Record()

Empty initialization.
"""
function GFF3Record()
    return GFF3Record(Vector{Any}(missing, 9)...)
end

"""
    GFF3Record(record::GFF3.Record)

Convert the provided GFF3.Record to GFF3Record. 
This enables easier access to data fields as well as validity checks on the data.
"""
function GFF3Record(record::GFF3.Record)
    result = Vector{Any}(missing, 9)
    functions = [GFF3.seqid, GFF3.source, GFF3.featuretype, GFF3.seqstart, 
                 GFF3.seqend, GFF3.score, GFF3.strand, GFF3.phase, GFF3Attributes]

    for i in eachindex(result)
        try
            result[i] = functions[i](record)
        catch e
            if !(isa(e, BioCore.Exceptions.MissingFieldException))
                println(e)
                return GFF3Record()
            end
        end
    end

    if isempty(result[9])
        result[9] = missing
    end

    return GFF3Record(result...)
end

# third constructor for building a GFF3Record from a String 
# makes use of the fancy GFF3.jl FSM/parser (thus the extra conversion step).
"""
    GFF3Record(record::AbstractString)

Convert the provided String to a GFF3Record. 
The String is parsed to ensure it is correctly formatted.
"""
function GFF3Record(record::AbstractString)
    return GFF3Record(GFF3.Record(record))
end

function BioCore.isfilled(record::GFF3Record)
    return !(ismissing(record.seqid) && ismissing(record.source) && ismissing(record.type) && 
            ismissing(record.start) && ismissing(record.end_) && ismissing(record.score) && 
            ismissing(record.strand) && ismissing(record.phase) && ismissing(record.attributes))

end

# convert GFF3Record back to GFF3.Record --> for GFF3.Writer
function Base.convert(::Type{GFF3.Record}, record::GFF3Record) 
    return GFF3.Record(String(record))
end

function GFF3.Record(record::GFF3Record)
    return convert(GFF3.Record, record)
end

function Base.convert(::Type{String}, record::GFF3Record)
    tmp = replace(join([record.seqid, record.source, record.type, record.start, record.end_, 
                        record.score, record.strand, record.phase], "\t"), "missing" => ".")
    
    if ismissing(record.attributes)
        return tmp * "\t"
    else
        return tmp * "\t" * String(record.attributes)
    end
end

function Base.String(record::GFF3Record)
    return convert(String, record)
end

function Base.show(io::IO, record::GFF3Record)
    if BioCore.isfilled(record)
        println(io, "seqid: ", record.seqid)
        println(io, "source: ", record.source)
        println(io, "type: ", record.type)
        println(io, "start: ", record.start)
        println(io, "end: ", record.end_)
        println(io, "score: ", record.score)
        println(io, "strand: ", record.strand)
        println(io, "phase: ", record.phase)
        println(io, "attributes: ", String(record.attributes))
    else
        println(io, "empty record")
    end
end

function Base.convert(::Type{DataFrames.DataFrame}, record::GFF3Record)
    data = [[record.seqid], [record.source], [record.type], [record.start], [record.end_], 
            [record.score], [record.strand], [record.phase], [record.attributes]]
    result = DataFrames.DataFrame(data, columns)
    return result
end

function DataFrames.DataFrame(gff3record::GFF3Record)
    return convert(DataFrames.DataFrame, gff3record)
end

"""
    checkstartend(record::GFF3Record)

Verify that `start` and `end` of this record are valid.
`start` must be before or same as `end`.
Both are given as positive 1-based integers relative to the landmark `seqid`.
"""
function checkstartend(record::GFF3Record)
    if !ismissing(record.start) || !ismissing(record.end_)
        @assert !ismissing(record.seqid) "Sequence start and end must be given in relation to the landmark seqid."
        if !ismissing(record.start) 
            @assert record.start >= 1 "Sequence start and end must be given as positive 1-based integers relativ to the landmark seqid."
            if !ismissing(record.end_)
                @assert record.start <= record.end_ "Sequence end must be greater than or equal to start."
            end
        end
    end
end

"""
    checkcdsphase(record::GFF3Record)

Verify that a record of `type`== "CDS" also specifies a `phase`. 
"""
function checkcdsphase(record::GFF3Record)
    if !ismissing(record.type)
        @assert !(record.type == "CDS" && ismissing(record.phase)) "Feature type \"CDS\" requires a phase."
    end
end


"""
    gap_hastarget(record::GFF3Record)

Returns false if this record specifies a `gap` attribute but does not specify a `target` attribute, and true otherwise.
"""
function gap_hastarget(record::GFF3Record)
    if !ismissing(record.attributes) && !ismissing(record.attributes.gap) && ismissing(record.attributes.target)
        return false
    else
        return true
    end
end
