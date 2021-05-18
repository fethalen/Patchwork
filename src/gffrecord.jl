# Provides I/O utilities and basic types for working with the general feature
# format (GFF)

import Base
import GFF3
import DataFrames
using GenomicFeatures
using BioCore

TESTDIR = "/home/clara/Desktop/SHK-Job_Bleidorn/Projects/BioFmtSpecimens/GFF3/"
TESTFILES = [TESTDIR * "au9_scaffold_subset.gff3", TESTDIR * "directives.gff3", TESTDIR * "TAIR10.part.gff.bgz"]
TESTRECORD = GFF3.Record("CCDS1.1\tCCDS\tgene\t801943\t802434\t.\t-\t.\tID=LINC00115;Name=LINC00115")

columns = [:seqid, :source, :type, :start, :end_, :score, :strand, :phase, :attributes]
attributelist = ["ID", "Name", "Alias", "Parent", "Target", "Gap", "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular"]

###############################################################################################################
################################################ GFF3Attributes ###############################################
###############################################################################################################

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

    # attention: case sensitivity in attribute keys ! 
    function GFF3Attributes(gff3record::GFF3.Record)
        allowsmultiplevalues = [3,4,9,10,12]
        result = Vector{Any}(missing, length(attributelist))
        result[12] = [] # other
        recordattributes = Dict(GFF3.attributes(gff3record))
        recordkeys = collect(keys(recordattributes))
        
        for key in recordkeys
            if key in attributelist
                index = findfirst(isequal(key), attributelist)
                if !(index in allowsmultiplevalues)
                    @assert length(recordattributes[key]) == 1 "Attribute $key only allows one value per record."
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

        return new(result...)
    end
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

    function GFF3Record()
        return new(Vector{Any}(missing, 9)...)
    end

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

        return new(result...)
    end

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
    print(io, summary(record), ":")
    if BioCore.isfilled(record)
        println(io)
        println(io, "        seqid: ", record.seqid)
        println(io, "       source: ", record.source)
        println(io, "         type: ", record.type)
        println(io, "        start: ", record.start)
        println(io, "          end: ", record.end_)
        println(io, "        score: ", record.score)
        println(io, "       strand: ", record.strand)
        println(io, "        phase: ", record.phase)
        println(io, "   attributes: ", String(record.attributes))
    else
        println(io, " empty record")
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

function checkstartend(record::GFF3Record)
    if record.start > record.end_ 
        error("Sequence start must be before or same as end.")
    end
end

function checkcdsphase(record::GFF3Record)
    if record.type == "CDS" && ismissing(record.phase)
        error("Feature type \"CDS\" requires a phase.")
    end
end
