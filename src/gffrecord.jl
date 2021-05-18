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

###############################################################################################################
######################################### GFF3RecordCollection ################################################
###############################################################################################################

mutable struct GFF3RecordCollection
    records::Vector{GFF3Record}
end

function GFF3RecordCollection()
    return GFF3RecordCollection([])
end

# read an entire gff3 file, omitting any non-feature records: 
function GFF3RecordCollection(reader::GFF3.Reader)
    result = GFF3RecordCollection()
    for record in reader
        if BioCore.isfilled(record) && GFF3.isfeature(record) # omit other record types
            push!(result, GFF3Record(record))
        end
    end
    checkgff3input(result)
    return result
end

function GFF3RecordCollection(file::AbstractString; index=:auto)
    reader = GFF3.Reader(file; index=index)
    result = GFF3RecordCollection(reader)
    close(reader)
    checkgff3input(result)
    return result
end

function GFF3RecordCollection(input::IO; index=nothing)
    reader = GFF3.Reader(input; index=index)
    result = GFF3RecordCollection(reader)
    close(reader)
    checkgff3input(result)
    return result
end

function Base.length(records::GFF3RecordCollection)
    return length(records.records)
end

function Base.push!(records::GFF3RecordCollection, record::GFF3Record)
    return push!(records.records, record)
end

function Base.getindex(records::GFF3RecordCollection, index::Integer)
    return getindex(records.records, index)
end

Base.firstindex(records::GFF3RecordCollection) = 1
Base.lastindex(records::GFF3RecordCollection) = length(records)
Base.eachindex(records::GFF3RecordCollection) = Base.OneTo(lastindex(records))

function Base.isempty(records::GFF3RecordCollection)
    return isempty(records.records)
end

function Base.iterate(records::GFF3RecordCollection)
    if isempty(records)
        return nothing
    else
        i = firstindex(records)
        return getindex(records, i), i + 1
    end
end

function Base.iterate(records::GFF3RecordCollection, i::Int)
    if i > lastindex(records)
        return nothing
    else
        return getindex(records, i), i + 1
    end
end

function Base.eltype(::Type{GFF3RecordCollection})
    return GFF3Record
end

function Base.deleteat!(records::GFF3RecordCollection, i::Int)
    return deleteat!(records.records, i)
end

function Base.popat!(records::GFF3RecordCollection, i::Int)
    return popat!(records.records, i)
end

function Base.pop!(records::GFF3RecordCollection)
    return popat!(records, lastindex(records))
end

function Base.filter(f, records::GFF3RecordCollection)
    return filter(f, records.records)
end

function Base.filter!(f, records::GFF3RecordCollection)
    return filter!(f, records.records)
end

function Base.convert(::Type{DataFrames.DataFrame}, records::GFF3RecordCollection)
    result = DataFrames.DataFrame(seqid = Union{String, Missing}[], source = Union{String, Missing}[], 
                                  type = Union{String, Missing}[], start = Union{Int64, Missing}[], 
                                  end_ = Union{Int64, Missing}[], score = Union{Float64, Missing}[], 
                                  strand = Union{GenomicFeatures.Strand, Missing}[], phase = Union{Char, Missing}[], 
                                  attributes = Union{GFF3Attributes, Missing}[])

    for record in records
        if BioCore.isfilled(record)
            push!(result, (record.seqid, record.source, record.type, record.start, record.end_, 
                           record.score, record.strand, record.phase, record.attributes))
        end
    end
    return result
end

function DataFrames.DataFrame(records::GFF3RecordCollection)
    return convert(DataFrames.DataFrame, records)
end

################################################################################################################
###################################### checking input data for errors ##########################################
################################################################################################################

function checkstartend(record::GFF3Record)
    if record.start > record.end_ 
        error("Sequence start must be before or same as end.")
    end
end

function checkstartend(records::GFF3RecordCollection)
    for record in records
       checkstartend(record)
    end
    return
end

function checkcdsphase(record::GFF3Record)
    if record.type == "CDS" && ismissing(record.phase)
        error("Feature type \"CDS\" requires a phase.")
    end
end

function checkcdsphase(records::GFF3RecordCollection)
    for record in records
        checkcdsphase(record)
    end
    return
end

######### check the attributes #################################################################################
# only alias, parent, dbxref and ontologyterm allow multiple values (checked in GFF3Attributes constructor)

function allids(records::GFF3RecordCollection)
    allids = []
    for record in records
        if !ismissing(record.attributes) && !ismissing(record.attributes.id)
            push!(allids, record.attributes.id)
        end
    end
    return allids
end

function allparents(records::GFF3RecordCollection)
    parents = []
    for record in records
        if !ismissing(record.attributes) && !ismissing(record.attributes.parent) && !(record.attributes.parent in parents)
            push!(parents, record.attributes.parent...)
        end
    end
    return parents
end

function checkparents(records::GFF3RecordCollection)
    ids = allids(records)
    parents = allparents(records)
    for parentid in parents
        if !(parentid in ids)
            error("Parent $parentid must be found within scope of file")
        end
    end
    return
end

function alltargets(records::GFF3RecordCollection)
    targets = []
    for record in records
        if !ismissing(record.attributes) && !ismissing(record.attributes.target) && !(record.attributes.target in targets)
            push!(targets, record.attributes.target)
        end
    end
    return targets
end

function checktargets(records::GFF3RecordCollection)
    ids = allids(records)
    targets = alltargets(records)
    for t in targets
        if isnothing(match(r"^.+ [0-9]+ [0-9]+( [+-])?$", t)) 
            error("Target $(t) has wrong format.")
        end
        target = split(t, ' ')
        if !(target[1] in ids) 
            error("Target $(target[1]) must be found within scope of file.")
        end
    end
    return
end

function gap_hastarget(record::GFF3Record)
    if !ismissing(record.attributes) && !ismissing(record.attributes.gap) && ismissing(record.attributes.target)
        return false
    else
        return true
    end
end

function checkgaps(records::GFF3RecordCollection)
    pattern = r"(D[0-9]+|F[0-9]+|I[0-9]+|M[0-9]+|R[0-9]+)( (D[0-9]+|F[0-9]+|I[0-9]+|M[0-9]+|R[0-9]+))*"
    for record in records
        if !gap_hastarget(record) 
            error("Gap attribute requires an alignment target.")
        end
        if !ismissing(record.attributes) && !ismissing(record.attributes.gap) && isnothing(match(pattern, record.attributes.gap))
            error("Gap $(record.attributes.gap) has wrong format.")
        end
    end
    return
end

function allderivesfrom(records::GFF3RecordCollection)
    derivesfrom = []
    for record in records
        if !ismissing(record.attributes) && !ismissing(record.attributes.derivesfrom)
            push!(derivesfrom, record.attributes.derivesfrom)
        end
    end
    return derivesfrom
end

function checkderivesfrom(records::GFF3RecordCollection)
    ids = allids(records)
    derivesfrom = allderivesfrom(records)
    for originid in derivesfrom
        if !(originid in ids) 
            error("Derives_from origin $originid must be found within scope of file")
        end
    end
    return
end

######### the whole check-up: ##################################################################################
function checkgff3input(records::GFF3RecordCollection)
    checkstartend(records)
    checkcdsphase(records)
    checkparents(records)
    checktargets(records)
    checkgaps(records)
    checkderivesfrom(records)
end

TESTCOLLECTION = GFF3RecordCollection(TESTFILES[1]; index=nothing)
TESTCOLLECTION_empty = GFF3RecordCollection(TESTFILES[2]; index=nothing)
TESTCOLLECTION_bgz = GFF3RecordCollection(TESTFILES[3])

DataFrames.DataFrame(TESTCOLLECTION)
