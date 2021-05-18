# Collection of zero or more GFF3Records

include("gffrecord.jl")

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

function checkstartend(records::GFF3RecordCollection)
    for record in records
       checkstartend(record)
    end
    return
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

