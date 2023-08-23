# Provides I/O utilites for working with entire GFF3-files.
# Collection of zero or more GFF3Records

"""
    struct GFF3RecordCollection

Datastructure that holds one or more `GFF3Record`s read from a GFF3-file or entered as a vector of `GFF3Record`s.
- `records::Vector{GFF3Record}`


    GFF3RecordCollection(records::AbstractVector{GFF3Record}; extendedchecks::Bool=false)

Construct a new collection from a vector of `GFF3Record`s.
The function verifies the validity of the records in `records`. By default, this includes only start, end and phase.
If `extendedchecks=true` is set, run all available tests. This may take some time.
"""
mutable struct GFF3RecordCollection
    records::Vector{GFF3Record}

    function GFF3RecordCollection(records::AbstractVector{GFF3Record}; extendedchecks::Bool=false)
        checkgff3input(records; extended=extendedchecks)
        return new(records)
    end
end

"""
    GFF3RecordCollection()

Empty initialization.
"""
function GFF3RecordCollection()
    return GFF3RecordCollection([])
end

# read an entire gff3 file, omitting any non-feature records:
"""
    GFF3RecordCollection(reader::GFF3.Reader; extendedchecks::Bool=false)

Read a GFF3-file and save all non-empty feature records to a new `GFF3RecordCollection`.
The function verifies the validity of the input data. By default, this includes only start, end and phase.
If `extendedchecks=true` is set, run all available tests. This may take some time.
"""
function GFF3RecordCollection(
    reader::GFF3.Reader;
    extendedchecks::Bool=false
)
    result = []
    for record in reader
        if BioCore.isfilled(record) && GFF3.isfeature(record) # omit other record types
            push!(result, GFF3Record(record))
        end
    end
    return GFF3RecordCollection(Vector{GFF3Record}(result); extendedchecks=extendedchecks)
end

"""
    GFF3RecordCollection(file::AbstractString; index=:auto, extendedchecks::Bool=false)

    GFF3RecordCollection(input::IO; index=nothing, extendedchecks::Bool=false)

Create a GFF3.Reader for a GFF3-file with the provided data source.
If the file path ends with ".bgz", the function will try to find a tabix index file ".tbi" and read it.
Save all non-empty feature records to a new `GFF3RecordCollection`.
The function verifies the validity of the input data. By default, this includes only start, end and phase.
If `extendedchecks=true` is set, run all available tests. This may take some time.
"""
function GFF3RecordCollection(
    file::AbstractString;
    index=:auto,
    extendedchecks::Bool=false
)
    reader = GFF3.Reader(file; index=index)
    result = GFF3RecordCollection(reader; extendedchecks=extendedchecks)
    close(reader)
    return result
end

function GFF3RecordCollection(
    input::IO;
    index=nothing,
    extendedchecks::Bool=false
)
    reader = GFF3.Reader(input; index=index)
    result = GFF3RecordCollection(reader; extendedchecks=extendedchecks)
    close(reader)
    return result
end

function Base.length(records::GFF3RecordCollection)
    return length(records.records)
end

function Base.push!(
    records::GFF3RecordCollection,
    record::GFF3Record
)
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
    filteredrecords = filter(f, records.records)
    return GFF3RecordCollection(filteredrecords)
end

function Base.filter!(f, records::GFF3RecordCollection)
    filter!(f, records.records)
    return records
end

function Base.convert(
    ::Type{DataFrames.DataFrame},
    records::GFF3RecordCollection
)
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

# these are used in the GFF3RecordCollection constructors.

"""
    checkstartend(records::AbstractVector{GFF3Record})

Verify that each record in `records` has valid `start` and `end`.
"""
function checkstartend(records::AbstractVector{GFF3Record})
    for record in records
       checkstartend(record)
    end
    return
end

"""
    checkcdsphase(records::AbstractVector{GFF3Record})

Verify that each record in `records` has a `phase` if required by their `type`.
"""
function checkcdsphase(records::AbstractVector{GFF3Record})
    for record in records
        checkcdsphase(record)
    end
    return
end

######### check the attributes #################################################################################
# only alias, parent, dbxref and ontologyterm allow multiple values (checked in GFF3Attributes constructor)

"""
    allids(records::AbstractVector{GFF3Record})

Returns a vector that contains all attribute IDs found in `records`.
IDs can appear multiple times.
"""
function allids(records::AbstractVector{GFF3Record})
    allids = []
    for record in records
        if !ismissing(record.attributes) && !ismissing(record.attributes.id)
            push!(allids, record.attributes.id)
        end
    end
    return allids
end

"""
    allids(records::GFF3RecordCollection)

Returns a vector that contains all attribute IDs found in `records`.
IDs can appear multiple times.
"""
function allids(records::GFF3RecordCollection)
    return allids(records.records)
end

function allparents(
    records::AbstractVector{GFF3Record};
    onlyunique::Bool=true
) # HERE!!
    parents = []
    for record in records
        if onlyunique
            if !ismissing(record.attributes) && !ismissing(record.attributes.parent) && !(record.attributes.parent in parents)
                push!(parents, record.attributes.parent...)
            end
        else
            push!(parents, record.attributes.parent...)
        end
    end
    return parents
end

function allparents(records::GFF3RecordCollection)
    return allparents(records.records)
end

"""
    checkparents(records::AbstractVector{GFF3Record})

Verify that the Parent attributes of every record in `records` correspond to ID attributes.
"""
function checkparents(records::AbstractVector{GFF3Record})
    ids = allids(records)
    parents = allparents(records)
    for parentid in parents
        if !(parentid in ids)
            error("Parent $parentid must be found within scope of file")
        end
    end
    return
end

function alltargets(records::AbstractVector{GFF3Record})
    targets = []
    for record in records
        if !ismissing(record.attributes) && !ismissing(record.attributes.target) && !(record.attributes.target in targets)
            push!(targets, record.attributes.target)
        end
    end
    return targets
end

function alltargets(records::GFF3RecordCollection)
    return alltargets(records.records)
end

"""
    checktargets(records::AbstractVector{GFF3Record})

Verify that the Target attributes of every record in `records` correspond to ID attributes.
"""
function checktargets(records::AbstractVector{GFF3Record})
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

"""
    checkgaps(records::AbstractVector{GFF3Record})

Verify that the Gap attributes of every record in `records` match the required pattern.
Verify that each record which specifies a `gap` also specifies a `target`.
"""
function checkgaps(records::AbstractVector{GFF3Record})
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

function allderivesfrom(records::AbstractVector{GFF3Record})
    derivesfrom = []
    for record in records
        if !ismissing(record.attributes) && !ismissing(record.attributes.derivesfrom)
            push!(derivesfrom, record.attributes.derivesfrom)
        end
    end
    return derivesfrom
end

function allderivesfrom(records::GFF3RecordCollection)
    return allderivesfrom(records.records)
end

"""
    checkderivesfrom(records::AbstractVector{GFF3Record})

Verify that the Derives_from attributes of every record in `records` correspond to ID attributes.
"""
function checkderivesfrom(records::AbstractVector{GFF3Record})
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
"""
    checkgff3input(records::AbstractVector{GFF3Record}; extended::Bool=false)

Check start, end and phase for each record.
If `extended=true` is set, run all available tests on the records in `records`.
This may take some time.
"""
function checkgff3input(
    records::AbstractVector{GFF3Record};
    extended::Bool=false
)
    for record in records
        checkstartend(record)
        checkcdsphase(record)
    end
    if extended
        checkparents(records)
        checktargets(records)
        checkgaps(records)
        checkderivesfrom(records)
    end
end
