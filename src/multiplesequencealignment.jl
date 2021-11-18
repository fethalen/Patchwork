# Utilities and types for working with multiple sequence alignments (MSAs)

#include("sequencerecord.jl")

"""
Datastructures which holds a collection of `SequenceRecord`s.
"""
mutable struct MultipleSequenceAlignment
    name::AbstractString
    sequences::Vector{SequenceRecord}

    function MultipleSequenceAlignment()::MultipleSequenceAlignment
        return new("", Vector{SequenceRecord}())
    end

    function MultipleSequenceAlignment(name::String)::MultipleSequenceAlignment
        return new(name, Vector{SequenceRecord}())
    end

    function MultipleSequenceAlignment(alignments::Vector{SequenceRecord}
        )::MultipleSequenceAlignment
        return new("", alignments)
    end

    function MultipleSequenceAlignment(name::String, alignments::Vector{SequenceRecord}
        )::MultipleSequenceAlignment
        return new(name, alignments)
    end
end

function addalignment!(
    msa::MultipleSequenceAlignment,
    alignment::SequenceRecord
)::MultipleSequenceAlignment
    push!(msa.sequences, alignment)
    return msa
end

function removealignment!(
    msa::MultipleSequenceAlignment,
    alignment::SequenceRecord
)::MultipleSequenceAlignment
    delete!(msa.sequences, alignment)
    return msa
end

Base.push!(msa::MultipleSequenceAlignment, alignment::SequenceRecord) = push!(msa.sequences, alignment)
Base.deleteat!(msa::MultipleSequenceAlignment, index::Int64) = deleteat!(msa.sequences, index)

Base.length(msa::MultipleSequenceAlignment) = length(msa.sequences)
Base.firstindex(msa::MultipleSequenceAlignment) = 1
Base.lastindex(msa::MultipleSequenceAlignment) = lastindex(msa.sequences)
Base.isempty(msa::MultipleSequenceAlignment) = isempty(msa.sequences)

function Base.getindex(
    msa::MultipleSequenceAlignment,
    index::Int
)::SequenceRecord
    return msa.sequences[ind]
end

function Base.setindex!(
    msa::MultipleSequenceAlignment,
    sequence::SequenceRecord,
    index::Int
)
    msa.sequences[index] = sequence
end

function Base.iterate(msa::MultipleSequenceAlignment)
    isempty(msa) && return nothing
    i = firstindex(msa)
    return getindex(msa, i), i + 1
end

function Base.iterate(
    msa::MultipleSequenceAlignment,
    i::Int
)
    i > lastindex(msa) && return nothing
    return getindex(msa, i), i + 1
end

function Base.sort(msa::MultipleSequenceAlignment, bysequence::Bool=true)
    isempty(msa) && return msa
    if bysequence
        order = sortperm(map(record -> (record.sequencedata), msa))
    else
        order = sortperm(map(record -> (record.id.id), msa))
    end
    sortedmsa = MultipleSequenceAlignment(msa.name)
    for index in order
        addalignment!(sortedmsa, msa[index])
    end
    return sortedmsa
end

function countgaps(msa::MultipleSequenceAlignment)::Int
    return map(alignment -> length(alignment) - length(ungap(alignment)), msa.sequences) |>
        gapcount -> reduce(+, gapcount)
end

"""
    hasgaps(msa::MultipleSequenceAlignment)

Returns `true` if there are gap characters in `msa`.
"""
hasgaps(msa::MultipleSequenceAlignment)::Bool = countgaps(msa) > 0

"""
    ungap(msa::MultipleSequenceAlignment)

Returns a new `MultipleSequenceAlignment` object where gap characters within
`msa` are removed.
"""
function BioSequences.ungap(msa::MultipleSequenceAlignment)::MultipleSequenceAlignment
    return MultipleSequenceAlignment(
        msa.name,
        map(alignment -> ungap(alignment), msa.sequences))
end

"""
    ungap!(msa::MultipleSequenceAlignment)

Removes all gap characters from `msa`.
"""
function BioSequences.ungap!(msa::MultipleSequenceAlignment
    )::MultipleSequenceAlignment
    map(alignment -> ungap!(alignment), msa.sequences)
    return msa
end

"""
    otus(msa::MultipleSequenceAlignment)

Returns an `Array` of operational taxonomical units (OTUs; typically species)
in `msa`.
"""
function otus(msa::MultipleSequenceAlignment)::Array{String,1}
    otulist = []
    for alignment in msa.sequences
        push!(otulist, alignment.otu)
    end
    return unique!(otulist)
end

"""
    otufrequencies(msa::MultipleSequenceAlignment)

Returns a `Dict` object, showing how often each operational taxonomical units
(OTU) within `msa` occurs.
"""
function otufrequencies(msa::MultipleSequenceAlignment)::Dict{String, Int}
    frequencies = Dict{String, Int}()
    for alignment in msa.sequences
        otu = alignment.otu
        if haskey(frequencies, otu)
            frequencies[otu] += 1
        else
            frequencies[otu] = 1
        end
    end
    return frequencies
end

"""
    countotus(msa::MultipleSequenceAlignment)

Returns how many unique operational taxonomical units (OTUs) there are within
`msa`.
"""
countotus(msa::MultipleSequenceAlignment)::Int = length(otus(msa))

function coverage(msa::MultipleSequenceAlignment)::Float64
    length(msa) == 0 && return 0.0
    return map(alignment -> coverage(alignment), msa.sequences) |>
    alignment -> reduce(+, alignment) / length(msa)
end

"""
    equal_length(msa::MultipleSequenceAlignment)::Bool

Returns True if all of the sequences within the provided alignmment are of
equal length.
"""
function equal_length(msa::MultipleSequenceAlignment)::Bool
    return map(sequence -> length(sequence), msa.sequences) |>
    sequence_length -> unique(sequence_length) |>
    unique_seqlengths -> length(unique_seqlengths) == 1
end

"""
    gapmatrix(msa::MultipleSequenceAlignment)::Array{Int,2}

Returns a 2-dimensional array where each element represents the presence of a
gap character (1 indicate presence of a gap and 0 absence of a gap). Each row
represents an OTU and each column represents a position in a sequence.
"""
function gapmatrix(msa::MultipleSequenceAlignment)::Array{Int,2}
    equal_length(msa) || error("unequal ")
    alignmentlen = length(msa.sequences[1])
    gap_positions = zeros(Int, length(msa.sequences), alignmentlen)
    for (i, record) in enumerate(msa.sequences)
        for (j, position) in enumerate(record.sequencedata)
            if isgap(position)
                gap_positions[i,j] = 1
            end
        end
    end
    return gap_positions
end

"""
    gapmatrix(msa::MultipleSequenceAlignment)::Array{Int,2}

Returns an array where each element represents the proportion of gaps found in
all sequences in that position.
"""
function gapfrequencies(msa::MultipleSequenceAlignment)::Array{Float64,2}
    msalen = length(msa)
    msalen == 0 && error("cannot generate gap frequencies for an empty alignment")
    return map(gap_count -> gap_count / msalen, sum(gapmatrix(msa), dims=1))
end

"""
    mktemp_fasta(alignment; removehyphens)

Writes `alignment` to a temporary file and returns the path to that file.
Gaps are removed when `removehyphens=true`.
"""
function mktemp_fasta(
    alignment::MultipleSequenceAlignment;
    removehyphens::Bool=false
)::AbstractString
    path, io = mktemp()
    removehyphens && ungap!(alignment)

    for record in alignment.sequences
        write(io, *('>', record.id.id,
                    '\n', String(record.sequencedata), '\n'))
    end
    close(io)
    return path
end

function mktemp_fasta(
    alignment::AbstractString;
    removehyphens::Bool=false
)::AbstractString
    path, io = mktemp()
    msa = readmsa(alignment)
    removehyphens && ungap!(msa)

    for record in msa.sequences
        write(io, *('>', record.id.id,
                    '\n', String(record.sequencedata), '\n'))
    end
    close(io)
    return path
end

function pool(
    files::AbstractString...;
    name::AbstractString="", 
    removeduplicates::Bool=true
)::MultipleSequenceAlignment
    allsequences = Vector{SequenceRecord}()
    for file in files
        (isfastafile(file) || isfastqfile(file)) || error("Can only pool fasta or fastq files.")
        tmp = readmsa(file)
        append!(allsequences, (ungap(tmp)).sequences)
    end
    result = MultipleSequenceAlignment(name, allsequences)
    removeduplicates && remove_duplicates!(result)
    return result
end
