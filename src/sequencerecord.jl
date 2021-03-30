using BioSequences

"""
A `SequenceRecord` represents a record in a multiple sequence alignment
(MSA). A sequence record has an operational taxonomic unit (OTU; usually a
species), a sequence identifier, and the sequence data itself.
"""
struct SequenceRecord
    otu::AbstractString
    identifier::AbstractString
    sequencedata::LongSequence
end

function SequenceRecord()
    return SequenceRecord("", "", LongSequence())
end

function SequenceRecord(otuid::String, sequencedata::LongSequence, separator::Char)
    if count(separator in otuid) > 1
        error(*("multiple separators (", separator, ") found in ", otuid))
    elseif count(separator in otuid) == 0
        error(*("couldn't find separator (", separator, ") in ", otuid))
    end

    otu, identifier = split(otuid, separator)
    return SequenceRecord(otu, identifier, sequencedata)
end

Base.length(alignment::SequenceRecord) = length(alignment.sequencedata)

function BioSequences.ungap(alignment::SequenceRecord)::SequenceRecord
    return SequenceRecord(alignment.otu, alignment.identifier,
                           ungap(alignment.sequencedata))
end

function BioSequences.ungap!(alignment::SequenceRecord)::SequenceRecord
    ungap!(alignment.sequencedata)
    return alignment
end

countgaps(alignment::SequenceRecord) = length(alignment) - length(ungap(alignment))
hasgaps(alignment::SequenceRecord) = countgaps(alignment) > 0

function Base.iterate(alignment::SequenceRecord)
    isempty(alignment.sequencedata) && return nothing
    return iterate(alignment, 1)
end

function Base.iterate(alignment::SequenceRecord, state=1::Int)
    state > lastindex(alignment.sequencedata) && return nothing
    state += 1
    return (alignment.sequencedata, state)
end

function missingdata(alignment::SequenceRecord)::Float64
    alignmentlength = length(alignment)
    alignmentlength == 0 && return 0
    return 1 - length(ungap(alignment)) / alignmentlength
end

function coverage(alignment::SequenceRecord)::Float64
    alignmentlength = length(alignment)
    alignmentlength == 0 && return 0
    return length(ungap(alignment)) / alignmentlength
end

function gappositions(sequence::SequenceRecord)::Array{Int,1}
    gaprow = zeros(Int, length(sequence))
    for (index, position) in enumerate(sequence.sequencedata)
        if isgap(position)
            gaprow[index] = 1
        end
    end
    return gaprow
end

"""
    nongap_range(sequence)

For an aligned `sequence`, determine where the gap positions in that
`sequence` starts and stops.

Example 1. Returns `(4,13)` since the other characters in the `sequence` are
gap characters.

    ---REPIGKFHIQ--
    |...|....|....|
    1   5    10   15
"""
function nongap_range(sequence::LongSequence)
    first = 0
    last = 0
    for (index, position) in enumerate(sequence)
        if first == 0 && ! isgap(position)
            first = index
        elseif last == 0 && first != 0 && isgap(position)
            last = index - 1
        end
    end
    return (first, last)
end

function nongap_range(record::SequenceRecord)
    return nongap_range(record.sequencedata)
end

"""
    has_compoundregions(sequence)

Check whether an aligned `sequence` has *compound regions*. An aligned
sequence is composed of multiple regions if there are multiple non-overlapping
non-gap regions contained within that aligned sequence.

**Example 1.** This aligned sequence is a compound region because there are
two non-overlapping regions: `REPIGKFHIQ` and `FHIIGKPR`.

    ---REPIGKFHIQ----FHIIGKPR-----
    |...|....|....|....|....|....|
    1   5    10   15   20   25   30
"""
function has_compoundregions(sequence::SequenceRecord)
    leftmost, rightmost = nongap_range(sequence)
    regionsize = rightmost - leftmost + 1
    return regionsize < length(sequence)
end

"""
    fillmissing(sequence, totalrange, coveredrange)

Add gap characters to `sequence` based on the `totalrange` of an alignment and
the `coveredrange` by that `sequence`.
"""
function fillmissing(sequence::LongSequence, totalrange::Integer,
    coveredrange::Tuple)::LongSequence
    leftmost, rightmost = coveredrange
    coveredpositions = rightmost - leftmost + 1
    if length(sequence) != coveredpositions
        error("""Number of positions in the range $coveredrange is not equal to
                 the sequence length (""", length(sequence), ')')
    end
    gapsbefore = ""
    gapsafter = ""

    if leftmost > 1
        gapsbefore = "-"^(leftmost - 1)
    end

    if rightmost < totalrange
        gapsafter = "-"^(totalrange - rightmost)
    end

    filledsequence = *(gapsbefore, string(sequence), gapsafter)

    if eltype(sequence) == AminoAcid
        return LongAminoAcidSeq(filledsequence)
    elseif eltype(sequence) == DNA
        return LongDNASeq(filledsequence)
    elseif eltype(sequence) == RNA
        return LongRNASeq(filledsequence)
    end
end

"""
    compoundregions(sequence)

Return all *compound regions* within `sequence`. An aligned sequence is
composed of multiple regions if there are multiple non-overlapping non-gap
regions contained within that aligned sequence.

**Example 1.** Returns `REPIGKFHIQ` and `FHIIGKPR`, because this sequence is
composed of these two compound regions.

    ---REPIGKFHIQ----FHIIGKPR-----
    |...|....|....|....|....|....|
    1   5    10   15   20   25   30

**Example 2.** Returns `NILLCTL`, `FHIIGKPR`, because this sequence is
composed of these two regions.

    NILLCTL-----YKCCGC--TEVECLGKCC
    |...|....|....|....|....|....|
    1   5    10   15   20   25   30
"""
function compoundregions(record::SequenceRecord)::Array{LongSequence{},1}
    regions = []
    seqlength = length(record)
    println(length(seqlength))
    for (leftmost, rightmost) in compoundranges(record)
        push!(regions,
              fillmissing(record.sequencedata[leftmost:rightmost],
                          seqlength,
                          (leftmost, rightmost)))
    end
    return regions
end

"""
    compoundranges(sequence)

Check whether an aligned `sequence` has *compound regions*. An aligned
sequence is composed of multiple regions if there are multiple non-overlapping
non-gap regions contained within that aligned sequence.

**Example 1.** Returns `[(4,13), (18,25)]` because `(4,13)` is where the first
sequence starts and ends and `(18,25)` is where the second sequence starts and
ends.

    ---REPIGKFHIQ----FHIIGKPR-----
    |...|....|....|....|....|....|
    1   5    10   15   20   25   30

**Example 2.** Returns `[(1,7), (13,18), (21,30)]` because that is where each
non-gap region within the aligned `sequence` starts and ends.

    NILLCTL-----YKCCGC--TEVECLGKCC
    |...|....|....|....|....|....|
    1   5    10   15   20   25   30
"""
function compoundranges(record::SequenceRecord)
    regionranges = []
    inregion = false
    leftmost = nothing
    for i in eachindex(record.sequencedata)
        if !inregion && !isgap(record.sequencedata[i])
            leftmost = i
            inregion = true
        elseif inregion && isgap(record.sequencedata[i])
            rightmost = i - 1
            push!(regionranges, (leftmost, rightmost))
            inregion = false
        end
    end
    return regionranges
end

"""
    nongap_ranges(sequence)

For an aligned `sequence`, determine all start and stop positions based on
each gap-region within this `sequence`. Unlike `non_gapranges(sequence)`,
this function determines the gaps for all possible ranges and doesn't assume
a singular range.

Example 1. Returns `[(4,6), (9,14)]` since the other characters in the
`sequence` are gap characters.

    ---REP--KFHIQK-
    |...|....|....|
    1   5    10   15
"""
function nongap_ranges(sequence::SequenceRecord)
    ranges = []
    tovisit = eachindex(sequence)

    for (index, position) in enumerate(sequence.sequencedata)
        if first == 0 && ! isgap(position)
            first = index
        elseif last == 0 && first != 0 && isgap(position)
            last = index - 1
        end
    end
    return (first, last)
end

"""
    identifier(sequence, separator='@')

Returns the identifier of `sequence`, including the OTU separated by
`separator`.
"""
function identifier(sequence::SequenceRecord, separator::AbstractChar='@')
    return *(sequence.otu, separator, sequence.identifier)
end
