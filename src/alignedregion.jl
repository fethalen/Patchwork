# An AlignedRegion is a continuous sequence which is part of a larger
# alignment. Provides methods and types for working with the range and
# length of such regions.

import BioSequences
import BioAlignments

include("diamond.jl")
include("sequencerecord.jl")

"""
    struct AlignedRegion

- `record::SequenceRecord`: the sequence record associated with the interval.
- `subjectfirst::Int64`: the leftmost position.
- `subjectlast::Int64`: the rightmost position.
- `queryframe::Int64: the queryframe of the sequence (-3, -2, -1, 1, 2, or 3;
  0 = undefined).
- `percentidentity::Float64: percent identity in Bsubjectlast search results`
"""
struct AlignedRegion
    pairwisealignment::BioAlignments.PairwiseAlignment
    subjectfirst::Int64
    subjectlast::Int64
    queryid::SequenceIdentifier
    queryfirst::Int64
    querylast::Int64
    queryframe::Int64

    function AlignedRegion(result::DiamondSearchResult)
        if result.queryframe > 0
            alignment = BioAlignments.PairwiseAlignment(BioSequences.translate(result.querysequence),
                result.subjectsequence, result.cigar)
        else
            alignment = BioAlignments.PairwiseAlignment(BioSequences.translate(BioSequences.reverse_complement(result.querysequence)),
                result.subjectsequence, result.cigar)
        end
        return new(alignment, result.subjectstart, result.subjectend, result.queryid,
                   result.querystart, result.queryend, result.queryframe)
    end
end

function Base.show(io::IO, region::AlignedRegion)
    compact = get(io, :compact, false)

    if compact
        print(io, "query name: ", region.queryid.id)
    else
        println(io, "query name: ", region.queryid.id)
        println(io, "length: ", length(region))
        println(io, "subject interval: ", region.subjectfirst, " -> ", region.subjectlast)
        println(io, "query interval: ", region.queryfirst, " -> ", region.querylast)
        println(io, "query frame: ", region.queryframe)
        println(io, "alignment: \n", region.pairwisealignment)
    end
end

function Base.show(io::IO, ::MIME"text/plain", region::AlignedRegion)
    println(io, length(region), "aa long Aligned Region:")
    print(io, "query name: ", region.queryid.id)
end

"""
    query_isreverse(region::AlignedRegion)::Bool

Returns `true` if the first query position is larger than the last query position.
"""
function query_isreverse(region::AlignedRegion)::Bool
    return region.queryfirst > region.querylast
end

query2subject(region::AlignedRegion, i::Integer)::Integer = first(BioAlignments.seq2ref(region.pairwisealignment, i))
subject2query(region::AlignedRegion, i::Integer)::Integer = first(BioAlignments.ref2seq(region.pairwisealignment, i))

"""
    subject2fullquery(region::AlignedRegion, i::Integer)::UnitRange

Given a position, `i`, in the aligned part of the provided reference sequence,
return the corresponding _positions_ in the full query sequence.
"""
function subject2fullquery(region::AlignedRegion, i::Integer)::UnitRange
    querypos = subject2query(region, i)
    return query2fullquery(region, querypos)
end

"""
    fullsubject2subject(region::AlignedRegion, i::Integer)::Integer

Given a position, `i`, in the full reference sequence, return the corresponding position in
the aligned part of the subject sequence.
"""
function fullsubject2subject(region::AlignedRegion, i::Integer)::Integer
    region.subjectfirst <= i <= region.subjectlast || error("the provided integer, ", i,
    ", is outside the range of the provided region, ", region)
    return i - region.subjectfirst + 1
end

"""
    fullsubject2fullquery(region::AlignedRegion, i::Integer)::Integer

Given a position, `i`, in the full reference sequence, return the corresponding position in
the aligned part of the query sequence.
"""
function fullsubject2query(region::AlignedRegion, i::Integer)::Integer
    subjectpos = fullsubject2subject(region, i)
    return subject2query(region, subjectpos)
end

"""
    fullsubject2fullquery(region::AlignedRegion, i::Integer)::UnitRange

Given a position, `i`, in the full reference sequence, return the corresponding _positions_
in the full query sequence.
"""
function fullsubject2fullquery(region::AlignedRegion, i::Integer)::UnitRange
    querypos = fullsubject2query(region, i)
    return query2fullquery(region, querypos)
end

"""
    query2fullquery(region::AlignedRegion, i::Integer)::Integer

Given a position, `i`, in the aligned part of the query sequence, return the corresponding
_positions_ in the full part of the query sequence.
"""
function query2fullquery(region::AlignedRegion, i::Integer)::UnitRange
    if query_isreverse(region)
        full_querypos = query_rightposition(region) - ((i - 1) * 3)
        return full_querypos - 2:full_querypos
    else
        full_querypos = query_leftposition(region) + ((i - 1) * 3)
        return full_querypos:full_querypos + 2
    end
end

"""
    fullsubject_queryboundaries(region::AlignedRegion, range::UnitRange)::Tuple

Given a `range` in the full reference sequence, return the first and last positions in the
full query sequence, which correspond to the provided range. Because the query sequence
consist of nucleotide sequences, each position in the reference (subject) sequence
correspond to a range of positions. This function's _raison d'etre_ is thus to calculate
the first and the last position of the ranges, given the provided range. The results are
returned as a tuple. Note that the returned range may be inverted.
"""
function fullsubject_queryboundaries(region::AlignedRegion, range::UnitRange)::Tuple
    if query_isreverse(region)
        return (last(fullsubject2fullquery(region, first(range))),
                first(fullsubject2fullquery(region, last(range))))
    else
        return (first(fullsubject2fullquery(region, first(range))),
                last(fullsubject2fullquery(region, last(range))))
    end
end

"""
    subject_queryboundaries(region::AlignedRegion, range::UnitRange)::Tuple

Given a `range` in the aligned part of the reference sequence, return the first and last positions in the
full query sequence, which correspond to the provided range.
"""
function subject_queryboundaries(region::AlignedRegion, range::UnitRange)::Tuple
    if query_isreverse(region)
        return (last(subject2fullquery(region, first(range))),
                first(subject2fullquery(region, last(range))))
    else
        return (first(subject2fullquery(region, first(range))),
                last(subject2fullquery(region, last(range))))
    end
end

"""
    length(region:AlignedRegion)

Returns the length of this sequence region (without gaps).
"""
function Base.length(region::AlignedRegion)
    return 1 + region.subjectlast - region.subjectfirst - 1
end

Base.isempty(region::AlignedRegion) = length(region) < 1
Base.firstindex(region::AlignedRegion) = 1
Base.lastindex(region::AlignedRegion) = length(region)
Base.eachindex(region::AlignedRegion) = Base.OneTo(lastindex(region))

# TODO: Finish slicing functions
function Base.getindex(region::AlignedRegion, index::Integer)
    return getindex(region.records, index)
end

@inline function Base.iterate(region::AlignedRegion, i::Int = firstindex(regions))
    if i > lastindex(region)
        return nothing
    else
        return getindex(regions, i), i + 1
    end
end

"""
    query_leftposition(region::AlignedRegion)::Integer

Return the leftmost position of the query sequence in the provided `region`.
"""
function query_leftposition(region::AlignedRegion)::Integer
    return min(region.queryfirst, region.querylast)
end

"""
    query_rightposition(region::AlignedRegion)::Integer

Return the rightmost position of `region`.
"""
function query_rightposition(region::AlignedRegion)::Integer
    return max(region.queryfirst, region.querylast)
end

"""
    subject_leftposition(region::AlignedRegion)::Integer

Return the leftmost position of the subject sequence in the provided `region`.
"""
function subject_leftposition(region::AlignedRegion)::Integer
    return region.subjectfirst
end

"""
    subject_rightposition(region::AlignedRegion)::Integer

Return the rightmost position of the subject sequence in the provided `region`.
"""
function subject_rightposition(region::AlignedRegion)::Integer
    return region.subjectlast
end

"""
    subjectinterval(region::AlignedRegion)::UnitRange

Return the subjectfirst and subjectlast positions of `region` as a tuple (`(subjectfirst, subjectlast)`).
"""
function subjectinterval(region::AlignedRegion)::UnitRange
    return subject_leftposition(region):subject_rightposition(region)
end

"""
    queryinterval(region::AlignedRegion)::UnitRange

Returns the range of the query sequence in the provided `region`. Note that this returns
the leftmost to the rightmost position in the query sequence, i.e., the direction is lost.
"""
function queryinterval(region::AlignedRegion)::UnitRange
    return query_leftposition(region):query_rightposition(region)
end

"""
    queryframe(region::AlignedRegion)

Return the queryframe of `region`.
"""
function queryframe(region::AlignedRegion)
    return region.queryframe
end

"""
    isordered(a::AlignedRegion, b::AlignedRegion)

Check if two aligned regions are well-ordered. Two regions, `a` and `b` are
considered well-ordered if `leftposition(a) <= leftposition(b)`.
"""
function isordered(a::AlignedRegion, b::AlignedRegion)
    return subject_leftposition(a) <= subject_leftposition(b)
end

"""
    precedes(a::AlignedRegion, b::AlignedRegion)

Check whether the interval in `a` entirely precedes that of `b`.
"""
function precedes(a::AlignedRegion, b::AlignedRegion)
    return subject_rightposition(a) < subject_leftposition(b)
end

"""
    isoverlapping(a::AlignedRegion, b::AlignedRegion)

Returns `true` if the interval of `a` overlaps that of `b`.
"""
function isoverlapping(a::AlignedRegion, b::AlignedRegion)
    return subject_leftposition(a) <= subject_rightposition(b) &&
           subject_leftposition(b) <= subject_rightposition(a)
end

"""
    samerange(a::AlignedRegion, b::AlignedRegion)

Returns `true` if the interval of `a` is the same as that of `b`.
"""
function samerange(a::AlignedRegion, b::AlignedRegion)
    return subject_leftposition(a) == subject_rightposition(b) &&
           subject_leftposition(b) == subject_rightposition(a)
end

"""
    identifier(region::AlignedRegion, separator='@')

Returns the sequence identifier of `region` with the OTU and sequence ID
separated by `separator`.
"""
function identifier(region::AlignedRegion, separator::AbstractChar='@')
    return identifier(region.record, separator)
end

"""
    sameid(a, b)

Returns `true` if the sequence identifer of `a` and `b` are identical.
"""
function sameid(a::AlignedRegion, b::AlignedRegion)
    return identifier(a) == identifier(b)
end

"""
    samesequence(a, b)

Returns `true` if the sequence data of `a` and `b` are identical.
"""
function samesequence(a::AlignedRegion, b::AlignedRegion)
    return a.record.sequencedata == b.record.sequencedata
end

"""
    consensus(a, b)

Returns
"""
function consensus(a::AlignedRegion, b::AlignedRegion, skipcheck::Bool=false)
    !skipcheck && !isoverlapping(a, b) && error("region $a and $b are not overlapping")

    if !isoverlapping(a, b)
        # TODO: Combine non-overlapping regions by filling w. gaps
    elseif shadows(a, b)
        return a
    elseif shadows(b, a)
        return b
    elseif samerange(a, b)
        return bestpick(a, b)
    end
    leftmost, rightmost = totalrange(a, b)
    overlapping = overlap(a, b)
    preceding = beforeoverlap(a,b)
    afteroverlap = afteroverlap(a,b)
    # A sequence that doesn't "shadow" another sequence cannot start and end
    # at an earlier and a later position than that of the other sequence.
    if subject_leftposition(a) < subject_leftposition(b)
        sequence = *(a.record.sequencedata[beforeoverlap(a, b)],
                     bestpick(a, b).record.sequencedata[overlap(a, b)],
                     b.record.sequencedata[afteroverlap(a, b)])
    elseif subject_leftposition(b) < subject_leftposition(a)
        sequence = *(a.record.sequencedata[beforeoverlap(a, b)],
                     bestpick(a, b).record.sequencedata[overlap(a, b)],
                     b.record.sequencedata[afteroverlap(a, b)])
    end
    record = SequenceRecord(a.record.otu, a.record.identifier, sequence)
    identity = (a.percentidentical + b.percentidentical) / 2
    return AlignedRegion(record, leftmost, rightmost, identity)
end

"""
    leftintersection(a, b)

Returns the leftmost position at which the two regions `a` and `b` meet.
"""
function leftintersection(a::AlignedRegion, b::AlignedRegion)
    return max(subject_leftposition(a), subject_leftposition(b))
end

"""
    rightintersection(a, b)

Returns the rightmost position at which the two regions `a` and `b` meet.
"""
function rightintersection(a::AlignedRegion, b::AlignedRegion)
    return min(subject_rightposition(a), subject_rightposition(b))
end

"""
    overlap(a, b)

Returns the leftmost and the rightmost position of `a` and `b`.

Example 1: Returns `(5,8)` since this is the range in which `a` and `b` are
    overlapping.
a:  [------]
b:      [-------]
    |...|....|....|
    1   5    10   15
"""
function overlap(a::AlignedRegion, b::AlignedRegion)::UnitRange
    return leftintersection(a, b):rightintersection(a, b)
end

"""
    beforeoverlap(a, b)

Returns the leftmost and the rightmost position of `a` and `b`.

Example 1: Returns `(1,4)` since this is the range in which `a` and `b` are
    overlapping.
a:  [------]
b:      [-------]
    |...|....|....|
    1   5    10   15
"""
function beforeoverlap(a::AlignedRegion, b::AlignedRegion)
    if subject_leftposition(a) == subject_leftposition(b)
        error("Region $a and $b starts at the same position")
    elseif subject_leftposition(a) < subject_leftposition(b)
        return (subject_leftposition(a):leftintersection(a,b) - 1)
    else
        return (subject_leftposition(b):leftintersection(a,b) - 1)
    end
end

"""
    afteroverlap(a, b)

Returns the leftmost and the rightmost position of `a` and `b`.

Example 1: Returns `(9,13)` since this is the range *after* the overlapping
    region of `a` and `b`.
a:  [------]
b:      [-------]
    |...|....|....|
    1   5    10   15
"""
function afteroverlap(a::AlignedRegion, b::AlignedRegion)
    if subject_rightposition(a) == subject_rightposition(b)
        error("Region $a and $b end at the same position")
    elseif subject_rightposition(a) > subject_rightposition(b)
        return (rightintersection(a,b) + 1:subject_rightposition(a))
    elseif subject_rightposition(b) > subject_rightposition(a)
        return (rightintersection(a,b) + 1:subject_rightposition(b))
    end
end

"""
    shadows(a, b)

An `AlignedRegion`, `a`, shadows another `AlignedRegion`, `b`, if and only if
`a` is longer than `b` and `a` starts at an earlier, or the same, position as
that of `b` and `b` ends before `a` does.

Example 1: `true` since `a` is shadowing `b`.
a: [--------------]
b: [---------]

Example 2: `false` since `b` is shadowing `a` and not the other way.
a: [---------]
b: [--------------]

Example 3: `true` since `a` is shadowing `b`.
a: [--------------]
b:   [----]

Example 4: `false` since `a` is of equal length to `b`.
a: [--------------]
b: [--------------]
"""
function shadows(a::AlignedRegion, b::AlignedRegion)
    return length(a) > length(b) &&
           subject_leftposition(a) <= subject_leftposition(b) &&
           subject_rightposition(a) >= subject_rightposition(b)
end

"""
    equalrange(a, b)

Check whether `a` and `b` starts and stops at the same position.

Example 1: `true` since `a` starts and stop at the same position as `b`.
a: [-------]
b: [-------]
    |...|....|....|
    1   5    10   15

Example 2: `false` since `b` begins after `a`.
a:  [-------]
b:    [-----]
    |...|....|....|
    1   5    10   15
"""
function equalrange(a::AlignedRegion, b::AlignedRegion)
    return subject_leftposition(a) == subject_leftposition(b) &&
           subject_rightposition(a) == subject_rightposition(b)
end

"""
    totalrange(a, b)

Returns the total range covered by the two regions `a` and `b`, if and only
if these two regions overlap one-another.

Example 1: `10` since there are `10` positions covered by regions `a` and `b`.
a:    [------]
b:  [----]
    |...|....|....|
    1   5    10   15
"""
function totalrange(a::AlignedRegion, b::AlignedRegion)
    if !isoverlapping(a, b)
        error("Cannot do total range of non-overlapping regions, $a and $b")
    elseif equalrange(a, b)
        return (subject_leftposition(a), subject_rightposition(a))
    elseif shadows(a, b)
        return (subject_leftposition(a), subject_rightposition(a))
    elseif shadows(b, a)
        return (subject_leftposition(b), subject_rightposition(b))
    else
        return (min(subject_leftposition(a), subject_leftposition(b)),
                max(subject_rightposition(a), subject_rightposition(b)))
    end
end

"""
    longest(a, b)

Example 1: Returns `a` because `a` is longer than `b`.
a: [-----------]
b: [-----]

Example 2: Returns `a` because `a` and `b` are of the same length.
a: [-----------]
b: [-----------]

Example 3: Returns `b` because `b` is longer than `a`.
a: [--------]
b: [-----------]
"""
function longest(a::AlignedRegion, b::AlignedRegion)
    if length(a) >= length(b)
        return a
    else
        return b
    end
end

"""
    isnucleotide(region)

Returns true if this region's `record` consists of nucleotides.
"""
function isnucleotide(region::AlignedRegion)
    return isnucleotide(region.record)
end

"""
    isaminoacid(region)

Returns true if this record's `record` consists of amino acids.
"""
function isaminoacid(region::AlignedRegion)
    return isaminoacid(region.record)
end

function BioSequences.translate(region::AlignedRegion)
    isnucleotide(region) || error("cannot translate non-nucleic region: ", region)
    return BioSequences.translate(region.record)
end
