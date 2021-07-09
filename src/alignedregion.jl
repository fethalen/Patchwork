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
        alignment = BioAlignments.PairwiseAlignment(BioSequences.translate(result.querysequence),
            result.subjectsequence, result.cigar)
        return new(alignment, result.subjectstart, result.subjectend, result.queryid,
                   result.querystart, result.queryend, result.queryframe)
    end
end

function Base.show(io::IO, region::AlignedRegion)
    compact = get(io, :compact, false)

    if compact
        print(io, "query name: ", region.queryid)
    else
        println(io, "query name: ", *(otupart(region.queryid), '@', sequencepart(region.queryid)))
        println(io, "length: ", length(region))
        println(io, "subject interval: ", region.subjectfirst, " -> ", region.subjectlast)
        println(io, "query interval: ", region.queryfirst, " -> ", region.querylast)
        println(io, "queryframe: ", region.queryframe)
        println(io, "alignment: \n", region.pairwisealignment)
    end
end

function Base.show(io::IO, ::MIME"text/plain", region::AlignedRegion)
    println(io, typeof(region))
end

showrange(region::AlignedRegion) = println(region.subjectfirst, ' ', region.subjectlast)

"""
    compare(a:AlignedRegion, b:AlignedRegion)

Display the (translated) regions side-by-side.
TODO: Complete unfinished comparison function.
"""
function compare(a::AlignedRegion, b::AlignedRegion)
    subjectfirst = min(a.subjectfirst, b.subjectfirst)
    subjectlast = max(a.subjectlast, b.subjectlast)
    aprefix = repeat('-', a.subjectfirst - subjectfirst)
    bprefix = repeat('-', b.subjectfirst - subjectfirst)
    apostfix = repeat('-', subjectlast - a.subjectlast)
    bpostfix = repeat('-', subjectlast - b.subjectlast)
    println(aprefix, translate(a.record.sequencedata), apostfix)
    println(bprefix, translate(b.record.sequencedata), bpostfix)
end

"""
    length(region:AlignedRegion)

Returns the length of this sequence region (without gaps).
"""
function Base.length(region::AlignedRegion)
    return 1 + region.subjectlast - region.subjectfirst
end

"""
    leftposition(region::AlignedRegion)

Return the leftmost position of `region`.
"""
function leftposition(region::AlignedRegion)
    return region.subjectfirst
end

"""
    rightposition(region::AlignedRegion)

Return the rightmost position of `region`.
"""
function rightposition(region::AlignedRegion)
    return region.subjectlast
end

"""
    interval(region::AlignedRegion)

Return the subjectfirst and subjectlast positions of `region` as a tuple (`(subjectfirst, subjectlast)`).
"""
function interval(region::AlignedRegion)::Tuple
    return (leftposition(region), rightposition(region))
end

"""
    interval(region::AlignedRegion)

Print the subjectfirst and subjectlast positions of `region` (like so: subjectfirst -> subjectlast).
"""
function showinterval(region::AlignedRegion)
    println(leftposition(region), " -> ", rightposition(region))
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
    return leftposition(a) <= leftposition(b)
end

"""
    precedes(a::AlignedRegion, b::AlignedRegion)

Check whether the interval in `a` entirely precedes that of `b`.
"""
function precedes(a::AlignedRegion, b::AlignedRegion)
    return rightposition(a) < leftposition(b)
end

"""
    isoverlapping(a::AlignedRegion, b::AlignedRegion)

Returns `true` if the interval of `a` overlaps that of `b`.
"""
function isoverlapping(a::AlignedRegion, b::AlignedRegion)
    return leftposition(a) <= rightposition(b) &&
           leftposition(b) <= rightposition(a)
end

"""
    samerange(a::AlignedRegion, b::AlignedRegion)

Returns `true` if the interval of `a` is the same as that of `b`.
"""
function samerange(a::AlignedRegion, b::AlignedRegion)
    return leftposition(a) == rightposition(b) &&
           leftposition(b) == rightposition(a)
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
    if leftposition(a) < leftposition(b)
        sequence = *(a.record.sequencedata[beforeoverlap(a, b)],
                     bestpick(a, b).record.sequencedata[overlap(a, b)],
                     b.record.sequencedata[afteroverlap(a, b)])
    elseif leftposition(b) < leftposition(a)
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
    return max(leftposition(a), leftposition(b))
end

"""
    rightintersection(a, b)

Returns the rightmost position at which the two regions `a` and `b` meet.
"""
function rightintersection(a::AlignedRegion, b::AlignedRegion)
    return min(rightposition(a), rightposition(b))
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
function overlap(a::AlignedRegion, b::AlignedRegion)
    return (leftintersection(a, b), rightintersection(a, b))
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
    if leftposition(a) == leftposition(b)
        error("Region $a and $b starts at the same position")
    elseif leftposition(a) < leftposition(b)
        return (leftposition(a), leftintersection(a,b) - 1)
    else
        return (leftposition(b), leftintersection(a,b) - 1)
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
    if rightposition(a) == rightposition(b)
        error("Region $a and $b end at the same position")
    elseif rightposition(a) > rightposition(b)
        return (rightintersection(a,b) + 1, rightposition(a))
    elseif rightposition(b) > rightposition(a)
        return (rightintersection(a,b) + 1, rightposition(b))
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
           leftposition(a) <= leftposition(b) &&
           rightposition(a) >= rightposition(b)
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
    return leftposition(a) == leftposition(b) &&
           rightposition(a) == rightposition(b)
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
        return (leftposition(a), rightposition(a))
    elseif shadows(a, b)
        return (leftposition(a), rightposition(a))
    elseif shadows(b, a)
        return (leftposition(b), rightposition(b))
    else
        return (min(leftposition(a), leftposition(b)),
                max(rightposition(a), rightposition(b)))
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
    bestpick(a, b)

Choose between two `AlignedRegion`s, `a` and `b`, and retain the region with
the longest region. Or, if both regions cover the same region, return the
region with the highest amount of percent identity from the Bsubjectlast results.
"""
function bestpick(a::AlignedRegion, b::AlignedRegion)
    if equalrange(a, b)
        if a.percentidentity >= b.percentidentical
            return a
        else
            return b
        end
    else
        return longest(a, b)
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
