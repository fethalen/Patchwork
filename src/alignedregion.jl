# An AlignedRegion is a continuous sequence which is part of a larger
# alignment. Provides methods and types for working with the range and
# length of such regions.

"""
    struct AlignedRegion

- `pairwisealignment::BioAlignments.PairwiseAlignment`: the alignment associated with the
  interval.
- `subjectfirst::Int64`: the first position in the reference sequence.
- `subjectlast::Int64`: the last position in the reference sequence.
- `queryid::SequenceIdentifier`: the identifier associated with the aligned region.
- `queryfirst::Int64`: the first position in the query sequence.
- `querylast::Int64`: the first position in the query sequence.
- `queryframe::Int64`: the queryframe of the sequence (-3, -2, -1, 1, 2, or 3;
  0 = undefined).
"""
struct AlignedRegion
    pairwisealignment::BioAlignments.PairwiseAlignment
    full_querysequence::LongDNA
    subjectfirst::Int64
    subjectlast::Int64
    queryid::SequenceIdentifier
    queryfirst::Int64
    querylast::Int64
    queryframe::Int64

    function AlignedRegion()
        emptyaln = BioAlignments.PairwiseAlignment(BioSequences.LongAA(),
            BioSequences.LongAA(), "")
        return new(emptyaln, BioSequences.LongDNA(), 0, -1, SequenceIdentifier(), 0, -1, 0)
    end

    function AlignedRegion(alignment::BioAlignments.PairwiseAlignment, fullqseq::LongDNA,
        subjectfirst::Int64, subjectlast::Int64, queryid::SequenceIdentifier,
        queryfirst::Int64, querylast::Int64, queryframe::Int64)
        return new(alignment, fullqseq, subjectfirst, subjectlast, queryid, queryfirst,
            querylast, queryframe)
    end

    function AlignedRegion(result::DiamondSearchResult)
        alignment = BioAlignments.PairwiseAlignment(result.translated_querysequence,
            result.subjectsequence, cleancigar(result.cigar)) # should be clean from frameshift() call but nvm
        return new(alignment, result.full_querysequence, result.subjectstart,
                    result.subjectend, result.queryid, result.querystart, result.queryend,
                    result.queryframe)
    end
end

function BioAlignments.cigar(anchors::AbstractVector{BioAlignments.AlignmentAnchor})::String
    out = IOBuffer()
    if isempty(anchors)
        return ""
    end
    @assert anchors[1].op == BioAlignments.OP_START "Alignments must start with OP_START."
    for i in 2:lastindex(anchors)
        # `alnpos` is required when the version of BioAlignments.jl >= 2.1.0
        # if-else not necessary since BioAlignments >= 3.0.0 required in Project/Manifest
        #if packageversion("BioAlignments") >= 210
		    positions = anchors[i].alnpos - anchors[i-1].alnpos
        #else
        #    positions = max(anchors[i].seqpos - anchors[i-1].seqpos,
        #                    anchors[i].refpos - anchors[i-1].refpos)
        #end
        print(out, positions, anchors[i].op)
    end
    return String(take!(out))
end

function BioAlignments.cigar(alignment::BioAlignments.PairwiseAlignment)::String
    return BioAlignments.cigar(alignment.a.aln.anchors)
end


"""
    cleancigar(cigar::AbstractString)

Clean a CIGAR String from any frameshiftig operations `/` and `\\`, to be used in a
`BioAlignments.PairwiseAlignment`.
"""
function cleancigar(cigar::AbstractString)::String
    parts = collect(eachmatch(r"\d+[MIDNSHP=X]", cigar))
    buffer = IOBuffer()

    for part in parts
        print(buffer, part.match)
    end

    return String(take!(buffer))
end

function BioSequences.translate(region::AlignedRegion)
    isnucleotide(region) || error("cannot translate non-nucleic region: ", region)
    return BioSequences.translate(region.record)
end

function BioSequences.translate(
    sequence::BioSequences.LongDNA,
    cigar::AbstractString,
    frame::Int64=0 # not specified
)::BioSequences.LongAA
    return BioSequences.translate(frameshift(sequence, cigar, frame)[1])
end

function Base.show(
    io::IO,
    region::AlignedRegion
)
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

function Base.show(
    io::IO,
    ::MIME"text/plain",
    region::AlignedRegion
)
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
    subject2fullsubject(region::AlignedRegion, i::Integer)::UnitRange

Given a position, `i`, in the aligned part of the provided reference sequence,
return the corresponding _positions_ in the full reference sequence.
"""
function subject2fullsubject(
    region::Patchwork.AlignedRegion,
    i::Integer
)::Integer
    i <= lastindex(region) || error("the provided integer, ", i,
    ", is outside the range of the aligned part of the region, ", region)
    return i + region.subjectfirst - 1
end

"""
    subject2fullquery(region::AlignedRegion, i::Integer)::UnitRange

Given a position, `i`, in the aligned part of the provided reference sequence,
return the corresponding _positions_ in the full query sequence.
"""
function subject2fullquery(
    region::AlignedRegion,
    i::Integer
)::UnitRange
    subjectpos = subject2fullsubject(region, i)
    return fullsubject2fullquery(region, subjectpos)
end

"""
    fullsubject2subject(region::AlignedRegion, i::Integer)::Integer

Given a position, `i`, in the full reference sequence, return the corresponding position in
the aligned part of the subject sequence.
"""
function fullsubject2subject(
    region::AlignedRegion,
    i::Integer
)::Integer
    region.subjectfirst <= i <= region.subjectlast || error("the provided integer, ", i,
    ", is outside the range of the provided region, ", region)
    return i - region.subjectfirst + 1
end

function fullsubject2subject(
    region::AlignedRegion,
    range::UnitRange
)::UnitRange
    return fullsubject2subject(region, first(range)):fullsubject2subject(region, last(range))
end

"""
    fullsubject2fullquery(region::AlignedRegion, i::Integer)::Integer

Given a position, `i`, in the full reference sequence, return the corresponding position in
the aligned part of the query sequence.
"""
function fullsubject2query(
    region::AlignedRegion,
    i::Integer
)::Integer
    subjectpos = fullsubject2subject(region, i)
    return subject2query(region, subjectpos)
end

"""
    fullsubject2fullquery(region::AlignedRegion, i::Integer)::UnitRange

Given a position, `i`, in the full reference sequence, return the corresponding _positions_
in the full query sequence.
"""
function fullsubject2fullquery(
    region::AlignedRegion,
    i::Integer
)::UnitRange
    querypos = fullsubject2query(region, i)
    return query2fullquery(region, querypos)
end

"""
    query2fullquery(region::AlignedRegion, i::Integer)::Integer

Given a position, `i`, in the aligned part of the query sequence, return the corresponding
_positions_ in the full part of the query sequence.
"""
function query2fullquery(
    region::AlignedRegion,
    i::Integer
)::UnitRange
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
function fullsubject_queryboundaries(
    region::AlignedRegion,
    range::UnitRange
)::Tuple
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

Given a `range` in the aligned part of the reference sequence, return the first and last
positions in the full query sequence, which correspond to the provided range.
"""
function subject_queryboundaries(
    region::AlignedRegion,
    range::UnitRange
)::Tuple
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

Return the number of amino acids or nucleotides in this `region`'s subject sequence
(without gaps).
"""
function Base.length(region::AlignedRegion)
    return region.subjectlast - region.subjectfirst + 1
end

Base.isempty(region::AlignedRegion) = isempty(region.pairwisealignment)
#Base.isempty(region::AlignedRegion) = length(region) < 1
Base.firstindex(region::AlignedRegion) = 1
Base.lastindex(region::AlignedRegion) = length(region)
Base.eachindex(region::AlignedRegion) = Base.OneTo(lastindex(region))

function Base.getindex(
    aln::BioAlignments.PairwiseAlignment,
    indices::UnitRange
)
    queryinterval = BioAlignments.ref2seq(aln, indices)
    subjectseq = aln.b[indices]
    isequal(queryinterval, 0:0) &&
        return pairalign_global(BioSequences.LongAA(), subjectseq)
    queryinterval = max(1, first(queryinterval)):last(queryinterval)
    queryseq = aln.a.seq[queryinterval]
    return pairalign_global(queryseq, subjectseq)
end

function Base.getindex(
    region::AlignedRegion,
    indices::UnitRange
)
    subjectfirst = first(indices)
    subjectlast = last(indices)
    queryfirst, querylast = Patchwork.subject_queryboundaries(region, indices)
    return AlignedRegion(region.pairwisealignment[indices].aln, region.full_querysequence,
        subject2fullsubject(region, subjectfirst), subject2fullsubject(region, subjectlast),
        region.queryid, queryfirst, querylast, region.queryframe)#, region.fqueryfist, region.fquerylast)
end

function slicealignment(region::AlignedRegion, indices::UnitRange)::AlignedRegion
    query = region.pairwisealignment.a.seq
    subject = region.pairwisealignment.b
    anchors = deepcopy(region.pairwisealignment.a.aln.anchors)
    @assert first(anchors).op == BioAlignments.OP_START
    @assert length(anchors) >= 2
    subjectstart = first(indices)
    subjectend = last(indices)
    while length(anchors) >= 2 && (isdeleteop(anchors[2].op) ||
            anchors[2].refpos < first(indices))
        removeanchor = popat!(anchors, 2)
        if (length(anchors) >= 2 && isdeleteop(removeanchor.op) &&
                anchors[2].refpos > first(indices))
            subjectstart = anchors[2].refpos
        end
    end
    @assert length(anchors) >= 2
    lastanchor = last(anchors)
    while length(anchors) >= 2 && (isdeleteop(lastanchor.op) ||
            lastanchor.refpos > last(indices))
        toomuch = pop!(anchors)
        lastanchor = last(anchors)
        if lastanchor.refpos < last(indices)
            if isdeleteop(toomuch.op)
                subjectend = lastanchor.refpos
            else
                #packageversion("BioAlignments") >= 210
                    seqposition = subject2query(region, subjectend)
                    alnposition = lastanchor.alnpos + max(seqposition - lastanchor.seqpos,
                        subjectend - lastanchor.refpos)
                    push!(anchors, BioAlignments.AlignmentAnchor(seqposition, subjectend,
                        alnposition, toomuch.op))
                #else
                #    push!(anchors, BioAlignments.AlignmentAnchor(subject2query(region, subjectend),
                #    subjectend, toomuch.op))
                #end
                lastanchor = last(anchors)
            end
        end
    end
    @assert length(anchors) >= 2 # not strictly speaking necessary
    while isinsertop(last(anchors).op)
        pop!(anchors)
    end
    queryint = BioAlignments.ref2seq(region.pairwisealignment, subjectstart:subjectend)

    # neither of these should happen...
    isequal(queryint, 0:0) && return AlignedRegion() # empty/gap-only alignment
    #queryint = max(1, first(queryint)):last(queryint)
    @assert first(queryint) != 0

    #latest_bioalignments = packageversion("BioAlignments") >= 210

    querystart = first(queryint)
    newanchors = [anchors[1]]
    #if latest_bioalignments
	    offset = anchors[2].alnpos -
	    	max(anchors[2].seqpos-querystart+1, anchors[2].refpos-subjectstart+1)
    #end
    for anchor in anchors[2:lastindex(anchors)]
        #if latest_bioalignments
            seqposition = anchor.seqpos - querystart + 1
            refposition = anchor.refpos - subjectstart + 1
            alnposition = anchor.alnpos - offset
            push!(newanchors, BioAlignments.AlignmentAnchor(seqposition, refposition,
                alnposition, anchor.op))
        #else
        #    push!(newanchors, BioAlignments.AlignmentAnchor(anchor.seqpos - querystart + 1,
        #        anchor.refpos - subjectstart + 1, anchor.op))
        #end
    end
    sliced = BioAlignments.PairwiseAlignment(query[queryint],
        subject[subjectstart:subjectend], cigar(newanchors))
    queryfirst, querylast = subject_queryboundaries(region, subjectstart:subjectend)
    return AlignedRegion(sliced, region.full_querysequence,
        subject2fullsubject(region, subjectstart), subject2fullsubject(region, subjectend),
        region.queryid, queryfirst, querylast, region.queryframe)
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
function isordered(
    a::AlignedRegion,
    b::AlignedRegion
)::Bool
    return subject_leftposition(a) <= subject_leftposition(b)
end

"""
    precedes(a::AlignedRegion, b::AlignedRegion)

Check whether the interval in `a` starts before `b`.
"""
function precedes(
    a::AlignedRegion,
    b::AlignedRegion
)::Bool
    return subject_leftposition(a) < subject_leftposition(b)
end

"""
    isoverlapping(a::AlignedRegion, b::AlignedRegion)

Returns `true` if the interval of `a` overlaps that of `b`.
"""
function BioGenerics.isoverlapping(
    a::AlignedRegion,
    b::AlignedRegion
)::Bool
    return subject_leftposition(a) <= subject_rightposition(b) &&
           subject_leftposition(b) <= subject_rightposition(a)
end

"""
    samerange(a::AlignedRegion, b::AlignedRegion)

Returns `true` if the interval of `a` is the same as that of `b`.
"""
function samerange(
    a::AlignedRegion,
    b::AlignedRegion
)::Bool
    return subject_leftposition(a) == subject_leftposition(b) &&
           subject_rightposition(b) == subject_rightposition(a)
end

"""
    identifier(region::AlignedRegion, separator::AbstractChar='@')

Returns the sequence identifier of `region` with the OTU and sequence ID
separated by `separator`.
"""
function identifier(
    region::AlignedRegion,
    separator::AbstractChar='@'
)
    return identifier(region.record, separator)
end

"""
    sameid(a, b)

Returns `true` if the sequence identifer of `a` and `b` are identical.
"""
function sameid(
    a::AlignedRegion,
    b::AlignedRegion
)::Bool
    return identifier(a) == identifier(b)
end

"""
    samesequence(a, b)

Returns `true` if the sequence data of `a` and `b` are identical.
"""
function samesequence(
    a::AlignedRegion,
    b::AlignedRegion
)::Bool
    return a.record.sequencedata == b.record.sequencedata
end

"""
    merge(a::AlignedRegion, b::AlignedRegion, skipcheck::Bool=false)

Given two overlapping `AlignedRegion`s, slice the regions such that the best-scoring region
covers the overlapping interval. Returns an array of `AlignedRegion` objects.
"""
function merge(a::AlignedRegion, b::AlignedRegion, skipcheck::Bool=false)
    !skipcheck && !BioGenerics.isoverlapping(a, b) && error("region $a and $b are not overlapping")

    overlappingregion = overlap(a, b)
    bestscore, lowestscore = order(a, b, overlappingregion)
    if shadows(bestscore, lowestscore) || samerange(bestscore, lowestscore)
        bestfirst = fullsubject2subject(bestscore, bestscore.subjectfirst)
        bestlast = fullsubject2subject(bestscore, bestscore.subjectlast)
        return [slicealignment(bestscore, bestfirst:bestlast)]
    elseif precedes(bestscore, lowestscore) || bestscore.subjectfirst == lowestscore.subjectfirst
        bestfirst = fullsubject2subject(bestscore, bestscore.subjectfirst)
        bestlast = fullsubject2subject(bestscore, last(overlappingregion))
        lowestfirst = fullsubject2subject(lowestscore, last(overlappingregion) + 1)
        lowestlast = fullsubject2subject(lowestscore, lowestscore.subjectlast)
        return [slicealignment(bestscore, bestfirst:bestlast),
                slicealignment(lowestscore, lowestfirst:lowestlast)]
    else#if precedes(lowestscore, bestscore)
        beforeoverlap_first = fullsubject2subject(lowestscore, lowestscore.subjectfirst)
        beforeoverlap_last = fullsubject2subject(lowestscore, first(overlappingregion) - 1)
        bestfirst = fullsubject2subject(bestscore, first(overlappingregion))
        bestlast = fullsubject2subject(bestscore, last(overlappingregion))
        if lowestscore.subjectlast > last(overlappingregion)
            afteroverlap_first = fullsubject2subject(lowestscore, last(overlappingregion) + 1)
            afteroverlap_last = fullsubject2subject(lowestscore, lowestscore.subjectlast)
            return [slicealignment(lowestscore, beforeoverlap_first:beforeoverlap_last),
                slicealignment(bestscore, bestfirst:bestlast),
                slicealignment(lowestscore, afteroverlap_first:afteroverlap_last)]
        else
            return [slicealignment(lowestscore, beforeoverlap_first:beforeoverlap_last),
                slicealignment(bestscore, bestfirst:bestlast)]
        end
    end
end

"""
    leftintersection(a, b)

Returns the leftmost position at which the two regions `a` and `b` meet.
"""
function leftintersection(
    a::AlignedRegion,
    b::AlignedRegion
)
    return max(subject_leftposition(a), subject_leftposition(b))
end

"""
    rightintersection(a, b)

Returns the rightmost position at which the two regions `a` and `b` meet.
"""
function rightintersection(
    a::AlignedRegion,
    b::AlignedRegion
)
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
function overlap(
    a::AlignedRegion,
    b::AlignedRegion
)::UnitRange
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
function beforeoverlap(
    a::AlignedRegion,
    b::AlignedRegion
)
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
function afteroverlap(
    a::AlignedRegion,
    b::AlignedRegion
)
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
function shadows(
    a::AlignedRegion,
    b::AlignedRegion
)::Bool
    return length(a) > length(b) &&
           subject_leftposition(a) <= subject_leftposition(b) &&
           subject_rightposition(a) >= subject_rightposition(b) # b ends BEFORE a, see docstring?!
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
function equalrange(
    a::AlignedRegion,
    b::AlignedRegion
)::Bool
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
function totalrange(
    a::AlignedRegion,
    b::AlignedRegion
)
    if !BioGenerics.isoverlapping(a, b)
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
function longest(
    a::AlignedRegion,
    b::AlignedRegion
)
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

"""
    packageversion(name)

Returns the version number depedency with the provided `name` in the form of an integer.
For example, a package with version number 2.0.1 will be returned as 201.
"""
function packageversion(name::AbstractString)
    dependencies = [entry.second for entry in Pkg.dependencies()]
    match = only(filter(dependency -> dependency.name == name, dependencies))
    versionstring = string(match.version)
    versionnumber = parse(Int64, replace(versionstring, '.' => ""))
    return versionnumber
end
