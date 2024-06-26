"""
    pairalign_global(seq, ref, scoremodel)

Globally align the `seq` (sequence) against a `ref` (reference) sequence using the
provided `scoremodel`.
"""
function pairalign_global(
    seq::SequenceRecord,
    ref::SequenceRecord,
    scoremodel = DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel
)
    return BioAlignments.pairalign(BioAlignments.GlobalAlignment(), seq.sequencedata,
        ref.sequencedata, scoremodel)
end

function pairalign_global(
    seq::BioSequences.LongSequence,
    ref::BioSequences.LongSequence,
    scoremodel = DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel
)
    return BioAlignments.pairalign(BioAlignments.GlobalAlignment(), seq, ref, scoremodel)
end

"""
    pairalign_local(seq, ref, scoremodel)

Align the `seq` (sequence) against a `ref` (reference) sequence using the
provided `scoremodel`.
"""
function pairalign_local(
    seq::SequenceRecord,
    ref::SequenceRecord,
    scoremodel = DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel
)
    return BioAlignments.pairalign(BioAlignments.LocalAlignment(), seq.sequencedata,
        ref.sequencedata, scoremodel)
end


function pairalign_local(
    seq::BioSequences.LongSequence,
    ref::BioSequences.LongSequence,
    scoremodel = DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel
)
    return BioAlignments.pairalign(BioAlignments.LocalAlignment(), seq, ref, scoremodel)
end

"""
    bestscore(a, b, interval)

Realign the provided `AlignedRegion` objects at the provided interval and a tuple with the
best-scoring alignment first and the other alignment second.
"""
function order(
    a::AlignedRegion,
    b::AlignedRegion,
    interval::UnitRange
)::Tuple
    scorea = realign(a, fullsubject2subject(a, interval)).value
    scoreb = realign(b, fullsubject2subject(b, interval)).value
    if scorea > scoreb
        return (a, b)
    elseif scorea == scoreb && Random.rand(0:1) == 1
        return (a, b)
    else
        return (b, a)
    end
end

#"""
#    realign(region, interval)
#
#Realigned the `region` at the specific `interval`.
#"""
#function realign(region::AlignedRegion, interval::UnitRange)
#    alignment = region.pairwisealignment
#    queryinterval = BioAlignments.ref2seq(alignment, interval)
#    queryseq = alignment.a.seq[queryinterval]
#    subjectseq = alignment.b[interval]
#    pairalign_local(queryseq, subjectseq, DEFAULT_SCOREMODEL)
#end

"""
    realign(region, interval)

Realigned the `region` at the specific `interval`.
"""
function realign(
    region::AlignedRegion,
    interval::UnitRange
)
    alignment = region.pairwisealignment
    queryinterval = BioAlignments.ref2seq(alignment, interval)
    subjectseq = alignment.b[interval]
    isequal(queryinterval, 0:0) &&
        return pairalign_local(BioSequences.LongAA(), subjectseq[interval])
    queryinterval = max(1, first(queryinterval)):min(lastindex(alignment.a.seq), last(queryinterval))
    queryseq = alignment.a.seq[queryinterval]
    pairalign_local(queryseq, subjectseq, DEFAULT_SCOREMODEL)
end

"""
    PairwiseAlignment(seq, ref, cigar)

Construct a `PairwiseAlignment` from two sequences, `seq` and `ref`, and a `cigar` string.
"""
function BioAlignments.PairwiseAlignment(
    seq::BioSequences.LongSequence,
    ref::BioSequences.LongSequence,
    cigar::String
)
    alignment = BioAlignments.Alignment(cigar)
    alignedseq = BioAlignments.AlignedSequence(seq, alignment)
    return BioAlignments.PairwiseAlignment(alignedseq, ref)
end

function BioAlignments.seq2ref(
    aln::BioAlignments.PairwiseAlignment,
    interval::UnitRange
)
    # for "negative" frames, the first < stop, so we grab the leftmost and rightmost
    # positions
    firstposition = first(interval)
    lastposition = last(interval)
    leftmost = min(firstposition, lastposition)
    rightmost = max(firstposition, lastposition)
    return UnitRange(first(BioAlignments.seq2ref(aln, leftmost)),
        first(BioAlignments.seq2ref(aln, rightmost)))
end

function BioAlignments.ref2seq(
    aln::BioAlignments.PairwiseAlignment,
    interval::UnitRange
)
    # for "negative" frames, the first < stop, so we grab the leftmost and rightmost
    # positions
    firstposition = first(interval)
    lastposition = last(interval)
    leftmost = min(firstposition, lastposition)
    rightmost = max(firstposition, lastposition)
    return UnitRange(first(BioAlignments.ref2seq(aln, leftmost)),
        first(BioAlignments.ref2seq(aln, rightmost)))
end

Base.firstindex(aln::BioAlignments.PairwiseAlignment) = 1
Base.lastindex(aln::BioAlignments.PairwiseAlignment) = last(aln.a.aln.anchors).refpos
Base.eachindex(aln::BioAlignments.PairwiseAlignment) = Base.OneTo(lastindex(aln))

function Base.getindex(
    aln::BioAlignments.PairwiseAlignment,
    index::Integer
)
    return aln[index:index]
end

# TODO: Lose global alignment dependency, retrieve alignment without realigning
#function Base.getindex(aln::BioAlignments.PairwiseAlignment, indices::UnitRange)
#    queryinterval = BioAlignments.ref2seq(aln, indices)
#    queryseq = aln.a.seq[queryinterval]
#    subjectseq = aln.b[indices]
#    return pairalign_global(queryseq, subjectseq)
#end

# @inline function Base.iterate(aln::BioAlignments.PairwiseAlignment, i::Int = firstindex(aln))
#     if i > lastindex(aln)
#         return nothing
#     else
#         return getindex(aln, i), i + 1
#     end
# end

"""
    gapexcluded_identity(matches, mismatches, precision)

Calculate the gap-excluded identity by calculating #matches / (#matches + #mismatches).
The result is reported as a percentage, so 45% identity is reported as 45.0.
"""
function gapexcluded_identity(
    matches::Int,
    mismatches::Int,
    precision::Int = 2
)
    if matches + mismatches == 0
        return 0.0
    end
    percentidentity = round((matches / (matches + mismatches)) * 100.0, digits = precision)
    return percentidentity
end
