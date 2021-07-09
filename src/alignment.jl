import BioAlignments
import BioSequences
import Random

include("alignedregion.jl")
include("sequencerecord.jl")

# Default scoremodel taken from DIAMOND's defaults for BLOSUM62
const DEFAULT_SCOREMODEL = BioAlignments.AffineGapScoreModel(BioAlignments.BLOSUM62,
    gap_open=-11, gap_extend=-1)

"""
    pairalign_global(seq, ref, scoremodel)

Globally align the `seq` (sequence) against a `ref` (reference) sequence using the
provided `scoremodel`.
"""
function pairalign_global(seq::SequenceRecord, ref::SequenceRecord,
                          scoremodel=DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel)
    return BioAlignments.pairalign(BioAlignments.GlobalAlignment(), seq.sequencedata,
                                   ref.sequencedata, scoremodel)
end

function pairalign_global(seq::BioSequences.LongSequence, ref::BioSequences.LongSequence,
                         scoremodel=DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel)
    return BioAlignments.pairalign(BioAlignments.GlobalAlignment(), seq, ref, scoremodel)
end

"""
    pairalign_local(seq, ref, scoremodel)

Align the `seq` (sequence) against a `ref` (reference) sequence using the
provided `scoremodel`.
"""
function pairalign_local(seq::SequenceRecord, ref::SequenceRecord,
                         scoremodel=DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel)
    return BioAlignments.pairalign(BioAlignments.LocalAlignment(), seq.sequencedata,
                                   ref.sequencedata, scoremodel)
end


function pairalign_local(seq::BioSequences.LongSequence, ref::BioSequences.LongSequence,
                         scoremodel=DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel)
    return BioAlignments.pairalign(BioAlignments.LocalAlignment(), seq, ref, scoremodel)
end

"""
    slice(region::AlignedRegion, interval::UnitRange)

Return the subpart of the provided `region` as given by the `interval`.
"""
function slice(region::AlignedRegion, interval::UnitRange)
    alignment = region.pairwisealignment
    queryinterval = BioAlignments.ref2seq(alignment, ref2reg(region, interval))
    subjectinterval = BioAlignments.seq2ref(alignment, queryinterval)
    queryportion = region.pairwisealignment.a.seq[queryinterval]
    subjectportion = region.pairwisealignment.b[subjectinterval]
    println(region)
    println("new subject interval: ", interval)
    # println("new query interval: ", BioAlignments.ref2seq(alignment, interval))
    println("sequence queryinterval: ", queryinterval)
    println("sequence subjectinterval: ", subjectinterval)
    slicedaln = pairalign_global(queryportion, subjectportion, DEFAULT_SCOREMODEL)
end

"""
    bestscore(a, b, interval)

Realign the provided `AlignedRegion` objects at the provided interval and return the
name of the best-scoring sequence. If the alignment scores are equal, select a region by
random.
"""
function bestscore(a::AlignedRegion, b::AlignedRegion, interval::UnitRange)
    scorea = realign(a, ref2reg(a, interval)).value
    scoreb = realign(b, ref2reg(b, interval)).value

    if scorea > scoreb
        return a
    elseif scorea == scoreb && Random.rand(0:1) == 1
        return a
    else
        return b
    end
end

"""
    realign(region, interval)

Realigned the `region` at the specific `interval`.
"""
function realign(region::AlignedRegion, interval::UnitRange)
    alignment = region.pairwisealignment
    queryinterval = BioAlignments.ref2seq(alignment, interval)
    queryseq = alignment.a.seq[queryinterval]
    subjectseq = alignment.b[interval]
    pairalign_local(queryseq, subjectseq, DEFAULT_SCOREMODEL)
end

"""
    realign(region, interval)

Realigned the `region` at the specific `interval`.
"""
function realign(region::AlignedRegion, interval::UnitRange)
    alignment = region.pairwisealignment
    queryinterval = BioAlignments.ref2seq(alignment, interval)
    queryseq = alignment.a.seq[queryinterval]
    subjectseq = alignment.b[interval]
    pairalign_local(queryseq, subjectseq, DEFAULT_SCOREMODEL)
end

"""
    PairwiseAlignment(seq, ref, cigar)

Construct a `PairwiseAlignment` from two sequences, `seq` and `ref`, and a `cigar` string.
"""
function BioAlignments.PairwiseAlignment(
    seq::BioSequences.LongSequence, ref::BioSequences.LongSequence, cigar::String)
    alignment = BioAlignments.Alignment(cigar)
    alignedseq = BioAlignments.AlignedSequence(seq, alignment)
    return BioAlignments.PairwiseAlignment(alignedseq, ref)
end

function BioAlignments.PairwiseAlignment(aln::BioAlignments.PairwiseAlignment,
                                         interval::UnitRange)
    queryinterval = BioAlignments.ref2seq(aln, interval)
    queryseq = aln.a.seq[queryinterval]
    subjectseq = aln.b[interval]
    println(length(interval))
    println(queryseq)
    println(subjectseq)
end

function BioAlignments.seq2ref(aln::BioAlignments.PairwiseAlignment, interval::UnitRange)
    # for "negative" frames, the first < stop, so we grab the leftmost and rightmost
    # positions
    firstposition = first(interval)
    lastposition = last(interval)
    leftmost = min(firstposition, lastposition)
    rightmost = max(firstposition, lastposition)
    return UnitRange(first(BioAlignments.seq2ref(aln, leftmost)),
                     first(BioAlignments.seq2ref(aln, rightmost)))
end

function BioAlignments.ref2seq(aln::BioAlignments.PairwiseAlignment, interval::UnitRange)
    # for "negative" frames, the first < stop, so we grab the leftmost and rightmost
    # positions
    firstposition = first(interval)
    lastposition = last(interval)
    leftmost = min(firstposition, lastposition)
    rightmost = max(firstposition, lastposition)
    return UnitRange(first(BioAlignments.ref2seq(aln, leftmost)),
                     first(BioAlignments.ref2seq(aln, rightmost)))
end
