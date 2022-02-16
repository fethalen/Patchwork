import BioAlignments
import BioSequences
import DataFrames
import Random

include("sequencerecord.jl")

# Default scoremodel taken from DIAMOND's defaults for BLOSUM62
const DEFAULT_SCOREMODEL = BioAlignments.AffineGapScoreModel(BioAlignments.BLOSUM62,
    gap_open=-11, gap_extend=-1)

"""
    pairalign_global(seq, ref, scoremodel)

Globally align the `seq` (sequence) against a `ref` (reference) sequence using the
provided `scoremodel`.
"""
function pairalign_global(
    seq::SequenceRecord,
    ref::SequenceRecord,
    scoremodel=DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel
)
    return BioAlignments.pairalign(BioAlignments.GlobalAlignment(), seq.sequencedata,
                                   ref.sequencedata, scoremodel)
end

function pairalign_global(
    seq::BioSequences.LongSequence,
    ref::BioSequences.LongSequence,
    scoremodel=DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel
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
    scoremodel=DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel
)
    return BioAlignments.pairalign(BioAlignments.LocalAlignment(), seq.sequencedata,
                                   ref.sequencedata, scoremodel)
end


function pairalign_local(
    seq::BioSequences.LongSequence,
    ref::BioSequences.LongSequence,
    scoremodel=DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel
)
    return BioAlignments.pairalign(BioAlignments.LocalAlignment(), seq, ref, scoremodel)
end

function scorealignment(
    seq::BioSequences.LongSequence,
    ref::BioSequences.LongSequence
)
    reflen = length(ref)
    seqlen = length(seq)
    # alignment = pairalign_global(seq, ref)
    alignment = pairalign_local(seq, ref)

    if seqlen >= reflen
        occupancy = 100.
        percentidentity = round((BioAlignments.count_matches(alignment.aln) / length(alignment.aln)) * 100., digits=2)
    elseif reflen > 0
        occupancy = (seqlen / reflen) * 100.
        percentidentity = round(
            BioAlignments.count_matches(alignment.aln) /
            (BioAlignments.count_matches(alignment.aln) +
            BioAlignments.count_mismatches(alignment.aln)) * 100.,
            digits=2
        )
    else
        occupancy = 0.
        percentidentity = 0.
    end

    querycoverage = round(occupancy, digits=2)
    row = [seqlen, reflen, BioAlignments.count_matches(alignment.aln),
           BioAlignments.count_mismatches(alignment.aln),
           BioAlignments.count_insertions(alignment.aln),
           BioAlignments.count_deletions(alignment.aln),
           BioAlignments.distance(alignment), querycoverage, percentidentity]
    return (alignment, row)
end
