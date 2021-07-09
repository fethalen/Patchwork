import BioAlignments
import BioSequences

include("sequencerecord.jl")

# Default scoremodel taken from DIAMOND's defaults for BLOSUM62
DEFAULT_SCOREMODEL = BioAlignments.AffineGapScoreModel(BioAlignments.BLOSUM62, gap_open=-11,
                                         gap_extend=-1)

function pairalign_local(seq::SequenceRecord, ref::SequenceRecord,
                         scoremodel=DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel)
    return BioAlignments.pairalign(BioAlignments.LocalAlignment(), seq.sequencedata,
                                   seq.sequencedata, SCOREMODEL)
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
