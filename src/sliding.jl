# Provides a "sliding window" approach for removing spurious alignments

import BioAlignments
import BioSequences

function cutalignment(
    alignment::BioAlignments.PairwiseAlignment{BioSequences.LongAminoAcidSeq,BioSequences.LongAminoAcidSeq},
    cuttingpoints::Vector{Tuple{Int64,Int64}}
)
    keep = []
    keepopen = false
    cutbegin = first(first(cuttingpoints))
    cutend = last(last(cuttingpoints))
    keepbegin = 1
    anchors = alignment.a.aln.anchors
    alignmentend = last(anchors).refpos
    sequence = LongAminoAcidSeq()

    # If the first cut starts after the first position, include everything up until that
    # point
    if cutbegin > 1
        push!(keep, (keepbegin, cutbegin - 1))
    end

    # Retrieve positions in between the cutting points
    for (cutbegin, cutend) in cuttingpoints
        if !keepopen
            keepbegin = cutend + 1
            keepopen = true
        else
            keepend = cutbegin - 1
            push!(keep, (keepbegin, keepend))
            keepopen = false
        end
    end

    # If the last cut comes before the last position in the alignment, include everything
    # after that point
    if cutend < alignmentend
        push!(keep, (cutend + 1, alignmentend))
    end

    for (keepbegin, keepend) in keep
        sequence *= alignment[keepbegin:keepend].aln.a.seq
    end

    return sequence
end

"""
    slidingwindow(alignment, windowsize, distancethreshold)

Performs sliding window alignment trimming using a window of size `windowsize` on the
provided `alignment` where the average distance within the window falls below the
`distancethreshold`. Returns a vector of starts and stops where the
"""
function slidingwindow(
    alignment::BioAlignments.PairwiseAlignment{BioSequences.LongAminoAcidSeq,BioSequences.LongAminoAcidSeq},
    windowsize::Int64,
    distancethreshold::Float64,
    scoremodel = DEFAULT_SCOREMODEL::BioAlignments.AbstractScoreModel
)
    # If the alignment length is smaller or equal to that of the window size, then make
    # the window size the length of the alignment
    if length(alignment) <= windowsize
        windowsize = length(alignment)
    end

    anchors = alignment.a.aln.anchors
    alignmentend = last(anchors).refpos

    windowstart = 1
    windowstop = alignmentend - windowsize + windowstart
    flagged = Tuple{Int64,Int64}[]
    cutbegin = 0
    cutend = 0
    cutopen = false
    lastpos = 0

    for i in windowstart:windowstop
        lastpos = i + windowsize - 1
        windowtotal = BioAlignments.score(alignment[i:lastpos])
        windowaverage = windowtotal / windowsize

        if windowaverage < distancethreshold
            if cutopen
                cutend = lastpos
            else
                cutbegin = i
                cutend = lastpos
                cutopen = true
            end
        elseif cutopen && i > cutend + 1
            push!(flagged, (cutbegin, cutend))
            cutopen = false
        end
    end
    if cutopen
        push!(flagged, (cutbegin, lastpos))
    end

    if length(flagged) == 0 # no need for trimming
        return alignment
    else
        return pairalign_global(cutalignment(alignment, flagged), alignment.b, scoremodel).aln
    end
end
