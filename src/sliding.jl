# Provides a "sliding window" approach for removing spurious alignments

import BioAlignments
import BioSequences

function cutalignment(
    alignment::BioAlignments.PairwiseAlignment{BioSequences.LongAA,BioSequences.LongAA},
    cuttingpoints::Vector{Tuple{Int64,Int64}}
)
    keep = []
    keepopen = false
    cutbegin = first(first(cuttingpoints))
    cutend = last(last(cuttingpoints))
    keepbegin = 1
    anchors = alignment.a.aln.anchors
    alignmentend = last(anchors).refpos
    sequence = LongAA()

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

function cutsequence(
    seq::Union{BioSequences.LongAA, BioSequences.LongDNA},
    cuttingpoints::Vector{Tuple{Int64,Int64}},
    seqend::Int64
)
    keep = []
    keepopen = false
    cutbegin = first(first(cuttingpoints))
    cutend = last(last(cuttingpoints))
    keepbegin = 1
    sequence = typeof(seq)()

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
    if cutend < seqend
        push!(keep, (cutend + 1, seqend))
    end

    for (keepbegin, keepend) in keep
        sequence *= seq[keepbegin:keepend]
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
    alignment::BioAlignments.PairwiseAlignment{BioSequences.LongAA,BioSequences.LongAA},
    dna::BioSequences.LongDNA,
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
    dnaflagged = Tuple{Int64,Int64}[]
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
            cutbegin = ref2seq(alignment, cutbegin)[1] # get position from (position, operation)
            cutend = ref2seq(alignment, cutend)[1]
            push!(flagged, (cutbegin, cutend))
            push!(dnaflagged, (3*(cutbegin-1)+1, 3*cutend))
            cutopen = false
        end
    end
    if cutopen
        cutbegin = ref2seq(alignment, cutbegin)[1]
        cutend = ref2seq(alignment, lastpos)[1]
        push!(flagged, (cutbegin, cutend)) # lastpos
        push!(dnaflagged, (3*(cutbegin-1)+1, lastindex(dna)))
    end

    if length(flagged) == 0 # no need for trimming
        return alignment, dna
    else
        anchors = alignment.a.aln.anchors
        alignmentend = last(anchors).seqpos # last(anchors).refpos
        return pairalign_global(cutsequence(alignment.a.seq, flagged, alignmentend), 
            alignment.b, scoremodel).aln, cutsequence(dna, dnaflagged, lastindex(dna))
    end
end
