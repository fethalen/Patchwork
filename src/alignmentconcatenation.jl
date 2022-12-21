import BioAlignments
using BioSequences

#include("alignment.jl")
#include("alignedregion.jl")
#include("alignedregioncollection.jl")

#########################################################################################################################################################

"""
    createbridgealignment(reference::LongSequence)

Create a `PairwiseAlignment` object with an empty (gap-only) query and the provided
`reference`.
"""
function createbridgealignment(
    reference::LongSequence,
    range::UnitRange
)::BioAlignments.PairwiseAlignment
    first(range) < 1 || last(range) > lastindex(reference) &&
        error("The provided indices $range are out ouf bounds for this reference sequence.")
    gapquery = typeof(reference)()
    gapcigar = string(length(range)) * "D"
    return BioAlignments.PairwiseAlignment(gapquery, reference[range], gapcigar)
end

#function createbridge_dna(length::Int64)::BioSequences.LongDNA
#    return BioSequences.LongDNA{4}(repeat("N", length))
#end

"""
    concatenate(first::PairwiseAlignment, second::PairwiseAlignment)

Construct a new `PairwiseAlignment` object by joining the first and second alignment object
in the order in which they where provided. The sequences will not be realigned.
"""
function concatenate(
    first::BioAlignments.PairwiseAlignment,
    second::BioAlignments.PairwiseAlignment
)::BioAlignments.PairwiseAlignment
    if isempty(first) && isempty(second)
        return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
    elseif isempty(first)
        return second
    elseif isempty(second)
        return first
    end

    firstquery = first.a.seq
    secondquery = second.a.seq
    firstreference = first.b
    secondreference = second.b

    @assert Alphabet(firstquery) == Alphabet(secondquery) """Can only concatenate
        alignments of same type (i.e. two protein-protein alignments)."""

    # joinedquery = typeof(firstquery)(firstquery, secondquery)
    # joinedreference = typeof(firstquery)(firstreference, secondreference)
    joinedquery = firstquery * typeof(firstquery)(secondquery)
    joinedreference = typeof(firstquery)(firstreference) * typeof(firstquery)(secondreference)
    cigar = BioAlignments.cigar(first.a.aln) * BioAlignments.cigar(second.a.aln)
    return BioAlignments.PairwiseAlignment(joinedquery, joinedreference, cigar)
end

"""
    concatenate(alignments::AbstractVector{PairwiseAlignment})

Construct a new `PairwiseAlignment` object by joining the vector of alignment objects in
the order in which they where provided. The sequences will not be realigned.
"""
function concatenate(
    alignments::AbstractVector{<:BioAlignments.PairwiseAlignment}
)::BioAlignments.PairwiseAlignment
    if isempty(alignments)
        return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
    elseif length(alignments) == 1
        return alignments[1]
    end

    queries = [alignment.a.seq for alignment in alignments]
    references = [alignment.b for alignment in alignments]

    @assert eltype(queries) == typeof(queries[1]) """Can only concatenate alignments of
        same type (i.e. a collection of protein-protein alignments)."""

    # joinedquery = typeof(queries[1])(queries...)
    # joinedreference = typeof(queries[1])(references...)
    joinedquery = typeof(queries[1])(join(queries))
    joinedreference = typeof(queries[1])(join(references))
    out = IOBuffer()
    for alignment in alignments
        print(out, BioAlignments.cigar(alignment))
    end

    return BioAlignments.PairwiseAlignment(joinedquery, joinedreference, String(take!(out)))
end

# AlignedRegionCollection ###############################################################################################################################

"""
    concatenate(regions::AlignedRegionCollection)

Construct a new `PairwiseAlignment` object by joining the collection of alignment objects
in the order in which they where provided. The sequences will not be realigned.
This function expects `regions` to be sorted and free of aligned region overlaps.
The `referencesequence` of the collection may not be empty. Regions in the reference that
are not covered by alignments to query sequences will be aligned to gaps (empty queries).
"""
function concatenate(
    regions::Patchwork.AlignedRegionCollection,
    delimiter::Char = '@'
)::Tuple{BioAlignments.PairwiseAlignment, BioSequences.LongDNA}
    @assert !hasoverlaps(regions) "Detected overlaps in regions."
    @assert !isempty(regions.referencesequence) "Reference sequence must not be empty."
    reference = regions.referencesequence.sequencedata
    isempty(regions) &&
        return createbridgealignment(reference, 1:lastindex(regions.referencesequence))#, DNA???

    if regions[1].subjectfirst > 1
        alignments = [createbridgealignment(reference, 1:regions[1].subjectfirst-1),
            regions[1].pairwisealignment]
#        dnalength = 3 * (regions[1].subjectfirst-1)
#        dnaqueries = [createbridge_dna(dnalength), regions[1].full_querysequence[regions[1].queryfirst:regions[1].querylast]]
    else
        alignments = [regions[1].pairwisealignment]
#        dnaqueries = [regions[1].full_querysequence[regions[1].queryfirst:regions[1].querylast]]
    end
    dnaqueries = [regions[1].full_querysequence[regions[1].queryfirst:regions[1].querylast]] # no need for bridge DNA!

    if lastindex(regions) == 1
        if regions[1].subjectlast < lastindex(regions.referencesequence)
            push!(alignments, createbridgealignment(reference,
                regions[1].subjectlast+1:lastindex(regions.referencesequence)))
#            dnalength = 3 * (lastindex(regions.referencesequence)-regions[1].subjectlast)
#            push!(dnaqueries, createbridge_dna(dnalength))
        end
        return concatenate(alignments), LongDNA{4}(join(dnaqueries))
    end

    otu = otupart(regions[1].queryid, delimiter) # empty if no OTU part found
    # fill the gaps between end and next start)
    for i in 2:lastindex(regions)
        if !isempty(otu) # missing OTUs are always okay
            currentotu = otupart(regions[i].queryid, delimiter)
            !isempty(currentotu) && @assert isequal(currentotu, otu) """Can only concatenate
                contigs from same species."""
        else
            otu = otupart(regions[i].queryid, delimiter)
        end
        @assert regions[i-1].subjectlast < regions[i].subjectfirst """Regions incorrectly
            sorted."""
        if regions[i].subjectfirst > regions[i-1].subjectlast + 1
            push!(alignments, createbridgealignment(reference,
                    regions[i-1].subjectlast+1:regions[i].subjectfirst-1),
                regions[i].pairwisealignment)
#            dnalength = 3 * (regions[i].subjectfirst - regions[i-1].subjectlast + 1)
#            push!(dnaqueries, createbridge_dna(dnalength), regions[i].full_querysequence[regions[i].queryfirst:regions[i].querylast])
        else
            push!(alignments, regions[i].pairwisealignment)
#            push!(dnaqueries, regions[i].full_querysequence[regions[i].queryfirst:regions[i].querylast])
        end
        push!(dnaqueries, regions[i].full_querysequence[regions[i].queryfirst:regions[i].querylast])

        if (i == lastindex(regions)
            && regions[i].subjectlast < lastindex(regions.referencesequence))
            push!(alignments, createbridgealignment(reference,
                regions[i].subjectlast+1:length(regions.referencesequence)))
#            dnalength = 3 * (length(regions.referencesequence) - regions[i].subjectlast + 1)
#            push!(dnaqueries, createbridge_dna(dnalength))
        end
    end
    return concatenate(alignments), LongDNA{4}(join(dnaqueries))
end

"""
    countmatches(alignment::PairwiseAlignment)

Compute the absolute number of residues in the query sequence that align to residues in the
reference sequence.
"""
function countmatches(
    alignment::BioAlignments.PairwiseAlignment
)::Int64
    isempty(alignment) && return 0

    covered = 0
    anchors = alignment.a.aln.anchors
    @assert anchors[1].op == BioAlignments.OP_START "Alignments must start with OP_START."
    for i in 2:lastindex(anchors)
        if BioAlignments.ismatchop(anchors[i].op)
            covered += (anchors[i].seqpos - anchors[i-1].seqpos)
        end
    end
    return covered
end

"""
    countmatches(alignment::AlignedRegion)

Compute the absolute number of residues in the query sequence that align to residues in the
reference sequence.
"""
function countmatches(region::AlignedRegion)::Int64
    return countmatches(region.pairwisealignment)
end

#"""
#    countgaps(alignment::PairwiseAlignment)
#    countgaps(region::AlignedRegion)
#
#Compute the absolute number of gaps in the query sequence that align to residues in the
#reference sequence.
#"""
function countgaps(alignment::BioAlignments.PairwiseAlignment)::Int64
    isempty(alignment) && return 0

    gaps = 0
    anchors = alignment.a.aln.anchors
    @assert anchors[1].op == BioAlignments.OP_START "Alignments must start with OP_START."
    for i in 2:lastindex(anchors)
        if BioAlignments.isdeleteop(anchors[i].op)
            gaps += (anchors[i].refpos - anchors[i-1].refpos)
        end
    end
    return gaps
end

function countgaps(region::AlignedRegion)::Int64
    return countgaps(region.pairwisealignment)
end

"""
    occupancy(alignment::PairwiseAlignment)

Compute the relative amount of residues in the query sequence that align to residues in the
reference sequence.
"""
function occupancy(alignment::BioAlignments.PairwiseAlignment)::Float64
    return countmatches(alignment) / length(alignment.b)
end

"""
    occupancy(region::AlignedRegion)

Compute the relative amount of residues in the query sequence that align to residues in the
reference sequence.
"""
function occupancy(region::AlignedRegion)::Float64
    return occupancy(region.pairwisealignment)
end

"""
    occupancy(regions::AlignedRegionCollection)

Compute the relative amount of residues in the query sequence that align to residues in the
reference sequence. This function concatenates `regions` before computing occupancy, returning
the occupancy score based on the entire reference sequence of the collection.
"""
function occupancy(regions::AlignedRegionCollection)::Float64
    return occupancy(concatenate(regions)[1]) # use the pw aln, not the concatenated dna seq
end


"""
    maskgaps(alignment::BioAlignments.PairwiseAlignment)
Remove unaligned columns (i.e., insertions in the reference sequence) _in the query
sequence_.
"""
function maskgaps(
    alignment::BioAlignments.PairwiseAlignment, 
    dna::BioSequences.LongDNA
)::Tuple{BioAlignments.PairwiseAlignmentResult, BioSequences.LongDNA}
    anchors = alignment.a.aln.anchors
    maskedseq = BioSequences.LongAA()
    maskeddna = BioSequences.LongDNA{4}()
    from = 1
    dnafrom = 1
    gapcount = 0
    for i in 2:lastindex(anchors)
        if isinsertop(anchors[i].op) # === BioAlignments.OP_INSERT
            firstinsertion = anchors[i-1].seqpos + 1
            lastinsertion = anchors[i].seqpos
            to = anchors[i-1].seqpos
            gapcount += lastinsertion - firstinsertion + 1
            maskedseq *= alignment.a.seq[from:to]
            from = anchors[i].seqpos + 1
            # convert AA coords to DNA coords for masking DNA seq
            dnato = 3 * anchors[i-1].seqpos 
            maskeddna *= dna[dnafrom:dnato]
            dnafrom = 3 * anchors[i].seqpos + 1
        elseif i == lastindex(anchors)
            to = anchors[i].seqpos
            maskedseq *= alignment.a.seq[from:to]
            # convert AA coords to DNA coords for masking DNA seq
            dnato = 3 * anchors[i].seqpos 
            maskeddna *= dna[dnafrom:dnato]
        end
    end
    return pairalign_global(maskedseq, alignment.b), maskeddna
end

"""
    maskalignment(alignment, scoremodel, retainstops, retainambiguous)

Remove unaligned residues from the provided `alignment`. Also remove stop codons (`*`) and
or ambigious characters if `retainstops` and or `retainambigiuous` are not set. Realign
the resulting query sequences to the reference using the provided `scoremodel`.
"""
function maskalignment(
    alignment::BioAlignments.PairwiseAlignment,
    dnaseq::BioSequences.LongDNA,
    scoremodel::BioAlignments.AbstractScoreModel,
    retainstops::Bool = false,
    retainambiguous::Bool = false
)::Tuple{BioAlignments.PairwiseAlignmentResult, BioSequences.LongDNA}
    flagged = Int64[]
    dnaflagged = Int64[]
    queryseq = alignment.a.seq
    anchors = alignment.a.aln.anchors
    alignmentend = last(anchors).seqpos

    for (position, letters) in enumerate(alignment)
        position > alignmentend && break
        (queryletter, refletter) = letters
        if ((!retainstops && queryletter == AA_Term) ||
            (!retainambiguous && isambiguous(queryletter)) ||
            (refletter == AA_Gap)
        )
            push!(flagged, position)
            push!(dnaflagged, 3*(position-1)+1:3*position...)
        end
    end
    # We reverse the list of indices as deleting from left shifts the next deletion
    map(seqpos -> deleteat!(queryseq, seqpos), reverse(flagged))
    map(seqpos -> deleteat!(dnaseq, seqpos), reverse(dnaflagged))
    return pairalign_global(queryseq, alignment.b, scoremodel), dnaseq
end
