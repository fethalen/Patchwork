import BioAlignments
using BioSequences

include("alignedregion.jl")
include("alignedregioncollection.jl") 

#########################################################################################################################################################

function BioAlignments.cigar(alignment::BioAlignments.PairwiseAlignment)
	anchors = alignment.a.aln.anchors
	out = IOBuffer()
	if isempty(anchors)
		return ""
	end
	@assert anchors[1].op == BioAlignments.OP_START
	for i in 2:lastindex(anchors)
		positions = max(anchors[i].seqpos - anchors[i-1].seqpos,
						anchors[i].refpos - anchors[i-1].refpos)
		print(out, positions, anchors[i].op)
	end
	return String(take!(out))
end

"""
	createbridgealignment(reference::LongSequence)

Create a `PairwiseAlignment` object with an empty (gap-only) query and the provided 
`reference`.
"""
function createbridgealignment(reference::LongSequence)
	gapquery = typeof(reference)()
	gapcigar = string(length(reference)) * "D"
	return BioAlignments.PairwiseAlignment(gapquery, reference, gapcigar)
end

# CONCATENTATION ########################################################################################################################################

# PairwiseAlignment #####################################################################################################################################

"""
	concatenate(first::PairwiseAlignment, second::PairwiseAlignment)

Construct a new `PairwiseAlignment` object by joining the first and second alignment object 
in the order in which they where provided. The sequences will not be realigned.
"""
function concatenate(first::BioAlignments.PairwiseAlignment, 
					 second::BioAlignments.PairwiseAlignment)
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

	@assert Alphabet(firstquery) == Alphabet(secondquery)

	joinedquery = typeof(firstquery)(firstquery, secondquery)
	joinedreference = typeof(firstquery)(firstreference, secondreference)
	cigar = BioAlignments.cigar(first.a.aln) * BioAlignments.cigar(second.a.aln)]
	return BioAlignments.PairwiseAlignment(joinedquery, joinedreference, cigar)
end

"""
	concatenate(alignments::AbstractVector{PairwiseAlignment})

Construct a new `PairwiseAlignment` object by joining the vector of alignment objects in 
the order in which they where provided. The sequences will not be realigned.
"""
function concatenate(alignments::AbstractVector{<:BioAlignments.PairwiseAlignment})
	if isempty(alignments)
		return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
	elseif length(alignments) == 1
		return alignments[1]
	end
	
	queries = [alignment.a.seq for alignment in alignments]
	references = [alignment.b for alignment in alignments]
	
	@assert eltype(queries) == typeof(queries[1])

	joinedquery = typeof(queries[1])(queries...)
	joinedreference = typeof(queries[1])(references...)
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
function concatenate(regions::Patchwork.AlignedRegionCollection)
	if isempty(regions)
		return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
	elseif length(regions) == 1
		return regions[1].pairwisealignment
	end
	@assert !isempty(regions.referencesequence)
	
	if regions[1].subjectfirst > 1
		reference = regions.referencesequence.sequencedata[1:regions[1].subjectfirst - 1]
		alignments = [createbridgealignment(reference), regions[1].pairwisealignment]
	else
		alignments = [regions[1].pairwisealignment]
	end
	# fill the gaps between end and next start
	for i in 2:lastindex(regions)
		firstotu = Patchwork.splitdescription(regions[i-1].queryid)[1]
		secondotu = Patchwork.splitdescription(regions[i].queryid)[1]
		@assert firstotu == secondotu
		@assert regions[i-1].subjectlast < regions[i].subjectfirst		 	# check sorted
		if regions[i].subjectfirst > regions[i-1].subjectlast + 1
			reference = regions.referencesequence.sequencedata[
								regions[i-1].subjectlast + 1:regions[i].subjectfirst - 1]
			push!(alignments, createbridgealignment(reference), regions[i].pairwisealignment)
		else
			push!(alignments, regions[i].pairwisealignment)
		end
		if i == lastindex(regions) 
		   && regions[i].subjectlast < lastindex(regions.referencesequence)
			reference = regions.referencesequence.sequencedata[
						regions[i].subjectlast + 1:lastindex(regions.referencesequence)]
			push!(alignments, createbridgealignment(reference))
		end
	end
	return concatenate(alignments)
end

# OCCUPANCY #############################################################################################################################################

"""
	countmatches(alignment::PairwiseAlignment)

Compute the absolute number of residues in the query sequence that align to residues in the 
reference sequence.
"""
function countmatches(alignment::BioAlignments.PairwiseAlignment)
	if isempty(alignment)
		return 0
	end

	covered = 0
	anchors = alignment.a.aln.anchors

	@assert anchors[1].op == BioAlignments.OP_START

	for i in 2:lastindex(anchors)
		if anchors[i].op == BioAlignments.OP_MATCH
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
function countmatches(region::AlignedRegion)
	return countmatches(region.pairwisealignment)
end

"""
	countgaps(alignment::PairwiseAlignment)

Compute the absolute number of gaps in the query sequence that align to residues in the 
reference sequence.
"""
function countgaps(alignment::BioAlignments.PairwiseAlignment)
	if isempty(alignment)
		return 0
	end

	gaps = 0
	anchors = alignment.a.aln.anchors

	@assert anchors[1].op == BioAlignments.OP_START

	for i in 2:lastindex(anchors)
		if anchors[i].op == BioAlignments.OP_DELETE
			gaps += (anchors[i].refpos - anchors[i-1].refpos)
		end
	end
	return gaps
end

"""
	countgaps(region::AlignedRegion)

Compute the absolute number of gaps in the query sequence that align to residues in the 
reference sequence.
"""
function countgaps(region::AlignedRegion)
	return countgaps(region.pairwisealignment)
end

"""
	occupancy(alignment::PairwiseAlignment)

Compute the relative amount of residues in the query sequence that align to residues in the 
reference sequence.
"""
function occupancy(alignment::BioAlignments.PairwiseAlignment)
	return countmatches(alignment) / length(alignment.b)
end

"""
	occupancy(region::AlignedRegion)

Compute the relative amount of residues in the query sequence that align to residues in the 
reference sequence.
"""
function occupancy(region::AlignedRegion)
	return occupancy(region.pairwisealignment)
end