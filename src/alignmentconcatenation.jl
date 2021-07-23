using BioSequences: alignback!
import BioAlignments
import BioSequences
using FASTX

pwd()
cd("src")

include("alignedregion.jl")
include("alignedregioncollection.jl") 

subject = "../test/07673_Alitta_succinea.fa"
diamondresults = "../test/c_australis_x_07673.tsv"
hits = Patchwork.readblastTSV(diamondresults)
full_subjectseq = Patchwork.get_fullseq(subject)
regions = Patchwork.AlignedRegionCollection(full_subjectseq, hits)
first_testregion = regions[12]
second_testregion = regions[13]

show(first_testregion)
show(first_testregion.pairwisealignment.a)
show(first_testregion.pairwisealignment.a.aln)
show(first_testregion.pairwisealignment.a.seq)
show(first_testregion.pairwisealignment.b)
alignments = [region.pairwisealignment for region in regions]

# Get the PairwiseAlignment object's cigar string.
function BioAlignments.cigar(alignment::BioAlignments.PairwiseAlignment)
	anchors = alignment.a.aln.anchors
	cigar = repeat([""], 2 * length(anchors) - 2)
	if isempty(anchors)
		return ""
	end
	@assert anchors[1].op == BioAlignments.OP_START
	for i in 1:length(anchors)-1
		positions = max(anchors[i+1].seqpos - anchors[i].seqpos,
				anchors[i+1].refpos - anchors[i].refpos)
		#cigar = join([cigar, string(n), string(anchors[i].op)], "")
		cigar[2*i-1] = string(positions)
		cigar[2*i] = string(anchors[i+1].op)
	end
	return join(cigar, "")
end

# Concatenating 2 PairwiseAlignment objects: 
# Join the 2 queries and the 2 references together; 
# No need to add anything in between, right?
function concatenate(first::BioAlignments.PairwiseAlignment, 
					 second::BioAlignments.PairwiseAlignment)
	firstquery = first.a.seq
	secondquery = second.a.seq
	firstreference = first.b
	secondreference = second.b 

	@assert Alphabet(firstquery) == Alphabet(secondquery)

	joinedquery = typeof(firstquery)(firstquery, secondquery)
	joinedreference = typeof(firstquery)(firstreference, secondreference)
	cigar = join([BioAlignments.cigar(first.a.aln), BioAlignments.cigar(second.a.aln)], "")
	return BioAlignments.PairwiseAlignment(joinedquery, joinedreference, cigar)
end

function concatenate(alignments::AbstractVector{<:BioAlignments.PairwiseAlignment})
	if length(alignments) == 0
		return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
		#return nothing
	elseif length(alignments) == 1
		return alignments[1]
	end
	
	queries = [alignment.a.seq for alignment in alignments]
	references = [alignment.b for alignment in alignments]
	
	@assert eltype(queries) == typeof(queries[1])

	joinedquery = typeof(queries[1])(queries...)
	joinedreference = typeof(queries[1])(references...)
	cigar = join([BioAlignments.cigar(alignment) for alignment in alignments], "")
	
	return BioAlignments.PairwiseAlignment(joinedquery, joinedreference, cigar)
end

# make this 2 alignedregions
function concatenate(regions::Vector{AlignedRegion})
	# regions is sorted, no overlaps remaining.
	if isempty(regions)
		empty = BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
		return AlignedRegion(empty, 0, -1, SequenceIdentifier(), 0, -1, 0) # empty AlignedRegion
		#return nothing
	elseif length(regions) == 1
		return regions[1]
	end

	# make sure that the entire subject sequence is covered? 
	# or insert gaps in reference and query for missing reference positions?
	# this should happen in the concatenate method for AlignedRegionCollections...
	superalignment = concatenate([region.pairwisealignment for region in regions])
	
end

function concatenate(regions::AlignedRegionCollection)
	# regions is sorted, no overlaps remaining.
	if isempty(regions)
		return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
		#return nothing
	elseif length(regions) == 1
		return regions[1]
	end

	# fill the gaps between end and next start
	for i in 2:length(regions)
		unaligned = regions[i].subjectfirst - regions[i-1].subjectlast
		
	end
	
end

# working with the anchors
function occupancy(alignment::BioAlignments.PairwiseAlignment)
	# relative amount of nucleotides in the query aligning to nucleotides in the reference 
	if isempty(alignment)
		return 0
	end

	anchors = alignment.a.aln.anchors
	referencelength = length(alignment.b)
	covered = 0

	@assert anchors[1].op == BioAlignments.OP_START

	for i in 2:length(anchors)
		if anchors[i].op == BioAlignments.OP_MATCH
			covered += anchors[i].seqpos - anchors[i-1].seqpos
		end
	end

	return covered / referencelength
end

# working with the cigar string
function occupancy(alignment::BioAlignments.PairwiseAlignment)
	if isempty(alignment)
		return 0
	end

	cigarstring = BioAlignments.cigar(alignment)

end

# working with the sequences
function occupancy(alignment::BioAlignments.PairwiseAlignment)
	referencelength = length(alignment.b)
	covered = length(ungap(alignment.a.seq)) - (length(alignment.a.seq) - length(alignment.b))
	println(string(covered) * " " * string(referencelength))
	return covered / referencelength
end