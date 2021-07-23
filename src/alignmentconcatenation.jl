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

# first comes before second!
function concatenate(first::AlignedRegion, second::AlignedRegion)
	# regions is sorted, no overlaps remaining.
	if isempty(first) && isempty(second)
		empty = BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
		return AlignedRegion(empty, 0, -1, SequenceIdentifier(), 0, -1, 0) # empty AlignedRegion
	elseif isempty(first)
		return second
	elseif isempty(second)
		return first
	end

	# make sure that the entire subject sequence is covered? 
	# or insert gaps in reference and query for missing reference positions?
	# this should happen in the concatenate method for AlignedRegionCollections...

	@assert first.subjectlast == second.subjectfirst - 1
	@assert first.querylast == second.queryfirst - 1 # Else you would have to add gaps in between...
	# But that's supposed to happen in the concatenate(AligńedRegionCollection)
	@assert first.queryid == second.queryid # Is that so?
	# What about the frame?
	# It's hard to concatenate 2 regions if the queries are in different frames but the 
	# reference is in amino acids...

	joinedalignments = concatenate(first.pairwisealignment, second.pairwisealignment)
	# TODO: Think about this, what about the frame?
	return AlignedRegion(joinedalignments, first.subjectfirst, second.subjectlast, 
						 first.queryid, first.queryfirst, second.querylast, 
						 first.queryframe)	
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
	return countmatches(alignment) / length(alignment.b)
end

function(region::AlignedRegion)
	return occupancy(region.pairwisealignment)
end

# number of positions where query residues align to reference residues
function countmatches(alignment::BioAlignments.PairwiseAlignment)
	if isempty(alignment)
		return 0
	end

	covered = 0
	anchors = alignment.a.aln.anchors

	@assert anchors[1].op == BioAlignments.OP_START

	for i in 2:length(anchors)
		if anchors[i].op == BioAlignments.OP_MATCH
			covered += anchors[i].seqpos - anchors[i-1].seqpos
		end
	end
	
	return covered
end

# REDUNDANT
# working with the sequences (== working with the gap anchors instead of matches)
function occupancy(alignment::BioAlignments.PairwiseAlignment)
	return (length(alignment.b) - countgaps(alignment)) / length(alignment.b)
end

# count gaps in query
function countgaps(alignment::BioAlignments.PairwiseAlignment)
	if isempty(alignment)
		return 0
	end

	anchors = alignment.a.aln.anchors
	gaps = 0

	@assert anchors[1].op == BioAlignments.OP_START

	#gaps = sum([anchors[i].refpos - anchors[i-1].refpos 
	#			for i in 2:length(anchors) if anchors[i].op == BioAlignments.OP_DELETE])
	for i in 2:length(anchors)
		if anchors[i].op == BioAlignments.OP_DELETE
			gaps += (anchors[i].refpos - anchors[i-1].refpos)
		end
	end
	return gaps
end