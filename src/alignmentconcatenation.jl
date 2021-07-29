import BioAlignments
import BioSequences

pwd()
cd("src")

include("alignedregion.jl")
include("alignedregioncollection.jl") 

subject = "../test/07673_Alitta_succinea.fa"
diamondresults = "../test/c_australis_x_07673.tsv"
hits = Patchwork.readblastTSV(diamondresults)
full_subjectseq = Patchwork.get_fullseq(subject)
regions = Patchwork.AlignedRegionCollection(full_subjectseq, hits)
regions = sort(regions)
nonoverlapping = Patchwork.AlignedRegionCollection(regions.referencesequence, [regions[1]])
alignments = [region.pairwisealignment for region in regions]

for i in 2:length(regions)
	overlapping = false
	for j in 1:length(nonoverlapping)
		if Patchwork.isoverlapping(nonoverlapping[j], regions[i])
			overlapping = true
		end
	end
	if !overlapping && regions[i].subjectfirst <= 354 && regions[i].subjectlast <= 354
		push!(nonoverlapping, regions[i])
	end
end

sortedregions = Patchwork.AlignedRegionCollection(nonoverlapping.referencesequence)
order = sortperm(map(region -> (subject_leftposition(region), subject_rightposition(region)), nonoverlapping))
for i in order
    push!(sortedregions, regions[i])
end
sortedregions

superalignment = concatenate(sortedregions)

function BioAlignments.cigar(alignment::BioAlignments.PairwiseAlignment)
	anchors = alignment.a.aln.anchors
	out = IOBuffer()
	if isempty(anchors)
		return ""
	end
	@assert anchors[1].op == BioAlignments.OP_START
	for i in 1:length(anchors)-1
		positions = max(anchors[i+1].seqpos - anchors[i].seqpos,
						anchors[i+1].refpos - anchors[i].refpos)
		print(out, positions, anchors[i+1].op)
	end
	return String(take!(out))
end

# Concatenating 2 PairwiseAlignment objects: 
# Join the 2 queries and the 2 references together; 
# No need to add anything in between, right?
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
	cigar = join([BioAlignments.cigar(first.a.aln), BioAlignments.cigar(second.a.aln)], "")
	return BioAlignments.PairwiseAlignment(joinedquery, joinedreference, cigar)
end

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
	#cigar = join([BioAlignments.cigar(alignment) for alignment in alignments], "")
	out = IOBuffer()
	for alignment in alignments
		print(out, BioAlignments.cigar(alignment))
	end
	
	return BioAlignments.PairwiseAlignment(joinedquery, joinedreference, String(take!(out)))
end

# AlignedRegionCollection ###############################################################################################################################

# TODO: Benchmark!

function concatenate(regions::Patchwork.AlignedRegionCollection)
	# regions is sorted, no overlaps remaining.
	if isempty(regions)
		return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
	elseif length(regions) == 1
		return regions[1].pairwisealignment
	end
	
	alignments = [regions[1].pairwisealignment]
	# fill the gaps between end and next start
	# this requires reallocating memory every time for the alignments vector
	# would it be better to reallocate memory every time for the concatenation step 
	# and omit the vector by incrementally building the alignmnent? 
	for i in 2:length(regions)
		@assert regions[i-1].subjectlast < regions[i].subjectfirst		 	# check sorted
		#firstotu = splitdescription(regions[i-1].queryid)
		#secondotu = splitdescription(regions[i].queryid)
		#@assert firstotu == secondotu 										# Is that so?
		if regions[i].subjectfirst > regions[1].subjectlast + 1
			#gapquery = repeat("-", regions[i].subjectfirst - regions[i-1].subjectlast)
			gapcigar = string(regions[i].subjectfirst - regions[i-1].subjectlast) * "D"
			reference = regions.referencesequence.sequencedata[(regions[i-1].subjectlast + 1):(regions[i].subjectfirst)]
			gapquery = typeof(reference)()
			bridgealignment = BioAlignments.PairwiseAlignment(gapquery, reference, gapcigar)
			push!(alignments, bridgealignment, regions[i].pairwisealignment)
		else
			push!(alignments, regions[i].pairwisealignment)
		end
	end
	return concatenate(alignments)
end

function concatenate(regions::Patchwork.AlignedRegionCollection)
	# regions is sorted, no overlaps remaining.
	if isempty(regions)
		return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
	elseif length(regions) == 1
		return regions[1].pairwisealignment
	end
	
	superalignment = regions[1].pairwisealignment
	# fill the gaps between end and next start
	# this requires reallocating memory every time for the superalignment
	for i in 2:length(regions)
		@assert !(regions[i-1].subjectlast < regions[i].subjectfirst - 1) 	# check sorted
		firstotu = splitdescription(regions[i-1].queryid)
		secondotu = splitdescription(regions[i].queryid)
		@assert firstotu == secondotu 										# Is that so?
		if regions[i].subjectfirst > regions[1].subjectlast + 1
			gapquery = repeat("-", regions[i].subjectfirst - regions[i-1].subjectlast)
			gapcigar = string(regions[i].subjectfirst - regions[i-1].subjectlast) * "D"
			reference = regions.referencesequence[(regions[i-1].subjectlast + 1)
												   :regions[i].subjectfirst]
			bridgealignment = BioAlignments.PairwiseAlignment(gapquery, reference, gapcigar)
			superalignment = concatenate([superalignment, bridgealignment, 
										  regions[i].pairwisealignment])
		else
			superalignment = concatenate(superalignment, regions[i].pairwisealignment)
		end
	end
	return superalignment
end

# UNUSED RIGHT NOW ######################################################################################################################################

# first comes before second!
function concatenate(first::Main.Patchwork.AlignedRegion, 
					 second::Main.Patchwork.AlignedRegion)
	# regions is sorted, no overlaps remaining.
	if isempty(first) && isempty(second)
		return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
	elseif isempty(first)
		return second
	elseif isempty(second)
		return first
	end

	# make sure that the entire subject sequence is covered? 
	# or insert gaps in reference and query for missing reference positions?
	# this should happen in the concatenate method for AlignedRegionCollections...
	@assert first.subjectlast == second.subjectfirst - 1
	@assert first.querylast == second.queryfirst - 1 
	# Else you would have to add gaps in between...
	# But that's supposed to happen in the concatenate(AligńedRegionCollection)
	firstotu = splitdescription(regions[i-1].queryid)
	secondotu = splitdescription(regions[i].queryid)
	@assert firstotu == secondotu 										# Is that so?

	return concatenate(first.pairwisealignment, second.pairwisealignment)
end

function concatenate(regions::AbstractVector{Main.Patchwork.AlignedRegion})
	if isempty(regions)
		return BioAlignments.PairwiseAlignment(LongSequence(), LongSequence(), "")
		#return nothing
	elseif length(regions) == 1
		return regions[1].pairwisealignment
	end

	for i in 2:length(regions)
		#@assert regions[i-1].subjectlast == regions[i].subjectfirst - 1
		#@assert regions[i-1].querylast == regions[i].queryfirst - 1 
		firstotu = splitdescription(regions[i-1].queryid)
		secondotu = splitdescription(regions[i].queryid)
		@assert firstotu == secondotu 										# Is that so?
	end
	alignments = [region.pairwisealignment for region in regions]

	return concatenate(alignments)
end

# OCCUPANCY #############################################################################################################################################
# relative amount of nucleotides in the query aligning to nucleotides in the reference 

# working with the anchors
function occupancy(alignment::BioAlignments.PairwiseAlignment)
	return countmatches(alignment) / length(alignment.b)
end

function occupancy(region::AlignedRegion)
	return occupancy(region.pairwisealignment)
end

# REDUNDANT
# working with the sequences (== working with the gap anchors instead of matches)
function occupancy(alignment::BioAlignments.PairwiseAlignment)
	return (length(alignment.b) - countgaps(alignment)) / length(alignment.b)
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

# count gaps in query
function countgaps(alignment::BioAlignments.PairwiseAlignment)
	if isempty(alignment)
		return 0
	end

	gaps = 0
	anchors = alignment.a.aln.anchors

	@assert anchors[1].op == BioAlignments.OP_START

	for i in 2:length(anchors)
		if anchors[i].op == BioAlignments.OP_DELETE
			gaps += (anchors[i].refpos - anchors[i-1].refpos)
		end
	end
	return gaps
end