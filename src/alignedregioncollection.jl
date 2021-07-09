# Collection of zero or more AlignedRegions

include("alignedregion.jl")
include("diamond.jl")
include("mafft.jl")

"""
    struct AlignedRegionCollection

- `records::Vector{AlignedRegion}`: the sequence record associated with the interval.
"""
mutable struct AlignedRegionCollection
    referencesequence::SequenceRecord
    records::Vector{AlignedRegion}

    "Empty initialization."
    function AlignedRegionCollection()
        return new(SequenceRecord(), [])
    end

    function AlignedRegionCollection(reference::SequenceRecord)
        return new(reference, [])
    end

    function AlignedRegionCollection(records::Vector{AlignedRegion})
        return new(SequenceRecord(), records)
    end

    function AlignedRegionCollection(reference::SequenceRecord, records::Vector{AlignedRegion})
        return new(reference, records)
    end

    function AlignedRegionCollection(reference::SequenceRecord, results::Vector{DiamondSearchResult})
        regions = []
        for result in results
            push!(regions, AlignedRegion(result))
        end
        return new(reference, regions)
    end
end

function AlignedRegionCollection(msa::MultipleSequenceAlignment,
                                 searchresults::Array{DiamondSearchResult})
    querycount = length(msa)
    merge!(msa, searchresults)
    hitcount = length(msa) - querycount
    alignment = mafft_linsi(msa, ["--thread", Sys.CPU_THREADS])
    regions = AlignedRegionCollection()

    for index in (querycount + 1):(querycount + hitcount)
        record = alignment.sequences[index]
        # search results don't contain alignment sequences
        result = searchresults[index - querycount]

        for compoundregion in compoundregions(record)
            (leftmost, rightmost) = nongap_range(compoundregion)
            compound_seqrecord = SequenceRecord(record.otu,
                *(record.identifier, '_', string(leftmost), '-',
                    string(rightmost)), compoundregion)
            region = AlignedRegion(compound_seqrecord, leftmost, rightmost,
                result.subjectframe, result.percentidentical)
            push!(regions, region)
        end
    end
    return regions
end

function Base.length(regions::AlignedRegionCollection)
    return length(regions.records)
end

function Base.push!(regions::AlignedRegionCollection, region::AlignedRegion)
    return push!(regions.records, region)
end

function Base.getindex(regions::AlignedRegionCollection, index::Integer)
    return getindex(regions.records, index)
end

Base.isempty(regions::AlignedRegionCollection) = length(regions) < 1
Base.firstindex(regions::AlignedRegionCollection) = 1
Base.lastindex(regions::AlignedRegionCollection) = length(regions)
Base.eachindex(regions::AlignedRegionCollection) = Base.OneTo(lastindex(regions))

@inline function Base.iterate(regions::AlignedRegionCollection, i::Int = firstindex(regions))
    if i > lastindex(regions)
        return nothing
    else
        return getindex(regions, i), i + 1
    end
end

"""
    sameids(regions)

Returns a dictionary where (unique) sequence identifiers of each region
contained within `regions` are key(s) and the value is a vector with elements
of the type `AlignedRegion`.
"""
function sameids(regions::AlignedRegionCollection)
    idgroups = Dict()
    for region in regions
        id = identifier(region)
        if !haskey(idgroups, id)
            idgroups[id] = [region]
        else
           push!(idgroups[id], region)
        end
    end
    return idgroups
end

"""
    uniquesequences(regions)

Returns a subset of `regions`, containing only unique sequences. If two or
more regions exists such that they have identical sequences, then keep the
`AlignedRegion` with the highest percentidentical.
"""
# TODO: uniquersequences doesn't work correctly, test with has_uniquesequences()
function uniquesequences(regions::AlignedRegionCollection)
    tovisit = collect(1:lastindex(regions))
    uniqueseqs = AlignedRegionCollection()
    while !isempty(tovisit)
        keep = popfirst!(tovisit)
        remove = []
        for i in tovisit
            regiona = regions[keep]
            regionb = regions[i]
            samesequence(regiona, regionb) || continue
            if regiona.percentidentical < regionb.percentidentical
                push!(remove, keep)
                keep = i
            else
                push!(remove, i)
            end
        end
        filter!(e -> e ∉ remove, tovisit)
        push!(uniqueseqs, regions[keep])
    end
    return uniqueseqs
end

"""
    has_uniquesequences(regions)

Check whether all of the sequences within `regions` are unique.
"""
function has_uniquesequences(regions::AlignedRegionCollection)
    for (regiona, regionb) in Iterators.product(regions, regions)
        regiona == regionb && continue
        if samesequence(regiona, regionb)
            return false
        end
    end
    return true
end

"""
    eachoverlap(regions)

Returns a `Tuple` of each pairs of `regions` which sequence's are overlapping.
"""
function eachoverlap(regions::AlignedRegionCollection)
    overlaps = []
    for (regiona, regionb) in Iterators.product(regions, regions)
        if isoverlapping(regiona, regionb)
            push!(overlaps, (regiona, regionb))
        end
    end
    return overlaps
end

"""
    hasoverlaps(regions)

Returns `true` if there are overlapping sequences within `regions`.
"""
function hasoverlaps(regions::AlignedRegionCollection)
    return length(eachoverlap(regions)) > 0
end

"""
    mergeoverlapping(regions)

Returns an `AlignedRegionCollection` which is a subset of `regions` where
overlapping sequences has been merged into longer stretches.
"""
function mergeoverlapping(regions::AlignedRegionCollection, sorted=false)
    if !sorted
        regions = sort(regions)
    end

    length(regions) < 2 && return regions
    mergedregions = AlignedRegionCollection()
    push!(mergedregions, first(regions))

    for i in 2:lastindex(regions)
        if isoverlapping(regions[i], last(mergedregions))
            println("A IS overlapping with B")
            showinterval(regions[i])
            showinterval(last(mergedregions))
            println()
        else
            println("A is NOT overlapping with B")
            showinterval(regions[i])
            showinterval(last(mergedregions))
            println()
            push!(mergedregions, regions[i])
        end
    end

    return mergedregions
end

"""
    sort(regions)

Sort all `AlignedRegion`s within `regions` based on their position, starting
with the leftmost region and ending with the rightmost region.
"""
function Base.sort(regions::AlignedRegionCollection)
    sortedregions = AlignedRegionCollection()
    order = sortperm(map(region -> (leftposition(region), rightposition(region)), regions))
    for i in order
        push!(sortedregions, regions[i])
    end
    return sortedregions
end

"""
    removeredundant(regions)

Delete
"""
# function removeredundant(regions::AlignedRegionCollection)
# end

# TODO: function layout
# function layout(regions::AlignedRegionCollection)
# end

# TODO: function consensus
# function consensus(regions::AlignedRegionCollection)
#     order = sortperm(map(region -> (leftposition(region), rightposition(region)), uniqueregions))
#     for i in order
#         println(count)
#         regiona = uniqueregions[i]
#         regionb = uniqueregions[i + 1]
#         println((leftposition(regiona), rightposition(regiona)))
#         println(precedes(regiona, regionb))
#     end
# end

"""
    isnucleotide(region)

Returns true if this region collection's first `record` consists of nucleotides.
"""
function isnucleotide(regions::AlignedRegionCollection)
    return isnucleotide(regions[1])
end

"""
    isaminoacid(region)

Returns true if this region collection's first `record` consists of amino acids.
"""
function isaminoacid(regions::AlignedRegionCollection)
    return isaminoacid(regions[1])
end

function BioSequences.translate(regions::AlignedRegionCollection)
    translatedregions = []
    for region in regions
        push!(translatedregions, translate(region))
    end
    return translatedregions
end
