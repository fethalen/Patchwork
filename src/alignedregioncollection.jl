# Collection of zero or more AlignedRegions

#include("alignedregion.jl")
#include("diamond.jl")

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

    function AlignedRegionCollection(
        reference::SequenceRecord,
        records::Vector{AlignedRegion}
    )
        return new(reference, records)
    end

    function AlignedRegionCollection(
        reference::SequenceRecord,
        results::Vector{DiamondSearchResult}
    )
        regions = []
        for result in results
            push!(regions, AlignedRegion(result))
        end
        return new(reference, regions)
    end
end

function AlignedRegionCollection(
    msa::MultipleSequenceAlignment,
    searchresults::Array{DiamondSearchResult}
)
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

function Base.push!(
    regions::AlignedRegionCollection,
    region::AlignedRegion
)
    return push!(regions.records, region)
end

function Base.pop!(
    regions::AlignedRegionCollection
)
    return pop!(regions.records)
end

function Base.getindex(
    regions::AlignedRegionCollection,
    index::Integer
)
    return getindex(regions.records, index)
end

Base.isempty(regions::AlignedRegionCollection) = length(regions) < 1
Base.firstindex(regions::AlignedRegionCollection) = 1
Base.lastindex(regions::AlignedRegionCollection) = length(regions)
Base.eachindex(regions::AlignedRegionCollection) = Base.OneTo(lastindex(regions))

function Base.iterate(regions::AlignedRegionCollection)
    isempty(regions) && return nothing
    i = firstindex(regions)
    return getindex(regions, i), i + 1
end

function Base.iterate(regions::AlignedRegionCollection, i::Int)
    i > lastindex(regions) && return nothing
    return getindex(regions, i), i + 1
end

# @inline function Base.iterate(
#     regions::AlignedRegionCollection,
#     i::Int = firstindex(regions)
# )
#     if i > lastindex(regions)
#         return nothing
#     else
#         return getindex(regions, i), i + 1
#     end
# end

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
    hasoverlaps(regions)

Returns `true` if there are overlapping sequences within `regions`.
"""
function hasoverlaps(
    regions::AlignedRegionCollection,
    sorted=false
)
    if !sorted
        regions = sort(regions)
    end

    for i in 2:lastindex(regions)
        if isoverlapping(regions[i], regions[i - 1])
            return true
        end
    end
    return false
end

function show_subjectintervals(regions::AlignedRegionCollection)
    for region in regions
        println(region.subjectfirst, " -> ", region.subjectlast)
    end
    return
end

"""
    mergeoverlaps(regions::AlignedRegionCollection, sorted::Bool=false,
                  iteration::Int64=1, maxiterate::Int64=100)

Returns an `AlignedRegionCollection` where overlapping regions have been merged into
their non-overlapping subset and are selected based on how well they align to the reference
sequence. Because sometimes multiple rounds of merging are needed, there is a limit to
how many times the function can run before stopping.
"""
function mergeoverlaps(
    regions::AlignedRegionCollection,
    sorted::Bool=false,
    iteration::Int64=1,
    maxiterate::Int64=100
)
    if !sorted
        regions = sort(regions)
    end

    length(regions) < 2 && return regions
    mergedregions = AlignedRegionCollection(regions.referencesequence)
    push!(mergedregions, first(regions))

    for i in 2:lastindex(regions)
        currentregion = regions[i]
        lastregion = last(mergedregions)
        if isoverlapping(currentregion, lastregion)
            pop!(mergedregions)
            for region in merge(currentregion, lastregion)
                if !isempty(region)
                    push!(mergedregions, region)
                end
            end
        else
            push!(mergedregions, currentregion)
        end
    end
    if hasoverlaps(mergedregions, false) && iteration <= maxiterate
        return mergeoverlaps(mergedregions, false, iteration + 1)
    end
    return sort(mergedregions)
end

"""
    sort(regions)

Sort all `AlignedRegion`s within `regions` based on their position, starting
with the leftmost region and ending with the rightmost region.
"""
function Base.sort(regions::AlignedRegionCollection)
    sortedregions = AlignedRegionCollection(regions.referencesequence)
    order = sortperm(map(region -> (region.subjectfirst, region.subjectlast), regions))
    # order = sortperm(map(region -> region.subjectlast, regions))
    for i in order
        push!(sortedregions, regions[i])
    end
    return sortedregions
end

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

function queryids(
    regions::AlignedRegionCollection
)::Vector{String}
    return map(region -> region.queryid.id, regions)
end