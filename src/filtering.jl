using BioSequences
using DataFrames

include("multiplesequencealignment.jl")
include("alignedregion.jl")

function remove_duplicates!(searchresults::DataFrame)::DataFrame
    unique!(searchresults, [:qstart, :qend])
    return searchresults
end

function remove_duplicates(msa::MultipleSequenceAlignment)
    collection = []
    for record in msa.sequences
        record.otu == "Ceratonereis_australis" || continue
        first = 0
        last = 0
        for (index, position) in enumerate(record.sequencedata)
            if first == 0 && ! isgap(position)
                first = index
            elseif last == 0 && first != 0 && isgap(position)
                last = index
            end
        end
        println("first: ", first, " last: ", last)
        push!(collection, AlignedRegion(record, first, last, 1, 100.0))
    end
    return collection
end

function remove_duplicates(
    msa::MultipleSequenceAlignment,
    searchresults::BLA
)
    sotu = "Ceratonereis_australis"
    sequences = queryalignment

    for result in searchresults
        addalignment!(queryalignment, SequenceRecord(sotu, result.subjectid, result.subjectsequence))
    end
    return sequences::MultipleSequenceAlignment
end
