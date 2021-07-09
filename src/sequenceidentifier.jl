# Provides a type for sequence identifiers

mutable struct SequenceIdentifier
    id::String

    function SequenceIdentifier()
        return new("")
    end

    function SequenceIdentifier(id::String)
        return new(id)
    end
end

Base.length(id::SequenceIdentifier) = length(id.id)
Base.isempty(id::SequenceIdentifier) = length(id) >= 1
Base.firstindex(id::SequenceIdentifier) = 1
Base.lastindex(id::SequenceIdentifier) = length(id)
Base.eachindex(id::SequenceIdentifier) = Base.OneTo(lastindex(id))
hasspaces(id::SequenceIdentifier) = return ' ' in id.id

"""
    splitdescription(id; delimiter)

Split the provided `SequenceIdentifier` into two separate parts at the `delimiter`.

Example 1: `Drosophila_melanogaster@16S` becomes `Drosophila` and `16S`
"""
function splitdescription(id::SequenceIdentifier; delimiter='@')::Vector{String}
    if ! (delimiter in id.id)
        error("Missing species delimiter (\'$delimiter\') in description $description")
    end
    parts = split.(id.id, delimiter)
    return map(part -> string(part), parts)
end

"""
    otupart(id)

Returns the OTU part of the provided `SequenceIdentifier`

Example 1: `Drosophila_melanogaster@16S` becomes `Drosophila`
"""
function otupart(id::SequenceIdentifier)
    return first(splitdescription(id))
end

"""
    sequencepart(id)

Returns the sequence part (i.e., the identifier without the OTU) of the provided
`SequenceIdentifier`.

Example 1: `Drosophila_melanogaster@16S` becomes `16S`
"""
function sequencepart(id::SequenceIdentifier)
    return last(splitdescription(id))
end