# Wrapper for running BLAST+ from the command-line.

using CSV
using DataFrames
using BioSequences

include("multiplesequencealignment.jl")
include("sequencerecord.jl")
include("sequenceidentifier.jl")

const FIELDS = ["qseqid", "qseq", "full_qseq", "qstart", "qend", "qframe", "sseqid",
                "sseq", "sstart", "send", "cigar", "pident", "bitscore"]
const OUTPUT_FORMAT = [6; FIELDS]

mutable struct DiamondSearchResult
    queryid::SequenceIdentifier
    querysequence::LongDNASeq
    full_querysequence::LongDNASeq
    querystart::Int64
    queryend::Int64
    queryframe::Int64
    subjectid::SequenceIdentifier
    subjectsequence::LongAminoAcidSeq
    subjectstart::Int64
    subjectend::Int64
    cigar::AbstractString
    percentidentical::Float64
    bitscore::Float64
end

"""
    readblastTSV(path, min_percentidentity)

Read the contents of a tabular BLAST output into an array of `DiamondSearchResult`s.
Filters non-unique results and results with less percent identity than
`min_percentidentity`. Adhear to the following BLAST `-outfmt`:

`6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
 qframe sseq seq`
"""
function readblastTSV(path::AbstractString; delimiter='@')::Array{DiamondSearchResult,1}
    results = CSV.File(path; header=FIELDS, delim='\t') |> DataFrame
    unique!(results)

    # qsplits = [(qotu=first, qid=last)
    #     for (first, last) in splitdescription(Vector{String}(results.qseqid); delimiter)] |>
    #     DataFrames.DataFrame
    # ssplits = [(sotu=first, sid=last)
    #     for (first, last) in splitdescription(Vector{String}(results.sseqid); delimiter)] |>
    #     DataFrames.DataFrame
    # select!(results, Not([:qseqid, :sseqid]))
    # results = hcat(qsplits, ssplits, results)

    diamondsearchresults = []
    for row in eachrow(results)
        queryid = SequenceIdentifier(row.qseqid)
        subjectid = SequenceIdentifier(row.sseqid)
        result = DiamondSearchResult(
            queryid, BioSequences.LongDNASeq(row.qseq),
            BioSequences.LongDNASeq(row.full_qseq), row.qstart, row.qend, row.qframe,
            subjectid, BioSequences.LongAminoAcidSeq(row.sseq), row.sstart, row.send,
            row.cigar, row.pident, row.bitscore)
        push!(diamondsearchresults, result)
    end
    return diamondsearchresults
end

"""
    writeblastTSV(path, results; delimiter)

Write `results` to a TSV file to `path`, using the provided `delimiter` to separate
columns (default: `'\\t'`). No header is added when `header` is set to `false`.
"""
function writeblastTSV(path::AbstractString, results::Array{DiamondSearchResult,1};
                       delimiter='\t', header=false)::AbstractString
    CSV.write(path, select!(DataFrames.DataFrame(results), Not(subjectotu)),
              delim=delimiter, writeheader=header)
    return path
end

"""
    merge!(queryalignment::MultipleSequenceAlignment,
           searchresults::Array{DiamondSearchResult})

Insert the alignments within `searchresults` into `queryalignment`.
"""
function merge!(queryalignment::MultipleSequenceAlignment,
                searchresults::Array{DiamondSearchResult})::MultipleSequenceAlignment
    sotu = "Ceratonereis_australis"
    sequences = queryalignment

    for result in searchresults
        addalignment!(queryalignment, SequenceRecord(sotu, result.subjectid, result.subjectsequence))
    end
    return sequences::MultipleSequenceAlignment
end

"""
    diamond_blastx(query, subject, flags; db)

Runs DIAMOND BLASTX on `query` against `subject`. Subjects and queries may be
file names (as strings), or `MultipleSequenceAlignment`. May include optional
`flags` such as `["-num_threads", 4,]`. Hyphens ('-') are removed
automatically before the BLAST search.
"""
function diamond_blastx(query::AbstractString, subject::AbstractString,
                        flags=[])::AbstractString
    results_file, results_io = mktemp()
    write(results_file, join(FIELDS, '\t') * '\n')
    diamond_cmd = pipeline(`diamond blastx --query $query --db $subject $flags
                            --outfmt $OUTPUT_FORMAT --out $results_file`)
    run(diamond_cmd)
    close(results_io)
    return results_file
end

function diamond_blastx(query::MultipleSequenceAlignment, subject::AbstractString,
                        flags=[])::AbstractString
    results_file, results_io = mktemp()
    write(results_file, join(FIELDS, '\t') * '\n')
    diamond_cmd = pipeline(`diamond blastx --query $query --db $subject $flags
                            --outfmt $OUTPUT_FORMAT --out $results_file`)
    run(diamond_cmd)
    close(results_io)
    return results_file
end

"""
Runs `diamond makedb` on the provided `reference` FASTA sequence. May include
optional `flags` such as `["--taxonnodes"]`.
"""
function diamond_makeblastdb(reference::AbstractString, flags=[])
    if isdiamonddatabase(reference)
        return 
    end
    db_file, db_io = mktemp()
    if isfastafile(reference)
        makedb_cmd = pipeline(`diamond makedb --in $reference -d $db_file $flags`)
    else # BLAST DB
        makedb_cmd = pipeline(`diamond prepdb -d $reference`)
    end
    run(makedb_cmd)
    close(db_io)
    return db_file
end

function queryid(result::DiamondSearchResult, speciesdelimiter='@')::String
    return *(result.queryotu, speciesdelimiter, result.queryid)
end

function subjectid(result::DiamondSearchResult, speciesdelimiter='@')::String
    return *(result.subjectotu, speciesdelimiter, result.subjectid)
end

function queryids(results::Vector{DiamondSearchResult}, speciesdelimiter='@')::Vector{String}
    return map(result -> queryid(result), results)
end

function subjectids(results::Vector{DiamondSearchResult}, speciesdelimiter='@')::Vector{String}
    return map(result -> subjectid(result), results)
end

function isfastafile(path::AbstractString)::Bool
    splits = split(path, ".")
    if length(splits) > 1
        extension = last(splits)
        if extension in FASTAEXTENSIONS
            return true
        end
    end
    return false
end

function isdiamonddatabase(path::AbstractString)::Bool
    splits = split(path, ".")
    if length(splits) > 1 && isequal(last(splits), DIAMONDDB)
            return true
        end
    end
    return false
end