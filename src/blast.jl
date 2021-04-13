# Wrapper for running BLAST+ from the command-line.

using CSV
using DataFrames
using BioSequences

include("multiplesequencealignment.jl")
include("sequencerecord.jl")

HEADER = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
          "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
          "qseq", "qframe", "sseq"]
OUTPUT_FORMAT = [6; HEADER]

struct BLASTSearchResult
    queryid::String
    queryotu::String
    subjectid::String
    subjectotu::String
    percentidentical::Float64
    alignmentlength::Int64
    mismatches::Int64
    gapopenings::Int64
    querystart::Int64
    queryend::Int64
    subjectstart::Int64
    subjectend::Int64
    evalue::Float64
    bitscore::Float64
    subjectsequence::LongAminoAcidSeq
    queryframe::Int64
    querysequence::LongDNASeq
end

"""
    splitdescription(description; delimiter)

Split the provided `description` into two separate parts at the `delimiter`.

Example 1: `Drosophila_melanogaster@16S` becomes `Drosophila` and `16S`
"""
function splitdescription(description::AbstractString; delimiter='@')::Vector{String}
    if ! (delimiter in description)
        error("Missing species delimiter (\'$delimiter\') in description $description")
    end
    parts = split.(description, delimiter)
    return map(part -> string(part), parts)
end

function splitdescription(descriptions::Vector{String}; delimiter='@')::Vector{Vector{String}}
    return [splitdescription(description; delimiter) for description in descriptions]
end

"""
    readblastTSV(path, min_percentidentity)

Read the contents of a tabular BLAST output into an array of `BLASTSearchResult`s.
Filters non-unique results and results with less percent identity than 
`min_percentidentity`. Adhear to the following BLAST `-outfmt`:

`6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
 qframe sseq seq`
"""
function readblastTSV(path::AbstractString; delimiter='@')::Array{BLASTSearchResult,1}
    results = CSV.File(path; header=HEADER, delim='\t') |> DataFrame
    unique!(results)

    qsplits = [(qotu=first, qid=last) 
        for (first, last) in splitdescription(Vector{String}(results.qseqid); delimiter)] |>
        DataFrames.DataFrame
    ssplits = [(sotu=first, sid=last) 
        for (first, last) in splitdescription(Vector{String}(results.sseqid); delimiter)] |> 
        DataFrames.DataFrame
    select!(results, Not([:qseqid, :sseqid]))
    results = hcat(qsplits, ssplits, results)

    blastsearchresults = []
    for row in eachrow(results)
        result = BLASTSearchResult(
            row.qid, row.qotu, row.sid, row.sotu, row.pident, row.length, row.mismatch, 
            row.gapopen, row.qstart, row.qend, row.sstart, row.send, row.evalue, 
            row.bitscore, LongAminoAcidSeq(row.sseq), row.qframe, LongDNASeq(row.qseq))
        push!(blastsearchresults, result)
    end
    return blastsearchresults
end

"""
    writeblastTSV(path, results; delimiter)

Write `results` to a TSV file to `path`, using the provided `delimiter` to separate
columns (default: `'\\t'`). No header is added when `header` is set to `false`.
"""
function writeblastTSV(path::AbstractString, results::Array{BLASTSearchResult,1}; 
                       delimiter='\t', header=false)::AbstractString
    CSV.write(path, select!(DataFrames.DataFrame(results), Not(:subjectotu)), 
              delim=delimiter, writeheader=header)
    return path
end

"""
    merge!(queryalignment::MultipleSequenceAlignment,
           searchresults::Array{BLASTSearchResult})

Insert the alignments within `searchresults` into `queryalignment`.
"""
function merge!(queryalignment::MultipleSequenceAlignment,
                searchresults::Array{BLASTSearchResult})::MultipleSequenceAlignment
    sotu = "Ceratonereis_australis"
    sequences = queryalignment

    for result in searchresults
        addalignment!(queryalignment, SequenceRecord(sotu, result.subjectid, result.subjectsequence))
    end
    return sequences::MultipleSequenceAlignment
end

"""
    tblastn(query, subject, flags; db)

Runs tblastn on `query` against `subject`. Subjects and queries may be file
names (as strings), or `MultipleSequenceAlignment`. May include optional
`flags` such as `["-num_threads", 4,]`. Output `format` must include names
of custom columns such as `[10, "qseqid", "sseqid", "length", "evalue"]`.
Hyphens ('-') are removed automatically before the BLAST search.
"""
function tblastn(query::AbstractString, subject::AbstractString,
                 flags=[], format=[]; db::Bool=false)::Array{BLASTSearchResult}
    output_format = join(format, " ")
    csv_file, csv_io = mktemp()
    write(csv_io, join(format[2:end], ',') * '\n')
    query_fasta = mktemp_fasta(query, removehyphens=true)

    if db
        tblastn_cmd = pipeline(`tblastn -query $query_fasta -db $subject $flags
                                -outfmt $output_format`, stdout=csv_file)
    else
        tblastn_cmd = pipeline(`tblastn -query $query_fasta -subject $subject
                                $flags -outfmt $output_format`,
                                stdout=csv_file)
    end
    run(tblastn_cmd)
    close(csv_io)
    searchresults = readblastCSV(csv_file)
    return searchresults
end

function tblastn(query::MultipleSequenceAlignment, subject::AbstractString,
                 flags=[], format=[]; db::Bool=false)::Array{BLASTSearchResult}
    output_format = join(format, " ")
    csv_file, csv_io = mktemp()

    write(csv_io, join(format[2:end], ',') * '\n')
    query_fasta = mktemp_fasta(query, removehyphens=true)

    if db
        tblastn_cmd = pipeline(`tblastn -query $query_fasta -db $subject $flags
                                -outfmt $output_format`, stdout=csv_file)
    else
        tblastn_cmd = pipeline(`tblastn -query $query_fasta -subject $subject $flags
                                -outfmt $output_format`, stdout=csv_file)
    end
    run(tblastn_cmd)
    close(csv_io)
    searchresults = readblastCSV(csv_file)
    return searchresults
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
    results_file, results_io= mktemp()
    write(results_file, join(HEADER, '\t') * '\n')
    # query_fasta = mktemp_fasta(query, removehyphens=true)

    diamond_cmd = pipeline(`diamond blastx --query $query --db $subject $flags
                            --outfmt $OUTPUT_FORMAT --out $results_file`)
    run(diamond_cmd)
    close(results_io)
    return results_file
end

function diamond_blastx(query::MultipleSequenceAlignment, subject::AbstractString,
                        flags=[])::AbstractString
    results_file, results_io = mktemp()
    write(results_file, join(HEADER, '\t') * '\n')
    # query_fasta = mktemp_fasta(query, removehyphens=true)
    diamond_cmd = pipeline(`diamond blastx --query $query --db $subject $flags
                            --outfmt $OUTPUT_FORMAT --out $results_file`)
    run(diamond_cmd)
    close(results_io)
    return results_file
end

"""
Uses DIAMOND to create a BLAST database from the provided `reference` FASTA
sequence. May include optional `flags` such as `["--taxonnodes"]`.
"""
function diamond_makeblastdb(reference::AbstractString, flags=[])
    # TODO
    db_file, db_io = mktemp()
    makedb_cmd = pipeline(`diamond makedb --in $reference -d $db_file $flags`)
    run(makedb_cmd)
    close(db_io)
    return db_file
end