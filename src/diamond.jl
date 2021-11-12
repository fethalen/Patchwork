# Wrapper for running BLAST+ from the command-line.

using CSV
using DataFrames
using BioSequences

include("multiplesequencealignment.jl")
include("sequencerecord.jl")
include("sequenceidentifier.jl")

const DATABASE = "database.dmnd"
const DIAMONDDB_EXT = "dmnd"
const FASTAEXTENSIONS = ["aln", "fa", "fn", "fna", "faa", "fasta", "FASTA"]
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
`min_percentidentity`. Adhere to the following BLAST `-outfmt`:

`6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
 qframe sseq seq`
"""
function readblastTSV(
    path::AbstractString;
    delimiter='@'
)::Array{DiamondSearchResult,1}
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
        queryid = SequenceIdentifier(String(row.qseqid))
        subjectid = SequenceIdentifier(String(row.sseqid))
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
function writeblastTSV(
    path::AbstractString,
    results::Array{DiamondSearchResult,1};
    delimiter='\t',
    header=false
)::AbstractString
    dataframe = select!(DataFrames.DataFrame(results), Not(:subjectid))
    dataframe[!, :queryid] = map(result -> result.queryid.id, results)
    CSV.write(path, dataframe, delim=delimiter, writeheader=header)
    return path
end

"""
    merge!(queryalignment::MultipleSequenceAlignment,
           searchresults::Array{DiamondSearchResult})

Insert the alignments within `searchresults` into `queryalignment`.
"""
function merge!(
    queryalignment::MultipleSequenceAlignment,
    searchresults::Array{DiamondSearchResult}
)::MultipleSequenceAlignment
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
function diamond_blastx(
    query::AbstractString,
    subject::AbstractString,
    outdir::AbstractString,
    flags=[]
)::AbstractString
    results_file, results_io = mktemp()
    logfile = outdir * "/diamond_blastx.log"
    diamond_cmd = pipeline(`diamond blastx --query $query --db $subject $flags
                            --outfmt $OUTPUT_FORMAT --out $results_file`, stdout=logfile,
                            stderr=logfile)
    run(diamond_cmd)
    close(results_io)
    return results_file
end

function diamond_blastx(
    query::MultipleSequenceAlignment,
    subject::AbstractString,
    outdir::AbstractString,
    flags=[]
)::AbstractString
    querypath = mktemp_fasta(query)
    return diamond_blastx(querypath, subject, outdir, flags)
end

function diamond_blastx(
    query::AbstractString,
    subject::MultipleSequenceAlignment,
    outdir::AbstractString,
    flags=[]
)::AbstractString
    subjectpath = mktemp_fasta(subject)
    return diamond_blastx(query, subjectpath, outdir, flags)
end

function diamond_blastx(
    query::MultipleSequenceAlignment,
    subject::MultipleSequenceAlignment,
    outdir::AbstractString,
    flags=[]
)::AbstractString
    querypath = mktemp_fasta(query)
    subjectpath = mktemp_fasta(subject)
    return diamond_blastx(querypath, subjectpath, outdir, flags)
end

"""
Runs `diamond makedb` on the provided `reference` FASTA sequence. May include
optional `flags` such as `["--taxonnodes"]`.
"""
function diamond_makeblastdb(
    reference::AbstractString,
    outdir::AbstractString,
    flags=[]
)::AbstractString
    if isdiamonddatabase(reference)
        return reference
    elseif isfastafile(reference)
        db_file = outdir * "/" * DATABASE
        logfile = outdir * "/diamond_makedb.log"
        makedb_cmd = pipeline(`diamond makedb --in $reference $flags -d $db_file`,
                              stdout=logfile, stderr=logfile)
        run(makedb_cmd)
        return db_file
    else # BLAST DB
        logfile = outdir * "/diamond_prepdb.log"
        makedb_cmd = pipeline(`diamond prepdb -d $reference`, stdout=logfile,
                               stderr=logfile)
        return reference
    end
end

function diamond_makeblastdb(
    alignment::MultipleSequenceAlignment,
    outdir::AbstractString,
    flags=[]
)::AbstractString
    reference = mktemp_fasta(alignment)
    db_file = outdir * "/" * DATABASE
    logfile = outdir * "/diamond_makedb.log"
    makedb_cmd = pipeline(`diamond makedb --in $reference $flags -d $db_file`,
                          stdout=logfile, stderr=logfile)
    run(makedb_cmd)
    return db_file
end

function queryid(
    result::DiamondSearchResult,
    speciesdelimiter='@'
)::String
    return *(result.queryotu, speciesdelimiter, result.queryid)
end

function subjectid(
    result::DiamondSearchResult,
    speciesdelimiter='@'
)::String
    return *(result.subjectotu, speciesdelimiter, result.subjectid)
end

function queryids(
    results::Vector{DiamondSearchResult},
    speciesdelimiter='@'
)::Vector{String}
    return map(result -> queryid(result), results)
end

function subjectids(
    results::Vector{DiamondSearchResult},
    speciesdelimiter='@'
)::Vector{String}
    return map(result -> subjectid(result), results)
end

function isfastafile(path::AbstractString)::Bool
    splits = split(path, ".")
    length(splits) > 1 && last(splits) in FASTAEXTENSIONS && return true
    return false
end

function isdiamonddatabase(path::AbstractString)::Bool
    splits = split(path, ".")
    length(splits) > 1 && isequal(last(splits), DIAMONDDB_EXT) && return true
    return false
end
