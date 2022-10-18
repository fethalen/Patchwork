using CSV
using DataFrames
using BioSequences

include("sequenceidentifier.jl")

const DATABASE = "database.dmnd"
const DIAMONDDB_EXT = "dmnd"
const FASTAEXTENSIONS = ["aln", "fa", "fn", "fna", "faa", "fasta", "FASTA"]
const FIELDS = ["qseqid", "qseq", "full_qseq", "sseqid", "sseq", "full_sseq",
                "cigar", "pident", "bitscore"]
const OUTPUT_FORMAT = [6; FIELDS]

mutable struct DiamondSearchResult
    queryid::SequenceIdentifier
    querysequence::LongAminoAcidSeq
    full_querysequence::LongAminoAcidSeq
    subjectid::SequenceIdentifier
    subjectsequence::LongAminoAcidSeq
    full_subjectsequence::LongAminoAcidSeq
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
function readblastTSV(path::AbstractString)::Array{DiamondSearchResult,1}
    results = CSV.File(path; header=FIELDS, delim='\t') |> DataFrame
    unique!(results)

    diamondsearchresults = []
    for row in eachrow(results)
        queryid = SequenceIdentifier(String(row.qseqid))
        subjectid = SequenceIdentifier(String(row.sseqid))
        result = DiamondSearchResult(
            queryid, LongAminoAcidSeq(row.qseq), LongAminoAcidSeq(row.full_qseq),
            subjectid, LongAminoAcidSeq(row.sseq), LongAminoAcidSeq(row.full_sseq),
            row.cigar, row.pident, row.bitscore
        )
        push!(diamondsearchresults, result)
    end
    return diamondsearchresults
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
    db_file = outdir * "/" * DATABASE
    logfile = outdir * "/diamond_makedb.log"
    makedb_cmd = pipeline(`diamond makedb --in $reference $flags -d $db_file`,
                          stdout=logfile, stderr=logfile)
    run(makedb_cmd)
    return db_file
end

"""
Runs `diamond blastp` on the provided `query` and `subject` paths. May include optional
`flags` such as `["--threads", 6]`.
"""
function diamond_blastp(
    query::AbstractString,
    subject::AbstractString,
    outdir::AbstractString,
    flags=[]
)::Array{DiamondSearchResult,1}
    logfile = outdir * "/diamond_blastp.log"
    blastout = outdir * "/blastp_out.tsv"
    diamond_cmd = pipeline(`diamond blastp --query $query --db $subject $flags
                            --outfmt $OUTPUT_FORMAT --out $blastout`, stdout=logfile,
                            stderr=logfile)
    run(diamond_cmd)
    return readblastTSV(blastout)
end

"""
Runs `diamond blastx` on the provided `query` and `subject` paths. May include optional
`flags` such as `["--threads", 6]`.
"""
function diamond_blastx(
    query::AbstractString,
    subject::AbstractString,
    outdir::AbstractString,
    flags=[]
)::Array{DiamondSearchResult,1}
    logfile = outdir * "/diamond_blastx.log"
    blastout = outdir * "/blastx_out.tsv"
    diamond_cmd = pipeline(`diamond blastx --query $query --db $subject $flags
                            --outfmt $OUTPUT_FORMAT --out $blastout`, stdout=logfile,
                            stderr=logfile)
    run(diamond_cmd)
    return readblastTSV(blastout)
end