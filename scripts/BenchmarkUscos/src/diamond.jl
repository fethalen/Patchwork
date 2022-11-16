using CSV
using DataFrames
using BioSequences

# include("sequenceidentifier.jl")

const DATABASE = "database.dmnd"
const DIAMONDDB_EXT = "dmnd"
const FASTAEXTENSIONS = ["aln", "fa", "fn", "fna", "faa", "fasta", "FASTA"]
const FIELDS = ["qseqid", "qseq", "full_qseq", "sseqid", "sseq", "full_sseq",
                "cigar", "pident", "bitscore"]
const FIELDS_DNAQUERY = ["qseqid", "qseq", "qseq_translated", "full_qseq", "sseqid", "sseq", "full_sseq",
                "cigar", "pident", "bitscore"]
const OUTPUT_FORMAT = [6; FIELDS]
const OUTPUT_FORMAT_DNAQUERY = [6; FIELDS_DNAQUERY]

mutable struct DiamondSearchResult
    queryid::SequenceIdentifier
    querysequence::LongAA
    full_querysequence::LongAA
    subjectid::SequenceIdentifier
    subjectsequence::LongAA
    full_subjectsequence::LongAA
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
function readblastTSV(path::AbstractString; dnaquery::Bool=false)::Array{DiamondSearchResult,1}
    if !dnaquery 
        results = CSV.File(path; header=FIELDS, delim='\t') |> DataFrame
    else
        results = CSV.File(path; header=FIELDS_DNAQUERY, delim='\t') |> DataFrame
    end
    unique!(results)

    diamondsearchresults = []
    for row in eachrow(results)
        queryid = SequenceIdentifier(String(row.qseqid))
        subjectid = SequenceIdentifier(String(row.sseqid))
        if !dnaquery
            result = DiamondSearchResult(
                queryid, LongAA(row.qseq), LongAA(row.full_qseq),
                subjectid, LongAA(row.sseq), LongAA(row.full_sseq),
                row.cigar, row.pident, row.bitscore
            )
        else
            # if no frameshift: 
            # @assert length(row.qseq) % 3 == 0 """
            #     qseq:\n $(row.qseq)\nlength:\n $(length(row.qseq)).\n 
            #     Not divisible by three. 
            #     """

            # ASK FELIX: QUERYSEQUENCE TRANSLATION BY JULIA AND DIAMOND ARE DIFFERENT; 
            # WHICH ONE TO USE? 
            # IF DIAMOND, NEED TO "PATCH" THE FULL QUERY SEQUENCE TOGETHER FROM THE 
            # DIAMOND-TRANSLATED PART + JULIA-TRANSLATED UNALIGNED REGIONS 

            # # querysequence = BioSequences.translate(LongDNA{4}(row.qseq))    # this should work, providing no frameshifting (default) happened in the DIAMOND BLASTX search

            # # @assert isequal(String(querysequence), String(row.qseq_translated)) """
            # #     julia translation:\n $querysequence\nDIAMOND translation:\n $(row.qseq_translated).\n
            # #     DNA query:\n $(row.qseq).\n
            # #     """

            # # if length(row.full_qseq) % 3 != 0
            # #     trimmed = trim(row.full_qseq, row.qseq)
            # #     @assert length(trimmed) % 3 == 0 """
            # #         trimmed full_qseq:\n $(trimmed)\nlength:\n $(length(trimmed)).\n 
            # #         Not divisible by three. 
            # #         """
            # #     fullquery = BioSequences.translate(LongDNA{4}(trimmed))
            # # else
            # #     fullquery = BioSequences.translate(LongDNA{4}(row.full_qseq))   # (default) happened in the DIAMOND BLASTX search
            # # end

            # querysequence = translate_fromcigar(LongDNA{4}(row.qseq), row.cigar)  # julia-translation if frameshifting is allowed in DIAMOND
            # fullquery = patch_translate_fullquery(row.full_qseq, row.qseq, querysequences)
            querysequence = LongAA(row.qseq_translated)
            fullquery = patch_translate_fullquery(row.full_qseq, row.qseq, row.qseq_translated)

            result = DiamondSearchResult(
                queryid, querysequence, fullquery,
                subjectid, LongAA(row.sseq), LongAA(row.full_sseq),
                row.cigar, row.pident, row.bitscore
            )
        end
        push!(diamondsearchresults, result)
    end
    return diamondsearchresults
end

"""
For DNA queries, trim the full query sequence to a length that is divisible by three.
Ensure that the frame corresponds to the frame of the aligned part of the sequence. 
"""
function trim(full_qseq::AbstractString, qseq::AbstractString)::AbstractString
    subseqidx = findfirst(qseq, full_qseq)
    if isnothing(subseqidx)
        println("Query:")
        println(qseq)
        println("Full Query:")
        println(full_qseq)
        error("Query is not a subsequence of Full Query.")
    end
    startidx = (first(subseqidx) % 3) == 0 ? 3 : (first(subseqidx) % 3)
    
    #newlength = length(full_qseq) - startidx + 1
    #endidx = (last(subseqidx) == lastindex(full_qseq)) ? last(subseqidx) : length(full_qseq) - (newlength%3)

    lengthafterqseq = lastindex(full_qseq) - last(subseqidx)
    @assert lengthafterqseq >= 0 """
        query:\n $(qseq)\n 
        full query:\n $(full_qseq)\n 
        Query ends after Full Query. 
        """

    endidx = last(subseqidx) + (lengthafterqseq - lengthafterqseq%3)
    #println("first: ", startidx, " last: ", endidx, " len -qseq: ", length(full_qseq[startidx:endidx])-length(qseq))
    return full_qseq[startidx:endidx]
end

"""
Translate the full query sequence. This function "patches together" the translation using the 
DIAMOND-translated `qseq_translated`, which corresponds to the aligned part of the sequence, 
and the julia-translated "rest" of the full query sequence (i.e. the parts that precede and 
succeed the part of the full sequence that corresponds to `qseq_translated`).
"""
function patch_translate_fullquery(full_qseq::AbstractString, qseq::AbstractString, qseq_translated::AbstractString)::LongAA
    translatedseq = LongAA(qseq_translated)
    return patch_translate_fullquery(full_qseq, qseq, translatedseq)
end

function patch_translate_fullquery(full_qseq::AbstractString, qseq::AbstractString, qseq_translated::LongAA)::LongAA
    if isequal(full_qseq, qseq) #length(qseq) >= length(full_qseq)
        return qseq_translated
    elseif length(full_qseq) <= length(qseq)
        println("Query:")
        println(qseq)
        println("Full Query:")
        println(full_qseq)
        error("lengths Query <= Full Query; and sequences not equal.")
    end
    
    old_seq = full_qseq
    full_qseq = trim(full_qseq, qseq)

    # substract length(qseq) from length(full_qseq) 
    # bc with frameshift, length(qseq) not necessarily divisible by 3
    # but length of rest of full_qseq should be divisible by 3
    @assert (length(full_qseq)-length(qseq)) % 3 == 0 """
        trimmed full_qseq:\n $(full_qseq)\nlength:\n $(length(full_qseq)).\n
        untrimmed full_qseq:\n $(old_seq)\nlength:\n $(length(old_seq)).\n 
        New length not divisible by three. 
        """

    translatedpart = findfirst(qseq, full_qseq)
    @assert !isnothing(translatedpart) """
        query:\n $(qseq)\n 
        full query:\n $(full_qseq)\n 
        Query is not a subsequence of full query. 
        """

    if first(translatedpart) > firstindex(full_qseq)
        firstpart = firstindex(full_qseq):(first(translatedpart)-1)
        firstseq = translate(LongDNA{4}(full_qseq[firstpart]))
    else
        firstseq = LongAA()
    end
    if last(translatedpart) < lastindex(full_qseq)
        lastpart = (last(translatedpart)+1):lastindex(full_qseq)
        lastseq = translate(LongDNA{4}(full_qseq[lastpart]))
    else
        lastseq = LongAA()
    end

    return firstseq * qseq_translated * lastseq
end

# in case frameshifting operations are allowed in DIAMOND (--framshift 15)
# not needed if patch_translate_fullquery() (i.e., using DIAMOND translation
# qseq_translated of aligned part qseq)
# would only be needed for julia-translating qseq
# function translate_fromcigar(
#     sequence::BioSequences.LongDNA,
#     cigar::AbstractString
# )::BioSequences.LongAA
#     anchors = collect(eachmatch(r"[/\\MIDNSHP=X]", cigar))
#     positions = collect(eachmatch(r"\d+", cigar))
#     index = 1
#     buffer = IOBuffer()

#     for (p, anchor) in zip(positions, anchors)
#         pos = parse(Int64, p.match)
#         # anchor corresponding to pos * 1 nucleotide
#         if isequal(anchor.match, "\\")
#             index += pos
#         elseif isequal(anchor.match, "/")
#             index -= pos
#         # anchor corresponding to pos * nucleotide triplet (if DELETE, don't change index):
#         elseif !BioAlignments.isdeleteop(BioAlignments.Operation(anchor.match[1]))
#             print(buffer, sequence[index:index + 3 * pos - 1])
#             index += 3 * pos
#         end
#     end

#     return BioSequences.translate(BioSequences.LongDNA(String(take!(buffer))))
# end

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
                            --outfmt $OUTPUT_FORMAT_DNAQUERY --out $blastout`, stdout=logfile,
                            stderr=logfile)         # in add. to qseq and full_qseq, save qseq_translated
    run(diamond_cmd)
    return readblastTSV(blastout, dnaquery=true)    # query sequence needs translating, compare to qseq_translated
end