# Wrapper for running BLAST+ from the command-line.

using CSV
using DataFrames
using BioSequences

const DATABASE = "database.dmnd"
const DIAMONDDB_EXT = "dmnd"
const FASTAEXTENSIONS = ["aln", "fa", "fn", "fna", "faa", "fasta", "FASTA"]
const FIELDS = ["qseqid", "qseq", "full_qseq", "qseq_translated", "qstart",
                "qend", "qframe", "sseqid", "sseq", "sstart", "send", "cigar", "pident",
                "bitscore"]
const OUTPUT_FORMAT = [6; FIELDS]

mutable struct DiamondSearchResult
    queryid::SequenceIdentifier
    translated_querysequence::LongAA
    full_querysequence::LongDNA # ATTENTION: modified to include frameshifts, rev. comp.
    querystart::Int64
    queryend::Int64
    queryframe::Int64
    subjectid::SequenceIdentifier
    subjectsequence::LongAA
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
    file, io = mktemp()
    run(pipeline(`sed 's/\\/f/g' $path`, file)) # b for backwards frameshift
    run(pipeline(`sed 's/\//b/g' $path`, file)) # f for forwards frameshift

    results = CSV.File(file; header=FIELDS, delim='\t') |> DataFrame
#    results = CSV.File(path; header=FIELDS, delim='\t') |> DataFrame
    unique!(results)

    diamondsearchresults = []
    for row in eachrow(results)
        queryid = SequenceIdentifier(String(row.qseqid))
        subjectid = SequenceIdentifier(String(row.sseqid))
        frameshifted_fullqseq, newcigar, newqstart, newqstop, frame = frameshift_fullseq(
            row.full_qseq, row.cigar, row.qstart, row.qend, row.qframe)
        result = DiamondSearchResult(
            queryid, BioSequences.LongAA(row.qseq_translated), frameshifted_fullqseq, 
            newqstart, newqstop, frame, subjectid, BioSequences.LongAA(row.sseq),
            row.sstart, row.send, newcigar, row.pident, row.bitscore
        )
        push!(diamondsearchresults, result)
    end

    close(io)
    return diamondsearchresults
end

function frameshift(
    sequence::Union{BioSequences.LongDNA, AbstractString},
    cigar::AbstractString,
    frame::Int64
)
    if typeof(sequence) <: AbstractString
        sequence = LongDNA{4}(sequence)
    end
    if frame < 0
        return frameshift(reverse_complement(sequence), cigar)
    else
        return frameshift(sequence, cigar)
    end
end

function frameshift(
    sequence::Union{BioSequences.LongDNA, AbstractString},
    cigar::AbstractString
)::Tuple{BioSequences.LongDNA, AbstractString}
    cigar = replace(cigar, "\\" => "f", "/" => "b") # see DIAMOND Wiki on GitHub, \ = +1, / = -1 in translation direction
    anchors = collect(eachmatch(r"[fbMIDNSHP=X]", cigar)) # b for backwards, f for forwards frameshift
    positions = collect(eachmatch(r"\d+", cigar))
    index = 1
    seqbuffer = IOBuffer()

    for (p, anchor) in zip(positions, anchors)
        pos = parse(Int64, p.match)
        # anchor corresponding to pos * 1 nucleotide
        if isequal(anchor.match, "f")
            index += pos
        elseif isequal(anchor.match, "b")
            index -= pos
        # anchor corresponding to pos * nucleotide triplet (if DELETE, don't change index): (delete = gap in ref)
        elseif !BioAlignments.isdeleteop(BioAlignments.Operation(anchor.match[1]))
            print(seqbuffer, sequence[index:index + 3 * pos - 1])
            index += 3 * pos
        end
    end

    tmpcigar = replace(cigar, r"\d+b" => "", r"\d+f" => "") # delete frameshifts
    tmpanchors = collect(eachmatch(r"[MIDNSHP=X]", tmpcigar))
    tmppositions = collect(eachmatch(r"\d+", tmpcigar))
    newcigar = ""
    poscount = 0
    # if there are consecutive same operations after rm frameshifts, bundle them together
    for i in 1:lastindex(tmpanchors)-1
        poscount += parse(Int64, tmppositions[i].match)
        if !isequal(tmpanchors[i].match, tmpanchors[i+1].match)
            newcigar = *(newcigar, string(poscount), tmpanchors[i].match)
            poscount = 0
        end
    end
    poscount += parse(Int64, last(tmppositions).match)
    newcigar = *(newcigar, string(poscount), last(tmpanchors).match)

    return (LongDNA{4}(String(take!(seqbuffer))), newcigar)
end

function frameshift_fullseq(fullqseq::AbstractString, 
    cigar::AbstractString, tstart::Int64, tstop::Int64, frame::Int64
)::Tuple{LongDNA, AbstractString, Int64, Int64, Int64}
    if frame >= 0 # 0 = unspecified frame, assume forward strand...
        frameshifted, newcigar = frameshift(fullqseq[tstart:tstop], cigar)
        patched = fullqseq[firstindex(fullqseq):tstart-1] * 
            String(frameshifted) * fullqseq[tstop+1:lastindex(fullqseq)]
        newtstop = tstart + length(frameshifted) - 1
        return (BioSequences.LongDNA{4}(patched), newcigar, tstart, newtstop, frame)
    else # in case of negative reading frames, apply cigar ops to reverse complement
        revtstart = length(fullqseq) - tstart + 1
        revtstop = length(fullqseq) - tstop + 1
        revcomp = reverse_complement(LongDNA{4}(fullqseq))
        frameshifted, newcigar = frameshift(revcomp[revtstart:revtstop], cigar)
        revcomp_patched = LongDNA{4}(join([revcomp[firstindex(fullqseq):revtstart-1], 
            frameshifted, revcomp[revtstop+1:lastindex(fullqseq)]]))
        newrevtstop = revtstart + length(frameshifted) - 1
        return (revcomp_patched, newcigar, revtstart, newrevtstop, -frame) # return reverse complement
        # --> after this, assume positive frame!
    end
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
    header=false,
    omit::Vector{Symbol}=Symbol[]
)::AbstractString
    # ATTENTION: this function operates on the DiamondSearchResult Array, i.e. 
    # the full_querysequence is modified (includes frameshifts and is potentially 
    # reverse-complement of "original" DNA seq)
    # this results in frame being always positive
    # cigar should not contain frameshifting ops
    dataframe = select!(DataFrames.DataFrame(results), Not(omit))
    if !in(:queryid, omit)
        dataframe[!, :queryid] = map(result -> result.queryid.id, results)
    end
    if !in(:subjectid, omit)
        dataframe[!, :subjectid] = map(result -> result.subjectid.id, results)
    end
    # not really necessary; frameshifts should be absent from modified cigar string
    if !in(:cigar, omit) 
        dataframe[!, :cigar] = map(c -> replace(c, "b" => "/", "f" => "\\"), 
            dataframe[!, :cigar])
    end
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
        push!(queryalignment, SequenceRecord(sotu, result.subjectid, result.subjectsequence))
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

function queryids(
    results::Vector{DiamondSearchResult}
)::Vector{SequenceIdentifier}
    return map(result -> result.queryid, results)
end

function subjectids(
    results::Vector{DiamondSearchResult}
)::Vector{SequenceIdentifier}
    return map(result -> result.subjectid, results)
end

function isdiamonddatabase(path::AbstractString)::Bool
    splits = split(path, ".")
    length(splits) > 1 && isequal(last(splits), DIAMONDDB_EXT) && return true
    return false
end
