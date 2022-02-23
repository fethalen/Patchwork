using FASTX
using CodecZlib
using BioSequences
using BioAlignments
using Patchwork

# Subsample reads: 
# julia --project=/home/clara/Patchwork/scripts/Subsample /home/clara/Patchwork/scripts/Subsample/src/Subsample.jl --r1 ~/Data/felix2_EKDL190133942-1a_HVN23DSXX_L1_1.fq.gz --r2 ~/Data/felix2_EKDL190133942-1a_HVN23DSXX_L1_2.fq.gz --count 40000000 --random --outdir ~/Data 

# no internal gaps allowed --> look only at diagonals in alignment matrix
# only diagonals that are long enough to meet the threshold t are considered
# if inside such a diagonal, it is found that the sequences do not align well enough to meet t, 
# the alignment for that diagonal is stopped and we continue with the next diagonal.
# alignment scoring is very basic bc the sequences are evaluated according to identity 
# --> exact match +1, mismatch +-0 (for now, can be adjusted or added to fct signature later)
# score(x,y) = score(x-1, y-1) + m or - mm depending 
# any(score(x,y) / length(shorter sequence) >= t) ? success : fail
function gapfreealign_t(s1::T, s2::T, t::Float64)::Bool where {T<:BioSequence}
    matrix = zeros(Int64, length(s1), length(s2)) # s1 ^= rows, s2 ^= columns; 0st row/col ^= gap
    m = 1
    mm = 0
    # pos is rownumber for s1, colnumber for s2
    shorter, pos = length(s1) <= length(s2) ? (length(s1), 1) : (length(s2), 2)
    thresholdposition = ceil(Int64, shorter * (1 - t))
    # iterate: thresholdpos row -> 1 -> thresholdpos col
    diagonals = [collect(zip(1:shorter, 1:shorter))] # init: middle diagonal
    for diagstart in 2:thresholdposition
        diagrow = collect(zip(diagstart:shorter, 1:shorter-diagstart))
        diagcol = collect(zip(1:shorter-diagstart, diagstart:shorter))
        push!(diagonals, diagrow, diagcol)
    end
    for diag in diagonals
        for coords in diag
            x = coords[1]
            y = coords[2]
            if x == 1 || y == 1 # no internal gaps allowed
                matrix[x, y] = isequal(s1[x], s2[y]) ? m : -mm 
            else # no internal gaps allowed
                matrix[x, y] = isequal(s1[x], s2[y]) ? matrix[x-1, y-1] + m : matrix[x-1, y-1] - mm 
            end
            if matrix[x, y] / shorter >= t # threshold met! 
                return true
            elseif (matrix[x, y] + shorter - coords[pos]) / shorter < t 
                break # mark this diagonal as done/look at next diagonal
            end # else continue in this diagonal because t can still be met
        end
    end
    return false
end

function gapfreealign_t(first::SequenceRecord, second::SequenceRecord, t::Float64)
    return gapfreealign_t(first.sequencedata, second.sequencedata, t)
end

###########################################################################################
# simple (gapfree-align) LINCLUST: 

# (MORE OR LESS) COPIED FROM LINCLUST SUPPLEMENTARY MATERIAL ------------------------------
rotateleft(val, bits) = xor(val << bits, val >> (16 - bits))

function circularhash(kmer::LongDNASeq, k::Int64)
    # fixed 16bit randoms
    rands = Dict(DNA_A => 0x4567, DNA_C => 0x23c6, DNA_G => 0x9869, DNA_T => 0x4873, 
        DNA_M => 0xdc51, DNA_R => 0x5cff, DNA_W => 0x944a, DNA_S => 0x58ec, DNA_Y => 0x1f29, 
        DNA_K => 0x7ccd, DNA_V => 0x58ba, DNA_H => 0xd7ab, DNA_D => 0x41f2, DNA_B => 0x1efb, 
        DNA_N => 0xa9e3, DNA_Gap => 0xe146)
    h = 0x0
    h = xor(h, rands[kmer[1]])
    for i in 2:k
        h = rotateleft(h, 5)
        h = xor(h, rands[kmer[i]])
    end
    return h
end

# rolling hash variant for previous hash function
function circularhash_next(seq::LongDNASeq, k::Int64, prevbase::DNA, h::UInt16)
    rands = Dict(DNA_A => 0x4567, DNA_C => 0x23c6, DNA_G => 0x9869, DNA_T => 0x4873, 
        DNA_M => 0xdc51, DNA_R => 0x5cff, DNA_W => 0x944a, DNA_S => 0x58ec, DNA_Y => 0x1f29, 
        DNA_K => 0x7ccd, DNA_V => 0x58ba, DNA_H => 0xd7ab, DNA_D => 0x41f2, DNA_B => 0x1efb, 
        DNA_N => 0xa9e3, DNA_Gap => 0xe146)
    h = xor(h, rotateleft(rands[prevbase], (5 * k) % 16))
    h = rotateleft(h, 5)
    h = xor(h, rands[seq[k]])
    return h
end
#------------------------------------------------------------------------------------------

mutable struct KmerIndex
    kmer::LongDNASeq
    seqlen::Int64
    kmerpos::Int64
    msaindex::Int64 # reference to id and sequence in corresponding msa

    function KmerIndex(kmer::LongDNASeq, seqlen::Int64, kmerpos::Int64, msaindex::Int64)
        return new(kmer, seqlen, kmerpos, msaindex)
    end
end

function KmerIndex()
    return KmerIndex(LongDNASeq(), 0, 0, 0)
end

function KmerIndex(kmer::LongDNASeq, seq::LongDNASeq, kmerpos::Int64, msaindex::Int64)
    return KmerIndex(kmer, length(seq), kmerpos, msaindex)
end

Base.isempty(kmerindex::KmerIndex) = isempty(kmerindex.kmer)

# split the file first (splitfile in fastq.jl) in filelength reads per tempfile!!!
# then you can use the following functions sequentially for doing the clustering: 

function buildkmertable(
    name::AbstractString, 
    reader::Union{FASTA.Reader, FASTQ.Reader}, 
    k::Int64, m::Int64, 
    filelength::Int64=10000000
)
    record = typeof(reader) == FASTA.Reader ? FASTA.Record() : FASTQ.Record()
    all_kmertable = repeat([Vector{KmerIndex}()], filelength)
    sequences = repeat([SequenceRecord()], filelength)
    eof(reader) && return kmertable
    msaindex = 1

    while !eof(reader)
        read!(reader, record)
        seq = FASTX.sequence(record)
        seqlen = length(seq)

        # no need to relocate the array for storing new records: array was preallocated
        sequences[msaindex] = SequenceRecord(FASTX.identifier(record), seq) 
        # looking for a max of m lowest-hashing kmers
        num_kmers = m + k - 1 <= seqlen ? m : seqlen - k + 1
        #kmertable = repeat([KmerIndex()], num_kmers)
        kmertable = Vector{KmerIndex}(undef, num_kmers)
        # [(kmer, hash, kmerpos),...]
        seq_kmers = Vector{Tuple{LongDNASeq, UInt16, Int64}}(undef, seqlen - k + 1)
        seq_kmers[1] = (seq[1:k], circularhash(seq[1:k], k), 1)

        for i in 2:seqlen - k + 1 # compute hashes for all kmers
            seq_kmers[i] = (seq[i:i+k-1],
                circularhash_next(seq, i+k-1, seq[i-1], seq_kmers[i-1][2]), i)
        end

        if length(seq_kmers) > m # store the m lowest-hashing kmers in the index table
            for i in 1:m
                minhash, idx = findmin(map(tup -> tup[2], seq_kmers))
                kmertable[i] = KmerIndex(seq_kmers[idx][1], seqlen, seq_kmers[idx][3], 
                    msaindex)
                if minhash == 0xffff
                    deleteat!(seq_kmers, idx)
                else # avoid relocating m * ~l array entries where possible
                    seq_kmers[idx] = (seq_kmers[idx][1], 0xffff, seq_kmers[idx][3]) 
                end
            end
        else # only <= m kmers in total
            for (i, tup) in enumerate(seq_kmers)
                kmertable[i] = KmerIndex(tup[1], seqlen, tup[3], msaindex)
            end
        end
        all_kmertable[msaindex] = kmertable
        msaindex += 1
    end
    sequences = MultipleSequenceAlignment(name, sequences)
    all_kmertable = collect(Iterators.flatten(all_kmertable)) # [...] instead of [[...],...]
    return sequences, all_kmertable
end

function buildkmertable(
    file::AbstractString, 
    k::Int64, m::Int64, 
    filelength::Int64=10000000
)
	if isfastafile(file) 
		reader = isgzipcompressed(file) ? FASTA.Reader(GzipDecompressorStream(open(file))) : FASTA.Reader(open(file))
	elseif isfastqfile(file)
		reader = isgzipcompressed(file) ? FASTQ.Reader(GzipDecompressorStream(open(file))) : FASTQ.Reader(open(file))
	else
		println("File must be either in FASTA or FASTQ format (gzipped input is supported).")
		return MultipleSequenceAlignment(), Vector{KmerIndex()}()
	end
	return buildkmertable(file, reader, k, m, filelength)
end

function mktemp_kmertable(kmertable::Vector{KmerIndex})
    # sort by kmer block, sequence length, sequence id
    # --> get kmertable with blocks of kmers, inside of which the longest/center seq is at the top
    # one kmer block corresponds to exactly one center sequence, so the whole thing is also ordered 
    # by center sequence right from the beginning
    sort!(kmertable, by = kmerindex -> (kmerindex.kmer, -kmerindex.seqlen, kmerindex.msaindex))
    tmpfile, tmpio = mktemp()
    if isempty(kmertable)
        close(tmpio)
        return tmpfile
    end
    writestart = 1
    currentcenter = kmertable[writestart] # ===
    i = 2
    countlines = 0
    while i <= length(kmertable)
        if !isequal(kmertable[i].kmer, currentcenter.kmer) # entered new kmer block in sorted kmer table
            if writestart < i-1 # don't write trivial/1-sequence blocks
                write(tmpio, string(currentcenter.msaindex)) # write center'\t's1'\t's2'\t's3'\t'...'\t'sn'\n'
                for j in writestart+1:i-1 # write the entire previous group to file
                    if !isempty(kmertable[j]) # reduntant entries were marked empty instead of deleted
                        write(tmpio, "\t" * string(kmertable[j].msaindex))
                    end
                end
                write(tmpio, "\n") # group end --> new line
                countlines += 1
            end
            writestart = i # current position is new kmer group start
            currentcenter = kmertable[writestart] # current position is new center sequence; ===
            i += 1
        else # still in same kmer block, sequence is not the current center
            kmertable[i].kmerpos = currentcenter.kmerpos - kmertable[i].kmerpos # overwrite absolute with relative position
            if kmertable[i].msaindex == kmertable[i-1].msaindex # same sequence has more than one matching kmer with center
                # keep only entry with smaller relative kmer pos
                # make sure to keep currentcenter (writestart), if current center seq has same kmer in diff pos
                if kmertable[i].kmerpos >= kmertable[i-1].kmerpos || writestart == i-1 
                    # deleteat!(kmertable, i)
                    kmertable[i] = KmerIndex() # avoid relocating the array
                else
                    # deleteat!(kmertable, i - 1)
                    kmertable[i-1] = KmerIndex()
                end
            # else
            #     i += 1 # no entry deleted --> next seq == next index
            end
            i += 1 # entries aren't deleted but just emptied, so it's always next seq == next index
        end
    end
    if writestart < i-1 # don't write trivial/1-sequence blocks
        write(tmpio, string(currentcenter.msaindex)) # write center'\t's1'\t's2'\t's3'\t'...'\t'sn'\n'
        for j in writestart+1:i-1 # write the entire previous group to file
            if !isempty(kmertable[j]) # reduntant entries were marked empty instead of deleted
                write(tmpio, "\t" * string(kmertable[j].msaindex))
            end
        end
        write(tmpio, "\n") # group end --> new line
        countlines += 1
    end
    close(tmpio)
    return tmpfile, countlines # number of centers == number of lines
end

# ungapped alignment filter
function gapfreealign_tocenter(kmertable_file::AbstractString, filelength::Int64, msa::MultipleSequenceAlignment, t::Float64)
    # edges = Vector{Int64}[] # [[center idx, s1 idx, s2 idx, ...],...]
    edges = repeat([Vector{Int64}()], filelength) # [[center idx, s1 idx, s2 idx, ...],...]
    i = 1
    open(kmertable_file) do reader
        while !eof(reader)
            @assert i <= filelength "Some error in counting lines"
            newcluster = true
            positions = map(x -> parse(Int64, x), split(readline(reader), "\t")) # center msaindex, sequence msaindex
            for seqindex in positions[2:end]
                aln = gapfreealign_t(msa[positions[1]], msa[seqindex], t)
                if aln
                    if newcluster 
                        # push!(edges, [positions[1], seqindex])
                        edges[i] = [positions[1], seqindex]
                        newcluster = false
                    else
                        # push!(edges[end], seqindex)
                        push!(edges[i], seqindex)
                    end
                end
            end
            if !newcluster # at least one alignment was found for this index --> index++
                i += 1
            end
        end
    end
    sort!(edges, by = vec -> length(vec), rev = true) # sort by decreasing cluster size
    tmpfile, tmpio = mktemp()
    # for edge in edges
    #     write(tmpio, join(string.(edge), "\t") * "\n")
    # end
    i = 1
    while !isempty(edges[i])
        write(tmpio, join(string.(edges[i]), "\t") * "\n")
        i += 1
    end
    close(tmpio)
    return tmpfile
end

function greedyincremental_clustering(edgesfile::AbstractString, msa::MultipleSequenceAlignment)
    keepers = repeat([true], length(msa))
    open(edgesfile) do reader
        while !eof(reader)
            edge = map(x -> parse(Int64, x), split(readline(reader), "\t"))
            # if this center doesn't belong to a different, larger cluster: remove all its neighbours
            if keepers[edge[1]] 
                for seqindex in edge[2:end]
                    keepers[seqindex] = false
                end
            end # else, keep them: they're orphaned, the cluster center was already removed
        end
    end
    msa = MultipleSequenceAlignment(msa.name, msa[findall(keepers)])
    return msa
end

# THIS WOULD BE THE WHOLE WORKFLOW: #######################################################

#file = "/scratch/users/koehne19/patchwork_benchmarking/data/ERR4013119_Dimorphilus_gyrociliatus_reads/subsampled_reads/ERR4013119.sra_1.fastq.gz_40000000.fq.gz"
#println("count records:")
#@time filelength = round(Int, countlines(GzipDecompressorStream(open(file))) / 4)
#k = 14
#m = 20
#threshold = 0.95

#println("build kmer table:")
#@time msa, evolvingclusters = buildkmertable(file, k, m, filelength) 
## 202.322562 seconds (1.69 G allocations: 106.001 GiB, 28.16% gc time)
#println("length of input msa: ", length(msa)) # 1000000
#println("length of kmertable: ", length(evolvingclusters)) # m * length(msa) = 20000000
#if !isempty(msa) # if isempty, print errormessage and return.
#	println("save kmer table:")
#	@time evolvingclusters, numlines = mktemp_kmertable(evolvingclusters) # file that contains all possible center/seq edges 
#	# 107.737633 seconds (20.89 M allocations: 960.595 MiB, 14.66% gc time)
#	println("gapfree alignment filter:")
#	@time evolvingclusters = gapfreealign_tocenter(evolvingclusters, numlines, msa, threshold) # file that contains all clusters center->s1,s2,s3..., sorted by cluster size
#	# 226.887146 seconds (142.79 M allocations: 996.792 GiB, 28.36% gc time, 0.20% compilation time)
#	println("greedy incremental clustering:")
#	@time msa = greedyincremental_clustering(evolvingclusters, msa) 
#	# 1.302488 seconds (6.69 M allocations: 546.244 MiB, 10.39% gc time, 12.41% compilation time)
#	println("length of output msa: ", length(msa)) # 890927
#end

# THAT STILL TAKES TOO LONG, BUT HEY IT'S BETTER THAN BEFORE AT LEAST.
# --> ~=10 minutes for clustering 1000000 sequences
# --> 10 minutes * 182 files = 1820 minutes ~= 30 hrs
# --> + 45 minutes for splitting the primary input file ~= 31 hrs
# --> for both pe-files, 62 hrs...
# or maybe, can you keep only the _2 reads that belong to the surviving _1 reads?
# (instead of doing the clustering for the second set all over again...)
