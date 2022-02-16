# reduce coverage/subsample reads from .fastq.gz file

# coverage estimate: RESPECT or GenomeScope2
# for D. gyrociliatus, coverage estimate is known
# reads are distributed randomly in the file, but in same order in _R1 and _R2 files
# --> be consistent for forward and reverse reads! 
# --> user gives original coverage and percentage of reads they want 
# --> compute new coverage based on read percentage and genome size/old coverage
# --> input gzipped files if possible, else just fastq files

# julia --project=. src/compile.jl . src/subsample_precompiled.jl build
module Subsample

using ArgParse
using CodecZlib
using StatsBase

function parse_parameters()
    overview = """
    Subsampling of reads from a FASTQ file.
    """
    settings = ArgParseSettings(description = overview,
        version = "0.1.0",
        add_version = false)
    @add_arg_table! settings begin
        "--r1"
			help = "Path to one or more sequences in FASTQ format. If `--2` is also 
					specified, this file contains the forward reads of a paired-end 
					sequencing run."
			required = true
			arg_type = String
			metavar = "PATH"
        "--r2"
			help = "Only applicable if `--1` was also specified. This file contains the 
					reverse reads of a paired-end sequencing run."
			required = false
			default = ""
			arg_type = String
			metavar = "PATH"
        "--out"
			help = "Prefix for output file(s). File type extension (.fq or .fq.gz) will be 
				added automatically."
			required = false
			arg_type = String
			metavar = "PATH"
		"--outdir"
			help = "Output directory."
			required = false
			default = "."
			arg_type = String
			metavar = "PATH"
        "--p"
			help = "Percentage of reads you want to retain. The actual number of reads in 
					output file(s) will be rounded. Overrides both `--count` and `--cov-out`."
			required = false
			arg_type = Float64
			metavar = "NUMBER"
		"--cov-out"
			help = "The sequencing coverage for yout output file(s). Requires the coverage 
					of the input file(s), to be given with `--cov`. Overrides `--count`."
			required = false
			arg_type = Int64
			metavar = "NUMBER"
        "--count"
			help = "Exact number of reads you want to retain. Use this option instead of 
					`--p`."
			required = false
			arg_type = Int64
			metavar = "NUMBER"
        "--cov"
			help = "Estimated coverage of your sequencing run."
			required = false
			default = -1.0
			arg_type = Float64
			metavar = "NUMBER"
		"--genome-size"
			help = "Estimated genome size (in bases)."
			required = false
			arg_type = Int64
			metavar = "NUMBER"
		"--log"
			help = "Store information about this subsampling run."
			required = false
			arg_type = String
			default = "subsample.log"
		"--rstart"
			help = "Select reads from a random position in the file(s). If not specified, 
					reads will be selected from the beginning of the file(s)."
			action = :store_true
		"--random"
			help = "Select all reads from random positions. If two (paired-end) input files
				are given, corresponding mate pairs will be sampled from the files."
			action = :store_true
		"--gzip-out"
			help = "Compress output file(s)."
			action = :store_true
		"--fasta-out"
			help = "Output file(s) in FASTA format."
			action = :store_true
		"--quiet"
			help = "Suppress output on stdout."
			action = :store_true
    end

    return ArgParse.parse_args(settings)
end

function tofastaheader(str::String)
	newstr = str[nextind(str, 1):lastindex(str)]
	newstr = replace(newstr, " " => "_")
	return ">" * newstr
end

function subsample(
	file_r1::String, 
	file_r2::String, 
	outfile_r1::String, 
	outfile_r2::String,
	count::Int64, 
	miss::UnitRange, 
	fastaout::Bool=false,
	compress::Bool=true
)
	misslines = (4 * first(miss)):(4 * last(miss))
	endlines = last(misslines) + 1 + 4*count
	reader_r1 = isgzipcompressed(file_r1) ? GzipDecompressorStream(open(file_r1)) : open(file_r1)
    writer_r1 = compress ? GzipCompressorStream(open(outfile_r1, "w")) : open(outfile_r1, "w")
	reader_r2 = isgzipcompressed(file_r2) ? GzipDecompressorStream(open(file_r2)) : open(file_r2)
	writer_r2 = compress ? GzipCompressorStream(open(outfile_r2, "w")) : open(outfile_r2, "w")
    written = 0
	iter = enumerate(zip(eachline(reader_r1), eachline(reader_r2)))
	# @time begin
		for (line, (r1, r2)) in iter
			if !in(line, misslines)
				if !fastaout
					write(writer_r1, r1 * "\n")
					write(writer_r2, r2 * "\n")
				elseif 1 <= line%4 <= 2 # fastaout
					if line%4 == 1
						r1 = tofastaheader(r1)
						r2 = tofastaheader(r2)
					end
					write(writer_r1, r1 * "\n")
					write(writer_r2, r2 * "\n")
				end
				if line%4 == 0 # 1 record is 4 lines
					written += 1
				end
			end
			if line == endlines
				break
			end
		end
    # end
    close(reader_r1)
    close(writer_r1)
	close(reader_r2) 
	close(writer_r2)
	@assert written == count "written $written != $count"
	return written
end

function subsample(
	file::String, 
	outfile::String, 
	count::Int64, 
	miss::UnitRange, 
	fastaout::Bool=false,
	compress::Bool=true
)
	misslines = (4 * first(miss)):(4 * last(miss))
	endlines = last(misslines) + 1 + 4*count
	reader = isgzipcompressed(file) ? GzipDecompressorStream(open(file)) : open(file)
    writer = compress ? GzipCompressorStream(open(outfile, "w")) : open(outfile, "w")
    written = 0
	iter = enumerate(eachline(reader))
	# @time begin
		for (line, r) in iter
			if !in(line, misslines)
				if !fastaout
					write(writer, r * "\n")
				elseif 1 <= line%4 <= 2 # fastaout
					if line%4 == 1
						r = tofastaheader(r)
					end
					write(writer, r * "\n")
				end
				if line%4 == 0 # 1 record is 4 lines
					written += 1
				end
			end
			if line == endlines
				break
			end
		end
    # end
    close(reader)
    close(writer)
	@assert written == count "written $written != $count"
	return written
end

function subsample(
	file::String, 
	outfile::String,
	count::Int64, 
	positions::Vector{Int64}, 
	fastaout::Bool=true,
	compress::Bool=true
)
	#positions = sort(pos)
	reader = isgzipcompressed(file) ? GzipDecompressorStream(open(file)) : open(file)
    writer = compress ? GzipCompressorStream(open(outfile, "w")) : open(outfile, "w")
    written = 0
	iter = enumerate(eachline(reader))
	i = 1
	writing = false
	# @time begin
		for (line, r) in iter
			if !writing 
				if written == length(positions)
					break
				elseif line == 4 * positions[written+1] - 3 
					writing = true
				else
					continue
				end
			end
			if writing
				if !fastaout
					write(writer, r * "\n")
				elseif 1 <= line%4 <= 2
					if line%4 == 1
						r = tofastaheader(r)
					end
					write(writer, r * "\n")
				end
				if line%4 == 0 # 1 record is 4 lines
					written += 1
					#i += 1
					writing = false
				end
			end
		end
    # end
    close(reader)
    close(writer)
	@assert written == count "written $written != $count"
	return written
end

function subsample(
	file::String, 
	outfile::String, 
	records::Int64, 
	miss::UnitRange, 
	compress::Bool=true
)
	lines = records * 4 # 4 lines are 1 record
	if last(miss) + 1 > 1 # < 5 minutes
		@assert !isequal(miss, 1:0)
		beforemiss = (first(miss)-1) * 4
		remaining = lines - beforemiss # (records - (first(miss)-1)) * 4
		endposition = (last(miss)+1) * 4 + remaining
		run(pipeline(`head -n $beforemiss`, stdin=pipeline(`gunzip -c $file`), stdout=outfile))
		run(pipeline(`tail -n $remaining`, stdin=pipeline(`gunzip -c $file`, `head -n $endposition`), stdout=outfile, append=true))
	else # < 1 minute
		run(pipeline(`head -n $lines`, stdin=pipeline(`gunzip -c $file`) , stdout=outfile)) 
	end
	if compress
		run(`gzip -q $outfile`) # < 11 minutes
	end
end

function computecoverage(file::AbstractString, genomesize::Int64)
	io = GzipDecompressorStream(open(file))
	readcount = 0 #convert(Int64, countlines(io) / 4) # zcat ERR4013119.sra_1.fastq.gz | wc -l 
	avglength = 0
	for (i, line) in enumerate(eachline(io))
		i%4 != 2 && continue
		avglength += length(line)
		readcount += 1
	end
	close(io)
	avglength /= readcount
	return readcount * avglength / genomesize
end

function computecoverage(avglength::Int64, readcount::Int64, genomesize::Int64)
	return readcount * avglength / genomesize
end

function main()
	args = parse_parameters()

	file_r1 = args["r1"]
    file_r2 = args["r2"]
	paired = !isempty(file_r2)
    out = args["out"]
	outdir = args["outdir"]
	logwriter = open(outdir * "/" * args["log"], "w")
	p = args["p"]
	records = args["count"]
	covout = args["cov-out"]
	genomesize = args["genome-size"]
	keeprecords = records
	randstart = args["rstart"]
	random = args["random"]
	cov = args["cov"]
	compress = args["gzip-out"]
	fastaout = args["fasta-out"]
	loud = !args["quiet"]

	if count([isnothing(p), isnothing(records), isnothing(covout)]) != 2 #!xor(isnothing(p), isnothing(records))
		println("Please provide either a percentage, an output coverage, or an exact number of reads you ",
				"would like to retain.")
		return
	elseif !isnothing(p) && (p < 0.0 || p > 1.0)
			println("Read Percentage should be between 0 and 1. For an exact read number, ",
					"please use `--count`.")
			return
	elseif !isnothing(covout) && covout <= 0
			prinltn("Please provide a positive integer for `--cov-out`.")
			return
	elseif !isnothing(records) && records < 0
			println("Please provide a positive integer for `--count`.")
			return
	end

	if isnothing(out)
		if !isnothing(covout)
			suffix = "_" * string(covout) * "x"
		else
			suffix = "_" * (isnothing(p) ? string(records) : string(p))
		end
		prefix1 = outdir * "/" * last(split(file_r1, "/"))
		outfile_r1 = prefix1 * suffix * (fastaout ? ".fa" : ".fq")
		if paired
			prefix2 = outdir * "/" * last(split(file_r2, "/"))
			outfile_r2 = prefix2 * suffix * (fastaout ? ".fa" : ".fq")
		else
			outfile_r2 = ""
		end
	else
		prefix = outdir * "/" * out
		if paired
			outfile_r1 = prefix * (fastaout ? "_1.fa" : "_1.fq")
    		outfile_r2 = prefix * (fastaout ? "_2.fa" : "_2.fq")
		else
			outfile_r1 = prefix * (fastaout ? ".fa" : ".fq")
			outfile_r2 = ""
		end
	end
	# if compress # only julia-internal subsample() without calling external/bash functions
	# 	outfile_r1 = outfile_r1 * ".gz"
	# 	if paired
	# 		outfile_r2 = outfile_r2 * ".gz"
	# 	end
	# end

	println(logwriter, "Subsample.")
	if paired 
		println(logwriter, "Forward strand input file: $file_r1.")
		println(logwriter, "Reverse strand input file: $file_r2.")
	else
		println(logwriter, "Input file: $file_r1.")
	end
	println(logwriter)

	s = isgzipcompressed(file_r1) ? GzipDecompressorStream(open(file_r1)) : open(file_r1)
	if !isnothing(genomesize) && cov <= 0.0 # compute coverage from genome size
		#convert(Int64, countlines(io) / 4) # zcat ERR4013119.sra_1.fastq.gz | wc -l : 
		countrecords = 0 
		avglength = 0
		for (i, line) in enumerate(eachline(s))
			i%4 != 2 && continue
			avglength += length(line)
			countrecords += 1
		end
		avglength /= countrecords
		cov = computecoverage(avglength, countrecords, genomesize)
		println(logwriter, "Genome size is $genomesize.")
		println(logwriter, "Average read length in $file_r1 is $avglength.")
	else
		countrecords = convert(Int64, countlines(s) / 4) # < 5 minutes for ~182 Million reads
	end
	close(s)
	println(logwriter, "Your input file(s) contain(s) $countrecords reads (each).") 
	# Will always be printed when genome size is printed: 
	cov >= 0.0 && println(logwriter, "The estimated coverage of your sequencing run is $cov.")

	if !isnothing(p)
		loud && println("Computing absolute number of reads to be subsampled (from each ", 
			"input file) for percentage $p...")
		keeprecords = round(Int64, p * countrecords)
    	power = 10^(floor(Int64, log10(keeprecords))) # order of magnitude
    	records = ceil(Int64, keeprecords / power) * power # round up to next "step" 
	elseif !isnothing(covout)
		loud && println("Computing absolute number of reads to be subsampled (from each ", 
			"input file) for resulting coverage $covout...")
		if cov <= 0.0
			println("The option `--cov-out` requires either a coverage estimate of your ", 
				"sequencing run, to be provided with `--cov`, or a genome size estimate, ", 
				"to be provided with `--genome-size`.")
			return
		elseif covout >= cov
			println("Your output coverage is equal to or larger than the estimated/computed ", 
				"coverage of your sequencing run.")
			return
		end
		# cov / countrecords == covout / records
		records = round(Int64, covout * countrecords / cov)
		keeprecords = records # for the next check, but shouldn't happen if covout < cov
	end
	if records >= countrecords # possible side effect of always rounding up
		if keeprecords < countrecords # keep exact read number instead of rounded 
			records = keeprecords
		else
			str = *("The provided/computed read number $records is equal to or larger than ", 
				"the input file.")
			println(str)
			println(logwriter, str)
			close(logwriter)
			return
		end
	end
	println(logwriter, "Subsampling $records reads from (each of) the input file(s).")

	if random
		loud && println("Subsampling $records reads randomly from (each of) the input ", 
			"files(s)")
		println(logwriter, "Selecting reads randomly from the input file(s).")
		positions = sort(sample(1:countrecords, records, replace = false))
		subsample(file_r1, outfile_r1, records, positions, fastaout, compress)
		if paired
			subsample(file_r2, outfile_r2, records, positions, fastaout, compress)
		end
	else
		if randstart
			start = rand(1:countrecords)
			loud && println("Subsampling $records reads from a random starting position in ", 
				"the input file(s)...")
			println(logwriter, "The random starting position in the file(s) is $start.")
		else
			start = 1
			loud && 
				println("Subsampling $records reads from the beginning of the input file(s)...")
			println(logwriter, "The starting position in the file(s) is 1.")
		end
		# skip the records at these indices in the file:
		miss = start > countrecords - records + 1 ? ((start-records):(start-1)) : (1:start-1)
		randstart && println(logwriter, "Skipping reads ", first(miss), " to ", last(miss), 
			", both included.")
		# @time begin
			subsample(file_r1, outfile_r1, records, miss, fastaout, compress)
			if paired
				subsample(file_r2, outfile_r2, records, miss, fastaout, compress)
			end
		# end # ~ 30 minutes for paired, gzip being the most time-intensive step.
	end

	# if paired
	# 	@time records = subsample(file_r1, file_r2, outfile_r1, outfile_r2, records, miss, 
	#		compress) # ~ 40 minutes
	# else
	# 	@time records = subsample(file_r1, outfile_r1, records, miss, compress)
	# end

	println(logwriter)

	if compress
		str = "Subset of $file_r1 saved to $outfile_r1.gz."
		str2 = "Subset of $file_r2 saved to $outfile_r2.gz."
	else
		str = "Subset of $file_r1 saved to $outfile_r1."
		str2 = "Subset of $file_r2 saved to $outfile_r2."
	end
	if loud
		println(str)
		paired && println(str2)
		println()
	end
	println(logwriter, str)
	paired && println(logwriter, str2)
	println(logwriter)

	str = "Number of Retained Reads: $records"
    loud && println(str)
	println(logwriter, str)
	if cov > 0
        # cov / countrecords == newcovestimate / keeprecords
        newcovestimate = cov * records / countrecords
		str = "New Coverage Estimate: $newcovestimate"
		loud && println(str)
		println(logwriter, str)
    end
	close(logwriter)
end

function julia_main()::Cint
    try
        main()
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end

    return 0
end

if length(ARGS) >= 1
    julia_main()
end

# coverage = numreads * avgreadlength / genomesize

# output: num reads, estimate coverage for num reads

# 10 % reads for Christoph's beetle assembly
# to reduce NUMTs --> integration of mito sequences into nuclear genome in assemblies, 
# and that's shit you know.
# --> dropping the coverage gets rid of NUMTs (we expect at least)
end
