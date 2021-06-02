include("src/gffrecordcollection.jl")

TESTDIR = "/home/clara/Desktop/SHK-Job_Bleidorn/Projects/BioFmtSpecimens/GFF3/"
TESTFILES = [TESTDIR * "au9_scaffold_subset.gff3", TESTDIR * "directives.gff3", 
             TESTDIR * "TAIR10.part.gff.bgz", TESTDIR * "knownGene_out_of_order.gff3", 
             TESTDIR * "messy_protein_domains.gff3", TESTDIR * "Homo_sapiens.GRCh38.85.n1000.gff3"]

TESTRECORD = GFF3Record("CCDS1.1\tCCDS\tCDS\t802434\t801943\t.\t-\t.\tID=LINC00115.cds;Name=LINC00115.cds")
TESTRECORD_2 = GFF3Record("CCDS1.2\tCCDS\tgene\t801943\t802434\t.\t-\t.\tID=LINC00115;Name=LINC00115;Gap=M12 G4 M24.")

try 
    checkstartend(TESTRECORD)
catch e
    if isa(e, AssertionError)
        println("Start and End OK")
    else
        println(e)
    end
end

try 
    checkcdsphase(TESTRECORD)
catch e
    if isa(e, AssertionError)
        println("Type and Phase OK")
    else
        println(e)
    end
end
@assert !gap_hastarget(TESTRECORD_2)

TESTCOLLECTION = GFF3RecordCollection(TESTFILES[1]; index=nothing)
TESTCOLLECTION_empty = GFF3RecordCollection(TESTFILES[2]; index=nothing) # only directives
@assert isempty(TESTCOLLECTION_empty)
TESTCOLLECTION_bgz = GFF3RecordCollection(TESTFILES[3])
TESTCOLLECTION_4 = GFF3RecordCollection(TESTFILES[4]; index=nothing) # error: BLAST Target not found in file 
                                                                     # not sure if that should be an error though
TESTCOLLECTION_5 = GFF3RecordCollection(TESTFILES[5]; index=nothing)
TESTCOLLECTION_6 = GFF3RecordCollection(TESTFILES[6]; index=nothing)

DataFrames.DataFrame(TESTCOLLECTION)
DataFrames.DataFrame(TESTCOLLECTION_bgz)
DataFrames.DataFrame(TESTCOLLECTION_5)
DataFrames.DataFrame(TESTCOLLECTION_6)

