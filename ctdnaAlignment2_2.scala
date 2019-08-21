#!/usr/bin/env anduril
//$OPT --threads 10
//$OPT -d /mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

object ctdna{
  val reference = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.fasta")
  val list2_2 = INPUT(path = "/mnt/storage2/rawdata/ctDNA/CHICseq2_2/list.csv")
  val listPilot = INPUT(path = "/mnt/storage1/rawdata/ctDNA/list.csv")
  val list3 = INPUT(path = "/mnt/storage2/rawdata/ctDNA/CHICseq3/list.csv")
  val dbsnp = INPUT(path="/mnt/storage1/rawdata/resources/hg19/dbsnp_138.hg19.vcf")
  val indelsGold = INPUT(path="/mnt/storage1/rawdata/resources/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")
  val targets = INPUT(path = "/mnt/storage1/rawdata/ctDNA/metadata/targets.bed")
  val referenceDict = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.dict")
  
  val bwa = "/mnt/storage1/tools/bwa/bwa-0.7.17/bwa"
  val picard = "/mnt/storage1/tools/picard/picard-2.18.26.jar"
  val picardDir = "/mnt/storage1/tools/picard2"
  val fgbio = "/mnt/storage1/tools/fgbio/fgbio-0.7.0.jar"
  val gatk = "/mnt/storage1/tools/gatk/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"
  val java8 = "/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java"
  
  val list = CSVDplyr(csv1 = listPilot,
       csv2 = list2_2,
       csv3 = list3,
       function1 = """mutate(readGroup = paste0(sample, "_pilot"))""",
       function2 = """select(sample, readGroup, reads, UMI, mates)""",
       function3 = """bind_rows(csv2)""",
       function4 = """bind_rows(csv3)""")

  val targetsPadded = CSVDplyr(csv1 = targets,
        function1 = """mutate(Start = pmax(0,Start-200), End = End + 200)""")
  
  
  val targetsSorted = CSVDplyr(csv1 = targets,
        function1 = """arrange(factor(Chromosome, levels = paste0("chr",c(1:22,"X"))), Start)""")

  val targetsBed = BashEvaluate(var1 = targetsSorted,
       script = "cat @var1@ | awk 'NR > 1' > @out1@")
  targetsBed._filename("out1","targets.bed")

  val targetsIL = BashEvaluate(var1 = targets,
       var2 = referenceDict,
       param1 = picard,
       script = """
             cat @var1@ | awk 'NR>1' > @out2@
             java -jar @param1@ BedToIntervalList I=@out2@ O=@out1@ SD=@var2@
                """)
  targetsIL._filename("out1","targets.interval_list")
  targetsIL._filename("out2","out2.bed")

  val reads   = NamedMap[INPUT]("reads")
  val mates   = NamedMap[INPUT]("mates")
  val umi     = NamedMap[INPUT]("umi")
  val aligned = NamedMap[BashEvaluate]("aligned")
  val alignedUmi = NamedMap[BashEvaluate]("alignedUmi")
  val sorted  = NamedMap[BashEvaluate]("sorted")
  val sortedOut = NamedMap[Any]("sortedOut")
  val bamCSVIn = NamedMap[INPUT]("bamCSVIn")
  val bamArrayIn = NamedMap[CSV2Array]("bamArrayIn")
  val merged = NamedMap[BamCombiner]("merged")
  val marked  = NamedMap[BashEvaluate]("marked")
  val groupedUmi = NamedMap[BashEvaluate]("groupedUmi")
  val libraryComplexity = NamedMap[BashEvaluate]("libraryComplexity")
  val sortedTemplates = NamedMap[BashEvaluate]("sortedTemplates")
  val consensused = NamedMap[BashEvaluate]("consensused")
  val consensusedFilter = NamedMap[BashEvaluate]("consensusedFilter")
  val consensusedFastq = NamedMap[BashEvaluate]("consensusedFastq")
  val consensusAligned = NamedMap[BashEvaluate]("consensusAligned")
  val consensusSorted = NamedMap[BashEvaluate]("consensusSorted")
  val bamOut = NamedMap[Any]("bamOut")
  val bamCorrectedOut = NamedMap[Any]("bamCorrectedOut")
  val complexOut = NamedMap[Any]("complexOut")
  val consensusRecal = NamedMap[BashEvaluate]("consensusRecal")
  val groupedUmiOut = NamedMap[Any]("groupedUmiOut")
  val alignmentSummary = NamedMap[BashEvaluate]("alignmentSummary")
  val alignmentSummaryCombined = NamedMap[BashEvaluate]("alignmentSummaryCombined")
  val alignmentSummaryOut = NamedMap[Any]("alignmentSummaryOut")


  for ( rowMap <- iterCSV(list) ) {
	  
    val Sample  = rowMap("sample")
    val readGroup = rowMap("readGroup")

    reads(readGroup) = INPUT(path = rowMap("reads"))
    mates(readGroup) = INPUT(path = rowMap("mates"))
    umi(readGroup)   = INPUT(path = rowMap("UMI"))
    

    // alignment using bwa-mem
    aligned(readGroup) = BashEvaluate(var1 = reads(readGroup),
        var2 = mates(readGroup),
        var3 = reference,
        param1 = readGroup,
        param2 = Sample,
        script = s"""
                 $bwa mem -M -t 2 -R "@RG\\tID:@param1@\\tPL:ILLUMINA\\tLB:Library1\\tSM:@param2@" @var3@ @var1@ @var2@ > @out1@
                 """)
    aligned(readGroup)._keep = false
    aligned(readGroup)._filename("out1","aligned.bam")
   // aligned(readGroup)._execute = "once"

    // attach Umis to aligned bam file from a fastq
    alignedUmi(readGroup) = BashEvaluate(var1 = aligned(readGroup).out1,
        var2 = umi(readGroup),
        script = s"""
                  java -jar $fgbio AnnotateBamWithUmis -i @var1@ -f @var2@ -o @out1@ 
                  """)
    alignedUmi(readGroup)._keep = false
    alignedUmi(readGroup)._filename("out1","alignedUmi.bam")
  //  alignedUmi(readGroup)._execute = "once"

    // sort reads by coordinates
    sorted(readGroup) = BashEvaluate(var1 = alignedUmi(readGroup).out1,
        var2   = reference,
        param1 = picard,
        script = """java -Xmx4g -jar @param1@ SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate \
               CREATE_INDEX=true CREATE_MD5_FILE=true INPUT=@var1@ OUTPUT=@out1@ TMP_DIR=$( gettempdir )
               """)
    sorted(readGroup)._filename("out1","sorted.bam")
   // sorted(readGroup)._execute = "once"
  
    sortedOut(readGroup) = sorted(readGroup).out1
   }

  val bamCSV = Array2CSV(in = sortedOut)
  val bamCSVperSample = CSVDplyr(csv1 = bamCSV,
          csv2 = list,
          function1 = """merge(csv2[,1:2], by.x = 1, by.y = 2, sort  = F)""")
        

  val bySampleArray = Rpipe(csv1=bamCSVperSample,
           function1 = """arrayOut <- csv1 %>% select(Key, File)""",
           function2 = """split(csv1$sample)""")
        
  
  val bySampleCSV = Array2CSV(in = bySampleArray.arrayOut)

  for ( rowMap <- iterCSV(bySampleCSV) ) {
    val Sample = rowMap("Key")
  
    bamCSVIn(Sample) = INPUT(path = rowMap("File"))
    bamArrayIn(Sample) = CSV2Array(in = bamCSVIn(Sample))

    merged(Sample) = BamCombiner(in = bamArrayIn(Sample),
              picard = picardDir)
    merged(Sample)._keep = false


 
    // markd duplicates (UMI aware)
    marked(Sample) = BashEvaluate(var1 = merged(Sample).out,
    	param1 = picard,
        script = """
                  java -Xmx4G -jar @param1@ UmiAwareMarkDuplicatesWithMateCigar \
                  I=@var1@  O=@out1@ M=@out2@ UMI_METRICS=@out3@ \
                  MOLECULAR_IDENTIFIER_TAG=MI BARCODE_TAG=BC CREATE_INDEX=true CREATE_MD5_FILE=true \
                  ASSUME_SORT_ORDER=coordinate TMP_DIR=$( gettempdir )
                  """)
    marked(Sample)._filename("out1","out1.bam")
    marked(Sample)._filename("out2","duplicateMetrics.txt")
    marked(Sample)._filename("out3","umiMetrics.txt")
    
    // extract markDuplicates metrics
    libraryComplexity(Sample) = BashEvaluate(var1 = marked(Sample).out2,
        param1 = picard,
        script = """
                 head -n 8 @var1@ | tail -n 2 > @out1@   
                      """)   
    
    complexOut(Sample) = libraryComplexity(Sample).out1

    // enhance SortBam by --tmp-dir=/my/big/scratch/volume after fgbio.jar
    groupedUmi(Sample) = BashEvaluate(var1 = merged(Sample).out,
       var2 = reference,
       param1 = fgbio,
       script = s"""
            java -jar @param1@ SortBam -s Queryname -i @var1@ -o @out3@ 
            java -jar @param1@ SetMateInformation -r @var2@ -i @out3@ -o @out2@
            java -jar @param1@ GroupReadsByUmi -i @out2@ -o @out1@ -f @out2@ -s adjacency -l 6   
                 """)
    groupedUmi(Sample)._filename("out1",Sample + "_groupedUmi.bam")
    groupedUmi(Sample)._filename("out2",Sample + "_withmateset.bam")
    groupedUmi(Sample)._filename("out3",Sample + "_sorted.bam")
    groupedUmi(Sample)._keep = false

    groupedUmiOut(Sample) = groupedUmi(Sample).out1
    // call consensus bases based on duplicates
    // M: min-reads 1
    consensused(Sample) = BashEvaluate(var1 = groupedUmi(Sample).out1,
       script = s"""
                java -jar $fgbio CallMolecularConsensusReads -i @var1@ -o @out1@ -r @out2@ -M 1
                """)
    consensused(Sample)._filename("out1","out1.bam")
    consensused(Sample)._filename("out2","rejected.bam")
    consensused(Sample)._keep = false

    // filtering reads and masking bases
    // min-mean-base-quality q: 10
    // max-no-call-fraction  n: 0.3
    // min-base-quality N: 10
    // max-base-error-rate e: 0.3
    // max-read-error-rate E: 0.05
    // min-reads supporting a consensus M: 1
    consensusedFilter(Sample) = BashEvaluate(var1 = consensused(Sample).out1,
        var2 = reference,
        param1 = fgbio,
        script = """
                 java -jar @param1@ FilterConsensusReads \
                   -i @var1@ -o @out1@ -r @var2@ -M 1 -N 15 -q 20
                 """)
    consensusedFilter(Sample)._filename("out1", Sample + ".bam")
    consensusedFilter(Sample)._keep = false

    consensusedFastq(Sample) = BashEvaluate(var1 = consensusedFilter(Sample).out1,
        script = s"""
                 java -jar $picard SamToFastq I=@var1@ F=@out1@ F2=@out2@
                 """)
    consensusedFastq(Sample)._filename("out1",Sample + "_R1.fastq")
    consensusedFastq(Sample)._filename("out2",Sample + "_R2.fastq")
    consensusedFastq(Sample)._keep = false

    consensusAligned(Sample) = BashEvaluate(var1 = consensusedFastq(Sample).out1,
        var2 = consensusedFastq(Sample).out2,
        var3 = reference,
        param1 = Sample,
        script = s"""
                 $bwa mem -M -t 4 -R "@RG\\tID:ID1\\tPL:ILLUMINA\\tLB:Library1\\tSM:@param1@" @var3@ @var1@ @var2@ > @out1@
                 """)
    consensusAligned(Sample)._filename("out1",Sample + "_ConsensusAligned.bam")
    consensusAligned(Sample)._keep = false

    consensusSorted(Sample) = BashEvaluate(var1 = consensusAligned(Sample).out1,
        var2   = reference,
        param1 = picard,
        script = """java -Xmx4g -jar @param1@ SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate \
               CREATE_INDEX=true CREATE_MD5_FILE=true INPUT=@var1@ OUTPUT=@out1@ TMP_DIR=$( gettempdir )
               """)
    consensusSorted(Sample)._filename("out1",Sample + "_CorrectedSorted.bam")
    consensusSorted(Sample)._keep = false

    consensusRecal(Sample) = BashEvaluate(var1 = consensusSorted(Sample).out1,
         var2 = reference,
         var3 = dbsnp,
         var4 = indelsGold,
         param1 = gatk,
         param2 = java8,
         script = """
                  @param2@ -jar @param1@ BaseRecalibrator -I @var1@ -R @var2@ \
                  --known-sites @var3@ --known-sites @var4@ -O @out1@
                  @param2@ -jar @param1@ ApplyBQSR -R @var2@ -I @var1@ \
                  --bqsr-recal-file @out1@ -O @out2@
                  """)
    consensusRecal(Sample)._filename("out1","recalibration.table")
    consensusRecal(Sample)._filename("out2", Sample + "_consensusRecal.bam")

  

    bamOut(Sample) = marked(Sample).out1
    bamCorrectedOut(Sample) = consensusRecal(Sample).out2

     alignmentSummary(Sample) = BashEvaluate(var1 = consensusRecal(Sample).out2,
       var2 = reference,
       var3 = targetsIL.out1,
       param1 = picard,
       script = """
            java -jar @param1@ CollectAlignmentSummaryMetrics R=@var2@ I=@var1@ O=@out1@ 
            java -jar @param1@ CollectHsMetrics I=@var1@ O=@out2@ R=@var2@ TARGET_INTERVALS=@var3@ BAIT_INTERVALS=@var3@
                """)
  
    alignmentSummaryCombined(Sample) = BashEvaluate(var1 = alignmentSummary(Sample).out1,
      var2 = alignmentSummary(Sample).out2,
      script = """
            cat @var1@ | grep 'CATEGORY\|^PAIR' > @out1@
            cat @var2@ | grep 'BAIT_SET\|^targets' > @out2@
            paste @out1@ @out2@ > @out3@
               """)


   alignmentSummaryOut(Sample) = alignmentSummaryCombined(Sample).out3

    }
 
  val stats = CSVListJoin(in = alignmentSummaryOut,
       fileCol = "Sample")
  

  val dupStats = CSVListJoin(in = complexOut,
       fileCol = "Sample")
  

  val bamCorrectedOutCSV = Array2CSV(in = bamCorrectedOut)

}
 