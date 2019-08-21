#!/usr/bin/env anduril
//$OPT --threads 20
//$OPT -d /mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentNovaseq/

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

object ctdna{
  val reference = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.fasta")
  val list = INPUT(path = "/mnt/storage2/rawdata/ctDNA/list.csv")
  val dbsnp = INPUT(path="/mnt/storage1/rawdata/resources/hg19/dbsnp_138.hg19.vcf")
  val indelsGold = INPUT(path="/mnt/storage1/rawdata/resources/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf")
  val targets = INPUT(path = "/mnt/storage1/rawdata/ctDNA/metadata/targets.bed")
  val concentration = INPUT(path = "/mnt/storage1/rawdata/ctDNA/metadata/concentration.csv")
  val referenceDict = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.dict")
  
  val bwa = "/mnt/storage1/tools/bwa/bwa-0.7.17/bwa"
  val picard = "/mnt/storage1/tools/picard/picard-2.18.26.jar"
  val picardDir = "/mnt/storage1/tools/picard2"
  val fgbio = "/mnt/storage1/tools/fgbio/fgbio-0.7.0.jar"
  val gatk = "/mnt/storage1/tools/gatk/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"
  val java8 = "/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java"

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
  val trimmed = NamedMap[TrimGalore]("trimmed")
  val tileFiltered = NamedMap[BashEvaluate]("tileFiltered")
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
  val consensusedFilterStrict = NamedMap[BashEvaluate]("consensusedFilterStrict")
  val consensusedFastqStrict = NamedMap[BashEvaluate]("consensusedFastqStrict")
  val consensusAlignedStrict = NamedMap[BashEvaluate]("consensusAlignedStrict")
  val consensusSortedStrict = NamedMap[BashEvaluate]("consensusSortedStrict")
  val consensusRecalStrict = NamedMap[BashEvaluate]("consensusRecalStrict")
  val consensusRecal = NamedMap[BashEvaluate]("consensusRecal")
  val bamCorrectedStrictOut = NamedMap[Any]("bamCorrectedStrictOut")
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
    
    trimmed(readGroup) = TrimGalore(reads = reads(readGroup),
        mates = mates(readGroup),
        clipR1 = 5,
        clipR2 = 5,
        gzip = true,
        minQuality = 20)
    
    tileFiltered(readGroup) = BashEvaluate(var1 = trimmed(readGroup).trimmed("Reads"),
         var2 = trimmed(readGroup).trimmed("Mates"),
         script = "/mnt/storage1/tools/BBMap/bbmap/filterbytile.sh in1=@var1@ in2=@var2@ out1=@out1@ out2=@out2@ ud=0.7 qd=0.7 ed=0.7 ua=.4 qa=.4 ea=.4 overwrite=true ")
    tileFiltered(readGroup)._filename("out1","reads.fq.gz")
    tileFiltered(readGroup)._filename("out2","mates.fq.gz")


    // alignment using bwa-mem
    aligned(readGroup) = BashEvaluate(var1 = tileFiltered(readGroup).out1,
        var2 = tileFiltered(readGroup).out2,
        var3 = reference,
        param1 = readGroup,
        param2 = Sample,
        script = s"""
                 $bwa mem -M -t 4 -R "@RG\\tID:@param1@\\tPL:ILLUMINA\\tLB:Library1\\tSM:@param2@" @var3@ @var1@ @var2@ > @out1@
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
        

  val bySampleArray = REvaluate(table1=bamCSVperSample,
           script="""
                array.out <- split(table1[,c("Key","File")],table1$sample)
                table.out <- data.frame()
                  """)
  
  val bySampleCSV = Array2CSV(in = bySampleArray.outArray)

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

    
    groupedUmi(Sample) = BashEvaluate(var1 = merged(Sample).out,
       var2 = reference,
       param1 = fgbio,
       script = s"""
            java -jar @param1@ SortBam -s Queryname -i @var1@ -o @out3@ 
            java -jar @param1@ SetMateInformation -r @var2@ -i @out3@ -o @out2@
            java -jar @param1@ GroupReadsByUmi -i @out2@ -o @out1@ -f @out2@ -s adjacency -l 5   
                 """)
    groupedUmi(Sample)._filename("out1",Sample + "_groupedUmi.bam")
    groupedUmi(Sample)._filename("out2",Sample + "_withmateset.bam")
    groupedUmi(Sample)._filename("out3",Sample + "_sorted.bam")

    groupedUmiOut(Sample) = groupedUmi(Sample).out1
    // call consensus bases based on duplicates
    // M: min-reads 1
    consensused(Sample) = BashEvaluate(var1 = groupedUmi(Sample).out1,
       script = s"""
                java -jar $fgbio CallMolecularConsensusReads -i @var1@ -o @out1@ -r @out2@ -M 1
                """)
    consensused(Sample)._filename("out1","out1.bam")
    consensused(Sample)._filename("out2","rejected.bam")
 

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
                   -i @var1@ -o @out1@ -r @var2@ -M 1 -N 10 -q 10 -n 0.3 -e 0.3 -E 0.05
                 """)
    consensusedFilter(Sample)._filename("out1", Sample + ".bam")
    

    consensusedFastq(Sample) = BashEvaluate(var1 = consensusedFilter(Sample).out1,
        script = s"""
                 java -jar $picard SamToFastq I=@var1@ F=@out1@ F2=@out2@
                 """)
    consensusedFastq(Sample)._filename("out1",Sample + "_R1.fastq")
    consensusedFastq(Sample)._filename("out2",Sample + "_R2.fastq")

    consensusAligned(Sample) = BashEvaluate(var1 = consensusedFastq(Sample).out1,
        var2 = consensusedFastq(Sample).out2,
        var3 = reference,
        param1 = Sample,
        script = s"""
                 $bwa mem -M -t 4 -R "@RG\\tID:ID1\\tPL:ILLUMINA\\tLB:Library1\\tSM:@param1@" @var3@ @var1@ @var2@ > @out1@
                 """)
    consensusAligned(Sample)._filename("out1",Sample + "_ConsensusAligned.bam")


    consensusSorted(Sample) = BashEvaluate(var1 = consensusAligned(Sample).out1,
        var2   = reference,
        param1 = picard,
        script = """java -Xmx4g -jar @param1@ SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate \
               CREATE_INDEX=true CREATE_MD5_FILE=true INPUT=@var1@ OUTPUT=@out1@ TMP_DIR=$( gettempdir )
               """)
    consensusSorted(Sample)._filename("out1",Sample + "_CorrectedSorted.bam")


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

    consensusedFilterStrict(Sample) = BashEvaluate(var1 = consensused(Sample).out1,
        var2 = reference,
        param1 = fgbio,
        script = """
                 java -jar @param1@ FilterConsensusReads \
                   -i @var1@ -o @out1@ -r @var2@ -M 2 -N 10 -q 20
                 """)
    consensusedFilterStrict(Sample)._filename("out1", Sample + ".bam")
    
    consensusedFastqStrict(Sample) = BashEvaluate(var1 = consensusedFilterStrict(Sample).out1,
        script = s"""
                 java -jar $picard SamToFastq I=@var1@ F=@out1@ F2=@out2@
                 """)
    consensusedFastqStrict(Sample)._filename("out1",Sample + "_R1.fastq")
    consensusedFastqStrict(Sample)._filename("out2",Sample + "_R2.fastq")

    consensusAlignedStrict(Sample) = BashEvaluate(var1 = consensusedFastqStrict(Sample).out1,
        var2 = consensusedFastqStrict(Sample).out2,
        var3 = reference,
        param1 = Sample,
        script = s"""
                 $bwa mem -M -t 4 -R "@RG\\tID:ID1\\tPL:ILLUMINA\\tLB:Library1\\tSM:@param1@" @var3@ @var1@ @var2@ > @out1@
                 """)
    consensusAlignedStrict(Sample)._filename("out1",Sample + "_ConsensusAligned.bam")


    consensusSortedStrict(Sample) = BashEvaluate(var1 = consensusAlignedStrict(Sample).out1,
        var2   = reference,
        param1 = picard,
        script = """java -Xmx4g -jar @param1@ SortSam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate \
               CREATE_INDEX=true CREATE_MD5_FILE=true INPUT=@var1@ OUTPUT=@out1@ TMP_DIR=$( gettempdir )
               """)
    consensusSortedStrict(Sample)._filename("out1",Sample + "_CorrectedSorted.bam")

    consensusRecalStrict(Sample) = BashEvaluate(var1 = consensusSortedStrict(Sample).out1,
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
    consensusRecalStrict(Sample)._filename("out1","recalibration.table")
    consensusRecalStrict(Sample)._filename("out2", Sample + "_consensusRecal.bam")

    bamOut(Sample) = marked(Sample).out1
    bamCorrectedOut(Sample) = consensusRecal(Sample).out2
    bamCorrectedStrictOut(Sample) = consensusRecalStrict(Sample).out2

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
  /*
  Is it worth it to sequence more given the complexities of the current libraries?
  Estimate the number of unique reads when doubling/tripling the sequencing
  This uses the statistics from Picard's EstimateLibraryComplexity and
  utilizes the following formula:
  C/x = 1 - exp(-N/x)
  where:
  C is the number of distinct fragments
  N is the number of read pairs
  x is the number of unique molecules
  */
  /*
  val libComplexities = REvaluate(inArray = complexOut,
        table1 = concentration,
        script = """
           library(dplyr)
           library(ggplot2)
           library(viridis)
           library(reshape)
             EstimateNewUniqueReads <- function(librarySize, uniqueReads, totalReads, newTotalReads){
                leftSide <- uniqueReads/librarySize
                rightSide <- 1 - exp((-1*totalReads)/librarySize)
                if(round(leftSide,3) != round(rightSide,3))
                  stop("The numbers don't match")
   
               newUniqueReads <- round((1 - exp((-1*newTotalReads)/librarySize)) * librarySize)
               return(newUniqueReads)
		   }

           table.out <- data.frame(sample = names(array), bind_rows(array)) %>%
                   mutate(uniqueReadsDoubleSeq = EstimateNewUniqueReads(ESTIMATED_LIBRARY_SIZE, READ_PAIRS_EXAMINED-READ_PAIR_DUPLICATES, READ_PAIRS_EXAMINED, READ_PAIRS_EXAMINED *2),
                          doubleSeqChangePercent = 100 * round((uniqueReadsDoubleSeq - (READ_PAIRS_EXAMINED-READ_PAIR_DUPLICATES))/(READ_PAIRS_EXAMINED-READ_PAIR_DUPLICATES),3),
                          uniqueReadsTripleSeq = EstimateNewUniqueReads(ESTIMATED_LIBRARY_SIZE, READ_PAIRS_EXAMINED-READ_PAIR_DUPLICATES, READ_PAIRS_EXAMINED, READ_PAIRS_EXAMINED *3),
                          tripleSeqChangePercent = 100 * round((uniqueReadsTripleSeq - (READ_PAIRS_EXAMINED-READ_PAIR_DUPLICATES))/(READ_PAIRS_EXAMINED-READ_PAIR_DUPLICATES),3))                 
           x <- table.out %>% mutate(uniqueReads = READ_PAIRS_EXAMINED -READ_PAIR_DUPLICATES) %>% 
                   mutate(uniqueReadsTripleSeq = uniqueReadsTripleSeq - uniqueReadsDoubleSeq, uniqueReadsDoubleSeq = uniqueReadsDoubleSeq - uniqueReads) %>% 
                   select(sample, uniqueReads, uniqueReadsDoubleSeq,uniqueReadsTripleSeq)
           m <- melt(x) %>% mutate(variable = factor(variable,levels = rev(levels(variable))),
                                   sample = factor(sample, levels = rev(levels(sample))))
           g <- m %>% ggplot(aes(x = sample, y = value, fill = variable)) + 
                geom_bar(stat = "identity") + coord_flip() + 
                scale_fill_viridis(discrete=T, direction = -1, alpha = 0.8, end = 0.85,begin = 0.15, option = "B") + 
                theme_minimal() + labs(y = "Number of Unique reads")
          setwd(document.dir)
          ggsave(g, file = "barplot.pdf")
          
          table1[,1] <- gsub("-","_",table1[,1])
          xx <- merge(x, table1, by = 1)
          gg <- xx %>% ggplot(aes(x = uniqueReads, y = ngul)) + 
               geom_point(size = 3, color = magma(100)[20]) +
               geom_smooth(method = "lm", se = F, color = magma(100)[50],size=1.5 ) + #xlim(c(1000000,10000000)) +
               theme_minimal() + theme(axis.line = element_line(size=0.7)) + labs(x = "Number of Unique Reads", y = "Concentration (ng/ul)")
          gg2 <- xx %>% ggplot(aes(x = uniqueReads, y = `max input ng`)) + 
               geom_point(size = 3, color = magma(100)[20]) +
               geom_smooth(method = "lm", se = F, color = magma(100)[50],size=1.5 ) + #xlim(c(1000000,10000000)) +
               theme_minimal() + theme(axis.line = element_line(size=0.7)) + labs(x = "Number of Unique Reads", y = "starting amount")

          ggsave(gg, file = "concVsReads.pdf")
          ggsave(gg2, file = "startAmountVsReads.pdf")
                 """)


  val bamOutCSV = Array2CSV(in = bamOut)
  val bamCorrectedOutCSV = Array2CSV(in = bamCorrectedOut)
  val bamCorrectedStrictOutCSV = Array2CSV(in = bamCorrectedStrictOut)

  val bamOutMixed = CSVDplyr(csv1 = bamCorrectedOutCSV,
       csv2 = bamCorrectedStrictOutCSV,
       function1 = """filter(!grepl("FFPE",Key))""",
       function2 = """bind_rows(filter(csv2, grepl("FFPE",Key)))""")
}
*/

