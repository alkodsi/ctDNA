#!/usr/bin/env anduril
//$OPT --threads 10
//$OPT -d /mnt/storage2/work/amjad/ctdna/result_SVpipeline

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

object ctdna{
  
    val reference = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.fasta")
    val list = INPUT(path = "/mnt/storage2/work/amjad/ctdna/Bams/list.csv")
  
    val svaba = "/mnt/storage1/tools/svaba/svaba/bin/svaba"
    val snpeff = "/mnt/storage1/tools/snpEff/snpEff/snpEff.jar"
    val snpSift = "/mnt/storage1/tools/snpEff/snpEff/SnpSift.jar"

    val listTumors = CSVDplyr(csv1 = list,
       function1 = """mutate(Patient = sapply(strsplit(Key,"_"),function(x)paste(x[2],x[3],sep="_")))""",
       function2 = """filter(!grepl("WB",Key))""",
       function3 = """select(KeyTumor = Key, Tumor = File, Patient)""")
  
    val listNormals = CSVDplyr(csv1 = list,
       function1 = """mutate(Patient = sapply(strsplit(Key,"_"),function(x)paste(x[2],x[3],sep="_")))""",
       function2 = """filter(grepl("WB",Key))""",
       function3 = """select(KeyNormal = Key,Normal = File, Patient)""",
       function4 = """bind_rows(data.frame(KeyNormal = "ctDNA_ctrl_1", 
                        Normal = "/mnt/storage2/work/amjad/ctdna/Bams/ctDNA_ctrl_1_consensusRecal.bam",
                        Patient = "CHIC_138", stringsAsFactors=F))""")
  
    val listMatched = CSVDplyr(csv1=listTumors,
        csv2 = listNormals,
        script = """library(plyr)""",
        function1 = """plyr::join(csv2,by="Patient")""",
        function2 = """filter(!is.na(KeyNormal))""",
        function3 = """mutate(Nocontrol = ifelse(Patient %in% c("CHIC_143","CHIC_138", "CHIC_62"), "yes", "no"))""")

    val control      = NamedMap[INPUT]("control")
    val tumor        = NamedMap[INPUT]("tumor")
    val SVs          = NamedMap[BashEvaluate]("SVs")
    val SVsAnn       = NamedMap[BashEvaluate]("SVsAnn")
    val SVsAnnCSV    = NamedMap[BashEvaluate]("SVsAnnCSV")
    val SVsAnnCSVFixed = NamedMap[CSVDplyr]("SVsAnnCSVFixed")

    for ( rowMap <- iterCSV(listMatched) ) { 

        val key = rowMap("KeyTumor")
        val nocontrol = rowMap("Nocontrol")

        tumor(key) = INPUT(path = rowMap("Tumor"))
        control(key) = INPUT(path=rowMap("Normal"))

        if(nocontrol == "no"){
        
            SVs(key)   = BashEvaluate(var1 = tumor(key),
                var2   = control(key),
                var3   = reference,
                script = s"$svaba run -t @var1@ -n @var2@ -a @folder1@/$key -G @var3@ -C 100000")
        
            SVsAnn(key) = BashEvaluate(var1 = SVs(key).folder1,
                script  = s"java -jar $snpeff hg19 @var1@/$key.svaba.somatic.sv.vcf > @out1@")
            SVsAnn(key)._filename("out1", key + "_ann.vcf")

            SVsAnnCSV(key) = BashEvaluate(var1 = SVsAnn(key).out1,
                script = s"java -jar $snpSift extractFields @var1@ CHROM POS ID REF ALT FILTER " +
                     "ANN[0].EFFECT ANN[0].IMPACT ANN[0].GENE ANN[0].GENEID ANN[0].FEATURE ANN[0].FEATUREID " +
                     "ANN[0].BIOTYPE ANN[0].RANK ANN[0].HGVS_C ANN[0].HGVS_P " +
                     "SPAN INSERTION SCTG EVDNC BX SVTYPE NM SUBN DISC_MAPQ REPSEQ MATEID HOMSEQ " +
                     "READNAMES IMPRECISE MAPQ MATEMAPQ SECONDARY HOMLEN MATENM NUMPARTS " +
                     "GEN[0].SL GEN[0].LO GEN[0].DR GEN[0].DP GEN[0].AD GEN[0].SR " +
                     "GEN[1].SL GEN[1].LO GEN[1].DR GEN[1].DP GEN[1].AD GEN[1].SR > @out1@")

            SVsAnnCSVFixed(key) = CSVDplyr(csv1 = SVsAnnCSV(key).out1,
                function1 = """rename_all(~ gsub("ANN\\[0\\].", "", .x))""",
                function2 = """rename_all(~ gsub("GEN\\[0\\]", "Normal", .x))""",
                function3 = """rename_all(~ gsub("GEN\\[1\\]", "Tumor", .x))""",
                function4 = """rename_all(~ gsub(".LO", ".LogOdd", .x))""",
                function5 = """rename_all(~ gsub(".SL", ".AlignmentQuality", .x))""",
                function6 = """rename_all(~ gsub(".SR", ".SpanningReads", .x))""",
                function7 = """rename_all(~ gsub(".LR", ".SpanningReads", .x))""")

        } else {

        	SVs(key)   = BashEvaluate(var1 = tumor(key),
                var2   = control(key),
                var3   = reference,
                script = s"$svaba run -t @var1@ -a @folder1@/$key -G @var3@ -C 100000")
        
            SVsAnn(key) = BashEvaluate(var1 = SVs(key).folder1,
                script  = s"java -jar $snpeff hg19 @var1@/$key.svaba.sv.vcf > @out1@")
            SVsAnn(key)._filename("out1", key + "_ann.vcf")

            SVsAnnCSV(key) = BashEvaluate(var1 = SVsAnn(key).out1,
                script = s"java -jar $snpSift extractFields @var1@ CHROM POS ID REF ALT FILTER " +
                     "ANN[0].EFFECT ANN[0].IMPACT ANN[0].GENE ANN[0].GENEID ANN[0].FEATURE ANN[0].FEATUREID " +
                     "ANN[0].BIOTYPE ANN[0].RANK ANN[0].HGVS_C ANN[0].HGVS_P " +
                     "SPAN INSERTION SCTG EVDNC BX SVTYPE NM SUBN DISC_MAPQ REPSEQ MATEID HOMSEQ " +
                     "READNAMES IMPRECISE MAPQ MATEMAPQ SECONDARY HOMLEN MATENM NUMPARTS " +
                     "GEN[0].SL GEN[0].LO GEN[0].DR GEN[0].DP GEN[0].AD GEN[0].SR > @out1@")
        
            SVsAnnCSVFixed(key) = CSVDplyr(csv1 = SVsAnnCSV(key).out1,
                function1 = """rename_all(~ gsub("ANN\\[0\\].", "", .x))""",
                function2 = """rename_all(~ gsub("GEN\\[0\\]", "Tumor", .x))""",
                function4 = """rename_all(~ gsub(".LO", ".LogOdd", .x))""",
                function5 = """rename_all(~ gsub(".SL", ".AlignmentQuality", .x))""",
                function6 = """rename_all(~ gsub(".SR", ".SpanningReads", .x))""",
                function7 = """rename_all(~ gsub(".LR", ".SpanningReads", .x))""")
        }
    }

  val SVsAll = CSVListJoin(in = SVsAnnCSVFixed,
       fileCol = "Sample")
  
}