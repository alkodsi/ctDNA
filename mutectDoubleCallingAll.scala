#!/usr/bin/env anduril
//$OPT --threads 15
//$OPT -d /mnt/storage2/work/amjad/ctdna/result_mutectDoubleCallingFinal

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

object ctdna{

  val list = INPUT(path = "/mnt/storage2/work/amjad/ctdna/result_ctDNAlign/bamCorrectedOutCSV/out.csv")
  val gnomad = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/af-only-gnomad.raw.sites.b37.vcf.gz")
  val chain = INPUT(path = "/mnt/storage1/rawdata/resources/chains/b37tohg19.chain")
  val reference = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.fasta")
  val referenceDict = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.dict")
  val targets = INPUT(path = "/mnt/storage1/rawdata/ctDNA/metadata/targets.bed")
  val exac = INPUT(path = "/mnt/storage1/rawdata/resources/hg38/small_exac_common_3.hg38.vcf.gz")
  val hg38tohg19chain = INPUT(path = "/mnt/storage1/rawdata/resources/chains/hg38ToHg19.over.chain")
  val hg38reference = INPUT(path = "/mnt/storage1/rawdata/resources/hg38/Homo_sapiens_assembly38.fasta")
  val pon = INPUT(path = "/mnt/storage2/work/amjad/ctdna/result_newMutectAll4_5/pon/pon.vcf.gz")
  val gnomadhg19 = INPUT(path = "/mnt/storage2/work/amjad/ctdna/result_newMutectAll4_5/gnomadhg19/gnomAD_hg19.vcf.gz")

  val picard = "/mnt/storage1/tools/picard/picard-2.21.4.jar"
  val fgbio = "/mnt/storage1/tools/fgbio/fgbio-0.7.0.jar"
  val gatk = "/mnt/storage1/tools/gatk/gatk-4.1.4.0/gatk-package-4.1.4.0-local.jar"
  val java8 = "/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java"
  val annovar = "/mnt/storage1/tools/Annovar/annovar/"
  val annovardb = "/mnt/storage1/tools/Annovar/annovar/humandb/"
  val gatk3 = "/mnt/storage1/tools/gatk/gatk-3.8/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
  

  val targetsIL = BashEvaluate(var1 = targets,
       var2 = referenceDict,
       param1 = picard,
       script = """
             cat @var1@ | awk 'NR>1' > @out2@
             java -jar @param1@ BedToIntervalList I=@out2@ O=@out1@ SD=@var2@
                """)
  targetsIL._filename("out1","targets.interval_list")
  targetsIL._filename("out2","out2.bed")

  val exachg19 = BashEvaluate(var1 = exac,
       var2 = hg38tohg19chain,
       var3 = reference,
       param1 = picard,
       script = """
          java -jar @param1@ LiftoverVcf \
           INPUT=@var1@ OUTPUT=@out1@ CHAIN=@var2@ REJECT=@out2@ REFERENCE_SEQUENCE=@var3@
                """)
  exachg19._filename("out1","Exac_hg19.vcf.gz")
  exachg19._filename("out2","rejected.vcf.gz")
  

  val listTumors = CSVDplyr(csv1 = list,
       function1 = """mutate(Patient = sapply(strsplit(Key,"_"),function(x)paste(x[2],x[3],sep="_")))""",
       function2 = """filter(!grepl("WB",Key))""",
       function3 = """select(KeyTumor = Key, Tumor = File, Patient)""")
  
  val listNormals = CSVDplyr(csv1 = list,
       function1 = """mutate(Patient = sapply(strsplit(Key,"_"),function(x)paste(x[2],x[3],sep="_")))""",
       function2 = """filter(grepl("WB",Key))""",
       function3 = """select(KeyNormal = Key,Normal = File, Patient)""",
       function4 = """bind_rows(data.frame(KeyNormal = "ctDNA_ctrl_1", 
                        Normal = "/mnt/storage2/work/amjad/ctdna/result_ctDNAlign/consensusRecal_ctDNA_ctrl_1/ctDNA_ctrl_1_consensusRecal.bam",
                        Patient = "CHIC_138", stringsAsFactors=F))""")
  
  val listMatched = CSVDplyr(csv1=listTumors,
        csv2 = listNormals,
        script = """library(plyr)""",
        function1 = """plyr::join(csv2,by="Patient")""",
        function2 = """filter(!is.na(KeyNormal))""")
 
 val listMatchedSelSamples = CSVDplyr(csv1 = listMatched,
       function1 = """filter(grepl("_0$", KeyTumor) | grepl("FFPE", KeyTumor))""")

   
  val control   = NamedMap[INPUT]("control")
  val tumor     = NamedMap[INPUT]("tumor")
  val contam    = NamedMap[BashEvaluate]("contam")
  val outContam = NamedMap[Any]("outContam")
  val outSegs   = NamedMap[Any]("outSegs")

  for ( rowMap <- iterCSV(listMatchedSelSamples) ) { 

    val keyT = rowMap("KeyTumor")
    val keyN = rowMap("KeyNormal")
    tumor(keyT) = INPUT(path = rowMap("Tumor"))
    control(keyT) = INPUT(path=rowMap("Normal"))

    contam(keyT) = BashEvaluate(var1 = tumor(keyT),
     var2 = exachg19.out1,
     var3 = control(keyT),
     var4 = targetsIL.out1,
     param1 = gatk,
     param2 = java8,
     script = """
         @param2@ -jar @param1@  GetPileupSummaries -I @var1@ -L @var4@ -ip 300 \
        -V @var2@ -L @var2@ -O @folder1@/tumor.pilups 
        @param2@ -jar @param1@  GetPileupSummaries -I @var3@ -L @var4@ -ip 300 \
        -V @var2@ -L @var2@ -O @folder1@/normal.pileups
         @param2@ -jar @param1@ CalculateContamination -I @folder1@/tumor.pilups \
         -matched @folder1@/normal.pileups -O @out1@ -tumor-segmentation @out2@
              """)
    contam(keyT)._filename("out1", keyT + "_contam.table")
    contam(keyT)._filename("out2", keyT + "_segments.table")

    outContam(keyT) = contam(keyT).out1
    outSegs(keyT) = contam(keyT).out2
}

 val outContamCSV = Array2CSV(in = outContam)
 val outSegsCSV   = Array2CSV(in = outSegs)


val outContamArray = REvaluate(table1 = outContamCSV,
     table2 = listMatchedSelSamples,
     script = """
         table1 <- table1[table1$Key %in% table2$KeyTumor,]
         table1 <- table1[order(match(table1$Key, table2$KeyTumor)),]
         table.out <- data.frame()
         array.out <- split(table1, table2$Patient)
              """)

val outContamArrayCSV = Array2CSV(in = outContamArray.outArray)

val outSegsArray = REvaluate(table1 = outSegsCSV,
     table2 = listMatchedSelSamples,
     script = """
         table1 <- table1[table1$Key %in% table2$KeyTumor,]
         table1 <- table1[order(match(table1$Key, table2$KeyTumor)),]
         table.out <- data.frame()
         array.out <- split(table1, table2$Patient)
              """)

val outSegsArrayCSV = Array2CSV(in = outSegsArray.outArray)



val outBamArray = REvaluate(table1 = listMatchedSelSamples,
     script = """
         library(dplyr)
         table.out <- data.frame()
         table2 <- select(table1, Key = KeyTumor, File = Tumor)
         array.out <- split(table2, table1$Patient)
              """)
val outBamArrayCSV = Array2CSV(in = outBamArray.outArray)

val outBamArrayCSVmerged = CSVDplyr(csv1 = outBamArrayCSV,
     csv2 = listMatchedSelSamples,
     function1 = """merge(csv2[!duplicated(csv2$Patient),c("Patient","Normal","KeyNormal")], by = 1, sort = F)""")



/* 
   Loop over patients and
   call variants from multiple tumor samples
   learn a model for F1R2 bias
   filter variants
   annotate variants
*/

val F1R2Model              = NamedMap[BashEvaluate]("F1R2Model")
val BamsByPatient          = NamedMap[CSV2Array]("BamsByPatient")
val contamByPatient        = NamedMap[CSV2Array]("contamByPatient")
val segsByPatient          = NamedMap[CSV2Array]("segsByPatient")
val normalControl          = NamedMap[INPUT]("normalControl")
val vars                   = NamedMap[BashEvaluate]("vars")
val varsFiltered           = NamedMap[BashEvaluate]("varsFiltered")
val varsPass               = NamedMap[BashEvaluate]("varsPass")
val varsAnnot              = NamedMap[Annovar]("varsAnnot")
val varsAnnotCSV           = NamedMap[BashEvaluate]("varsAnnotCSV")
val varsAnnotCSVFixed      = NamedMap[CSVDplyr]("varsAnnotCSVFixed")
val varsPassOut            = NamedMap[Any]("varsPassOut")
val varsFilteredAlignments = NamedMap[BashEvaluate]("varsFilteredAlignments")
val varsAnnotFixed         = NamedMap[BashEvaluate]("varsAnnotFixed")

for ( rowMap <- iterCSV(outBamArrayCSVmerged) ) { 
  val patient = rowMap("Key")
  val keyNormal = rowMap("KeyNormal")
  normalControl(patient) = INPUT(path = rowMap("Normal"))
  
  BamsByPatient(patient) = CSV2Array(in = outBamArray.outArray(patient),
      keys = "column")

  contamByPatient(patient) = CSV2Array(in = outContamArray.outArray(patient),
      keys = "column")

  segsByPatient(patient) = CSV2Array(in = outSegsArray.outArray(patient),
      keys = "column")

  if(patient != "CHIC_143" & patient != "CHIC_138" & patient != "CHIC_62"){
  vars(patient) = BashEvaluate(array1 = BamsByPatient(patient),
        var1 = normalControl(patient),
        var2 = reference,
        var3 = pon,
        var4 = gnomadhg19,
        var5 = targetsIL.out1,
        param1 = keyNormal,
        script = s"$java8 -jar $gatk Mutect2 -R @var2@ -I @var1@ -O @out1@ --max-reads-per-alignment-start 0 --pcr-indel-model HOSTILE --bam-output @out3@ " +
                  " --germline-resource @var4@ --panel-of-normals @var3@ -L @var5@ -ip 300 --f1r2-tar-gz @out2@ -normal @param1@ --max-mnp-distance 0 " +
                  """ $( paste -d ' ' <(getarrayfiles array1)  | sed 's,^, -I ,' | tr -d '\\\n' ) """)        
  vars(patient)._filename("out1", patient + "_rawVariants.vcf.gz")
  vars(patient)._filename("out2", patient + "_f1r2.tar.gz")
  vars(patient)._filename("out3", patient + "_mutect2.bam")
  // --force-active true --tumor-lod-to-emit 0 --initial-tumor-lod 0 -bamout @out3@
  // --pcr-indel-model AGGRESSIVE
 } else {
    vars(patient) = BashEvaluate(array1 = BamsByPatient(patient),
        var1 = normalControl(patient),
        var2 = reference,
        var3 = pon,
        var4 = gnomadhg19,
        var5 = targetsIL.out1,
        param1 = keyNormal,
        script = s"$java8 -jar $gatk Mutect2 -R @var2@ -O @out1@ --max-reads-per-alignment-start 0 --pcr-indel-model HOSTILE --bam-output @out3@ " +
                  " --germline-resource @var4@ --panel-of-normals @var3@ -L @var5@ -ip 300 --f1r2-tar-gz @out2@ --max-mnp-distance 0 " +
                  """ $( paste -d ' ' <(getarrayfiles array1)  | sed 's,^, -I ,' | tr -d '\\\n' ) """)        
  vars(patient)._filename("out1", patient + "_rawVariants.vcf.gz")
  vars(patient)._filename("out2", patient + "_f1r2.tar.gz")
  vars(patient)._filename("out3", patient + ".bam")
  }


  F1R2Model(patient) = BashEvaluate(var1 = vars(patient).out2,
      script = s"$java8 -jar $gatk LearnReadOrientationModel -I @var1@ -O @out1@")
  F1R2Model(patient)._filename("out1", patient + "_F1R2Model.tar.gz")


  varsFiltered(patient) = BashEvaluate(var1 = vars(patient).out1,
     var2 = F1R2Model(patient).out1,
     var3 = reference, 
     var4 = targetsIL.out1,
     array1 = contamByPatient(patient),
     array2 = segsByPatient(patient),
     script = s"$java8 -jar $gatk FilterMutectCalls -R @var3@ -V @var1@  -O @out1@ --stats @var1@.stats --filtering-stats @out2@ " +
               " --max-events-in-region 50 --min-median-read-position 15  -L @var4@ -ip 300 -ob-priors @var2@ " +
              """ $( paste -d ' ' <(getarrayfiles array1)  | sed 's,^, --contamination-table ,' | tr -d '\\\n' ) """ +
              """ $( paste -d ' ' <(getarrayfiles array2)  | sed 's,^, --tumor-segmentation ,' | tr -d '\\\n' ) """)
// -ob-priors @var2@
  varsFiltered(patient)._filename("out1", patient + "_filteredVariants.vcf.gz")
  varsFiltered(patient)._filename("out2", patient + "_filteringStats.csv")
/*
  varsFilteredAlignments(patient) = BashEvaluate(var1 = varsFiltered(patient).out1,
     var2 = reference,
     var3 = hg38Image.out1,
     array1 = BamsByPatient(patient),
     script = s"$java8 -jar $gatk FilterAlignmentArtifacts -R @var2@ -V @var1@ --bwa-mem-index-image @var3@ -O @out1@ " +
              """ $( paste -d ' ' <(getarrayfiles array1)  | sed 's,^, -I ,' | tr -d '\\\n' ) """)
  varsFilteredAlignments(patient)._filename("out1", patient + "_alignFiltered.vcf.gz")
*/
  varsPass(patient) = BashEvaluate(var1 = varsFiltered(patient).out1,
       var2 = targetsIL.out1,
       script = s"$java8 -jar $gatk SelectVariants --exclude-filtered -L @var2@ -ip 10 -V @var1@ -O @out1@")
  varsPass(patient)._filename("out1", patient + "_passed.vcf.gz")

  varsPassOut(patient) = varsPass(patient).out1



  varsAnnot(patient) = Annovar(vcfIn = varsPass(patient).out1,
            annovarPath = annovar,
            annovardb = annovardb,
            buildver = "hg19",
            inputType = "vcf",
            protocol = "refGene,avsnp147,cosmic68,dbnsfp30a,icgc21",
            operation = "g,f,f,f,f")

  varsAnnotFixed(patient) = BashEvaluate(var1 = varsAnnot(patient).vcfOut,
       script = """
              cat @var1@ | sed '/^[[:space:]]*$/d' > @out1@
                """)
  varsAnnotFixed(patient)._filename("out1", patient + "_annot.vcf")

  varsAnnotCSV(patient) = BashEvaluate(var1 = reference,
            var2 = varsAnnotFixed(patient).out1,
            param1 = gatk,
            param2 = java8,
            script = """ 
                @param2@ -jar @param1@ VariantsToTable \
                    -R @var1@ -V @var2@ -O @out1@ \
                    -F CHROM -F POS -F REF -F ALT -F ID \
                    -F Func.refGene -F Gene.refGene -F GeneDetail.refGene -F ExonicFunc.refGene -F AAChange.refGene \
                    -F avsnp147 -F cosmic68 -F SIFT_score -F SIFT_pred \
                    -F Polyphen2_HDIV_score -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_score -F Polyphen2_HVAR_pred \
                    -F MutationTaster_score -F MutationTaster_pred -F MutationAssessor_score -F MutationAssessor_pred -F PopFreqMax \
                    -F CADD_phred -F DANN_score \
                    -F CONTQ -F GERMQ -F ROQ -F MBQ -F ECNT -F MMQ -F MPOS -F NALOD -F POPAF -F SEQQ  \
                    -GF AD -GF DP -GF AF -GF F1R2 -GF F2R1 -GF PGT -GF PID 
            """)

  varsAnnotCSVFixed(patient) = CSVDplyr(csv1 = varsAnnotCSV(patient).out1,
     function1 = """mutate_if(is.character, function(x)gsub("\\\\x3b",";",gsub("\\\\x3d","=",x)))""",
     function2 = s"""rename_all(~ gsub("_$patient", "", .x))""",
     function3 = """rename_all(~ gsub("PGT", "PhasingInfo", .x))""",
     function4 = """rename_all(~ gsub("PID", "PhasingID", .x))""",
     function5 = """rename_all(~ gsub("SB", "PerSampleStrandBias", .x))""",
     function6 = """select(1:38, starts_with("FFPE"), starts_with("ctDNA_0"), starts_with("ctDNA_1"), starts_with("ctDNA_2"), starts_with("ctDNA_3"), starts_with("WB"))""")

 
}

val allVars = CSVListJoin(in = varsAnnotCSVFixed,
      fileCol = "Patient")

val allVarsFixed = CSVDplyr(csv1 = allVars,
     function1 = """mutate(ID = paste(Patient, CHROM, POS, REF, ALT, sep = "_"))""",
     function2 = """select(Patient, ID, everything())""",
     function3 = """filter(!grepl(",",ALT))""")

val allVarsFixedFilIndels = CSVDplyr(csv1 = allVarsFixed,
    script = """
              library(purrr)
              RefCount <- function(AD)as.numeric(map_chr(strsplit(as.character(AD),","),1))
              AltCount <- function(AD)as.numeric(map_chr(strsplit(as.character(AD),","),2))
              VAF <- function(AD)round(AltCount(AD)/(RefCount(AD) + AltCount(AD)),4)
              """,
     function1 = """mutate_at(vars(ends_with("AD")), ~ifelse(is.na(.x),"NA,NA",.x))""",
     function2 = """mutate_at(vars(ends_with("AD")), list(RefCount = RefCount, AltCount = AltCount, VAF = VAF))""",         
     function3 = """mutate(maxVAF = pmax(FFPE.AD_VAF, ctDNA_0.AD_VAF,  na.rm=T),
                           maxAlt = pmax(FFPE.AD_AltCount, ctDNA_0.AD_AltCount,  na.rm=T), 
                           Type = ifelse(nchar(REF) == 1 & nchar(ALT) == 1, "SNV", ifelse(nchar(REF) > 1 & nchar(ALT) > 1, "MNP", "INDEL")))""",
     function4 = """filter(Type %in% c("SNV","MNP") | (maxVAF > 0.01 & maxAlt > 15))""")
     

val allVarsEasyFormat = CSVDplyr(csv1 = allVarsFixedFilIndels,
     function1 = """select(Patient, CHROM, POS, REF, ALT, Func.refGene, Gene.refGene,ExonicFunc.refGene, CONTQ, GERMQ, ROQ, MBQ, ECNT, MPOS, NALOD, SEQQ, ends_with("AD"))""")

val allVarsExcel = CSV2Excel(csv = allVarsFixedFilIndels)

}