#!/usr/bin/env anduril
//$OPT --threads 14
//$OPT -d /mnt/storage2/work/amjad/ctdna/result_newMutectAll

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

object ctdna{

  val list = INPUT(path = "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/bamCorrectedOutCSV/out.csv")
  val gnomad = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/af-only-gnomad.raw.sites.b37.vcf.gz")
  val chain = INPUT(path = "/mnt/storage1/rawdata/resources/chains/b37tohg19.chain")
  val reference = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.fasta")
  val referenceDict = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.dict")
  val targets = INPUT(path = "/mnt/storage1/rawdata/ctDNA/metadata/targets.bed")
  val exac = INPUT(path = "/mnt/storage1/rawdata/resources/hg38/small_exac_common_3.hg38.vcf.gz")
  val hg38tohg19chain = INPUT(path = "/mnt/storage1/rawdata/resources/chains/hg38ToHg19.over.chain")
  val hg38reference = INPUT(path = "/mnt/storage1/rawdata/resources/hg38/Homo_sapiens_assembly38.fasta")
 
  val picard = "/mnt/storage1/tools/picard/picard-2.18.26.jar"
  val fgbio = "/mnt/storage1/tools/fgbio/fgbio-0.7.0.jar"
  val gatk = "/mnt/storage1/tools/gatk/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar"
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
  
  val hg38Image = BashEvaluate(var1 = hg38reference,
        script = s"$java8 -jar $gatk BwaMemIndexImageCreator -I @var1@ -O @out1@")
  hg38Image._filename("out1","referencehg38.fasta.img")


  val listTumors = CSVDplyr(csv1 = list,
       function1 = """mutate(Patient = sapply(strsplit(Key,"_"),function(x)paste(x[2],x[3],sep="_")))""",
       function2 = """filter(!grepl("WB",Key))""",
       function3 = """select(KeyTumor = Key, Tumor = File, Patient)""")
  
  val listNormals = CSVDplyr(csv1 = list,
       function1 = """mutate(Patient = sapply(strsplit(Key,"_"),function(x)paste(x[2],x[3],sep="_")))""",
       function2 = """filter(grepl("WB",Key))""",
       function3 = """select(KeyNormal = Key,Normal = File, Patient)""")
  
  val listMatched = CSVDplyr(csv1=listTumors,
        csv2 = listNormals,
        script = """library(plyr)""",
        function1 = """plyr::join(csv2,by="Patient")""",
        function2 = """filter(!is.na(KeyNormal))""")
                                  
  val gnomadhg19 = BashEvaluate(var1 = gnomad,
       var2 = chain,
       var3 = reference,
       param1 = picard,
       script = """
        java -jar @param1@ LiftoverVcf \
        INPUT=@var1@ OUTPUT=@out1@ CHAIN=@var2@ REJECT=@out2@ REFERENCE_SEQUENCE=@var3@
        """)
  gnomadhg19._filename("out1","gnomAD_hg19.vcf.gz")
  gnomadhg19._filename("out2","rejected.vcf.gz")

  

  
  val normal = NamedMap[INPUT]("normal")
  val ponNormal = NamedMap[BashEvaluate]("ponNormal")
  val ponNormalOut = NamedMap[Any]("ponNormalOut")
  
  for ( rowMap <- iterCSV(listNormals) ) { 
    val key = rowMap("KeyNormal")
    normal(key) = INPUT(path=rowMap("Normal"))


    ponNormal(key) = BashEvaluate(var1 = normal(key),
       var2 = targets,
       var3 = reference,
       param2 = gatk,
       param3 = java8,
       script = """
            cat @var2@ | awk 'NR>1' > @out2@
            @param3@ -jar @param2@ Mutect2 \
              -R @var3@ \
              -I @var1@ \
              -L @out2@ -ip 300 \
              -O @out1@ --max-mnp-distance 0          
                """)
   ponNormal(key)._filename("out1","Normal_"+key+".vcf.gz")
   ponNormal(key)._filename("out2","intervals.bed")
   ponNormal(key)._filename("out3","Normal_filtered" + key + ".vcf.gz")


   ponNormalOut(key) = ponNormal(key).out1
  }

  val ponDB =  BashEvaluate(var1 = reference,
        var2 = targetsIL.out1,  
        array1 = ponNormalOut,
        script = s"$java8 -jar $gatk GenomicsDBImport -R @var1@ -L @var2@ -ip 300 --genomicsdb-workspace-path @folder1@/genDB" +
                       """ $( paste -d ' ' <(getarrayfiles array1)  | sed 's,^, -V ,' | tr -d '\\\n' ) """)
   ponDB._filename("out1","pon_db")             

  val pon = BashEvaluate(var1 = reference,
        var2 = ponDB.folder1,
        param1 = java8,
        param2 = gatk,
        script = """
         cd @var2@
         @param1@ -jar @param2@ CreateSomaticPanelOfNormals -R @var1@ -V gendb://genDB -O @out1@
          """)
  pon._filename("out1","pon.vcf.gz")



   
  val control   = NamedMap[INPUT]("control")
  val tumor     = NamedMap[INPUT]("tumor")
  val contam    = NamedMap[BashEvaluate]("contam")
  val outContam = NamedMap[Any]("outContam")
  val outSegs   = NamedMap[Any]("outSegs")

  for ( rowMap <- iterCSV(listMatched) ) { 

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
     table2 = listMatched,
     script = """
         table.out <- data.frame()
         array.out <- split(table1, table2$Patient)
              """)

val outContamArrayCSV = Array2CSV(in = outContamArray.outArray)

val outSegsArray = REvaluate(table1 = outSegsCSV,
     table2 = listMatched,
     script = """
         table.out <- data.frame()
         array.out <- split(table1, table2$Patient)
              """)
val outSegsArrayCSV = Array2CSV(in = outSegsArray.outArray)


val outBamArray = REvaluate(table1 = listMatched,
     script = """
         library(dplyr)
         table.out <- data.frame()
         table2 <- select(table1, Key = KeyTumor, File = Tumor)
         array.out <- split(table2, table1$Patient)
              """)
val outBamArrayCSV = Array2CSV(in = outBamArray.outArray)

val outBamArrayCSVmerged = CSVDplyr(csv1 = outBamArrayCSV,
     csv2 = listMatched,
     function1 = """merge(csv2[!duplicated(csv2$Patient),c("Patient","Normal","KeyNormal")], by = 1, sort = F)""")



/* 
   Loop over patients and
   call variants from multiple tumor samples
   learn a model for F1R2 bias
   filter variants
   annotate variants
*/

val F1R2Model             = NamedMap[BashEvaluate]("F1R2Model")
val BamsByPatient         = NamedMap[CSV2Array]("BamsByPatient")
val contamByPatient       = NamedMap[CSV2Array]("contamByPatient")
val segsByPatient         = NamedMap[CSV2Array]("segsByPatient")
val normalControl         = NamedMap[INPUT]("normalControl")
val vars                  = NamedMap[BashEvaluate]("vars")
val varsFiltered          = NamedMap[BashEvaluate]("varsFiltered")
val varsPass              = NamedMap[BashEvaluate]("varsPass")
val varsAnnot             = NamedMap[Annovar]("varsAnnot")
val varsAnnotCSV          = NamedMap[BashEvaluate]("varsAnnotCSV")
val varsAnnotCSVFixed     = NamedMap[CSVDplyr]("varsAnnotCSVFixed")
val varsPassOut           = NamedMap[Any]("varsPassOut")
val varsFilteredAlignments = NamedMap[BashEvaluate]("varsFilteredAlignments")

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

  if(patient != "CHIC_143"){
  vars(patient) = BashEvaluate(array1 = BamsByPatient(patient),
        var1 = normalControl(patient),
        var2 = reference,
        var3 = pon.out1,
        var4 = gnomadhg19.out1,
        var5 = targetsIL.out1,
        param1 = keyNormal,
        script = s"$java8 -jar $gatk Mutect2 -R @var2@ -I @var1@ -O @out1@ --max-reads-per-alignment-start 0 --pcr-indel-model HOSTILE  --bam-output @out3@ " +
                  " --germline-resource @var4@ --panel-of-normals @var3@ -L @var5@ -ip 300 --f1r2-tar-gz @out2@ -normal @param1@ --max-mnp-distance 0 " +
                  """ $( paste -d ' ' <(getarrayfiles array1)  | sed 's,^, -I ,' | tr -d '\\\n' ) """)        
  vars(patient)._filename("out1", patient + "_rawVariants.vcf.gz")
  vars(patient)._filename("out2", patient + "_f1r2.tar.gz")
  vars(patient)._filename("out3", patient + ".bam")
  // --force-active true --tumor-lod-to-emit 0 --initial-tumor-lod 0 -bamout @out3@
  // --pcr-indel-model AGGRESSIVE
  } else {
    vars(patient) = BashEvaluate(array1 = BamsByPatient(patient),
        var1 = normalControl(patient),
        var2 = reference,
        var3 = pon.out1,
        var4 = gnomadhg19.out1,
        var5 = targetsIL.out1,
        param1 = keyNormal,
        script = s"$java8 -jar $gatk Mutect2 -R @var2@ -O @out1@ --max-reads-per-alignment-start 0 --pcr-indel-model HOSTILE  --bam-output @out3@ " +
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
               " --max-events-in-region 50 --min-median-read-position 15 -L @var4@ -ip 300 " +
             // """ $( paste -d ' ' <(getarrayfiles array1)  | sed 's,^, --contamination-table ,' | tr -d '\\\n' ) """ +
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
       script = s"$java8 -jar $gatk SelectVariants --exclude-filtered  -L @var2@ -ip 10 -V @var1@ -O @out1@")
  varsPass(patient)._filename("out1", patient + "_passed.vcf.gz")

  varsPassOut(patient) = varsPass(patient).out1

  varsAnnot(patient) = Annovar(vcfIn = varsPass(patient).out1,
            annovarPath = annovar,
            annovardb = annovardb,
            buildver = "hg19",
            inputType = "vcf",
            protocol = "refGene,avsnp147,cosmic68,dbnsfp30a,icgc21",
            operation = "g,f,f,f,f")

// if(patient != "CHIC_143"){
  varsAnnotCSV(patient) = BashEvaluate(var1 = reference,
            var2 = varsAnnot(patient).vcfOut,
            param1 = gatk,
            param2 = java8,
            script = """ 
                @param2@ -jar @param1@ VariantsToTable \
                    -R @var1@ -V @var2@ -O @out1@ \
                    -F CHROM -F POS -F REF -F ALT -F ID \
                    -F Func.refGene -F Gene.refGene -F GeneDetail.refGene -F ExonicFunc.refGene -F AAChange.refGene \
                    -F avsnp147 -F cosmic68 -F SIFT_score -F SIFT_pred \
                    -F Polyphen2_HDIV_score -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_score -F Polyphen2_HVAR_pred \
                    -F MutationTaster_score -F MutationTaster_pred -F MutationAssessor_score -F MutationAssessor_pred \
                    -F CONTQ -F GERMQ -F ROQ -F MBQ -F ECNT -F MMQ -F MPOS -F NALOD -F POPAF -F SEQQ  \
                    -F CADD_phred -F DANN_score -GF AD -GF DP -GF AF -GF F1R2 -GF F2R1 -GF PGT -GF PID 
            """)

  varsAnnotCSVFixed(patient) = CSVDplyr(csv1 = varsAnnotCSV(patient).out1,
     script = """
              library(purrr)
              RefCount <- function(AD)as.numeric(map_chr(strsplit(AD,","),1))
              AltCount <- function(AD)as.numeric(map_chr(strsplit(AD,","),2))
              VAF <- function(AD)round(AltCount(AD)/(RefCount(AD) + AltCount(AD)),4)
              """,
     function1 = """mutate_if(is.character, function(x)gsub("\\\\x3b",";",gsub("\\\\x3d","=",x)))""",
     function2 = s"""rename_all(~ gsub("_$patient", "", .x))""",
     function3 = """rename_all(~ gsub("PGT", "PhasingInfo", .x))""",
     function4 = """rename_all(~ gsub("PID", "PhasingID", .x))""",
     function5 = """rename_all(~ gsub("SB", "PerSampleStrandBias", .x))""",
     function6 = """mutate_at(vars(ends_with("AD")), list(RefCount = RefCount, AltCount = AltCount, VAF = VAF))""",
     function7 = """select(1:38, starts_with("FFPE"), starts_with("ctDNA_0"), starts_with("ctDNA_1"), starts_with("ctDNA_2"), starts_with("ctDNA_3"), starts_with("WB"))""")

// }
}

val allVars = CSVListJoin(in = varsAnnotCSVFixed,
      fileCol = "Patient")

val allVarsFixed = CSVDplyr(csv1 = allVars,
     function1 = """mutate(ID = paste(Patient, CHROM, POS, REF, ALT, sep = "_"))""",
     function2 = """select(Patient, ID, everything())""",
     function3 = """filter(!grepl(",",ALT))""")

val allVarsFixedFilIndels = CSVDplyr(csv1 = allVarsFixed,
     function1 = """mutate(maxVAF = pmax(FFPE.AD_VAF, ctDNA_0.AD_VAF, ctDNA_1.AD_VAF, ctDNA_2.AD_VAF, ctDNA_3.AD_VAF, na.rm=T),
                           maxAlt = pmax(FFPE.AD_AltCount, ctDNA_0.AD_AltCount, ctDNA_1.AD_AltCount, ctDNA_2.AD_AltCount, ctDNA_3.AD_AltCount, na.rm=T), 
                           Type = ifelse(nchar(REF) == 1 & nchar(ALT) == 1, "SNV", ifelse(nchar(REF) > 1 & nchar(ALT) > 1, "MNP", "INDEL")))""",
     function2 = """filter(Type %in% c("SNV","MNP") | (maxVAF > 0.01 & maxAlt > 15))""",
     function3 = """filter(!Patient %in% c("CHIC_123","CHIC_143") | ctDNA_2.AD_AltCount < 100)""")

val allVarsEasyFormat = CSVDplyr(csv1 = allVarsFixedFilIndels,
     function1 = """select(Patient, CHROM, POS, REF, ALT, Func.refGene, Gene.refGene,ExonicFunc.refGene, ends_with("AD"))""")

val allVarsExcel = CSV2Excel(csv = allVarsFixedFilIndels)


// Estimate the background
// --------------------------------------------

/*
val oneVCF = BashEvaluate(var1 = reference,
     array1 = varsPassOut,
     script = s"$java8 -jar $gatk3 -T CombineVariants -R @var1@ -o @out1@ -genotypeMergeOptions UNIQUIFY " +
     """ $( paste -d ' ' <(getarraykeys array1) <(getarrayfiles array1)  | sed 's,^, --variant:,' | tr -d '\\\n' ) """)
oneVCF._filename("out1","variants.vcf")   
*/

val background = NamedMap[BashEvaluate]("background")
val bamIn = NamedMap[INPUT]("bamIn")
val backgroundCSV = NamedMap[BashEvaluate]("backgroundCSV")
val backgroundCSVFixed = NamedMap[CSVDplyr]("backgroundCSVFixed")

for ( rowMap <- iterCSV(list) ) { 
  
  val sample = rowMap("Key")
  bamIn(sample) = INPUT(path = rowMap("File"))

  background(sample) = BashEvaluate(
        var1 = bamIn(sample),
        var2 = reference,
        var3 = targetsIL.out1,
        script = s"$java8 -jar $gatk Mutect2 -R @var2@ -I @var1@ -O @out1@ --max-reads-per-alignment-start 0 -L @var3@ -ERC BP_RESOLUTION  --bam-output @out2@ ")
  background(sample)._filename("out1", sample + "_background.vcf.gz")
  background(sample)._filename("out2", sample + "_background.bam")
 
  backgroundCSV(sample) = BashEvaluate(var1 = reference,
        var2 = background(sample).out1,
        param1 = gatk,
        param2 = java8,
        script = """ 
                   @param2@ -jar @param1@ VariantsToTable \
                    -R @var1@ -V @var2@ -O @out1@ \
                    -F CHROM -F POS -F REF -F ALT -F ID \
                    -GF AD -GF DP -GF TLOD 
                     """)

 backgroundCSVFixed(sample) = CSVDplyr(csv1 = backgroundCSV(sample).out1,
       script = "library(tidyr)",
       function1 = """filter(ALT == "<NON_REF>")""",
       function2 = s"""separate($sample.AD, into = c("REFreads","ALTreads"), convert = T)""",
       function3 = s"""select(CHROM, POS, REF, ALT, REFreads, ALTreads, Depth = $sample.DP, TLOD = $sample.TLOD)""",
       function4 = """filter(TLOD < (0))""")
}


 val backgroundAll = Rpipe(array1 = backgroundCSVFixed,
      function1 = """csvOut <- data.frame(Sample = names(array1), Rate = map_dbl(array1, ~ sum(.x$ALTreads)/sum(as.numeric(.x$Depth))))""")

// ------------------------------------------------------------
// Background estimation ends here

}