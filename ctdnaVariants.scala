#!/usr/bin/env anduril
//$OPT --threads 10

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

object ctdna{

  val list = INPUT(path = "/mnt/storage1/work/amjad/ctdna/result_ctdnaAlignment/bamOutCSV/out.csv")
  val gnomad = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/af-only-gnomad.raw.sites.b37.vcf.gz")
  val chain = INPUT(path = "/mnt/storage1/rawdata/resources/chains/b37tohg19.chain")
  val reference = INPUT(path = "/mnt/storage1/rawdata/resources/hg19/ucsc.hg19.fasta")
  val targets = INPUT(path = "/mnt/storage1/rawdata/ctDNA/metadata/targets.bed")
  val exac = INPUT(path = "/mnt/storage1/rawdata/resources/hg38/small_exac_common_3.hg38.vcf.gz")
  val hg38tohg19chain = INPUT(path = "/mnt/storage1/rawdata/resources/chains/hg38ToHg19.over.chain")
  
  val picard = "/mnt/storage1/tools/picard/picard-2.18.26.jar"
  val fgbio = "/mnt/storage1/tools/fgbio/fgbio-0.7.0.jar"
  val gatk = "/mnt/storage1/tools/gatk/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"
  val java8 = "/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java"
  val annovar = "/mnt/storage1/tools/Annovar/annovar/"
  val annovardb = "/mnt/storage1/tools/Annovar/annovar/humandb/"

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
       function3 = """select(KeyNormal = Key,Normal = File, Patient)""")
  
   val listMatched = CSVDplyr(csv1=listTumors,
        csv2=listNormals,
        script="""library(plyr)""",
        function1="""plyr::join(csv2,by="Patient")""")
                                  
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
       param1 = key,
       param2 = gatk,
       param3 = java8,
       script = """
            cat @var2@ | awk 'NR>1' > @out2@
            @param3@ -jar @param2@ Mutect2 \
              -R @var3@ \
              -I @var1@ \
              -L @out2@ -ip 300 \
              -tumor @param1@ \
              -O @out1@
                """)
   ponNormal(key)._filename("out1","Normal_"+key+".vcf.gz")
   ponNormal(key)._filename("out2","intervals.bed")

   ponNormalOut(key) = ponNormal(key).out1
  }

  val pon =  BashEvaluate(var1=reference,
                array1=ponNormalOut,
                script=s"$java8 -jar $gatk CreateSomaticPanelOfNormals -O @out1@" +
                       """ $( paste -d ' ' <(getarrayfiles array1)  | sed 's,^, -vcfs ,' | tr -d '\\\n' ) """)
   pon._filename("out1","pon.vcf.gz")             

   
  val control = NamedMap[INPUT]("control")
  val tumor  = NamedMap[INPUT]("tumor")
  val vars = NamedMap[BashEvaluate]("vars")
  val contam = NamedMap[BashEvaluate]("contam")
  val varsFiltered = NamedMap[BashEvaluate]("varsFiltered")
  val artifactMetrics = NamedMap[BashEvaluate]("artifactMetrics")
  val varsFiltered2 = NamedMap[BashEvaluate]("varsFiltered2")
  val varsPass = NamedMap[BashEvaluate]("varsPass")
  val varsAnnot = NamedMap[Annovar]("varsAnnot")
  val varsAnnotCSV = NamedMap[BashEvaluate]("varsAnnotCSV")
  val varsOut = NamedMap[CSVTransformer]("varsOut")
  
  for ( rowMap <- iterCSV(listMatched) ) { 

    val keyT = rowMap("KeyTumor")
    val keyN = rowMap("KeyNormal")
    tumor(keyT) = INPUT(path = rowMap("Tumor"))
    control(keyT) = INPUT(path=rowMap("Normal"))
    
    vars(keyT) = BashEvaluate(var1 = reference,
              var2 = tumor(keyT),
              var3 = control(keyT),
              var4 = gnomadhg19.out1,
              var5 = pon.out1,
              var6 = ponNormal(keyN).out2,
              param1 = keyT,
              param2 = keyN,
              param3 = gatk,
              param4 = java8,
              script = """
                   @param4@ -jar @param3@  Mutect2 -R @var1@ -I @var2@ -I @var3@ -tumor @param1@ -normal @param2@ \
                   -pon @var5@ --germline-resource @var4@ --af-of-alleles-not-in-resource 0.0000025 \
                   -L @var6@ -ip 300 -O @out1@ -bamout @out2@ 
                       """)
   vars(keyT)._filename("out1", keyT + ".vcf.gz")
   vars(keyT)._filename("out2", keyT + ".bam")                   
   
   
   contam(keyT) = BashEvaluate(var1 = tumor(keyT),
     var2 = exachg19.out1,
     param1 = gatk,
     param2 = java8,
     script = """
         @param2@ -jar @param1@  GetPileupSummaries -I @var1@ \
        -V @var2@ -L @var2@ -O @out1@ 
         @param2@ -jar @param1@ CalculateContamination -I @out1@ -O @out2@     
              """)
   contam(keyT)._filename("out1", keyT + "_pileupSummary.table")
   contam(keyT)._filename("out2", keyT + "_contam.table")
   
   
   // basic filtering
   varsFiltered(keyT) = BashEvaluate(var1 = vars(keyT).out1,
      var2 = contam(keyT).out2,
      param1 = gatk,
      param2 = java8,
      script = """
          @param2@ -jar @param1@  FilterMutectCalls -V @var1@ \
          --max-events-in-region 7 --normal-artifact-lod -1 --tumor-lod 6.3  --min-median-read-position 10 --contamination-table @var2@  -O @out1@   
               """)
   varsFiltered(keyT)._filename("out1", keyT + "_filteredVariants.vcf.gz")
   
   
   artifactMetrics(keyT) = BashEvaluate(var1 = reference,
            var2 = tumor(keyT),
            param1 = gatk,
            param2 = java8,
            script = """
                @param2@ -Xmx4G -jar @param1@ CollectSequencingArtifactMetrics \
                    -R @var1@ -I @var2@ -O @out1@ 
            """)
   artifactMetrics(keyT)._filename("out1", keyT + "_artifactMetrics.txt")
   
   
   // artifact filtering
   varsFiltered2(keyT) = BashEvaluate(var1 = varsFiltered(keyT).out1,
            var2 = artifactMetrics(keyT).out1,
            param1 = gatk,
            param2 = java8,
            script = """
                @param2@ -jar @param1@ FilterByOrientationBias \
                    -AM 'C/T' -AM 'G/T' \
                    -V @var1@ \
                    -P @var2@.pre_adapter_detail_metrics \
                    -O @out1@
            """)
   varsFiltered2(keyT)._filename("out1",keyT + "_doubleFilteredVariants.vcf.gz")     
   
   
   varsPass(keyT) = BashEvaluate(var1 = varsFiltered2(keyT).out1,
       param1 = gatk,
       param2 = java8,
       script = """
            @param2@ -jar @param1@ SelectVariants --exclude-filtered \
              -V @var1@ -O @out1@     
                """)
   varsPass(keyT)._filename("out1", keyT + "_passed.vcf.gz")
   
   
   varsAnnot(keyT) = Annovar(vcfIn = varsPass(keyT).out1,
            annovarPath = annovar,
            annovardb = annovardb,
            buildver = "hg19",
            inputType = "vcf",
            protocol = "refGene,avsnp147,cosmic68,dbnsfp30a,icgc21",
            operation = "g,f,f,f,f")
   
   
   varsAnnotCSV(keyT) = BashEvaluate(var1 = reference,
            var2 = varsAnnot(keyT).vcfOut,
            param1 = gatk,
            param2 = java8,
            script = """ 
                cat @var2@ | sed -r '/^\s*$/d' > @out2@
                @param2@ -jar @param1@ VariantsToTable \
                    -R @var1@ -V @out2@ -O @out1@ \
                    -F CHROM -F POS -F REF -F ALT -F ID \
                    -F Func.refGene -F Gene.refGene -F GeneDetail.refGene -F ExonicFunc.refGene -F AAChange.refGene \
                    -F avsnp147 -F cosmic79 -F SIFT_score -F SIFT_pred \
                    -F Polyphen2_HDIV_score -F Polyphen2_HDIV_pred -F Polyphen2_HVAR_score -F Polyphen2_HVAR_pred \
                    -F MutationTaster_score -F MutationTaster_pred -F MutationAssessor_score -F MutationAssessor_pred \
                    -F CADD_phred -F DANN_score -GF AD -GF DP -GF AF -GF SAPP
            """)
   varsAnnotCSV(keyT)._filename("out2", keyT + "_annotated.vcf")
   
   varsOut(keyT) = CSVTransformer(csv1 = varsAnnotCSV(keyT).out1,
       transform1 = s"""
              colnames(csv1) = gsub("$keyT","Tumor",colnames(csv1))
              colnames(csv1) = gsub("$keyN","Normal",colnames(csv1))
              csv1      
                    """)
   
   }

   val varsAll = CSVListJoin(in = varsOut,
     fileCol = "sample") 

  val varsAllFixed = CSVDplyr(csv1 = varsAll,
     script = "library(tidyr)",
     function1 = "select(-Normal.SAPP)",
     function2 = """mutate_if(is.character, function(x)gsub("\\\\x3b",";",gsub("\\\\x3d","=",x)))""",
     function3 = """filter(!grepl(",",ALT))""",
     function4 = """separate(Normal.AD, into = c("Normal.Ref","Normal.Var"), convert = T)""",
     function5 = """separate(Tumor.AD, into = c("Tumor.Ref","Tumor.Var"), convert = T)""",
     function6 = """mutate(Tumor.ratioVAF = round(Tumor.Var/(Tumor.Var+Tumor.Ref),4), Normal.ratioVAF = round(Normal.Var/(Normal.Var + Normal.Ref),4))""")
        //            filter(Tumor.ratioVaf > 0.03
     
  val varsAllFixedExcel = CSV2Excel(csv = varsAllFixed)   
}
