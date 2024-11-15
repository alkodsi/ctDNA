library(tidyverse)
y <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutect2_2/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")
x <- read.table("/mnt/storage2/work/amjad/ctdna/VariantsVersions/variants_seq2.2_v1.csv",header=T,stringsAsFactors = F, sep = "\t")
back <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutect2_2/backgroundAll/csvOut.csv", header = T, stringsAsFactors = F, sep = "\t")

backgroundRates <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll4_5/backgroundAll/csvOut.csv", header = T, stringsAsFactors = F, sep = "\t")
y <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")
x <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectNova/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")
x <- read.table("/mnt/storage2/work/amjad/ctdna/result_forcedCallingPipeline/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")
x <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll4_5/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")

backgroundRates <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectNova/backgroundAll/csvOut.csv", header = T, stringsAsFactors = F, sep = "\t")

mutations <- read.table("/mnt/storage2/work/amjad/ctdna/code/ctdna/result_mutectDoubleCalling/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")
m <- read.table("/mnt/storage2/work/amjad/myoma/result_myomaPipeline2/allVarsFiltered/out.csv", header = T, stringsAsFactors = F, sep = "\t")
a = get_fragment_size("/mnt/storage2/work/amjad/myoma/Bams//Patient3/18X-0203_S7_recal_reads.bam", mutations = m[m$Patient == "Patient3" & nchar(m$REF) == 1 & nchar(m$ALT) == 1,])

convolve.binomial <- function(p) {
  # p is a vector of probabilities of Bernoulli distributions.
  # The convolution of these distributions is returned as a vector
  # `z` where z[i] is the probability of i-1, i=1, 2, ..., length(p)+1.
  n <- length(p) + 1
  z <- c(1, rep(0, n-1))
  for (p in p) z <- (1-p)*z + p*c(0, z[-n])
  return(z)
}

findProb <- function(depths, nVar, errorRate, thr = 0){
  depths <- depths[!is.na(depths)]
  if(length(depths) == 0){
    return(NA)
  } else {
    n <- length(depths)
    n_nonzero <- sum(nVar > thr, na.rm = T)
    prob_error <- rep(errorRate, n)
    thresh_reads <- rep(thr + 1, n)
  
    mismatch_probs <- pmap_dbl(list(prob_error, depths, thresh_reads), function(x, y, z)sum(dbinom(z:y, y, x)))
    prob <- sum(convolve.binomial(mismatch_probs)[(n_nonzero+1):(n+1)])
    return(prob)
  }
}

#x = stackStringsFromBam("result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_74_2/ctDNA_CHIC_74_2_consensusRecal.bam", use.names=T, param = gr)

samplesToTest <- backgroundRates %>%
  mutate(patient = map_chr(strsplit(Sample,"_"), ~paste(.x[2],.x[3], sep = "_")),
         series = map_chr(strsplit(Sample,"_"), ~paste(.x[1],.x[4], sep = "_"))) %>%
  filter((patient %in% c("CHIC_15","CHIC_32") & grepl("_1$", Sample)) |
           (patient %in% c("CHIC_2","CHIC_27","CHIC_45") & grepl("_3$", Sample)) |
           (!patient %in% c("CHIC_2","CHIC_27","CHIC_45") & grepl("_2$", Sample))) %>%
  filter(!patient %in% c("CHIC_143","CHIC_64"))

samplesToTest <- backgroundRates %>% 
  filter(grepl("ctDNA",Sample) & !grepl("_0", Sample) & !grepl("Normal",Sample) ) %>%
  mutate(patient = map_chr(strsplit(Sample,"_"), ~paste(.x[2],.x[3], sep = "_")),
         series = map_chr(strsplit(Sample,"_"), ~paste(.x[1],.x[4], sep = "_"))) %>%
  filter(!patient %in% c("CHIC_64")) 
  filter(!patient %in% c("CHIC_52","CHIC_64","CHIC_16","CHIC_38")) #c("CHIC_16","CHIC_38","CHIC_52"))

wbSamples <-  backgroundRates %>% 
  filter(grepl("WB",Sample)) %>%
  mutate(patient = map_chr(strsplit(Sample,"_"), ~paste(.x[2],.x[3], sep = "_"))) %>%
  filter(!patient %in% c("CHIC_16","CHIC_38","CHIC_52","CHIC_143","CHIC_28","CHIC_4"))

mutations =  x %>% 
  filter(ctDNA_0.AD_VAF > 0.001 | FFPE.AD_VAF > 0.01) %>%  
  filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%
  filter(ECNT <= 10) %>%
  mutate(WB.PhasingID = ifelse(is.na(WB.PhasingID), ID, WB.PhasingID)) %>%
  filter(!duplicated(WB.PhasingID))

library(BSgenome.Hsapiens.UCSC.hg19)
reference = BSgenome.Hsapiens.UCSC.hg19
targets <- read.table("/mnt/storage1/rawdata/ctDNA/metadata/targets.bed", header = T, stringsAsFactors = F, sep= "\t")
colnames(targets) <- c("chr","start","end","gene")

samples <- data.frame(Sample = c(#"ctDNA_CHIC_100_2", "ctDNA_CHIC_102_2","ctDNA_CHIC_118_2","ctDNA_CHIC_123_2","ctDNA_CHIC_136_2","ctDNA_CHIC_143_2",
                                 #"ctDNA_CHIC_2_2","ctDNA_CHIC_27_3", "ctDNA_CHIC_28_2", "ctDNA_CHIC_29_2","ctDNA_CHIC_12_2","ctDNA_CHIC_15_1","ctDNA_CHIC_16_2",
                                 #"ctDNA_CHIC_3_2","ctDNA_CHIC_32_1","ctDNA_CHIC_35_2","ctDNA_CHIC_38_2","ctDNA_CHIC_4_2","ctDNA_CHIC_45_3","ctDNA_CHIC_49_2",
                                 #"ctDNA_CHIC_52_2","ctDNA_CHIC_61_2","ctDNA_CHIC_63_2","ctDNA_CHIC_65_2","ctDNA_CHIC_70_2","ctDNA_CHIC_72_2","ctDNA_CHIC_73_2",
                                 #"ctDNA_CHIC_74_2","ctDNA_CHIC_79_2","ctDNA_CHIC_88_2","ctDNA_CHIC_91_2","ctDNA_CHIC_92_2", 
                                 #"ctDNA_CHIC_94_2", "ctDNA_CHIC_97_2","ctDNA_CHIC_99_2",
                                 "ctDNA_CHIC_10_2", "ctDNA_CHIC_112_2","ctDNA_CHIC_113_2","ctDNA_CHIC_120_2",
                                 "ctDNA_CHIC_122_2", "ctDNA_CHIC_126_2","ctDNA_CHIC_129_2", "ctDNA_CHIC_141_2", "ctDNA_CHIC_14_2","ctDNA_CHIC_18_2",
                                 "ctDNA_CHIC_1_2", "ctDNA_CHIC_25_2", "ctDNA_CHIC_33_2", "ctDNA_CHIC_34_2", "ctDNA_CHIC_36_2",
                                 "ctDNA_CHIC_53_2", "ctDNA_CHIC_54_2", "ctDNA_CHIC_57_2","ctDNA_CHIC_60_2", "ctDNA_CHIC_67_2",
                                 "ctDNA_CHIC_81_2", "ctDNA_CHIC_84_2", "ctDNA_CHIC_89_2")) %>%
           mutate(outcome = ifelse(Sample %in% c("ctDNA_CHIC_100_2", "ctDNA_CHIC_102_2","ctDNA_CHIC_136_2","ctDNA_CHIC_143_2","ctDNA_CHIC_27_3","ctDNA_CHIC_15_1",
                                                 "ctDNA_CHIC_3_2","ctDNA_CHIC_32_1","ctDNA_CHIC_4_2","ctDNA_CHIC_52_2","ctDNA_CHIC_72_2","ctDNA_CHIC_74_2",
                                                 "ctDNA_CHIC_97_2","ctDNA_CHIC_99_2"), "relapse","cured"))

samples <- data.frame(Sample = c("ctDNA_CHIC_100_2", "ctDNA_CHIC_102_2","ctDNA_CHIC_118_2","ctDNA_CHIC_123_2","ctDNA_CHIC_136_2","ctDNA_CHIC_143_2",
  "ctDNA_CHIC_2_2","ctDNA_CHIC_27_3", "ctDNA_CHIC_28_2", "ctDNA_CHIC_29_2","ctDNA_CHIC_12_2","ctDNA_CHIC_15_1","ctDNA_CHIC_16_2",
  "ctDNA_CHIC_3_2","ctDNA_CHIC_32_1","ctDNA_CHIC_35_2","ctDNA_CHIC_38_2","ctDNA_CHIC_4_2","ctDNA_CHIC_45_3","ctDNA_CHIC_49_2",
  "ctDNA_CHIC_52_2","ctDNA_CHIC_61_2","ctDNA_CHIC_63_2","ctDNA_CHIC_65_2","ctDNA_CHIC_70_2","ctDNA_CHIC_72_2","ctDNA_CHIC_73_2",
  "ctDNA_CHIC_74_2","ctDNA_CHIC_79_2","ctDNA_CHIC_88_2","ctDNA_CHIC_91_2","ctDNA_CHIC_92_2", 
  "ctDNA_CHIC_94_2", "ctDNA_CHIC_97_2","ctDNA_CHIC_99_2")) %>%
  mutate(outcome = ifelse(Sample %in% c("ctDNA_CHIC_100_2", "ctDNA_CHIC_102_2","ctDNA_CHIC_136_2","ctDNA_CHIC_143_2","ctDNA_CHIC_27_3","ctDNA_CHIC_15_1",
                                        "ctDNA_CHIC_3_2","ctDNA_CHIC_32_1","ctDNA_CHIC_4_2","ctDNA_CHIC_52_2","ctDNA_CHIC_72_2","ctDNA_CHIC_74_2",
                                        "ctDNA_CHIC_97_2","ctDNA_CHIC_99_2"), "relapse","cured"))
results = future_map2_dfr(samples$Sample, map_chr(strsplit(as.character(samples$Sample),"_"), ~paste(.x[2],.x[3], sep = "_")), 
                            ~test_ctDNA(mutations[mutations$Patient == .y,], reference = reference, targets = targets, 
                                        bam = paste0("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_",.x,"/",.x,"_consensusRecal.bam"),
                                        by_substitution = F, ID_column = "ctDNA_1.PhasingID",min_samples = 4, min_alt_reads = 1 ),
                            .progress = T)
colAUC(samples$p, samples$outcome)

test_ctDNA(mutations[mutations$Patient == "CHIC_97",], reference = reference, targets = targets, 
           bam = "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_97_2/ctDNA_CHIC_97_2_consensusRecal.bam",
           black_list = bl$black_list, ID_column = "ctDNA_1.PhasingID" )

test_ctDNA(mutations[mutations$Patient == "CHIC_91",], reference = reference, targets = targets, 
           bam = "/mnt/storage2/work/amjad/ctdna/result_newMutectAll/vars_CHIC_91/CHIC_91.bam",tag = "ID1.3", by_substitution = F, min_base_quality = 10)
           
mut <- filter_mutations(mutations[mutations$Patient == "CHIC_91",], 
                        bams = c("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_118_2/ctDNA_CHIC_118_2_consensusRecal.bam", 
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_2_3/ctDNA_CHIC_2_3_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_12_2/ctDNA_CHIC_12_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_123_2/ctDNA_CHIC_123_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_88_2/ctDNA_CHIC_88_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_92_2/ctDNA_CHIC_92_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_45_3/ctDNA_CHIC_45_3_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_49_2/ctDNA_CHIC_49_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_70_2/ctDNA_CHIC_70_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_79_2/ctDNA_CHIC_79_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_91_2/ctDNA_CHIC_91_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_94_2/ctDNA_CHIC_94_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_29_2/ctDNA_CHIC_29_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_73_2/ctDNA_CHIC_73_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_63_2/ctDNA_CHIC_63_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_65_2/ctDNA_CHIC_65_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_61_2/ctDNA_CHIC_61_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_42_2/ctDNA_CHIC_42_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_39_2/ctDNA_CHIC_39_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_38_2/ctDNA_CHIC_38_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_35_2/ctDNA_CHIC_35_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_28_2/ctDNA_CHIC_28_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_16_2/ctDNA_CHIC_16_2_consensusRecal.bam",
                                 "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_Normal_1//ctDNA_Normal_1_consensusRecal.bam"), 
                        min_alt_reads = 0, min_samples = 2)

bams = c("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_118_0/ctDNA_CHIC_118_0_consensusRecal.bam", 
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_2_0/ctDNA_CHIC_2_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_12_0/ctDNA_CHIC_12_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_123_0/ctDNA_CHIC_123_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_88_0/ctDNA_CHIC_88_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_92_0/ctDNA_CHIC_92_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_45_0/ctDNA_CHIC_45_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_49_0/ctDNA_CHIC_49_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_70_0/ctDNA_CHIC_70_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_79_0/ctDNA_CHIC_79_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_91_0/ctDNA_CHIC_91_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_94_0/ctDNA_CHIC_94_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_29_2/ctDNA_CHIC_29_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_73_2/ctDNA_CHIC_73_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_63_2/ctDNA_CHIC_63_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_65_2/ctDNA_CHIC_65_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_61_2/ctDNA_CHIC_61_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_42_2/ctDNA_CHIC_42_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_39_2/ctDNA_CHIC_39_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_38_2/ctDNA_CHIC_38_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_35_2/ctDNA_CHIC_35_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_28_2/ctDNA_CHIC_28_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_16_2/ctDNA_CHIC_16_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_Normal_1//ctDNA_Normal_1_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_100_0/ctDNA_CHIC_100_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_102_0/ctDNA_CHIC_102_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_97_0/ctDNA_CHIC_97_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_72_0/ctDNA_CHIC_72_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_74_0/ctDNA_CHIC_74_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_15_0/ctDNA_CHIC_15_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_99_0/ctDNA_CHIC_99_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_3_0/ctDNA_CHIC_3_0_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_100_2/ctDNA_CHIC_100_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_102_2/ctDNA_CHIC_102_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_97_2/ctDNA_CHIC_97_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_72_2/ctDNA_CHIC_72_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_74_2/ctDNA_CHIC_74_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_15_1/ctDNA_CHIC_15_1_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_99_2/ctDNA_CHIC_99_2_consensusRecal.bam",
         "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_3_2/ctDNA_CHIC_3_2_consensusRecal.bam")

out <- future_map_dfc(bams, function(x){
   h <- hist(get_fragment_size(x)$size, breaks = seq(0, 400, by = 2), plot = F)
   return(h$counts/sum(h$counts))
}, .progress = T)
  
results <- data.frame(Sample = samplesToTest$Sample,
                      Pvalue = pmap_dbl(list(samplesToTest$patient, samplesToTest$series, samplesToTest$Rate),
                                        function(x,y,z) positivity_test(depths = mutations[mutations$Patient == x, paste0(y,".DP")],
                                               rate = list(rate = z/3, CT = NA, CA = NA, CG = NA, TC = NA, TG = NA, TA = NA),  
                                               altReads = mutations[mutations$Patient == x, paste0(y,".AD_AltCount")],
                                               seed = 123, n_simulations = 10000))) %>%
           mutate( positivity = ifelse(Pvalue < 0.05, "Pos","Neg"))
write.table(resultsClass, file = "/mnt/storage2/work/amjad/ctdna/ctDNAprediction.csv", col.names= T, row.names = F, sep = "\t", quote = F)

wbResults <- data.frame(Sample = wbSamples$Sample,
                        Pvalue = pmap_dbl(list(wbSamples$patient, "WB", wbSamples$Rate),
                                          function(x,y,z) positivityTest(depths = mutations[mutations$Patient == x, paste0(y,".DP")],
                                                                         rate = z/3,  
                                                                         altReads = mutations[mutations$Patient == x, paste0(y,".AD_AltCount")],
                                                                         seed = 123, nPermutation = 10000)))
  
scores <- x %>% mutate(
  ctDNA_0.pval = binom.test.na(x = ctDNA_0.AD_AltCount, n = ctDNA_0.AD_RefCount + ctDNA_0.AD_AltCount, 
                            alternative = "greater", p = 0.0025/20),
  ctDNA_1.pval = binom.test.na(x = ctDNA_1.AD_AltCount, n = ctDNA_1.AD_RefCount + ctDNA_1.AD_AltCount, 
                            alternative = "greater", p = 0.0025/20),
  ctDNA_2.pval = binom.test.na(x = ctDNA_2.AD_AltCount, n = ctDNA_2.AD_RefCount + ctDNA_2.AD_AltCount, 
                            alternative = "greater", p = 0.0025/20),
  ctDNA_3.pval = binom.test.na(x = ctDNA_3.AD_AltCount, n = ctDNA_3.AD_RefCount + ctDNA_3.AD_AltCount, 
                            alternative = "greater", p = 0.0025/20)
) %>%
  filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%
  group_by(Patient) %>%
  mutate( 
    ctDNA_0.fdr = p.adjust(ctDNA_0.pval, method = "fdr"),
    ctDNA_1.fdr = p.adjust(ctDNA_1.pval, method = "fdr"),
    ctDNA_2.fdr = p.adjust(ctDNA_2.pval, method = "fdr"),
    ctDNA_3.fdr = p.adjust(ctDNA_3.pval, method = "fdr")
  ) %>% ungroup() %>% as.data.frame()

summary <- scores %>% 
  filter(ctDNA_0.fdr < 0.05) %>%
  group_by(Patient) %>%
  summarize(
    score1 = sum(ctDNA_1.fdr < 0.1)/n(),
    score2 = sum(ctDNA_2.fdr < 0.1)/n()
  )

scores <- x %>% 
  filter(ctDNA_0.AD_VAF > 0.005) %>%
  filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%
  group_by(Patient) %>%
  summarize(prob1 = findProb(ctDNA_1.DP, ctDNA_1.AD_AltCount, 0.0025/3),
         prob2 = findProb(ctDNA_2.DP, ctDNA_2.AD_AltCount, 0.0025/3),
         prob3 = findProb(ctDNA_3.DP, ctDNA_3.AD_AltCount, 0.0025/3))

summary <- x %>% group_by(Patient) %>% 
  filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%
  filter(ctDNA_0.AD_VAF > 0.005) %>%
  summarize(
  sumAlt1 = sum(ctDNA_1.AD_AltCount, na.rm = T),
  sumAlt2 = sum(ctDNA_2.AD_AltCount, na.rm = T),
 # meanVaf1 = mean(ctDNA_1.AD_VAF, na.rm = T),
#  meanVaf2 = mean(ctDNA_2.AD_VAF, na.rm = T),
 # meanVaf1Norm = mean(ctDNA_1.AD_VAF, na.rm = T)/mean(ctDNA_0.AD_VAF, na.rm = T),
#  meanVaf2Norm = mean(ctDNA_2.AD_VAF, na.rm = T)/mean(ctDNA_0.AD_VAF, na.rm = T),
  nonzeropercent1 = sum(ctDNA_1.AD_AltCount>0, na.rm = T)/n(),
  nonzeropercent2 = sum(ctDNA_2.AD_AltCount>0, na.rm = T)/n(),
  nonzeropercent3 = sum(ctDNA_3.AD_AltCount>0, na.rm = T)/n()
#  meanNonZero1 = sum(ctDNA_1.AD_VAF, na.rm = T)/sum(ctDNA_1.AD_VAF>0, na.rm = T),
#  meanNonZero2 = sum(ctDNA_2.AD_VAF, na.rm = T)/sum(ctDNA_2.AD_VAF>0, na.rm = T),
#  sumAltNorm1 = sum(ctDNA_1.AD_AltCount, na.rm = T)/(sum(ctDNA_1.AD_AltCount, na.rm = T) + sum(ctDNA_1.AD_RefCount, na.rm = T)),
#  sumAltNorm2 = sum(ctDNA_2.AD_AltCount, na.rm = T)/(sum(ctDNA_2.AD_AltCount, na.rm = T) + sum(ctDNA_2.AD_RefCount, na.rm = T))
) %>%
  
  
  mutate_if(is.numeric, ~round(.x,5))




#' A function to plot distribution of a feature between train and test stratified by quake
#' @param train training set data frame
#' @param test test set data frame
#' @param feature name of the feature
#' @param binsPerRange a histogram binning parameter
#' @return a ggplot
#' @export
plotDist <- function(sample1, color = "sample", binwidth = 2.5){
  require(ggthemes)
  require(ggplot2)

  ggplot() + geom_freqpoly(data = sample1,aes_string(x = "size", y = "..density..", color = color), binwidth = binwidth) +
    #scale_colour_wsj() +
    theme_minimal()
}    


list <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll4_5/list/out.csv", header = T, stringsAsFactors = F, sep = "\t")
targets <- read.table("/mnt/storage1/rawdata/ctDNA/metadata/targets.bed", header = T, stringsAsFactors = F, sep= "\t")
targets$start = targets$start - 100
colnames(targets) <- c("chr","start","end","gene")
summaries <- future_map(list[grepl("ctDNA",list$Key),"File"], ~ summarize_fragment_size(bam = .x, regions = targets), .progress = T)
frag_matrix <- summaries %>% purrr::reduce(inner_join, by = "Region")

bam <- "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_100_0/ctDNA_CHIC_100_0_consensusRecal.bam"
a = estimate_ctDNA_level(mutations = mutations[mutations$Patient == "CHIC_70",],
                      bam = "/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_ctDNA_CHIC_70_0/ctDNA_CHIC_70_0_consensusRecal.bam",
                      reference = reference ,targets = targets, vaf_column = "ctDNA_0.AD_VAF", ref_reads_column = "ctDNA_0.AD_RefCount", 
                      alt_reads_column = "ctDNA_0.AD_AltCount")

list <- read.table("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/bamCorrectedOutCSV/out.csv",header=T,stringsAsFactors=F,sep="\t")
samples <-  list[grepl("ctDNA",list[,1]) & grepl("_0$",list[,1]) ,1]
bams <- list[grepl("ctDNA",list[,1]) & grepl("_0$",list[,1]) ,2]
patients = map_chr(str_split(samples, "_"), ~paste(.x[2], .x[3], sep = "_"))

vafs <- future_map2_dbl(patients, bams, ~estimate_ctDNA_level(mutations = mutations[mutations$Patient == .x,], bam = .y, reference = reference,
                                                       targets = targets, vaf_column = "ctDNA_0.AD_VAF", ref_reads_column = "ctDNA_0.AD_RefCount", 
                                                       alt_reads_column = "ctDNA_0.AD_AltCount", only_SNVs = T, use_clustering = F)$VAF, .progress = T)


mutations <- mutations %>% mutate(motif = ifelse(substitution %in% c("CT","CG","CA") & substring(context,1,1) %in% c("A","G") & substring(context,3,3) %in% c("A","C","T"), "RCH",
       ifelse(substitution %in% c("TA","TC","TG")  & substring(context,3,3) %in% c("T","A"),"TW","Other")))
