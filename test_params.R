library(ctDNAtools)
library(tidyverse)
library(furrr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(caTools)

x <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll4_5/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")
list <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll4_5/list/out.csv",header=T,stringsAsFactors=F,sep="\t")
reference <- BSgenome.Hsapiens.UCSC.hg19
targets <- read.table("/mnt/storage1/rawdata/ctDNA/metadata/targets.bed", header = T, stringsAsFactors = F, sep= "\t")
colnames(targets) <- c("chr","start","end","gene")

samples <- data.frame(Sample = c("ctDNA_CHIC_100_2", "ctDNA_CHIC_102_2","ctDNA_CHIC_118_2","ctDNA_CHIC_123_2","ctDNA_CHIC_136_2","ctDNA_CHIC_143_2",
                                 "ctDNA_CHIC_2_2","ctDNA_CHIC_27_3", "ctDNA_CHIC_28_2", "ctDNA_CHIC_29_2","ctDNA_CHIC_12_2","ctDNA_CHIC_15_1","ctDNA_CHIC_16_2",
                                 "ctDNA_CHIC_3_2","ctDNA_CHIC_32_1","ctDNA_CHIC_35_2","ctDNA_CHIC_38_2","ctDNA_CHIC_4_2","ctDNA_CHIC_45_3","ctDNA_CHIC_49_2",
                                 "ctDNA_CHIC_52_2","ctDNA_CHIC_61_2","ctDNA_CHIC_63_2","ctDNA_CHIC_65_2","ctDNA_CHIC_70_2","ctDNA_CHIC_72_2","ctDNA_CHIC_73_2",
                                 "ctDNA_CHIC_74_2","ctDNA_CHIC_79_2","ctDNA_CHIC_88_2","ctDNA_CHIC_91_2","ctDNA_CHIC_92_2", 
                                 "ctDNA_CHIC_94_2", "ctDNA_CHIC_97_2","ctDNA_CHIC_99_2")) %>%
  mutate(outcome = ifelse(Sample %in% c("ctDNA_CHIC_100_2", "ctDNA_CHIC_102_2","ctDNA_CHIC_136_2","ctDNA_CHIC_143_2","ctDNA_CHIC_27_3","ctDNA_CHIC_15_1",
                                        "ctDNA_CHIC_3_2","ctDNA_CHIC_32_1","ctDNA_CHIC_4_2","ctDNA_CHIC_52_2","ctDNA_CHIC_72_2","ctDNA_CHIC_74_2",
                                        "ctDNA_CHIC_97_2","ctDNA_CHIC_99_2"), "relapse","cured"))


bams <- list[grepl("ctDNA",list[,1]) & grepl("_2$",list[,1]) & !list[,1] %in% c("ctDNA_CHIC_102_2","ctDNA_CHIC_72_2","ctDNA_CHIC_74_2"),2]
bamSamples <- list[grepl("ctDNA",list[,1]) & grepl("_2$",list[,1]) & !list[,1] %in% c("ctDNA_CHIC_102_2","ctDNA_CHIC_72_2","ctDNA_CHIC_74_2"),1]

plan(multiprocess)
auc_phasing <- c()
auc_noPhasing <- c()
ecnt <- seq(10, 50, by = 5)

for(i in ecnt){
  mutations =  x %>% 
    filter(ctDNA_0.AD_VAF > 0.001 | FFPE.AD_VAF > 0.01) %>%  
    filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%
    filter(ECNT <= i)
  
  samples$p = future_map2_dbl(samples$Sample, map_chr(strsplit(as.character(samples$Sample),"_"), ~paste(.x[2],.x[3], sep = "_")), 
                              ~test_ctDNA(mutations[mutations$Patient == .y,], reference = reference, targets = targets, 
                                          bam = paste0("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_",.x,"/",.x,"_consensusRecal.bam"),
                                          bam_list = bams[bamSamples!=.x], by_substitution = F, ID_column = "ctDNA_1.PhasingID", use_unique_molecules = F,min_samples = 4, min_alt_reads = 1 )$pvalue,
                              .progress = T)
  
  auc_phasing <- c(auc_phasing, colAUC(samples$p, samples$outcome)[1,1])
  
  samples$p = future_map2_dbl(samples$Sample, map_chr(strsplit(as.character(samples$Sample),"_"), ~paste(.x[2],.x[3], sep = "_")), 
                              ~test_ctDNA(mutations[mutations$Patient == .y,], reference = reference, targets = targets, 
                                          bam = paste0("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_",.x,"/",.x,"_consensusRecal.bam"),
                                          bam_list = bams[bamSamples!=.x], by_substitution = F,  use_unique_molecules = F,min_samples = 4, min_alt_reads = 1 )$pvalue,
                              .progress = T)
 
  auc_noPhasing <- c(auc_noPhasing, colAUC(samples$p, samples$outcome)[1,1])
  
}


filterComb <- expand.grid(min_samples = c(1:6), min_alt_reads = c(0:3))
filterComb$auc <- NA

mutations =  x %>% 
  filter(ctDNA_0.AD_VAF > 0.001 | FFPE.AD_VAF > 0.01) %>%  
  filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%
  filter(ECNT <= 10)

samples <- samples[samples$Sample != "ctDNA_CHIC_4_2",]

for(i in 1:nrow(filterComb)){
      samples$p = future_map2_dbl(samples$Sample, map_chr(strsplit(as.character(samples$Sample),"_"), ~paste(.x[2],.x[3], sep = "_")), 
                                ~test_ctDNA(mutations[mutations$Patient == .y,], reference = reference, targets = targets, 
                                            bam = paste0("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/consensusRecal_",.x,"/",.x,"_consensusRecal.bam"),
                                            bam_list = bams[bamSamples!=.x], by_substitution = F, ID_column = "ctDNA_1.PhasingID", use_unique_molecules = F,
                                            min_samples = filterComb$min_samples[i], min_alt_reads = filterComb$min_alt_reads[i])$pvalue,
                                .progress = T)
    
       colAUC(samples$p, samples$outcome)
       filterComb$auc[i] <- colAUC(samples$p, samples$outcome)[1,1]
}


