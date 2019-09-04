library(tidyverse)
library(reshape)
library(viridis)
library(ggthemes)
library(cowplot)
library(ggsci)

#x <- read.table("/mnt/storage1/work/amjad/ctdna/result_ctdnaVariants/varsAllFixed/out.csv", header = T, stringsAsFactors = F, sep= "\t")
#list <- read.table("/mnt/storage1/work/amjad/ctdna/result_ctdnaAlignment/bamOutCSV/out.csv", header = T, stringsAsFactors = F, sep= "\t")
#oldVars <- read.table("/mnt/storage1/work/amjad/ctdna/oldVars.csv", header = T, stringsAsFactors = F, sep = "\t")
#corrected <- read.table("/mnt/storage1/work/amjad/ctdna/result_ctdnaVariantsCorrected/varsAllFixed/out.csv", header = T, stringsAsFactors = F, sep = "\t")
newMutect <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll/allVarsFixedFilIndels/out.csv", header = T, stringsAsFactors = F, sep = "\t")
newList <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll/listMatched/out.csv", header = T, stringsAsFactors = F, sep = "\t")
stats <- read.table("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/stats/out.csv",header=T, stringsAsFactors = F, sep = "\t")
dupStats <- read.table("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentAll/dupStats/out.csv",header=T, stringsAsFactors = F, sep = "\t")

stats %>% 
  filter(!grepl("Normal", Sample)) %>%
  mutate(class = case_when(grepl("ctDNA",Sample) ~ "ctDNA",
                           grepl("FFPE", Sample) ~ "FFPE",
                           grepl("WB", Sample) ~ "WB")) %>%
  arrange(class, -MEAN_BAIT_COVERAGE) %>%
  mutate(Sample = factor(Sample, levels = Sample)) %>%
          ggplot(aes(x = Sample, y = MEAN_BAIT_COVERAGE, fill = class)) + geom_bar(stat = "identity") +
          theme_wsj() + scale_fill_wsj() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                                 legend.title = element_blank())

dupStats %>% 
  filter(!grepl("Normal", Sample)) %>%
  mutate(class = case_when(grepl("ctDNA",Sample) ~ "ctDNA",
                           grepl("FFPE", Sample) ~ "FFPE",
                           grepl("WB", Sample) ~ "WB")) %>%
  arrange(class, -PERCENT_DUPLICATION) %>%
  mutate(Sample = factor(Sample, levels = Sample)) %>%
  ggplot(aes(x = Sample, y = PERCENT_DUPLICATION, fill = class)) + geom_bar(stat = "identity") +
  theme_wsj() + scale_fill_wsj() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                         legend.title = element_blank())


sequencedSamples <- newList %>%
  mutate(sample = map_chr(strsplit(KeyTumor, "_"),
                          ~ paste0(.x[2], "_", .x[3], "_", .x[1], ifelse(length(.x) == 4, paste0("_", .x[4]),"")))) %>% .$sample
  
  
newMutectSelectColumns <- newMutect %>%
        select(Patient, ID, REF, ALT, Func.refGene, ExonicFunc.refGene, ends_with("VAF"), ends_with("RefCount"), ends_with("AltCount")) 

vars <- map_dfr(c("WB","FFPE", "ctDNA_0","ctDNA_1","ctDNA_2","ctDNA_3"), 
                ~ mutate(newMutectSelectColumns, key = .x) %>% 
                  select(Patient, key, ID, REF, ALT, Func.refGene, ExonicFunc.refGene, starts_with(.x)) %>% 
                  rename_all(function(y) gsub(paste0(.x,"."),"",y))) %>%
        filter(!grepl("WB", key)) %>%
        mutate(sample = paste(Patient, key, sep = "_")) %>%
        mutate(sample = map_chr(strsplit(sample, "_"),
                                ~ paste0(.x[3], "_", .x[1], "_", .x[2], ifelse(length(.x) == 4, paste0("_", .x[4]),"")))) %>%
        filter(AD_VAF > 0.001) %>%
        mutate(Tumor.ratioVAF = AD_VAF)


ctDNABarPlot <- function(mutations, type = c("mutation","FFPE","effect","overlap"), whatOverlap = c("FFPE","ctDNA_0","ctDNA_1","ctDNA_2","ctDNA_3"), customTheme = theme_wsj(),
                         customFill = scale_fill_viridis(discrete = T,option = "A",begin = 0.2, end  = 0.8, direction = -1), ncol = 7){
  require(dplyr)
  require(viridis)
  require(ggplot2)
  require(purrr)

  sampleLevels <- unlist(map(unique(mutations$Patient), ~ c(paste0(.x, "_FFPE"),paste0(.x, "_ctDNA_0"), 
                     paste0(.x, "_ctDNA_1"), paste0(.x, "_ctDNA_2"),
                     paste0(.x, "_ctDNA_3"))))
  
  type <- match.arg(type)
  whatOverlap <- match.arg(whatOverlap)
  if(type == "mutation"){
    mutations <- mutations %>% mutate(type = ifelse(nchar(REF) == 1 & nchar(ALT) == 1,"SNV",ifelse(nchar(REF)>1 & nchar(ALT)>1,"MNV","INDEL")))
  } else if(type == "FFPE"){
    mutations <- mutations %>% mutate(type = ifelse((REF == "C" & ALT == "T") | (REF == "G" & ALT == "A"), "FFPEmotif","Else"))
  } else if(type == "effect"){
    mutations <- mutations %>% mutate(type = ifelse(ExonicFunc.refGene == ".", Func.refGene,ExonicFunc.refGene))
  } else if(type == "overlap"){
    IDs <- mutations %>% filter(key == whatOverlap) %>% .$ID
    mutations <- mutations %>% mutate(type = ifelse(ID %in% IDs, paste("In", whatOverlap), paste("Not in", whatOverlap)))
  }
  mut <- mutations %>% 
      mutate(sample = paste(Patient, key, sep = "_")) %>%
      select(sample, type, Patient) %>%
      mutate(sample = factor(sample, levels = sampleLevels))
      

  mut %>% count(sample, type, .drop = F) %>%
    mutate(Patient = factor(unlist(map(strsplit(as.character(sample),"_"), ~paste(.x[1],.x[2],sep = "_"))))) %>%
    ggplot(aes(x = sample, y = n, fill = type)) + 
    geom_bar(stat = "identity", position = "stack") + 
    facet_wrap(~Patient, scales = "free_x", #space = "free_x",scales = "free_x",
               ncol = ncol) + 
    customTheme + 
    #theme(axis.text.x = element_text(angle = 90)) +
    theme(axis.text.x = element_blank()) +
    customFill
}

ctDNAVafPlot <- function(mutations, threshold = 2, patients = unique(mutations$Patient),
                        sequencedSamples, customTheme = theme_wsj(), logY = T,
                        customColor = scale_color_d3(), customFill = scale_fill_d3()){
  sampleLevels <- unlist(map(unique(mutations$Patient), ~ c(paste0(.x, "_FFPE"),paste0(.x, "_ctDNA_0"), 
                                                            paste0(.x, "_ctDNA_1"), paste0(.x, "_ctDNA_2"),
                                                            paste0(.x, "_ctDNA_3"))))
  

  vaf = mutations %>% 
      mutate(sample = paste(Patient, key, sep = "_")) %>%
      mutate(sample = factor(sample,levels = sampleLevels)) %>%
      group_by(sample, .drop = F) %>% 
      summarize(N = n(),
                VAFmean = ifelse(N > threshold, mean(Tumor.ratioVAF),0),
                VAFsd   = ifelse(N > threshold, sd(Tumor.ratioVAF),0), 
                VAFse   = ifelse(N > threshold, VAFsd/sqrt(N),0)) %>%
      ungroup() %>%
      mutate(VAFmean = ifelse(sample %in% sequencedSamples, VAFmean, NA),
             VAFsd = ifelse(sample %in% sequencedSamples, VAFsd, NA),
             VAFse = ifelse(sample %in% sequencedSamples, VAFse, NA)) %>%
      mutate(series = case_when(
                 grepl("FFPE", sample) ~ 0,
                 grepl("ctDNA_0",   sample) ~ 1,
                 grepl("ctDNA_1",   sample) ~ 2,
                 grepl("ctDNA_2",   sample) ~ 3,
                 grepl("ctDNA_3",   sample) ~ 4)
             ) %>%
    mutate(Patient = factor(unlist(map(strsplit(as.character(sample),"_"), ~paste(.x[1],.x[2],sep = "_")))))
    
  plot = vaf %>%
        filter(Patient %in% patients) %>%
        ggplot(aes(x = series, y = VAFmean, colour = Patient, fill = Patient)) + 
        geom_errorbar(aes(ymin = VAFmean - VAFse, ymax = VAFmean + VAFse), width = 0.05, position = position_dodge(0.1), size = 1) +
        geom_line(position=position_dodge(0.1), size = 1) +
        geom_point(position = position_dodge(0.1), size=3, shape=21) +
        customTheme + 
        labs(x = "", y = "Variant Allele Frequency") +
        scale_x_continuous(breaks = c(0:4),labels = c("FFPE","ctDNA 0","ctDNA 1","ctDNA 2","ctDNA 3"))# + 
        #customColor + customFill
  if(logY){
    plot <- plot + scale_y_log10() 
  }
  vaf <- vaf %>% filter(Patient %in% patients)
  return(list(vaf = vaf, plot = plot))
}

relapseOverlap <-  ctDNABarPlot(vars, type = "overlap", whatOverlap = "ctDNA_3", customTheme = theme_linedraw() + theme(panel.grid = element_line(colour = "grey90")))
baselineOverlap <-  ctDNABarPlot(vars, type = "overlap", whatOverlap = "ctDNA_0", customTheme = theme_linedraw() + theme(panel.grid = element_line(colour = "grey90")))
ffpeOverlap <- ctDNABarPlot(vars, type = "overlap", whatOverlap = "FFPE",customTheme = theme_linedraw() + theme(panel.grid = element_line(colour = "grey90")))

ggsave(plot = relapseOverlap, filename = "/mnt/storage2/work/amjad/ctdna/plots/relapseOverlap.png")
ggsave(plot = baselineOverlap, filename = "/mnt/storage2/work/amjad/ctdna/plots/baselineOverlap.png")
ggsave(plot = ffpeOverlap, filename = "/mnt/storage2/work/amjad/ctdna/plots/ffpeOverlap.png")

withRelapse <- c("CHIC_136", "CHIC_52", "CHIC_97")
ffpeBaseline <- c("CHIC_100","CHIC_12","CHIC_15","CHIC_16","CHIC_2","CHIC_27","CHIC_29","CHIC_35","CHIC_38",
                  "CHIC_42","CHIC_49","CHIC_63","CHIC_70","CHIC_73","CHIC_88","CHIC_91","CHIC_92","CHIC_99",
                  "CHIC_118","CHIC_28","CHIC_39","CHIC_61", "CHIC_123", "CHIC_4","CHIC_94")

baseLineRelapsePlots <- newMutect %>% 
    filter(Patient %in% c("CHIC_97","CHIC_136")) %>% 
  ggplot(aes(x = ctDNA_0.AD_VAF, y = ctDNA_3.AD_VAF)) + 
  geom_point(size = 3, color = "skyblue4") + theme_linedraw() + 
  theme(panel.grid = element_line(colour = "grey90")) +
  facet_grid(~Patient, scales = "free")

ffpeRelapsePlots <- newMutect %>% 
  filter(Patient %in% c("CHIC_52")) %>% 
  ggplot(aes(x = FFPE.AD_VAF, y = ctDNA_3.AD_VAF)) + 
  geom_point(size = 3, color = "skyblue4") + theme_linedraw() + 
  theme(panel.grid = element_line(colour = "grey90")) +
  facet_grid(~Patient, scales = "free")

ffpeBaselinePlots <-  newMutect %>% 
  filter(Patient %in% ffpeBaseline) %>% 
  ggplot(aes(x = FFPE.AD_VAF, y = ctDNA_0.AD_VAF)) + 
  geom_point(size = 1.5, color = "skyblue4") + theme_linedraw() + 
  theme(panel.grid = element_line(colour = "grey90")) +
  facet_wrap(~Patient, scales = "free", ncol = 5)

ggsave(plot = baseLineRelapsePlots, filename = "/mnt/storage2/work/amjad/ctdna/plots/baseLineRelapsePlots.png")
ggsave(plot = ffpeRelapsePlots, filename = "/mnt/storage2/work/amjad/ctdna/plots/ffpeRelapsePlots.png")
ggsave(plot = ffpeBaselinePlots, filename = "/mnt/storage2/work/amjad/ctdna/plots/ffpeBaselinePlots.png")

overlaps <- newMutect %>%
  filter(Patient %in% ffpeBaseline) %>%
  filter(FFPE.AD_VAF > 0 | ctDNA_0.AD_VAF > 0) %>%
  group_by(Patient) %>%
  summarize(Shared = sum(FFPE.AD_VAF > 0 & ctDNA_0.AD_VAF > 0, na.rm = T),
            FFPE = sum(FFPE.AD_VAF > 0 & ctDNA_0.AD_VAF == 0, na.rm = T),
            Baseline = sum(FFPE.AD_VAF == 0 & ctDNA_0.AD_VAF > 0, na.rm =T)) %>%
  as.data.frame() 

overlapsPlot <- overlaps %>% melt() %>% mutate(variable = factor(variable, levels = rev(levels(variable)))) %>%
  ggplot(aes(x = Patient, y = value, fill = variable)) + geom_bar(stat = "identity", position = "stack") +
      theme_pander() + theme(axis.text.x = element_text(angle = 60, hjust = 1.1), legend.title = element_blank()) + 
      scale_fill_viridis(discrete = T, option = "A", alpha = 0.8, direction = -1,begin = 0.15, end = 0.8) +
      ylab("# Mutations") + xlab("")

ggsave(plot = overlapsPlot, filename = "/mnt/storage2/work/amjad/ctdna/plots/overlapsPlot.png")
write.table(overlaps, file = "/mnt/storage2/work/amjad/ctdna/plots/FFPEbaselineOverlaps.csv", col.names = T, row.names = F, sep = "\t", quote = F)

vafPlot <- ctDNAVafPlot(vars, sequencedSamples = sequencedSamples, logY = F)$plot
ggsave(plot = vafPlot, filename = "/mnt/storage2/work/amjad/ctdna/plots/vafPlot.png")

boxplots <- newMutect %>% select(Patient, ends_with("VAF"), -maxVAF) %>%
     melt() %>% filter(variable != "WB.AD_VAF") %>% ggplot(aes(x = variable, y = value, colour = variable)) + 
       geom_boxplot() + ggbeeswarm::geom_quasirandom() + 
       facet_wrap(~Patient, ncol = 5, scales = "free") + theme_linedraw() + 
       theme(panel.grid = element_line(colour = "grey90"), axis.text.x = element_blank()) + 
       scale_color_d3() + xlab("") + ylab("VAF")

ggsave(plot = boxplots, filename = "/mnt/storage2/work/amjad/ctdna/plots/VAFboxplots.png")

vafs <- newMutect %>% 
  filter(ctDNA_0.AD_VAF > 0.005) %>% 
  group_by(Patient) %>% summarize(meanVAF = mean(ctDNA_0.AD_VAF), n = n()) %>% 
  as.data.frame()

write.table(vafs, file = "/mnt/storage2/work/amjad/ctdna/plots/VAFs_threshold0.005.csv", col.names=T,row.names=F, sep = "\t", quote = F)
