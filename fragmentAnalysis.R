library(ctDNAtools)
library(tidyverse)
library(furrr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(caTools)
library(viridis)
library(ggsignif)
library(ggbeeswarm)
library(ggpubr)
library(pheatmap)
library(NMF)
library(cowplot)
library(factoextra)
library(randomForest)
library(caTools)

x <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll4_5/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")
mutations <- x %>% 
  filter(ctDNA_0.AD_VAF > 0.001 ) %>%  
  filter(nchar(REF) == 1 & nchar(ALT) == 1) 

vafs <- mutations %>% 
  group_by(Patient) %>% 
  summarize(VAF = mean(ctDNA_0.AD_VAF))

list <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll4_5/list/out.csv",header=T,stringsAsFactors=F,sep="\t")
reference <- BSgenome.Hsapiens.UCSC.hg19
targets <- read.table("/mnt/storage1/rawdata/ctDNA/metadata/targets.bed", header = T, stringsAsFactors = F, sep= "\t")
colnames(targets) <- c("chr","start","end","gene")
bams <- list$File[grepl("ctDNA",list$Key)]

#fl_example <- get_fragment_size(bams[4], mutations = filter(mutations, Patient == "CHIC_102"))
#plot_density(fl_example)
#plan(multiprocess)

#fl <- future_map_dfc(bams, bin_fragment_size)
fl <- read.table("/mnt/storage2/work/amjad/ctdna/fragment_length_bins.csv", header = T, stringsAsFactors = F, sep = "\t")

ann <- data.frame(sample = colnames(fl), Type = ifelse(grepl("_0$",colnames(fl)),"before_treatment","mid-after_treatment")) %>%
  column_to_rownames("sample")
pheatmap(fl, cluster_rows=F, show_rownames=T, show_colnames=F, annotation = ann,  color = rev(magma(256, begin = 0.1, end = 0.9)),
         filename = "/mnt/storage2/work/amjad/ctdna/fragmentsPlots/frag_raw.pdf",
         width = 9, height = 4,labels_row = ifelse(c(1:nrow(fl))%%50 == 0,c(1:nrow(fl))*2,rep("",nrow(fl))), fontsize_row = 8)

set.seed(123)
nm <- nmf(as.matrix(fl), 4, nrun = 10)
bas <- as.data.frame(basis(nm)) %>%
   set_names(paste0("C",c(1:4)))
coe <- as.data.frame(coef(nm))
rownames(coe) <- paste0("C",c(1:4))
anncoe <- data.frame(sample = colnames(coe), Type = ifelse(grepl("_0$",colnames(coe)),"before_treatment","mid-after_treatment")) %>%
  column_to_rownames("sample")
pheatmap(coe, cluster_cols=T,cluster_rows=F, show_colnames=F, annotation = anncoe, filename = "/mnt/storage2/work/amjad/ctdna/fragmentsPlots/frag_coeff.pdf",
         width = 9, height = 4, color = rev(magma(256, begin = 0.1, end = 0.9)))

basis_labels <- ifelse(c(1:nrow(bas))%%50 == 0, c(1:nrow(bas))*2, rep("",nrow(bas)))
basis_labels[1] <- "1"
pheatmap(t(bas), cluster_cols=F,cluster_rows=F, filename = "/mnt/storage2/work/amjad/ctdna/fragmentsPlots/frag_basis.pdf",
         width = 9, height = 4, color = rev(magma(256, begin = 0.1, end = 0.9)), show_colnames = T,
         labels_col = basis_labels, fontsize_col = 8)

contribution <- as.data.frame(t(coe)) %>%
  set_names(paste0("C",c(1:4))) %>%
  mutate(Sample = rownames(anncoe),Type = anncoe$Type)

boxplots <- 
  purrr::map(paste0("C",c(1:4)),
   ~ contribution %>% ggplot(aes_string(x = "Type", y = .x, color = "Type")) + 
     geom_boxplot() + geom_quasirandom() + theme_minimal() + 
     geom_signif(comparisons = list(c("before_treatment", "mid-after_treatment")), color = "black") +
     theme(legend.title = element_blank(), legend.position = "none",  axis.line = element_line(size = 0.6)) + #axis.text.x = element_text(hjust = 0.1, angle = 45)) +
     labs(x = ""))

boxplots_comb <- plot_grid(plotlist = boxplots, ncol = 4)
ggsave(plot = boxplots_comb, filename = "/mnt/storage2/work/amjad/ctdna/fragmentsPlots/boxplots_prepostTr.pdf", width = 15, height = 7)

contributionVAF <- contribution %>% 
  dplyr::filter(Type == "before_treatment") %>% 
  mutate(patient = map_chr(strsplit(Sample,"_"), ~paste(.x[2],.x[3],sep="_"))) %>% 
  merge(vafs, by.x = "patient",by.y=1)

corplots <- purrr::map(paste0("C",c(1:4)),
                  ~ contributionVAF %>% 
                    ggscatter(x = .x, y = "VAF", add = "reg.line", cor.coef = T, 
                              cor.method = "pearson", color = magma(100)[50], conf.int = T))

corplots_comb <- plot_grid(plotlist = corplots, ncol = 4)
ggsave(plot = corplots_comb, filename = "/mnt/storage2/work/amjad/ctdna/fragmentsPlots/corplots.pdf", width = 15, height = 7)

samples <- data.frame(Sample = c("ctDNA_CHIC_100_2", "ctDNA_CHIC_102_2","ctDNA_CHIC_118_2","ctDNA_CHIC_123_2","ctDNA_CHIC_136_2","ctDNA_CHIC_143_2",
                                 "ctDNA_CHIC_2_2","ctDNA_CHIC_27_3", "ctDNA_CHIC_28_2", "ctDNA_CHIC_29_2","ctDNA_CHIC_12_2","ctDNA_CHIC_15_1","ctDNA_CHIC_16_2",
                                 "ctDNA_CHIC_3_2","ctDNA_CHIC_32_1","ctDNA_CHIC_35_2","ctDNA_CHIC_38_2","ctDNA_CHIC_4_2","ctDNA_CHIC_45_3","ctDNA_CHIC_49_2",
                                 "ctDNA_CHIC_52_2","ctDNA_CHIC_61_2","ctDNA_CHIC_63_2","ctDNA_CHIC_65_2","ctDNA_CHIC_70_2","ctDNA_CHIC_72_2","ctDNA_CHIC_73_2",
                                 "ctDNA_CHIC_74_2","ctDNA_CHIC_79_2","ctDNA_CHIC_88_2","ctDNA_CHIC_91_2","ctDNA_CHIC_92_2", 
                                 "ctDNA_CHIC_94_2", "ctDNA_CHIC_97_2","ctDNA_CHIC_99_2",
                                 "ctDNA_CHIC_10_2", "ctDNA_CHIC_112_2","ctDNA_CHIC_113_2","ctDNA_CHIC_120_2",
                                 "ctDNA_CHIC_122_2", "ctDNA_CHIC_126_2","ctDNA_CHIC_129_2", "ctDNA_CHIC_141_2", "ctDNA_CHIC_14_2","ctDNA_CHIC_18_2",
                                 "ctDNA_CHIC_1_2", "ctDNA_CHIC_25_2", "ctDNA_CHIC_33_2", "ctDNA_CHIC_34_2", "ctDNA_CHIC_36_2",
                                 "ctDNA_CHIC_53_2", "ctDNA_CHIC_54_2", "ctDNA_CHIC_57_2","ctDNA_CHIC_60_2", "ctDNA_CHIC_67_2",
                                 "ctDNA_CHIC_81_2", "ctDNA_CHIC_84_2", "ctDNA_CHIC_89_2")) %>%
  mutate(outcome = ifelse(Sample %in% c("ctDNA_CHIC_100_2", "ctDNA_CHIC_102_2","ctDNA_CHIC_136_2","ctDNA_CHIC_143_2","ctDNA_CHIC_27_3","ctDNA_CHIC_15_1",
                                        "ctDNA_CHIC_3_2","ctDNA_CHIC_32_1","ctDNA_CHIC_4_2","ctDNA_CHIC_52_2","ctDNA_CHIC_72_2","ctDNA_CHIC_74_2",
                                        "ctDNA_CHIC_97_2","ctDNA_CHIC_99_2", "ctDNA_CHIC_25_2", "ctDNA_CHIC_18_2"), "relapse","cured"))

contributionsRelCured <- merge(contribution, samples, by= "Sample")

boxplotsRelCured <- 
  purrr::map(paste0("C",c(1:4)),
             ~ contributionsRelCured %>% ggplot(aes_string(x = "outcome", y = .x, color = "outcome")) + 
               geom_boxplot() + geom_quasirandom() + theme_minimal() + 
               geom_signif(comparisons = list(c("relapse", "cured")), color = "black") +
               theme(legend.title = element_blank(), legend.position = "none", axis.line = element_line(size = 0.6)) +
               labs(x = ""))

boxplotsRelCured_comb <- plot_grid(plotlist = boxplotsRelCured, ncol = 4)
ggsave(plot = boxplotsRelCured_comb, filename = "/mnt/storage2/work/amjad/ctdna/fragmentsPlots/boxplots_cured_rel.pdf", width = 9, height = 4)


flpca <- prcomp(t(fl), center = T, scale = T)

pcaplot <- fviz_pca_ind(flpca, geom.ind = "point", pointshape = 21, pointsize = 2, 
                        fill.ind = anncoe$Type, col.ind = "black", palette = "jco", 
                        addEllipses = T, col.var = "black", repel = T)
ggsave(plot = pcaplot, filename = "/mnt/storage2/work/amjad/ctdna/fragmentsPlots/PCAplot.pdf", width = 9, height = 9)

pcs <- as.data.frame(flpca$x) %>%
  mutate(Sample = rownames(flpca$x), Type = anncoe$Type)

boxplotsPCA <- 
  purrr::map(paste0("PC",c(1:5)),
             ~ pcs %>% ggplot(aes_string(x = "Type", y = .x, color = "Type")) + 
               geom_boxplot() + geom_quasirandom() + theme_minimal() + 
               geom_signif(comparisons = list(c("before_treatment", "mid-after_treatment")), color = "black") +
               theme(legend.title = element_blank(), legend.position = "none",  axis.line = element_line(size = 0.6)) + #axis.text.x = element_text(hjust = 0.1, angle = 45)) +
               labs(x = ""))

boxplotsPCA_comb <- plot_grid(plotlist = boxplotsPCA, ncol = 5, nrow = 1)
ggsave(plot = boxplotsPCA_comb, filename = "/mnt/storage2/work/amjad/ctdna/fragmentsPlots/boxplotsPCA_prePostTr.pdf", width = 16, height = 7)

pcsRelCured <- merge(pcs, samples, by= "Sample")

boxplotsPCARelCured <- 
  purrr::map(paste0("PC",c(1:5)),
             ~ pcsRelCured %>% ggplot(aes_string(x = "outcome", y = .x, color = "outcome")) + 
               geom_boxplot() + geom_quasirandom() + theme_minimal() + 
               geom_signif(comparisons = list(c("relapse", "cured")), color = "black") +
               theme(legend.title = element_blank(), legend.position = "none", axis.line = element_line(size = 0.6)) +
               labs(x = ""))
boxplotsPCARelCured_comb <- plot_grid(plotlist = boxplotsPCARelCured, ncol = 5, nrow = 1)
ggsave(plot = boxplotsPCARelCured_comb, filename = "/mnt/storage2/work/amjad/ctdna/fragmentsPlots/boxplotsPCARelCured.pdf", width = 12, height = 6)

set.seed(123)
tr <- dplyr::select(pcsRelCured, PC1, PC2, PC3, PC4, PC5,  outcome)
tr$outcome = as.factor(tr$outcome)
rf = randomForest(outcome ~., data = tr, ntree= 10001)
rf
pred = data.frame(Sample = pcsRelCured$Sample, outcome = pcsRelCured$outcome, prediction = rf$predicted, rf$votes)
colAUC(pred$relapse, pred$outcome)
confusionMatrix(pred$prediction, pred$outcome, positive = "relapse")

pcsPre <- pcs %>%
  dplyr::filter(grepl("_0$",Sample)) %>%
  mutate(Patient = map_chr(strsplit(Sample, "_"), ~paste(.x[2],.x[3],sep = "_"))) %>%
  select(PC1pre = PC1, PC2pre = PC2, PC3pre = PC3, PC4pre = PC4, PC5pre = PC5, Patient)

pcsRelCuredWPre <- pcsRelCured %>%
  mutate(Patient = map_chr(strsplit(Sample, "_"), ~paste(.x[2],.x[3],sep = "_"))) %>%
  merge(pcsPre, by = "Patient") %>%
  select(PC1, PC2, PC3, PC4, PC5, PC1pre, PC2pre, PC3pre,  Sample, outcome) %>%
  mutate(outcome = factor(outcome))

rf2 <- randomForest(outcome ~., data = select(pcsRelCuredWPre, -Sample))
