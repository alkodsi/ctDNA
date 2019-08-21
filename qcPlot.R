library(tidyverse)
oldMutect <- read.table("/mnt/storage1/work/amjad/ctdna/result_ctdnaVariantsCorrected/varsAllFixed/out.csv", header = T, stringsAsFactors = F, sep = "\t")
newMutectStrictIndel <- read.table("/mnt/storage1/work/amjad/ctdna/result_newMutect/allVarsFixed/out.csv", header = T, stringsAsFactors = F, sep = "\t")

oldMutect$ID <- paste(map_chr(strsplit(oldMutect$sample,"_"), ~ paste(.x[2:3],collapse = "_")), oldMutect$CHROM, oldMutect$POS, oldMutect$REF, oldMutect$ALT,sep="_")
uid <- unique(corrected$ID)

commonVariants <- intersect(uid, newMutect$ID)

newCalls <- filter(newMutect, !ID %in% commonVariants) # 101
missedCalls <- filter(oldMutect, !ID %in% commonVariants) # 504


statsNovaseqNotrim <- read.table("/mnt/storage2/work/amjad/ctdna/statsNoTrimming.csv", header = T, stringsAsFactors = F, sep = "\t") %>% 
  mutate(category = "NOVAseq-No-QC", sampleType = case_when(
    grepl("FFPE", Sample) ~ "FFPE",
    grepl("ctDNA_",   Sample) ~ "ctDNA",
    grepl("WB",   Sample) ~ "WB"))

statsNovaseqTrim <- read.table("/mnt/storage2/work/amjad/ctdna/statsTrimming.csv", header = T, stringsAsFactors = F, sep = "\t") %>%
  mutate(category = "NOVAseq-Trim", sampleType = case_when(
    grepl("FFPE", Sample) ~ "FFPE",
    grepl("ctDNA_",   Sample) ~ "ctDNA",
    grepl("WB",   Sample) ~ "WB"))

statsNovaseqTile <- read.table("/mnt/storage2/work/amjad/ctdna/statsTileFiltering.csv", header = T, stringsAsFactors = F, sep = "\t") %>%
  mutate(category = "NOVAseq-TileFiltering", sampleType = case_when(
    grepl("FFPE", Sample) ~ "FFPE",
    grepl("ctDNA_",   Sample) ~ "ctDNA",
    grepl("WB",   Sample) ~ "WB"))

pilotStats <- read.table("/mnt/storage1/work/amjad/ctdna/result_ctdnaAlignment/stats/out.csv", header = T, stringsAsFactors = F, sep = "\t") %>%
  mutate(category = "Pilot", sampleType = case_when(
    grepl("FFPE", Sample) ~ "FFPE",
    grepl("ctDNA_",   Sample) ~ "ctDNA",
    grepl("WB",   Sample) ~ "WB"))

combined <- list(pilotStats, statsNovaseqNotrim, statsNovaseqTrim, statsNovaseqTile) %>%
             map_dfr( ~ select(.x, Mismatch_Rate = PF_MISMATCH_RATE, MEAN_TARGET_COVERAGE, PF_BASES_ALIGNED, 
                               PCT_OFF_BAIT, category, sampleType)) %>%
             mutate(category = factor(category, levels = c("Pilot","NOVAseq-No-QC","NOVAseq-Trim","NOVAseq-TileFiltering")))

coverage <- combined %>% 
      ggplot(aes(x = category, y = MEAN_TARGET_COVERAGE)) + 
             geom_boxplot(outlier.color = NA)  + 
             ggbeeswarm::geom_quasirandom(aes(color = sampleType))

coverageNoBlood <- combined %>% 
  filter(sampleType != "WB") %>%
  ggplot(aes(x = category, y = MEAN_TARGET_COVERAGE)) + 
  geom_boxplot(outlier.color = NA)  + 
  ggbeeswarm::geom_quasirandom(aes(color = sampleType))

offTarget <- combined %>% 
  ggplot(aes(x = category, y = PCT_OFF_BAIT)) + 
  geom_boxplot(outlier.color = NA)  + 
  ggbeeswarm::geom_quasirandom(aes(color = sampleType))

mismatches <- combined %>% 
  ggplot(aes(x = category, y = Mismatch_Rate)) + 
  geom_boxplot(outlier.color = NA)  + 
  ggbeeswarm::geom_quasirandom(aes(color = sampleType))


dupStatsNova <- read.table("/mnt/storage2/work/amjad/ctdna/result_ctdnaAlignmentNovaseq/dupStats/out.csv", header = T, stringsAsFactors = F, sep ="\t") %>%
mutate(category = "NOVAseq-TileFiltering", sampleType = case_when(
  grepl("FFPE", Sample) ~ "FFPE",
  grepl("ctDNA_",   Sample) ~ "ctDNA",
  grepl("WB",   Sample) ~ "WB"))
dupStatsPilot<- read.table("/mnt/storage1/work/amjad/ctdna/result_ctdnaAlignment/dupStats/out.csv", header = T, stringsAsFactors = F, sep ="\t") %>%
  mutate(category = "Pilot", sampleType = case_when(
    grepl("FFPE", Sample) ~ "FFPE",
    grepl("ctDNA_",   Sample) ~ "ctDNA",
    grepl("WB",   Sample) ~ "WB"))

combinedDup <- rbind(dupStatsNova, dupStatsPilot) %>% 
  mutate(category = factor(category, levels = c("Pilot","NOVAseq-TileFiltering")))


dups <- combinedDup %>% 
  filter(sampleType != "WB") %>%
  ggplot(aes(x = category, y = PERCENT_DUPLICATION)) + 
  geom_boxplot(outlier.color = NA)  + 
  ggbeeswarm::geom_quasirandom(aes(color = sampleType))