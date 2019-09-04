library(tidyverse)
y <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutect2_2/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")
x <- read.table("/mnt/storage2/work/amjad/ctdna/VariantsVersions/variants_seq2.2_v1.csv",header=T,stringsAsFactors = F, sep = "\t")
back <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutect2_2/backgroundAll/csvOut.csv", header = T, stringsAsFactors = F, sep = "\t")

backgroundRates <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll/backgroundAll/csvOut.csv", header = T, stringsAsFactors = F, sep = "\t")
x <- read.table("/mnt/storage2/work/amjad/ctdna/result_newMutectAll/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")
mutations <- read.table("/mnt/storage2/work/amjad/ctdna/code/ctdna/result_mutectDoubleCalling/allVarsFixedFilIndels/out.csv",header=T,stringsAsFactors = F, sep = "\t")


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

simulator <- function(nVariants, depths, rate, altReads, seed){
  set.seed(seed)
  sim <- rbinom(n = nVariants, size = depths, prob = rate)
  comparison <- map_dbl(c(0:max(altReads)), ~as.numeric(sum(sim>=.x) >= sum(altReads>=.x)))
  out <- ifelse(sum(comparison) == length(comparison), 1, 0)
  return(out)
}

positivityTest <- function(depths, altReads, rate, seed, nPermutation){
  set.seed(seed)
  seeds <- round(runif(n = nPermutation, min = 0, max = 100000000))
  pvalue <- sum(map_dbl(seeds, ~simulator(length(depths), depths, rate,  altReads, seed = .x)))/nPermutation
  return(pvalue)
}

getReadsHoldingMutation <- function(chr, pos, alt, bam, tag = ""){
  require(GenomicAlignments)
  gr <- GRanges(chr, IRanges(pos, pos))
  if(tag == ""){
    stackedStrings <- stackStringsFromBam(bam, use.names=T, param = gr)
  } else {
    stackedStrings <- stackStringsFromBam(bam, use.names=T, 
                          param = ScanBamParam(tagFilter = list("RG" = tag), which = gr)) 
  }
  out <- data.frame(ID = names(stackedStrings), seq = as.data.frame(stackedStrings)[,1])
  return(as.character(out[out$seq == alt, "ID"]))
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
  filter(grepl("ctDNA",Sample) & !grepl("_0", Sample) & !grepl("Normal",Sample) & !grepl("CHIC_143",Sample)) %>%
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
  filter(nchar(REF) == 1 & nchar(ALT) == 1)

results <- data.frame(Sample = samplesToTest$Sample,
                      Pvalue = pmap_dbl(list(samplesToTest$patient, samplesToTest$series, samplesToTest$Rate),
                                        function(x,y,z) positivityTest(depths = mutations[mutations$Patient == x, paste0(y,".DP")],
                                               rate = z/3,  
                                               altReads = mutations[mutations$Patient == x, paste0(y,".AD_AltCount")],
                                               seed = 123, nPermutation = 10000))) %>%
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

binom.test(6,800,0.0025, alternative = "greater")$p.value