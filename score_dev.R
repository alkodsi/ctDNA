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

getReadCounts <- function(chr, pos, base, bam, tag = "", min_base_quality = 20, max_depth = 100000, include_indels = F, min_mapq = 30){
  require(Rsamtools)
  gr <- GRanges(chr, IRanges(pos, pos))
  if(tag == ""){
    sbp <- ScanBamParam(which = gr)
  } else {
    sbp <- ScanBamParam(which = gr, tagFilter = list("RG" = tag))
  }
  pileupParam = PileupParam(max_depth = max_depth, min_base_quality = min_base_quality,
                            min_mapq = min_mapq, distinguish_strands = F, 
                            include_deletions = include_indels, include_insertions = include_indels)

  p <- pileup(bam, scanBamParam = sbp, pileupParam = pileupParam)
  depth <- ifelse(!base %in% p$nucleotide, 0, p[p$nucleotide == base, "count"])
  return(depth)
}

getBackgroundRate <- function(bam, targets, reference, vafThreshold = 0.1, tag = "", min_base_quality = 20, max_depth = 100000, include_indels = F, min_mapq = 30){
  require(Rsamtools)
  gr <- GRanges(targets$chr, IRanges(targets$start, targets$end))
  if(tag == ""){
    sbp <- ScanBamParam(which = gr)
  } else {
    sbp <- ScanBamParam(which = gr, tagFilter = list("RG" = tag))
  }
  pileupParam <- PileupParam(max_depth = max_depth, min_base_quality = min_base_quality,
                            min_mapq = min_mapq, distinguish_strands = F, 
                            include_deletions = include_indels, include_insertions = include_indels)

  p <- pileup(bam, scanBamParam = sbp, pileupParam = pileupParam) %>%
          spread(-nucleotide, count, fill = 0) %>%
          mutate(ref = as.character(getSeq(reference, GRanges(seqnames, IRanges(pos, pos))))) %>%
          mutate(depth = A + C + G + T)
  
  pAnn <- p %>% mutate(refCount = map2_dbl(c(1:nrow(p)), p$ref, ~ p[.x, .y]),
                       nonRefCount = depth - refCount) %>%
               filter(nonRefCount/depth < vafThreshold)
  rate <- sum(as.numeric(pAnn$nonRefCount))/sum(as.numeric(pAnn$depth))
  return(rate)
}

testSample <- function(mutations, bam, targets, reference, tag = "", vafThreshold = 0.1, 
                       min_base_quality = 20, max_depth = 10000, include_indels = F, min_mapq = 30,
                       nPermutation = 10000, seed = 123){
  cat("Estimating background rate ...\n")
  bg <- getBackgroundRate(bam = bam, targets = targets, reference = reference, tag = tag,
                          vafThreshold = vafThreshold, min_base_quality = min_base_quality, max_depth = max_depth, 
                          include_indels = include_indels, min_mapq = min_mapq)
  cat(paste("Background rate is", bg, "\n"))
  cat("Getting ref and alt Counts \n")
  altReads <- pmap_dbl(list(mutations$CHROM, mutations$POS, mutations$ALT),
                       function(x, y, z) getReadCounts(chr = x, pos = y, base = z, bam = bam,
                                                       tag = tag, min_base_quality = min_base_quality,
                                                       min_mapq = min_mapq, max_depth = max_depth, include_indels = include_indels))
  
  refReads <- pmap_dbl(list(mutations$CHROM, mutations$POS, mutations$REF),
                       function(x, y, z) getReadCounts(chr = x, pos = y, base = z, bam = bam,
                                                       tag = tag, min_base_quality = min_base_quality,
                                                       min_mapq = min_mapq, max_depth = max_depth, include_indels = include_indels))
  refAlt <- data.frame(Ref = refReads, Alt = altReads)
  cat("Running permutation test \n")
  posTest <- positivityTest(depths = refReads + altReads, altReads = altReads, rate = bg/3, seed = seed, nPermutation = nPermutation)
  cat(paste("Pvalue = ", posTest, "\n"))
  cat(paste("Sample is ctDNA", ifelse(posTest < 0.05, "positive\n", "negative\n")))
  return(list(counts = refAlt, pvalue = posTest))
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
  filter(ECNT <= 10)

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