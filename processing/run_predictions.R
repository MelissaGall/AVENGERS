#Load libraries
libsToLoad <- c("dplyr", "plyr", "tidyr", "mclust", "openxlsx", "matrixStats", "ggplot2")
lapply(libsToLoad, library, character.only = TRUE)
source("functions.R")

#Load files

##Generate all files list
expDir <- "../data/Seqprep/merged/no_Ns/sam/variant_counts_no_thresh/annotated/"
editsData <- read.table("../reference/BRCA2_editing_data.txt", header = TRUE, fill = TRUE)
filesList <- list.files(expDir, pattern = ".xlsx")

##Generate amplicon list + filter
ampList <- unique(sub("R.*", "", filesList))

##Version
version <- ("_V1")

#Calculate the function score for all experiments to get the list of variants to be filtered

##Define the empty objects
excludedReads <- data.frame()
processedData <- data.frame()

##Filter variants: get the list of variant with discrepancies
for (amplicon in ampList){
  print(paste("Working on", amplicon))
  ###Get data for each amplicon
  edit5 <- editsData$HDR_5_SNV[editsData$amplicon == amplicon]
  edit3 <- editsData$HDR_3_SNV[editsData$amplicon == amplicon]
  mutStart <- editsData$mut_start[editsData$amplicon == amplicon]
  mutStop <- editsData$mut_end[editsData$amplicon == amplicon]
  exonStart <-  editsData$exon_start[editsData$amplicon == amplicon] - editsData$genomic_pos[editsData$amplicon == amplicon] + 1 
  exonStop <- editsData$exon_end[editsData$amplicon == amplicon] - editsData$genomic_pos[editsData$amplicon == amplicon] + 1
  posToKeep <- seq(mutStart, mutStop)
  freqThr = 0.00001 #Minimum threshold for reads frequency
  ###Load all files
  CISr1raw <- read.xlsx(paste(expDir, amplicon, "RP1CIS_annotated.xlsx", sep = ""), fill = TRUE)
  CISr1 <- createVarSimpDF(CISr1raw, edit3, edit5, posToKeep, exonStart, exonStop, FALSE)
  CISr2raw <- read.xlsx(paste(expDir, amplicon, "RP2CIS_annotated.xlsx", sep = ""), fill = TRUE)
  CISr2 <- createVarSimpDF(CISr2raw, edit3, edit5, posToKeep, exonStart, exonStop, FALSE)
  DMSOr1raw <- read.xlsx(paste(expDir, amplicon, "RP1DMSO_annotated.xlsx", sep = ""), fill = TRUE)
  DMSOr1 <- createVarSimpDF(DMSOr1raw, edit3, edit5, posToKeep, exonStart, exonStop, FALSE)
  DMSOr2raw <- read.xlsx(paste(expDir, amplicon, "RP2DMSO_annotated.xlsx", sep = ""), fill = TRUE)
  DMSOr2 <- createVarSimpDF(DMSOr2raw, edit3, edit5, posToKeep, exonStart, exonStop, FALSE)
  OLAr1raw <- read.xlsx(paste(expDir, amplicon, "RP1OLA_annotated.xlsx", sep = ""), fill = TRUE)
  OLAr1 <- createVarSimpDF(OLAr1raw, edit3, edit5, posToKeep, exonStart, exonStop, FALSE)
  OLAr2raw <- read.xlsx(paste(expDir, amplicon, "RP2OLA_annotated.xlsx", sep = ""), fill = TRUE)
  OLAr2 <- createVarSimpDF(OLAr2raw, edit3, edit5, posToKeep, exonStart, exonStop, FALSE)
  ###Get variant excluded due to read counts too low
  replicatesData <- merge(DMSOr1raw[, c("variant", "edit_string", "codon", "pre", "pre_freq", "post", "post_freq", "pos")], 
                          OLAr2raw[, c("variant", "edit_string", "codon", "pre", "pre_freq", "post", "post_freq")], 
                          by = c("variant", "codon", "edit_string"), 
                          suffixes = c(".RP1", ".RP2"),
                          all = T)
  replicatesDataExcluded <- replicatesData[(replicatesData$pre.RP1 < 10 | replicatesData$pre.RP2 < 10 | 
                                              replicatesData$pre_freq.RP1 < freqThr | replicatesData$pre_freq.RP2 < freqThr),]
  replicatesDataExcluded <- replicatesDataExcluded[!(is.na(replicatesDataExcluded$codon)),]
  replicatesDataExcluded <- replicatesDataExcluded[replicatesDataExcluded$pos %in% posToKeep,]
  if (nrow(replicatesDataExcluded) != 0){
    replicatesDataExcluded <- data.frame("exon" = amplicon, replicatesDataExcluded)
    excludedReads <- rbind(excludedReads, replicatesDataExcluded)
  }
  ###Merge data for the 3 conditions
  mergingColumns <- c("variant", "edit_string", "pos", "codon", "nuc_change", "AA_change", "c_nom",
                      "g_nom", "Affect_splicing", "dbSNP.ID", "clinical_sign", "Protein.change", "expected_clinical")
  mergedDMSO <- merge(DMSOr1, DMSOr2, by = mergingColumns, suffixes = c(".dmso.1", ".dmso.2"), all = F)
  mergedDMSO$post.pre.ratio.dmso <- rowMeans(mergedDMSO[, c('post.pre.ratio.sns.dmso.1', 'post.pre.ratio.sns.dmso.2')], na.rm = F)
  mergedCIS <- merge(CISr1, CISr2, by = mergingColumns, suffixes = c(".cis.1", ".cis.2"), all = F)
  mergedCIS$post.pre.ratio.cis <- rowMeans(mergedCIS[, c('post.pre.ratio.sns.cis.1', 'post.pre.ratio.sns.cis.2')], na.rm = F)
  mergedOla <- merge(OLAr1, OLAr2, by =mergingColumns, suffixes = c(".ola.1", ".ola.2"), all = F)
  mergedOla$post.pre.ratio.ola <- rowMeans(mergedOla[, c('post.pre.ratio.sns.ola.1', 'post.pre.ratio.sns.ola.2')], na.rm = F)
  mergedAll <- merge(mergedDMSO, mergedCIS, by = mergingColumns, all = TRUE)
  mergedAll <- merge(mergedAll, mergedOla, by = mergingColumns, all = TRUE)
  mergedAll <- mergedAll[!(mergedAll$variant %in% excludedReads$variant),]
  mergedAll <- mergedAll[!(is.na(mergedAll$pos)),]
  mergedAll <- mergedAll[order(mergedAll$pos),]
  mergedAll$post.pre.ratio.mean <- apply(mergedAll[c("post.pre.ratio.dmso",
                                                     "post.pre.ratio.cis",
                                                     "post.pre.ratio.ola")],
                                         1, weighted.mean, 
                                         w = c(0.4, 0.25, 0.25), 
                                         na.rm = TRUE)
  mergedAll <- removeOutliers(mergedAll, "post.pre.ratio.mean")
  ###Calculate position bias for all conditions
  mergedAll$post.pre.loess <- removePositionBiais(mergedAll, "post.pre.ratio.mean")
  mergedAll$post.pre.loess.dmso <- removePositionBiais(mergedAll, "post.pre.ratio.dmso")
  mergedAll$post.pre.loess.dmso.rp1 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.dmso.1")
  mergedAll$post.pre.loess.dmso.rp2 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.dmso.2")
  mergedAll$post.pre.loess.ola <- removePositionBiais(mergedAll, "post.pre.ratio.ola")
  mergedAll$post.pre.loess.ola.rp1 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.ola.1")
  mergedAll$post.pre.loess.ola.rp2 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.ola.2")
  mergedAll$post.pre.loess.cis <- removePositionBiais(mergedAll, "post.pre.ratio.cis")
  mergedAll$post.pre.loess.cis.rp1 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.cis.1")
  mergedAll$post.pre.loess.cis.rp2 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.cis.2")
  ###Concatenate data
  mergedAll$exon <- amplicon
  mergedAll$intronic <- ifelse(mergedAll$codon == "Intronic", "yes", "no")
  processedData <- rbind(processedData, mergedAll)
}

##Per exon normalisation
processedData$function.score <- perExonNorm(processedData, ampList, "post.pre.loess")
processedData$function.score.dmso <- perExonNorm(processedData, ampList, "post.pre.loess.dmso")
processedData$function.score.dmso.rp1 <- perExonNorm(processedData, ampList, "post.pre.loess.dmso.rp1")
processedData$function.score.dmso.rp2 <- perExonNorm(processedData, ampList, "post.pre.loess.dmso.rp2")
processedData$function.score.ola <- perExonNorm(processedData, ampList, "post.pre.loess.ola")
processedData$function.score.ola.rp1 <- perExonNorm(processedData, ampList, "post.pre.loess.ola.rp1")
processedData$function.score.ola.rp2 <- perExonNorm(processedData, ampList, "post.pre.loess.ola.rp2")
processedData$function.score.cis <- perExonNorm(processedData, ampList, "post.pre.loess.cis")
processedData$function.score.cis.rp1 <- perExonNorm(processedData, ampList, "post.pre.loess.cis.rp1")
processedData$function.score.cis.rp2 <- perExonNorm(processedData, ampList, "post.pre.loess.cis.rp2")

##Global norm
processedData$function.score <- globalNorm(processedData, ampList, "function.score")
processedData$function.score.dmso <- globalNorm(processedData, ampList, "function.score.dmso")
processedData$function.score.dmso.rp1 <- globalNorm(processedData, ampList, "function.score.dmso.rp1")
processedData$function.score.dmso.rp2 <- globalNorm(processedData, ampList, "function.score.dmso.rp2")
processedData$function.score.ola <- globalNorm(processedData, ampList, "function.score.ola")
processedData$function.score.ola.rp1 <- globalNorm(processedData, ampList, "function.score.ola.rp1")
processedData$function.score.ola.rp2 <- globalNorm(processedData, ampList, "function.score.ola.rp2")
processedData$function.score.cis <- globalNorm(processedData, ampList, "function.score.cis")
processedData$function.score.cis.rp1 <- globalNorm(processedData, ampList, "function.score.cis.rp1")
processedData$function.score.cis.rp2 <- globalNorm(processedData, ampList, "function.score.cis.rp2")

##Create training data frame (remove unknown clinical outcome and intronic)
trainingData <- processedData[processedData$expected_clinical != "Unknown",]
trainingData <- trainingData[!(is.na(trainingData$expected_clinical)),]
trainingData <- trainingData[trainingData$codon != "Intronic",]
trainingData$lof <- sapply(trainingData$expected_clinical, assingLOF)

##Run the prediction
processedData <- runMclust(trainingData, processedData, "function.score", "proba")
processedData <- runMclust(trainingData, processedData, "function.score.dmso", "proba.dmso")
processedData <- runMclust(trainingData, processedData, "function.score.dmso.rp1", "proba.dmso.rp1")
processedData <- runMclust(trainingData, processedData, "function.score.dmso.rp2", "proba.dmso.rp2")
processedData <- runMclust(trainingData, processedData, "function.score.ola", "proba.ola")
processedData <- runMclust(trainingData, processedData, "function.score.ola.rp1", "proba.ola.rp1")
processedData <- runMclust(trainingData, processedData, "function.score.ola.rp2", "proba.ola.rp2")
processedData <- runMclust(trainingData, processedData, "function.score.cis", "proba.cis")
processedData <- runMclust(trainingData, processedData, "function.score.cis.rp1", "proba.cis.rp1")
processedData <- runMclust(trainingData, processedData, "function.score.cis.rp2", "proba.cis.rp2")

##Annotate class
processedData <- addAnno(processedData, "proba", "class", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.dmso", "class.dmso", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.dmso.rp1", "class.dmso.rp1", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.dmso.rp2", "class.dmso.rp2", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.ola", "class.ola", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.ola.rp1", "class.ola.rp1", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.ola.rp2", "class.ola.rp2", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.cis", "class.cis", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.cis.rp1", "class.cis.rp1", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.cis.rp2", "class.cis.rp2", 0.95, 0.01)

##Get difference between condition/replicates
processedData <- processedData[!(is.na(processedData$variant)),]
processedData$dmso.diff <- abs(processedData$post.pre.ratio.sns.dmso.1 - processedData$post.pre.ratio.sns.dmso.2)
processedData$cis.diff <- abs(processedData$post.pre.ratio.sns.cis.1 - processedData$post.pre.ratio.sns.cis.2)
processedData$ola.diff <- abs(processedData$post.pre.ratio.sns.ola.1 - processedData$post.pre.ratio.sns.ola.2)
dmsoSD <- sd(processedData$dmso.diff, na.rm = T)
cisSD <- sd(processedData$cis.diff, na.rm = T)
olaSD <- sd(processedData$ola.diff, na.rm = T)

##Get list of variant to exclude
excludedVarDMSO <- processedData$variant[(
  (processedData$class.dmso.rp1 == "Non-functional" & processedData$class.dmso.rp2 == "Functional") |
    (processedData$class.dmso.rp1 == "Functional" & processedData$class.dmso.rp2 == "Non-functional") 
) & processedData$dmso.diff > 2 * dmsoSD]
excludedVarCIS <- processedData$variant[(
  (processedData$class.cis.rp1 == "Non-functional" & processedData$class.cis.rp2 == "Functional") |
    (processedData$class.cis.rp1 == "Functional" & processedData$class.cis.rp2 == "Non-functional") 
) & processedData$cis.diff > 2 * cisSD]
excludedVarOLA <- processedData$variant[(
  (processedData$class.ola.rp1 == "Non-functional" & processedData$class.ola.rp2 == "Functional") |
    (processedData$class.ola.rp1 == "Functional" & processedData$class.ola.rp2 == "Non-functional") 
) & processedData$ola.diff > 2 * olaSD]

allExcludedVar <- c(excludedVarDMSO, excludedVarCIS, excludedVarOLA)
allExcludedVar2Cond <- names(Filter(function(x) x >= 2, table(allExcludedVar)))
allExcludedVar2Cond <- unique(allExcludedVar2Cond[!(allExcludedVar2Cond %in% excludedVarDMSO)])
excludedVar <- c(excludedReads$variant, excludedVarDMSO, allExcludedVar2Cond)

##Save excluded variants
processedData$c_nom <- gsub("&gt;", ">", processedData$c_nom)
processedData$g_nom <- gsub("&gt;", ">", processedData$g_nom)
processedData <- processedData %>%  dplyr::select(exon, everything())
excludedDMSO <- processedData[processedData$variant %in% excludedVarDMSO,]
excludedDMSO <- excludedDMSO[!(excludedDMSO$variant %in% excludedReads$variant),]
allExcludedVarDiff <- processedData[processedData$variant %in% allExcludedVar2Cond,]
if(nrow(excludedReads) == 0){
  agg1 <- data.frame("exon" = ampList, "excluded.read.count" = 0)
} else {
  agg1 <- aggregate(excludedReads$exon, by = list(excludedReads$exon), FUN = length)
  colnames(agg1) <- c("exon", "excluded.read.count")  
}
if(nrow(excludedDMSO) == 0){
  agg2 <- data.frame("exon" = ampList, "excludedDMSO.diff" = 0)
} else {
  agg2 <- aggregate(excludedDMSO$exon, by=list(excludedDMSO$exon), FUN=length)
  colnames(agg2) <- c("exon", "excludedDMSO.diff") 
}
if(nrow(allExcludedVarDiff) == 0){
  agg3 <- data.frame("exon" = ampList, "excluded.diff.more.than.2.conditions" = 0)
} else {
  agg3 <- aggregate(allExcludedVarDiff$exon, by = list(allExcludedVarDiff$exon), FUN = length)
  colnames(agg3) <- c("exon", "excluded.diff.more.than.2.conditions") 
} 
excludedSummary <- data.frame("exon" = ampList)
excludedSummary <- merge(excludedSummary, agg1, by = "exon", all = T)
excludedSummary <- merge(excludedSummary, agg2, by = "exon", all = T)
excludedSummary <- merge(excludedSummary, agg3, by = "exon", all = T)
excludedSummary[is.na(excludedSummary)] <- 0
excludedSummary$total <- excludedSummary$excluded.read.count + excludedSummary$excludedDMSO.diff + excludedSummary$excluded.diff.more.than.2.conditions
expectedNum <- read.xlsx("../reference/expected_var_num.xlsx")
excludedSummary <- merge(excludedSummary, expectedNum, by = "exon", all.x = T)
excludedSummary$expected_variants_number <- as.numeric(excludedSummary$expected_variants_number)
processedDataFilt <- processedData[!(processedData$variant %in% excludedVar),]
exonCount <- dplyr::count(processedDataFilt, exon)
excludedSummary$recovered.number <- exonCount$n
excludedSummary$percentage.recovered <- (excludedSummary$recovered.number / excludedSummary$expected_variants_number) * 100
excludedTables <- list("Summary" = excludedSummary, 
                       "Read.counts" = excludedReads, 
                       "DMSO.difference" = excludedDMSO,
                       "Difference.in.two.conditions" = allExcludedVarDiff)
##Save table of filtered variants
write.xlsx(excludedTables, file = paste0("../results/excluded_variants_", version, ".xlsx"))

#Run the model, with the variants filtered

##Define the empty objects
processedData <- data.frame()
trainingData <- data.frame()

##Run the pipeline by removing excluded variants
for (amplicon in ampList){
  print(paste("Working on", amplicon))
  ###Get data for each amplicon
  edit5 <- editsData$HDR_5_SNV[editsData$amplicon == amplicon]
  edit3 <- editsData$HDR_3_SNV[editsData$amplicon == amplicon]
  mutStart <- editsData$mut_start[editsData$amplicon == amplicon]
  mutStop <- editsData$mut_end[editsData$amplicon == amplicon]
  exonStart <-  editsData$exon_start[editsData$amplicon == amplicon] - editsData$genomic_pos[editsData$amplicon == amplicon] + 1 
  exonStop <- editsData$exon_end[editsData$amplicon == amplicon] - editsData$genomic_pos[editsData$amplicon == amplicon] + 1
  posToKeep <- seq(mutStart, mutStop)
  ###Load all files
  CISr1raw <- read.xlsx(paste(expDir, amplicon, "RP1CIS_annotated.xlsx", sep = ""), fill = TRUE)
  CISr1 <- createVarSimpDF(CISr1raw, edit3, edit5, posToKeep, exonStart, exonStop, TRUE)
  CISr2raw <- read.xlsx(paste(expDir, amplicon, "RP2CIS_annotated.xlsx", sep = ""), fill = TRUE)
  CISr2 <- createVarSimpDF(CISr2raw, edit3, edit5, posToKeep, exonStart, exonStop, TRUE)
  DMSOr1raw <- read.xlsx(paste(expDir, amplicon, "RP1DMSO_annotated.xlsx", sep = ""), fill = TRUE)
  DMSOr1 <- createVarSimpDF(DMSOr1raw, edit3, edit5, posToKeep, exonStart, exonStop, TRUE)
  DMSOr2raw <- read.xlsx(paste(expDir, amplicon, "RP2DMSO_annotated.xlsx", sep = ""), fill = TRUE)
  DMSOr2 <- createVarSimpDF(DMSOr2raw, edit3, edit5, posToKeep, exonStart, exonStop, TRUE)
  OLAr1 <- read.xlsx(paste(expDir, amplicon, "RP1OLA_annotated.xlsx", sep = ""), fill = TRUE)
  OLAr1 <- createVarSimpDF(OLAr1, edit3, edit5, posToKeep, exonStart, exonStop, TRUE)
  OLAr2 <- read.xlsx(paste(expDir, amplicon, "RP2OLA_annotated.xlsx", sep = ""), fill = TRUE)
  OLAr2 <- createVarSimpDF(OLAr2, edit3, edit5, posToKeep, exonStart, exonStop, TRUE)
  ###Merge data for the 3 conditions
  mergingColumns <- c("variant", "edit_string", "pos", "codon", "nuc_change", "AA_change", "c_nom", "g_nom", "dbSNP.ID", "clinical_sign", "Protein.change", "expected_clinical", "Affect_splicing")
  mergedDMSO <- merge(DMSOr1, DMSOr2, by = mergingColumns, suffixes = c(".dmso.1", ".dmso.2"), all = F)
  mergedDMSO$post.pre.ratio.dmso <- rowMeans(mergedDMSO[, c('post.pre.ratio.sns.dmso.1', 'post.pre.ratio.sns.dmso.2')], na.rm = F)
  mergedCIS <- merge(CISr1, CISr2, by = mergingColumns, suffixes = c(".cis.1", ".cis.2"), all = F)
  mergedCIS$post.pre.ratio.cis <- rowMeans(mergedCIS[, c('post.pre.ratio.sns.cis.1', 'post.pre.ratio.sns.cis.2')], na.rm = F)
  mergedOla <- merge(OLAr1, OLAr2, by =mergingColumns, suffixes = c(".ola.1", ".ola.2"), all = F)
  mergedOla$post.pre.ratio.ola <- rowMeans(mergedOla[, c('post.pre.ratio.sns.ola.1', 'post.pre.ratio.sns.ola.2')], na.rm = F)
  mergedAll <- merge(mergedDMSO, mergedCIS, by = mergingColumns, all = TRUE)
  mergedAll <- merge(mergedAll, mergedOla, by = mergingColumns, all = TRUE)
  mergedAll <- mergedAll[!(is.na(mergedAll$pos)),]
  mergedAll <- mergedAll[order(mergedAll$pos),]
  ###Calculate post pre ratio
  mergedAll$post.pre.ratio.mean <- apply(mergedAll[c("post.pre.ratio.dmso", 
                                                     "post.pre.ratio.cis", 
                                                     "post.pre.ratio.ola")],
                                         1, weighted.mean, 
                                         w = c(0.4, 0.25, 0.25), 
                                         na.rm = TRUE)
  mergedAll$post.pre.ratio.mean.rp1 <- apply(mergedAll[c("post.pre.ratio.sns.dmso.1", 
                                                         "post.pre.ratio.sns.cis.1", 
                                                         "post.pre.ratio.sns.ola.1")],
                                             1, weighted.mean, 
                                             w = c(0.4, 0.25, 0.25), 
                                             na.rm = TRUE)
  mergedAll$post.pre.ratio.mean.rp2 <- apply(mergedAll[c("post.pre.ratio.sns.dmso.2", 
                                                         "post.pre.ratio.sns.cis.2", 
                                                         "post.pre.ratio.sns.ola.2")],
                                             1, weighted.mean, 
                                             w = c(0.4, 0.25, 0.25), 
                                             na.rm = TRUE)
  mergedAll <- removeOutliers(mergedAll, "post.pre.ratio.mean")
  ###Calculate position bias for all conditions
  mergedAll$post.pre.loess <- removePositionBiais(mergedAll, "post.pre.ratio.mean")
  mergedAll$post.pre.loess.rp1 <- removePositionBiais(mergedAll, "post.pre.ratio.mean.rp1")
  mergedAll$post.pre.loess.rp2 <- removePositionBiais(mergedAll, "post.pre.ratio.mean.rp2")
  mergedAll$post.pre.loess.dmso <- removePositionBiais(mergedAll, "post.pre.ratio.dmso")
  mergedAll$post.pre.loess.dmso.rp1 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.dmso.1")
  mergedAll$post.pre.loess.dmso.rp2 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.dmso.2")
  mergedAll$post.pre.loess.ola <- removePositionBiais(mergedAll, "post.pre.ratio.ola")
  mergedAll$post.pre.loess.ola.rp1 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.ola.1")
  mergedAll$post.pre.loess.ola.rp2 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.ola.2")
  mergedAll$post.pre.loess.cis <- removePositionBiais(mergedAll, "post.pre.ratio.cis")
  mergedAll$post.pre.loess.cis.rp1 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.cis.1")
  mergedAll$post.pre.loess.cis.rp2 <- removePositionBiais(mergedAll, "post.pre.ratio.sns.cis.2")
  mergedAll$exon <- amplicon
  mergedAll$intronic <- ifelse(mergedAll$codon == "Intronic", "yes", "no")
  processedData <- rbind(processedData, mergedAll)
}

##Save just in case
processedDataBeforeNorm <- processedData

##Normalize per exon
processedData$function.score <- perExonNorm(processedData, ampList, "post.pre.loess")
processedData$function.score.rp1 <- perExonNorm(processedData, ampList, "post.pre.loess.rp1")
processedData$function.score.rp2 <- perExonNorm(processedData, ampList, "post.pre.loess.rp2")
processedData$function.score.dmso <- perExonNorm(processedData, ampList, "post.pre.loess.dmso")
processedData$function.score.dmso.rp1 <- perExonNorm(processedData, ampList, "post.pre.loess.dmso.rp1")
processedData$function.score.dmso.rp2 <- perExonNorm(processedData, ampList, "post.pre.loess.dmso.rp2")
processedData$function.score.ola <- perExonNorm(processedData, ampList, "post.pre.loess.ola")
processedData$function.score.ola.rp1 <- perExonNorm(processedData, ampList, "post.pre.loess.ola.rp1")
processedData$function.score.ola.rp2 <- perExonNorm(processedData, ampList, "post.pre.loess.ola.rp2")
processedData$function.score.cis <- perExonNorm(processedData, ampList, "post.pre.loess.cis")
processedData$function.score.cis.rp1 <- perExonNorm(processedData, ampList, "post.pre.loess.cis.rp1")
processedData$function.score.cis.rp2 <- perExonNorm(processedData, ampList, "post.pre.loess.cis.rp2")

##Global norm
processedData$function.score <- globalNorm(processedData, ampList, "function.score")
processedData$function.score.rp1 <- globalNorm(processedData, ampList, "function.score.rp1")
processedData$function.score.rp2 <- globalNorm(processedData, ampList, "function.score.rp2")
processedData$function.score.dmso <- globalNorm(processedData, ampList, "function.score.dmso")
processedData$function.score.dmso.rp1 <- globalNorm(processedData, ampList, "function.score.dmso.rp1")
processedData$function.score.dmso.rp2 <- globalNorm(processedData, ampList, "function.score.dmso.rp2")
processedData$function.score.ola <- globalNorm(processedData, ampList, "function.score.ola")
processedData$function.score.ola.rp1 <- globalNorm(processedData, ampList, "function.score.ola.rp1")
processedData$function.score.ola.rp2 <- globalNorm(processedData, ampList, "function.score.ola.rp2")
processedData$function.score.cis <- globalNorm(processedData, ampList, "function.score.cis")
processedData$function.score.cis.rp1 <- globalNorm(processedData, ampList, "function.score.cis.rp1")
processedData$function.score.cis.rp2 <- globalNorm(processedData, ampList, "function.score.cis.rp2")

##Get training data
trainingData <- processedData[processedData$expected_clinical != "Unknown",]
trainingData <- trainingData[!(trainingData$codon == "Intronic"),]
trainingData <- trainingData[!(is.na(trainingData$expected_clinical)),]
trainingData$lof <- sapply(trainingData$expected_clinical, assingLOF)

##Run model
processedData <- runMclust(trainingData, processedData, "function.score", "proba")
processedData <- runMclust(trainingData, processedData, "function.score.dmso", "proba.dmso")
processedData <- runMclust(trainingData, processedData, "function.score.cis", "proba.cis")
processedData <- runMclust(trainingData, processedData, "function.score.ola", "proba.ola")
processedData <- addAnno(processedData, "proba", "class", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.dmso", "class.dmso", 0.95, 0.01)
processedData <- addAnno(processedData, "proba.ola", "class.ola", 0.95,  0.01)
processedData <- addAnno(processedData, "proba.cis", "class.cis", 0.95,  0.01)

##Merge with spliceAI results and save
processedData$c_nom <- gsub("&gt;", ">", processedData$c_nom)
processedData$g_nom <- gsub("&gt;", ">", processedData$g_nom)
spliceAI <- read.xlsx("../reference/splicing/all_exons_summary_spliceAI.xlsx")
processedData <- merge(processedData, spliceAI[, c("Exon", "g.nom")], by.x = c("exon", "g_nom"), by.y = c("Exon", "g.nom"), all.x = T)

#Classification

##If DMSO different from global, consider as "Uncertain"
processedData$processedData.class <- ifelse(processedData$class != processedData$class.dmso, "Uncertain", processedData$class.dmso)
processedData$processedData.FS <- processedData$function.score.dmso
processedData$processedData.proba <- processedData$proba.dmso
processedData$processedData.class <- ifelse(processedData$class.dmso == "Intermediate", "Uncertain", processedData$processedData.class)

#Stats

##Stats dataframe
processedDataStats <- processedData[!(processedData$Affect_splicing == "Yes" & processedData$codon == "Synonymous"),]
processedDataStats <- processedDataStats[processedDataStats$processedData.class != "Uncertain",]

##Calculate OddsPath
benign <- processedDataStats[processedDataStats$clinical_sign == "Benign",]
path <- processedDataStats[processedDataStats$clinical_sign == "Pathogenic",]
benign <- benign[!(is.na(benign$exon)),]
path <- path[!(is.na(path$exon)),]

num.benign.ctrl <- nrow(benign)
num.path.ctrl <- nrow(path)
num.benign.pred <- nrow(benign[benign$processedData.class == "Functional",])
num.path.pred <- nrow(path[path$processedData.class == "Non-functional",])
P1 = num.path.ctrl / (num.benign.ctrl + num.path.ctrl)
P2path = num.path.pred / (num.path.pred + 1)
P2benign = 1 / (num.benign.pred  + 1)
OP.path = (P2path * (1 - P1)) / ((1 - P2path) * P1)
OP.benign = (P2benign * (1 - P1)) / ((1 - P2benign) * P1)
print(paste0("OP.path: ", OP.path))
print(paste0("OP.benign: ", OP.benign))

##Calculate likelihood ratio based on ClinVar
processedDataNormPath <- processedDataStats[processedDataStats$clinical_sign %in% c("Pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic"),]
processedDataNormBenign <- processedDataStats[processedDataStats$clinical_sign %in% c("Benign", "Benign/Likely benign", "Likely benign"),]
a <- nrow(processedDataNormPath[processedDataNormPath$processedData.class == "Non-functional",])
b <- nrow(processedDataNormBenign[processedDataNormBenign$processedData.class == "Non-functional",])
c <- nrow(processedDataNormPath[processedDataNormPath$processedData.class == "Functional",])
d <- nrow(processedDataNormBenign[processedDataNormBenign$processedData.class == "Functional",])
Sensitivity <- (a / (a + c)) * 100
Specificity <- (d / (b + d)) * 100
print(paste0("Sensitivity based on ClinVar: ", Sensitivity))
print(paste0("Specificity based on ClinVar: ", Specificity))
###Positive LR = sensitivity / (100 – specificity).
###Negative LR = (100 – sensitivity) / specificity.
LR.plus = Sensitivity / (100 - Specificity)
LR.neg = (100 - Sensitivity) / Specificity
print(paste0("Positive LR based on ClinVar: ", LR.plus))
print(paste0("Negative LR based on ClinVar: ", LR.neg))

##Calculate likelihood ratio based on codon
processedDataNormPath <- processedDataStats[processedDataStats$codon == "Non-sense",]
processedDataNormBenign <- processedDataStats[processedDataStats$codon == "Synonymous",]
a <- nrow(processedDataNormPath[processedDataNormPath$processedData.class == "Non-functional",])
b <- nrow(processedDataNormBenign[processedDataNormBenign$processedData.class == "Non-functional",])
c <- nrow(processedDataNormPath[processedDataNormPath$processedData.class == "Functional",])
d <- nrow(processedDataNormBenign[processedDataNormBenign$processedData.class == "Functional",])
Sensitivity <- (a / (a + c)) * 100
Specificity <- (d / (b + d)) * 100
print(paste0("Sensitivity based on codon type: ", Sensitivity))
print(paste0("Specificity based on codon type: ", Specificity))
LR.plus = Sensitivity / (100 - Specificity)
LR.neg = (100 - Sensitivity) / Specificity
print(paste0("Positive LR based on codon type: ", LR.plus))
print(paste0("Negative LR based on codon type: ", LR.neg))

#Extract summary data
processedDataSummary <- processedData[, c("exon", "variant", "edit_string", "pos", "codon", "nuc_change", "AA_change", "c_nom", "g_nom", "dbSNP.ID", "clinical_sign", "Protein.change", "Affect_splicing", "processedData.class",  "processedData.FS", "processedData.proba", "class.dmso", "function.score.dmso", "proba.dmso",  "class", "function.score", "proba", "class.ola", "function.score.ola", "proba.ola", "class.cis", "function.score.cis", "proba.cis")]
colnames(processedDataSummary) <- c("Exon", "Variant", "edit_string", "Position", "Codon.change", "Nucleotide.change", "AA_change", "c_nom", "g_nom", "dbSNP.ID", "ClinVar.annotation", "Protein.change", "Affect.splicing", "Classification", "Function.score",  "Probability", "DMSO.class", "DMSO.FS", "DMSO.proba", "Global.class", "Global.FS", "Global.proba", "Olaparib.class", "Olaparib.FS", "Olaparib.proba", "Cisplatin.class", "Cisplatin.FS", "Cisplatin.proba")

processedDataSummary <- processedDataSummary[!(duplicated(processedDataSummary)),]

#Save file
write.xlsx(file = paste0("../results/all_exons_results_summary", version, ".xlsx"), processedDataSummary)
write.xlsx(file = paste0("../results/all_exons_results_full", version, ".xlsx"), processedData)

#Explore results per exon
ggplot() + 
  geom_point(data = processedData[processedData$codon %in% c("Intronic", "Non-synonymous"),], aes(x = function.score.dmso, y = proba.dmso), col = "grey") +
  geom_point(data = processedData[processedData$codon %in% c("Non-sense", "Synonymous"),], aes(x = function.score.dmso, y = proba.dmso, col = codon)) +
  scale_color_manual(values = c("Intronic" = "darkgrey", "Non-sense" = "#c11600", "Non-synonymous" = "grey", "Synonymous" = "#0c4ee5")) +
  facet_wrap(. ~ exon) +
  labs(col = "SNV type") +
  xlab("Function score DMSO") +
  ylab("Probability of pathogenicity")+
  theme(legend.position = "bottom")