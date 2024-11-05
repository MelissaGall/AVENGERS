# Functions for model prediction

# Create data frame with all variants and perform first pre processing steps
createVarSimpDF <- function(varData, edit3, edit5, posToKeep, exonStart, exonStop, remove){
  # Prepare data
  colToKeep <- c("variant", "edit_string", "pos", "pre",  "pre_freq", "post","post_freq", "codon", "post.pre.ratio.sns", "nuc_change", "AA_change", "c_nom", "g_nom", "dbSNP.ID", "clinical_sign", "Protein.change", "expected_clinical", "Affect_splicing")
  # Extract position of the mutation as a number
  varData$pos <- varData$edit_string
  varData$pos <- sub(edit3, "", varData$pos)
  varData$pos <- sub(edit5, "", varData$pos)
  varData$pos <- sub(",,", "", varData$pos)
  varData$pos <- sub(",", "", varData$pos)
  varData$pos <- sub("-X-.*", "", varData$pos)
  varData$pos <- as.numeric(varData$pos)
  # Add info on intronic variants
  varData$codon <- ifelse(varData$AA_change == "Intronic", "Intronic", varData$codon)
  varData$codon[(varData$pos < exonStart | varData$pos > exonStop)] <- "Intronic" 
  varData$AA_change[varData$codon == "Intronic"] <- "Intronic"
  # Filter to keep only SNV (remove WT and PAM)
  varData <- varData[!is.na(varData$pos),]
  varData <- varData[varData$pos %in% posToKeep,]
  varData <- varData[!duplicated(varData),]
  # If post freq is 0, convert to very low frequency to allow calculations
  varData$post_freq <- ifelse(varData$post == 0, 0.00001, varData$post_freq)
  # Filter SNV with low count
  varDataThr <- varData[(varData$pre >= 10 & varData$pre_freq > freqThr),]
  # If remove == TRUE (second run), remove filtered variants
  if (remove == TRUE){
    varDataThr <- varDataThr[!(varDataThr$variant %in% excludedVar),]
  }
  # Compute ratio
  varDataThr$post.pre.ratio <- varDataThr$post_freq / varDataThr$pre_freq
  # Assign expected behavior
  varDataThr <- varDataThr %>%
    mutate(expected_clinical = case_when((codon == "Synonymous") & !(clinical_sign %in% c("Pathogenic", "Likely pathogenic")) ~ "Functional",
                     (codon == "Non-sense") & !(clinical_sign %in% c("Benign", "Likely benign")) ~ "Non-functional",
                     TRUE ~ "Unknown"))
  # Set synonymous that affect splicing as "unknown"
  varDataThr$expected_clinical[varDataThr$Affect_splicing == "Yes" & varDataThr$codon != "Non-sense"] <- "Unknown"
  # Normalize
  medianSyn <- median(varDataThr$post.pre.ratio[varDataThr$expected_clinical == "Functional"])
  medianNS <- median(varDataThr$post.pre.ratio[varDataThr$expected_clinical == "Non-functional"])
  # If no non-functional data available, take lower data as reference
  if (is.na(medianNS) | nrow(varDataThr[varDataThr$expected_clinical == "Non-functional",]) < 1){
    medianNS <- min(varDataThr$post.pre.ratio)
    Q5.low <- quantile(varDataThr$post.pre.ratio, probs = 0.05, na.rm = T)
    medianNS <- median(varDataThr$post.pre.ratio[varDataThr$post.pre.ratio < Q5.low], na.rm = T) 
    }
  varDataThr$post.pre.ratio.sns <- varDataThr$post.pre.ratio/(medianSyn - medianNS)
  # Save
  varDataThr <- varDataThr[, colToKeep]
  varDataThr <- varDataThr[!duplicated(varDataThr),]  
  return(varDataThr)
}

# Remove outliers
removeOutliers <- function(data, col){
  lowerBoundSyn <- quantile(data[data$expected_clinical == "Functional", c(col)], 0.1, na.rm = T)
  upperBoundSyn <- quantile(data[data$expected_clinical == "Functional", c(col)], 0.9, na.rm = T)
  lowerBoundNS <- quantile(data[data$expected_clinical == "Non-functional", c(col)], 0.1, na.rm = T)
  upperBoundNS <- quantile(data[data$expected_clinical == "Non-functional", c(col)], 0.9, na.rm = T)
  data$outlier <- "No"
  data$outlier <- ifelse(data$codon == "Synonymous" & (data[, c(col)] > upperBoundSyn | data[, c(col)] < lowerBoundSyn), "Yes", data$outlier)
  data$outlier <- ifelse(data$codon == "Non-sense" & (data[, c(col)] > upperBoundNS | data[, c(col)] < lowerBoundNS), "Yes", data$outlier)
  data$expected_clinical[data$outlier == "Yes"] <- "Unknown"
  return(data)
}

# Extend fit for LOESS
extendFit <- function(fitVector) {
  if (length(which(is.na(fitVector))) == 0) {
    outVector <- fitVector
  }
  else {
    minNA <- min(which(is.na(fitVector) == TRUE))
    minPoint <- min(which(is.na(fitVector) == FALSE))
    minPointVal <- fitVector[minPoint]
    maxNA <- max(which(is.na(fitVector) == TRUE))
    maxPoint <- max(which(is.na(fitVector) == FALSE))
    maxPointVal <- fitVector[maxPoint]
    outVector <- fitVector
    if (minNA == 1) {
      outVector[minNA:minPoint] <- minPointVal
    }
    if (maxNA == length(fitVector)) {
      outVector[(maxPoint+1):length(fitVector)] <- maxPointVal 
    }
  }
  return(outVector)
}

# Get positional biais
removePositionBiais <- function(data, col){
  posRange <- range(as.numeric(as.character(data$pos)))
  posseq <- seq(from = posRange[1], to = posRange[2])
  Q3 <- quantile(data[, c(col)], probs = 0.25, na.rm = T) #Exclude lower values to avoid bias
  postPreRatioLoess <- loess(log2(data[, c(col)]) ~ as.numeric(as.character(pos)), 
                             data, 
                             subset = which(data[, c(col)] > Q3), 
                             span = 0.75,
                             model = TRUE)
  postPrePred <- predict(postPreRatioLoess, newdata = posseq, se = TRUE)
  postPrePosEffects <- postPrePred$fit
  postPrePosEffects.x <- extendFit(postPrePosEffects)
  postPrePosEffects.df <- data.frame(pos = posseq, 
                                     postPrePosEffects.x)
  postPrePosEffectsData <- merge(data, postPrePosEffects.df[, c("pos", "postPrePosEffects.x")], by = "pos", all.x = TRUE)
  return(log2(postPrePosEffectsData[, c(col)]) - postPrePosEffectsData$postPrePosEffects.x)
}

# Per exon normalization
perExonNorm <- function(data, ampList, col){
  # Calculate median for syn and ns variants 
  medianSyn <- median(data[data$expected_clinical == "Functional", c(col)], na.rm = T)
  medianNS <- median(data[data$expected_clinical == "Non-functional", c(col)], na.rm = T)
  normData <- data.frame()
  # For each experiment, normalize per exon
  for (amplicon in ampList){
    exon <- parse_number(amplicon)
    exonData <- data[data$exon == amplicon,]
    subsetMedianSyn <- median(exonData[exonData$codon == "Synonymous", c(col)], na.rm = T)
    subsetMedianNS <- median(exonData[exonData$codon == "Non-sense", c(col)], na.rm = T)
    # If not enough data available for calculation, use minimal values
    if (is.na(subsetMedianNS) | nrow(exonData[exonData$expected_clinical == "Non-functional",]) < 1){
      Q5.low <- quantile(exonData[, c(col)] , probs = 0.05, na.rm = T)
      subsetMedianNS <- median(exonData[exonData[,c(col)] < Q5.low, c(col)], na.rm = T) 
    }
  exonData$function.score <- (exonData[, c(col)]/subsetMedianNS) * (medianSyn - medianNS)
  normData <- rbind(normData, exonData)
  }
  return(as.vector(normData$function.score))
}

# Global normalization
globalNorm <- function(data, ampList, col){
  # Calculate median for syn and ns variants 
  medianSyn <- median(data[data$codon == "Synonymous", c(col)], na.rm = T)
  medianNS <- median(data[data$codon == "Non-sense", c(col)], na.rm = T)
  data$normData <- (data[, c(col)] / medianNS) * (medianSyn - medianNS)
  return(as.vector(data$normData))
}

# Assign class
assingLOF <- function(expected_clinical){
  if (expected_clinical == "Non-functional") return(1) else return(0)
}

# Run the model
runMclust <- function(trainingData, data, col, colNew){
  # Remove NA
  trainingData <- trainingData[!(is.na(trainingData[, c(col)])),]
  dataSubset <- data[!(is.na(data[, c(col)])),]
  Q10.syn <- quantile(trainingData[(trainingData$expected_clinical == "Functional"), c(col)] , probs = 0.2, na.rm = T)
  Q10.ns <- quantile(trainingData[(trainingData$expected_clinical == "Non-functional"), c(col)] , probs = 0.8, na.rm = T)
  # Remove outliers
  trainingData <- trainingData[!(trainingData$expected_clinical == "Functional" & trainingData[, c(col)] < Q10.syn),]
  trainingData <- trainingData[!(trainingData$expected_clinical == "Non-functional" & trainingData[, c(col)] > Q10.ns),]
  # Run the prediction
  x_mclust.train <- MclustDA(data = trainingData[, c(col)], class = trainingData$lof, G = 1:20, modelNames = c("V"), verbose = F)
  x_mclust.predict <- predict.MclustDA(x_mclust.train, newdata = dataSubset[, c(col)], newclass = trainingData$lof, prop = c(0.9, 0.1))
  # Merge data
  dataSubset[, c(colNew)] <- x_mclust.predict$z[,2]
  colToMerge <- c("variant", "edit_string", "pos")
  data <- merge(data, dataSubset[, c(colToMerge, colNew)], by = colToMerge, all.x = TRUE)
  return(data)
}
