library(plyr)
library(RJSONIO)
# library(rpart)
library(randomForest)
# library(survMisc)
library(stringr)

options(stringsAsFactors=F)

selectCoverage <- function(dataSet){
  dataSet$bestCoverage <- ifelse(dataSet$X.percent_b.ion_coverage > dataSet$X.percent_y.ion_coverage,
                                 as.numeric(dataSet$X.percent_b.ion_coverage),
                                 as.numeric(dataSet$X.percent_y.ion_coverage))
  return(dataSet)
}

decideCoverage <- function(dataSet){
  dataSet <- selectCoverage(dataSet)
  dataSet$call <- apply(cbind(dataSet$meanCoverage > 0.4, dataSet$numStubs > 0, dataSet$percentUncovered < 0.2),
                        1, all)
  return(dataSet)
}

shimPercentCalcs <- function(dataSet){
  print(names(dataSet))
  results <- adply(dataSet, 1, function(row){
    data.frame(
         percent_b_ion_coverage = ifelse(!("total_b_ions_possible" %in% names(row)) || (row[, "total_b_ions_possible"] == 0), 0,length(fromJSON(row[, "b_ion_coverage"])) / row[, "total_b_ions"]),
         percent_y_ion_coverage = ifelse(!("total_y_ions_possible" %in% names(row)) || (row[, "total_y_ions_possible"] == 0), 0, length(fromJSON(row[, "y_ion_coverage"])) / row[, "total_y_ions"]),
         percent_b_ion_with_HexNAc_coverage = ifelse(!("possible_b_ions_HexNAc" %in% names(row)) || (row[, "possible_b_ions_HexNAc"] == 0), 0, length(fromJSON(row[, "b_ions_with_HexNAc"])) / row[, "possible_b_ions_HexNAc"]),
         percent_y_ion_with_HexNAc_coverage = ifelse(!("possible_y_ions_HexNAc" %in% names(row)) || (row[, "possible_y_ions_HexNAc"] == 0), 0, length(fromJSON(row[, "y_ions_with_HexNAc"])) / row[, "possible_y_ions_HexNAc"])
         )
  })
  results
}

resolveCall <- function(potentialCalls){
  if((nrow(potentialCalls) == 1) || !any(potentialCalls$call == "Yes")){
    return(potentialCalls)
  }


  # Order by any of the columns that are available
  signs <- c(1,1,1,-1,1,-1,1,1)
  fields <- c("MS2_Score","numStubs","meanCoverage",
              "percentUncovered","peptideLens",
              "abs_ppm_error",
              "percent_b_ion_with_HexNAc_coverage",
              "percent_y_ion_with_HexNAc_coverage")

  mask <- fields %in% colnames(potentialCalls)
  signs <- signs[mask]
  fields <- fields[mask]
  orderings <- lapply(seq_along(fields), function(i){
    signs[i] * potentialCalls[,fields[i]]
  })

  orderings <- c(orderings, decreasing = T)
  o <- do.call(order, orderings)

  sortedCalls <- potentialCalls[o,]
  sortedCalls$ambiguity <- T
  return(sortedCalls)
}

labelAmbiguity <- function(dataSet){
  if(!("ambiguity" %in% colnames(dataSet))){
    dataSet$ambiguity <- F
  }
  # Group by unique identifier pairing of MS1 Score and Mass
  entities <- ddply(dataSet, c("MS1_Score", "Calc_mass"), resolveCall)
  return(entities)
}

numAmbiguous <- function(dataSet){
  #if(!("ambiguity" %in% colnames(dataSet))) stop("Data must have been marked ambiguous first!")
  res <- dataSet[dataSet$ambiguity,]
  res <- ddply(res,
               c("MS1_Score", "Calc_mass"),
               function(ent){
                 count <- nrow(ent)
                 bestSeq <- ent[1, "Glycopeptide_identifier"]

                 # Try to give a mapping from sequence to score in JSON?
                 # Not reasonable in current function
                 #competativeRows <- ent[ent$MS2_Score >= 0.5,]

                 allSeqs <- toJSON(as.list(ent[ent$MS2_Score >= 0.5, "Glycopeptide_identifier"]))
                 worstScore <- min(ent[,"MS2_Score"])
                 bestScore <- max(ent[,"MS2_Score"])
                 res <- data.frame(count = count, bestSeq = bestSeq, allSeqs = allSeqs,
                                   worstScore = worstScore, bestScore = bestScore)
                 res
               })
  return(list(n = nrow(res), counts = res))
}

numOxoniumIons <- function(dataSet){
  oxIons <- dataSet$Oxonium_ions
  oxIonsList <- lapply(as.character(oxIons), fromJSON)
  dataSet$numOxIons <- sapply(oxIonsList, length)
  return(dataSet)
}

getSequenceLength <- function(peptideSequence){
  sapply(peptideSequence, function(seq){
    aminoAcids <- paste(strsplit(seq, split = "\\([^\\(\\)]*\\)")[[1]],
                        collapse = '')
    nchar(aminoAcids)
  })
}

numStubs <- function(dataSet, thresholdLength = 0, thresholdPenalty = 1){
  stubs <- dataSet$Stub_ions
  dataSet$peptideLens <- getSequenceLength(as.character(dataSet$Peptide))
  modifiers <- sapply(dataSet$peptideLens, function(ln){
    if(ln > thresholdLength) 1
    else thresholdPenalty
  })
  stubList <- lapply(as.character(stubs), fromJSON)
  dataSet$numStubs <- sapply(stubList, length) * modifiers
  return(dataSet)
}

labelCall <- function(call){
  factor(call, labels=c("No", "Yes"), levels = c(F, T))
}

loadDir <- function(dir, includeAnno = F){
  dataFiles = list.files(dir, pattern = ".csv", full.names=T)
  preparedDataSets <- lapply(dataFiles, prepareModelHandler, includeAnno)
  names(preparedDataSets) <- basename(dataFiles)
  check <- sapply(preparedDataSets, function(x)class(x)=="data.frame")
  return(preparedDataSets[check])
}

prepareModelHandler <- function(modelFile, includeAnno = F, ...){
  if(grepl("annotated", modelFile)){
    if(includeAnno) return(prepareAnnotatedModel(modelFile, ...))
    else return(NA)
  } else {
    return(prepareModel(modelFile, ...))
  }
}

prepareModel <- function(modelFile, readFile = T, call = T, classifierFn = NULL){
  if(readFile){
    rawDataSet <- read.csv(modelFile, stringsAsFactors = F)
    if(!("percent_b_ion_coverage" %in% names(rawDataSet))){
      rawDataSet <- shimPercentCalcs(rawDataSet)
    }
  } else {
    rawDataSet <- modelFile
  }
  rawDataSet$X.y_ion_with_HexNAc_coverage <- as.numeric(rawDataSet$X.y_ion_with_HexNAc_coverage)
  rawDataSet$X.b_ion_with_HexNAc_coverage <- as.numeric(rawDataSet$X.b_ion_with_HexNAc_coverage)
  nOxDataSet <- numOxoniumIons(rawDataSet)
  nOxDataSet$abs_ppm_error <- abs(nOxDataSet$ppm_error)
  nStubsDataSet <- numStubs(nOxDataSet)
  preparedDataSet <- getIonCoverageScore(nStubsDataSet)
  if(call){
    if(!is.null(classifierFn)){
      calledDataSet <- classifierFn(preparedDataSet)
    } else{
      calledDataSet <- decideCoverage(preparedDataSet)
      calledDataSet$call <- labelCall(calledDataSet$call)
    }
    return(calledDataSet)
  }
  return(preparedDataSet)
}

#' Prepare a model frame from a manually annotated data file.
prepareAnnotatedModel <- function(modelFile, readFile = T, drop = T){
  if(readFile){
    rawDataSet <- read.csv(modelFile, stringsAsFactors = F)
    if(!("percent_b_ion_coverage" %in% names(rawDataSet))){
      rawDataSet <- shimPercentCalcs(rawDataSet)
    }
  } else {
    rawDataSet <- modelFile
  }
  if(drop){
    cleanDataSet <- rawDataSet[!apply(cbind(is.na(rawDataSet$`TRUE.FALSE`),
                                            rawDataSet$`TRUE.FALSE` == ""),1,any) ,]
  } else {
    cleanDataSet <- rawDataSet
  }
  calledDataSet <- selectCoverage(cleanDataSet)
  calledDataSet$abs_ppm_error <- abs(calledDataSet$ppm_error)
  preparedDataSet <- getIonCoverageScore(
    numOxoniumIons(
      numStubs(
        calledDataSet
      )
    )
  )
  preparedDataSet$call <- labelCall(preparedDataSet$`TRUE.FALSE`)
  return(preparedDataSet)
}

# complexModelFormula <- call ~ numOxIons + X.percent_b.ion_coverage +
#   X.percent_y.ion_coverage + X.b_ion_with_HexNAc_coverage +
#   X.y_ion_with_HexNAc_coverage +  numStubs + peptideLens

# complexModelFormula <- call ~ X.percent_b.ion_coverage +
#   X.percent_y.ion_coverage + X.b_ion_with_HexNAc_coverage +
#   X.y_ion_with_HexNAc_coverage +  numStubs + peptideLens +
#   ppm_error + meanCoverage + percentUncovered

complexModelFormula <- call ~ meanCoverage + percentUncovered + abs_ppm_error + X.b_ion_with_HexNAc_coverage +
  X.y_ion_with_HexNAc_coverage + numStubs + peptideLens

# complexModelFormulaWithInteractions <- call ~ X.percent_b.ion_coverage *
#   X.percent_y.ion_coverage * X.b_ion_with_HexNAc_coverage *
#   X.y_ion_with_HexNAc_coverage *  numStubs * peptideLens *
#   ppm_error * meanCoverage * percentUncovered

complexModelFormulaWithInteractions <- call ~ meanCoverage * percentUncovered *
  X.b_ion_with_HexNAc_coverage * X.y_ion_with_HexNAc_coverage *
  numStubs * peptideLens * ppm_error

#' Fits a model on the dataset
fitModel <- function(dataSet, modelFormula = complexModelFormula){
  if(all(dataSet$call == "No")) return(NA)
  model <- randomForest(modelFormula, data = dataSet, mtry = 6, ntree = 5000, importance = T)
  return(model)
}

#' Given a dataset with known information and a model's predictions about
#' class labels, describe the accuracy of the model's predictions.
checkModel <- function(dataSet, modelPred, threshold = 0.5, percent = F){
  predCall <- apply(modelPred, 1, function(row){
    (row["No"] < row["Yes"]) && (row["Yes"] >= threshold)
  })
  div = ifelse(percent, nrow(dataSet), 1)
  predCall <- labelCall(predCall)
  cmprCalls <- data.frame(predCall = predCall,
                          trueCall = dataSet$call)
  errorTable <- dlply(cmprCalls, c("predCall", "trueCall"),
                      function(responseType){
                        nrow(responseType)/div
                      })

  # Guarantee all rows are the same length
  if(is.null(errorTable$Yes.No)){
    errorTable$Yes.No <- 0
  }
  if(is.null(errorTable$No.No)){
    errorTable$No.No <- 0
  }
  if(is.null(errorTable$No.Yes)){
    errorTable$No.Yes <- 0
  }
  if(is.null(errorTable$Yes.Yes)){
    errorTable$Yes.Yes <- 0
  }

  errorTable[[".numNo"]]  <- sum(dataSet$call == "No")
  errorTable[[".numYes"]] <- sum(dataSet$call == "Yes")
  #errorTable[["TPR"]]  <- errorTable$Yes.Yes / (errorTable$Yes.Yes + errorTable$No.Yes)
  errorTable[["FPR"]]  <- errorTable$Yes.No / (errorTable$Yes.No + errorTable$No.No)
  errorTable[["FDR"]]  <- errorTable[["FPR"]] / (errorTable[["FPR"]] + errorTable$Yes.Yes)
  errorTable[["ErrorRate"]] <- (errorTable$No.Yes + errorTable$Yes.No) / nrow(dataSet)
  # Ensure that all names are in the same order, with match types first, statistics second
  errorTable <- errorTable[rev(sort(names(errorTable)))]
  names(errorTable)[1:4] <- c("TruePositive", "FalsePositive",
                              "FalseNegative", "TrueNegative")

  return(errorTable)
}

#' Given a dataset with known information and a model's prediction about
#' class labels, label the false positives
getFalsePositives <- function(dataSet, modelPred, threshold = 0.5){
  predCall <- apply(modelPred, 1, function(row){
    (row["No"] < row["Yes"]) && (row["Yes"] >= threshold)
  })
  predCall <- labelCall(predCall)
  cmprCalls <- data.frame(predCall = predCall,
                          trueCall = dataSet$call)

  fPos <- apply(cmprCalls, 1, function(row){
    (row[[1]] == "Yes") && (row[[2]] == "No")
  })

  res <- cbind(dataSet[fPos,], modelPred[fPos,])
  names(res)[length(names(res))] <- "fPos"

  return(res)
}

#' Given a dataset with known information and a model's prediction about
#' class labels, label the true positives
getTruePositives <- function(dataSet, modelPred, threshold = 0.5){
  predCall <- apply(modelPred, 1, function(row){
    (row["No"] < row["Yes"]) && (row["Yes"] >= threshold)
  })
  predCall <- labelCall(predCall)
  cmprCalls <- data.frame(predCall = predCall,
                          trueCall = dataSet$call)

  fPos <- apply(cmprCalls, 1, function(row){
    (row[[1]] == "Yes") && (row[[2]] == "Yes")
  })

  res <- cbind(dataSet[fPos,], modelPred[fPos,])
  names(res)[length(names(res))] <- "tPos"

  return(res)
}

#' Given a dataset with known information and a model's prediction about
#' class labels, label the false negatives
getFalseNegatives <- function(dataSet, modelPred, threshold = 0.5){
  predCall <- apply(modelPred, 1, function(row){
    (row["No"] < row["Yes"]) && (row["Yes"] >= threshold)
  })
  predCall <- labelCall(predCall)
  cmprCalls <- data.frame(predCall = predCall,
                          trueCall = dataSet$call)

  fNeg <- apply(cmprCalls, 1, function(row){
    (row[[1]] == "No") && (row[[2]] == "Yes")
  })

  res <- cbind(dataSet[fNeg,], modelPred[fNeg,])
  names(res)[length(names(res))] <- "fNeg"

  return(res)
}


#' Given a dataset with known information and a model's prediction about
#' class labels, label the true negatives
getTrueNegatives <- function(dataSet, modelPred, threshold = 0.5){
  predCall <- apply(modelPred, 1, function(row){
    (row["No"] > row["Yes"]) && (row["No"] >= threshold)
  })
  predCall <- labelCall(predCall)
  cmprCalls <- data.frame(predCall = predCall,
                          trueCall = dataSet$call)

  tNeg <- apply(cmprCalls, 1, function(row){
    (row[[1]] == "No") && (row[[2]] == "No")
  })

  res <- cbind(dataSet[tNeg,], modelPred[tNeg,])
  names(res)[length(names(res))] <- "tNeg"

  return(res)
}

computeIonCoverageMap <- function(row){
  bIons <- fromJSON(row[,"bare_b_ions"])
  bIonsWithHexNac <- fromJSON(row[,"b_ions_with_HexNAc"])

  yIons <- fromJSON(row[,"bare_y_ions"])
  yIonsWithHexNac <- fromJSON(row[,"y_ions_with_HexNAc"])

  peptideLen <- as.numeric(row[,"peptideLens"])

  totalCoverage <- numeric(peptideLen)

  bIonCov <- lapply(bIons, function(b){
    pos <- as.numeric(strsplit(x=b$key, split="B")[[1]][2])
    coverageInformation <- numeric(pos) + 1
    length(coverageInformation) <- peptideLen
    coverageInformation[is.na(coverageInformation)] <- 0
    coverageInformation
  })

  bIonCovTot <- Reduce(f=function(a,b) a + b, x=bIonCov, init=totalCoverage)

  bIonHexCov <- lapply(bIonsWithHexNac, function(b){
    pos <- as.numeric(gsub(pattern="B(\\d+).*", "\\1", b$key))
    coverageInformation <- numeric(pos) + 1
    length(coverageInformation) <- peptideLen
    coverageInformation[is.na(coverageInformation)] <- 0
    coverageInformation
  })

  bIonHexCovTot <- Reduce(f=function(a,b) a + b, x=bIonHexCov, init=totalCoverage)

  yIonCov <- lapply(yIons, function(y){
    pos <- as.numeric(strsplit(y$key, "Y")[[1]][2])
    coverageInformation <- numeric(pos) + 1
    length(coverageInformation) <- peptideLen
    coverageInformation[is.na(coverageInformation)] <- 0
    rev(coverageInformation)
  })

  yIonCovTot <- Reduce(f=function(a,b) a + b, x=yIonCov, init=totalCoverage)

  yIonHexCov <- lapply(yIonsWithHexNac, function(y){
    pos <- as.numeric(gsub(pattern="Y(\\d+).*", "\\1", y$key))
    coverageInformation <- numeric(pos) + 1
    length(coverageInformation) <- peptideLen
    coverageInformation[is.na(coverageInformation)] <- 0
    rev(coverageInformation)
  })

  yIonHexCovTot <- Reduce(f=function(a,b) a + b, x=yIonHexCov, init=totalCoverage)

  totalCoverage <- (bIonCovTot + bIonHexCovTot + yIonCovTot + yIonHexCovTot)

  percentUncovered <- sum(totalCoverage == 0)/peptideLen

  return(list(totalCoverage = totalCoverage, percentUncovered = percentUncovered))
}

getIonCoverageScore <- function(dataSet){
  mapData <- adply(dataSet, 1, function(row){
    covMap <- computeIonCoverageMap(row)
    meanCov <- mean(covMap$totalCoverage)
    return(data.frame(meanCoverage = meanCov, percentUncovered = covMap$percentUncovered))
  })

  return(mapData)
}


collapseLocalizeFrameRows <- function(locFrame){
  ddply(locFrame, "Position", function(pos){
    which.max(pos$nHexNac)
  })
}

localizeHexNac <- function(row){
  bIonsWithHexNac <- fromJSON(row[,"b_ions_with_HexNAc"])
  yIonsWithHexNac <- fromJSON(row[,"y_ions_with_HexNAc"])

  #print(row[,"Glycopeptide_identifier"])

  peptideLen <- as.numeric(row[,"peptideLens"])
  totalCoverage <- numeric(peptideLen)

  bIonHexCov <- ldply(bIonsWithHexNac, function(b){
    #print(b$key)
    pos <- as.numeric(str_match(b$key, "B(\\d+)\\+(\\d*).*")[,c(2,3)])
    pos[is.na(pos)] <- 1
    pos
  })
  if(nrow(bIonHexCov) > 0){
    names(bIonHexCov) <- c("Position", "nHexNac")
    bIonHexCov <- bIonHexCov[order(bIonHexCov$Position),]
  }

  yIonHexCov <- ldply(yIonsWithHexNac, function(y){
    #print(y$key)
    pos <- as.numeric(str_match(y$key, "Y(\\d+)\\+(\\d*).*")[,c(2,3)])
    pos[is.na(pos)] <- 1
    #print(pos)
    pos
  })
  if(nrow(yIonHexCov) > 0){
    names(yIonHexCov) <- c("Position", "nHexNac")
    yIonHexCov <- yIonHexCov[order(yIonHexCov$Position),]
  }

  return(list(B = bIonHexCov, Y = yIonHexCov))

}
