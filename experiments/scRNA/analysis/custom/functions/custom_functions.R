library(rprojroot)
library(dplyr)
library(tidyr)

# Append number to name of file to save if it already exists

appendNumberToPath <- function(proot, folder, name, extension) {
    counter <- 0
    exists <- TRUE
    path <- NULL

    while (exists) {
        nameNew <- paste0(
            name,
            ifelse(counter == 0, "", paste0("_", counter)),
            ".", extension
        )

        path <- file.path(folder, nameNew)
        exists <- file.exists(proot$find_file(path))
        counter <- counter + 1
    }

    return(path)
}

# Get tables with fractions of clusters by group and batch

getClusterFractionTables <- function(SObj, clusterCol, groupCol, batchCol) {
    selectionCols <- list()

    selectionCols[["batch"]] <- batchCol
    selectionCols[["group"]] <- groupCol
    selectionCols[["cluster"]] <- clusterCol

    dataSelection <- SObj@meta.data %>%
        select(
            .data[[selectionCols[["batch"]]]],
            .data[[selectionCols[["group"]]]],
            .data[[selectionCols[["cluster"]]]],
        )

    sumTables <- list()

    for (selCol in names(selectionCols)) {
        newColName <- paste0("count_", selCol)
        sumTables[[selCol]] <- dataSelection %>%
            group_by(.data[[selectionCols[[selCol]]]]) %>%
            summarise({{ newColName }} := n())
    }

    fractionCols <- c("group", "batch")
    clusterFractionTables <- list()
    clusterFractionPivotTables <- list()

    for (fractionCol in fractionCols) {
        clusterFractionTables[[fractionCol]] <- dataSelection %>%
            group_by(
                .data[[selectionCols[["cluster"]]]],
                .data[[selectionCols[[fractionCol]]]]
            ) %>% summarise(count = n())

        clusterFractionTables[[fractionCol]] <- inner_join(
            clusterFractionTables[[fractionCol]],
            sumTables[["cluster"]],
            by = selectionCols[["cluster"]]  
        )

        clusterFractionTables[[fractionCol]] <- inner_join(
            clusterFractionTables[[fractionCol]],
            sumTables[[fractionCol]],
            by = selectionCols[[fractionCol]]
        )

        sumCol <- paste0("count_", fractionCol)
        newCol <- paste0("fraction_", fractionCol)

        clusterFractionTables[[fractionCol]] <- clusterFractionTables[[fractionCol]] %>% 
        mutate(
            fraction_cluster = count/count_cluster,
            {{ newCol }} := count/.data[[sumCol]]
        )

        clusterFractionPivotTables[[fractionCol]] <- list()

        clusterFractionPivotTables[[fractionCol]][["cluster"]] <- pivot_wider(
            clusterFractionTables[[fractionCol]],
            id_cols = selectionCols[[fractionCol]],
            names_from = selectionCols[["cluster"]],
            values_from = "fraction_cluster"
        )
        clusterFractionPivotTables[[fractionCol]][[fractionCol]] <- pivot_wider(
            clusterFractionTables[[fractionCol]],
            id_cols = selectionCols[["cluster"]],
            names_from = selectionCols[[fractionCol]],
            values_from = newCol
        )
    }

    resultTables <- list(
        "long" = clusterFractionTables,
        "pivot" = clusterFractionPivotTables
    )

    return(resultTables)
}

# From the list of cluster marker genes derived from Seurat FindMarkers()
# and a list of celltypes with their known markers as character vectors,
# get information on known markers found in each cluster

getMarkersFoundInfo <- function(markers, knownMarkers, pAdjThreshold=0.05, logFCThreshold=0.25) {
    markers <- markers[markers$p_val_adj <= pAdjThreshold, ]
    markers <- markers[markers$avg_log2FC >= logFCThreshold, ]
    clusters <- unique(markers$cluster)
    clusterNames <- paste0("cluster_", clusters)
    names(clusters) <- clusterNames

    markersFoundInfo <- list()

    for (clusterName in names(clusters)) {
        markersCluster <- markers[markers[["cluster"]] == clusters[[clusterName]], ]
        markersFoundInfo[[clusterName]] <- list()
        
        for (celltype in names(knownMarkers)) {
            markersFoundInfo[[clusterName]][[celltype]] <- list()

            foundBool <- knownMarkers[[celltype]] %in% markersCluster[["gene"]]
            names(foundBool) <- knownMarkers[[celltype]]
            foundNames <- knownMarkers[[celltype]][foundBool]
            
            pAdj <- rep(NA, length = length(foundBool))
            names(pAdj) <- knownMarkers[[celltype]]
            pAdj[foundNames] <- markersCluster[match(foundNames, markersCluster[["gene"]]), ] %>% pull(p_val_adj)
            pAdj[pAdj == 0] <- .Machine$double.xmin

            avgLog2FC <- rep(NA, length = length(foundBool))
            names(avgLog2FC) <- knownMarkers[[celltype]]
            avgLog2FC[foundNames] <- markersCluster[match(foundNames, markersCluster[["gene"]]), ] %>% pull(avg_log2FC)

            meanLogFC <- mean(avgLog2FC, na.rm = TRUE)
            meanNegLogP <- mean((-1)*log10(pAdj), na.rm = TRUE)

            markersFoundInfo[[clusterName]][[celltype]][["gene_found"]] <- foundBool
            markersFoundInfo[[clusterName]][[celltype]][["p_val_adj"]] <- pAdj
            markersFoundInfo[[clusterName]][[celltype]][["avg_log2FC"]] <- avgLog2FC
            markersFoundInfo[[clusterName]][[celltype]][["fraction_found"]] <-
                sum(foundBool)/length(foundBool)
            markersFoundInfo[[clusterName]][[celltype]][["mean_log_FC"]] <-
                ifelse(is.nan(meanLogFC), 0, meanLogFC)
            markersFoundInfo[[clusterName]][[celltype]][["mean_neg_log_p"]] <-
                ifelse(is.nan(meanNegLogP), 0, meanNegLogP)
        }
    }

    return(markersFoundInfo)
}

getMarkersFound <- function(markersFoundInfo) {
    markersFound <- list()

    for (clusterName in names(markersFoundInfo)) {
        markersFound[[clusterName]] <- list()
        for (celltype in names(markersFoundInfo[[clusterName]])) {
            markersFound[[clusterName]][[celltype]] <-
                markersFoundInfo[[clusterName]][[celltype]][["gene_found"]]
        }
    }
    return(markersFound)
}

getMarkersFoundTable <- function(markersFoundInfo, onlyFractions = TRUE) {
    markersFoundSummary <- list()

    celltypes <- names(markersFoundInfo[[1]])

    if (onlyFractions) {
        summaryStats <- c("fraction_found")
    } else {
        summaryStats <- c("fraction_found", "mean_log_FC", "mean_neg_log_p")
    }

    rowNamesGrid <- expand.grid(summaryStats, celltypes)
    rowNames <- paste0(rowNamesGrid[,2], "_", rowNamesGrid[,1])

    for (clusterName in names(markersFoundInfo)) {
        markersFoundSummary[[clusterName]] <-
            rep(NA, length = length(rowNames))
        names(markersFoundSummary[[clusterName]]) <- rowNames

        for (celltype in names(markersFoundInfo[[clusterName]])) {
            for (stat in summaryStats) {
                rowName <- paste0(celltype, "_", stat)
                markersFoundSummary[[clusterName]][rowName] <-
                    markersFoundInfo[[clusterName]][[celltype]][[stat]]
            }
        }
    }

    markersFoundTable <- data.frame(markersFoundSummary)

    return(markersFoundTable)
}

# Regress out variables from feature matrix (features as rows, samples as columns)

# Metadata rows must correspond to samples

regressOutVars <- function(
    mat,
    metadata,
    varsToRegress,
    features = NULL,
    returnModel = FALSE,
    verbose = FALSE
) {
    if (is.null(features)) {
        features <- rownames(mat)
    }

    if (returnModel) {
        linearModels <- list()
    }

    fitResiduals <- list()

    modelFormula <- formula(paste0("y ~ ", paste0(varsToRegress, collapse = " + ")))

    featureNo <- 1

    for (feature in features) {
        if (verbose && ((featureNo %% 1000) == 0)) {
            message(paste0("Variables regressed out for ", featureNo, " features"))
        }

        currDf <- data.frame("y" = mat[feature, ])
        currDf[, varsToRegress] <- metadata[, varsToRegress]

        linearModel <- lm(modelFormula, data = currDf)
        fitResiduals[[feature]] <- residuals(linearModel)

        if (returnModel) {
            linearModels[[feature]] <- linearModel
        }

        featureNo <- featureNo + 1
    }

    fitResiduals <- t(as.matrix(data.frame(fitResiduals, check.names = FALSE)))

    results <- list(
        "fit_residuals" = fitResiduals
    )

    if (returnModel) {
        results[["linear_models"]] <- linearModels
    }

    if (verbose) {
        message(paste0("Variables regressed out for ", featureNo - 1, " features in total"))
    }

    return(results)
}

# Run VST according to Seurat FindVariableFeatures()

# Feature matrix with features as rows, samples as columns
# Possible to supply custom means if data are e.g. residuals after regressing out

runVst <- function(
    mat,
    logTransform = TRUE,
    means = NULL,
    loess_span = 0.3,
    clip_max = "auto",
    center = TRUE,
    verbose = FALSE
) {
    mat <- t(mat)

    if (is.null(means)) {
        means <- apply(mat, 2, mean)
    }
    variances <- apply(mat, 2, var)

    vstTable <- data.frame(
        mean = means,
        var = variances
    )

    nonConstantFeatures <- rownames(vstTable[vstTable$var > 0, ])
    vstTableFit <- vstTable[nonConstantFeatures, ]

    if (logTransform) {
        vstTableFit$mean <- log10(vstTableFit$mean)
        vstTableFit$var <- log10(vstTableFit$var)
    }

    loessFit <- loess(
        var ~ mean,
        data = vstTableFit,
        span = loess_span
    )
    variancesLoess <- predict(loessFit)

    if (logTransform) {
        variancesLoess <- 10 ^ variancesLoess
    }

    vstTable[["var_loess"]] <- 0
    vstTable[nonConstantFeatures, "var_loess"] <- variancesLoess

    matVst <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    rownames(matVst) <- rownames(mat)
    colnames(matVst) <- colnames(mat)

    featureNo <- 1

    for (feature in nonConstantFeatures) {
        featureMean <- vstTable[feature, "mean"]
        centerValue <- ifelse(center, featureMean, 0)
        featureValues <- mat[, feature]
        varLoess <- vstTable[feature, "var_loess"]

        matVst[, feature] <- (featureValues - centerValue)/sqrt(varLoess)

        if(verbose && ((featureNo %% 1000) == 0)) {
            message(paste0("VST computed for ", featureNo, " features"))
        }

        featureNo <- featureNo + 1
    }
    
    if (is.null(clip_max)) {
        clip_max <- Inf
    } else if (clip_max == "auto") {
        clip_max <- sqrt(nrow(mat))
    }

    matVst[matVst > clip_max] <- clip_max

    vstTable[["var_vst"]] <- apply(matVst, 2, var)

    matVst <- t(matVst)

    results <- list(
        "vsd" = matVst,
        "vst_table" = vstTable
    )

    if (verbose) {
        message(paste0("VST computed for ", featureNo - 1, " non-zero features in total"))
    }

    return(results)
}

# Find variable features in VST data after variables have been regressed out from VST data

findVarFeatAfterRegressOut <- function(
    SObj,
    varsToRegress,
    nfeatures = NULL,
    features = NULL,
    loess_span = 0.3,
    clip_max = "auto",
    verbose = FALSE
) {
    counts <- as.matrix(GetAssayData(SObj, slot = "counts"))

    vstResults <- runVst(
        mat = counts,
        logTransform = TRUE,
        means = NULL,
        loess_span = loess_span,
        clip_max = clip_max,
        center = TRUE,
        verbose = verbose
    )

    regressOutFit <- regressOutVars(
        mat = vstResults$vsd,
        metadata = SObj@meta.data,
        varsToRegress = varsToRegress,
        features = features,
        verbose = verbose
    )

    variances <- apply(t(regressOutFit$fit_residuals), 2, var)
    variancesSorted <- sort(variances, decreasing = TRUE)

    if (is.null(nfeatures)) {
        nfeatures <- length(variancesSorted)
    }

    return(names(variancesSorted)[1:nfeatures])
}