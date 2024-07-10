### Libraries

# Using rbioc3-14

library(rprojroot)
library(yaml)
library(Seurat)
library(SingleR)
library(monocle)
library(stringr)
library(ggplot2)
library(dplyr)
library(WriteXLS)
library(future)

### Parameters

annotCol <- "annot_SingleR_mod"
analysisLabel <- "Stem-TA_ref_regress-cc"

lfcThreshold <- 0.25
pAdjThreshold <- 0.05

### Paths

root_criterion <- ".projroot"
proot <- rprojroot::has_file(root_criterion)

source(proot$find_file("experiments/scRNA/analysis/custom/functions/custom_functions.R"))

params <- list()
params$config <- proot$find_file("experiments/scRNA/analysis/custom/hep62703_2023-03-24_A_harmony/config.yaml")
params$analparams <- proot$find_file("experiments/scRNA/analysis/custom/hep62703_2023-03-24_A_harmony/analparams_preprocessing_filters_2023-03-24.yaml")
params$date <- Sys.Date()

set.seed(7)

config <- yaml.load_file(proot$find_file(params$config))
analparams <- yaml.load_file(proot$find_file(params$analparams))
today <- params$date

plotPath <- proot$find_file(config$outputpath, "fig", "monocle", analysisLabel)

### Parallelization

# Increase max size, otherwise it is exceeded for these data
options(future.globals.maxSize = 1200 * 1024^2)

future.seed=TRUE
plan("multicore", workers = 6)
plan()

### Load data

SObj <- readRDS(proot$find_file(analparams$Seurat_clustered))

# Differential expression testing results between Stem and TA labels from Haber et al. 2017

diffExprRefStemTA <- read.csv(
    proot$find_file(config$outputpath, "table", "diff_expr_Stem-TA_ref-Haber.csv"),
    row.names = 1
)

# Get Stem/TA subset of Seurat object

Idents(SObj) <- annotCol
SObjSubset <- subset(SObj, idents = c("Stem", "TA"))

### Run Monocle

## Create Monocle object

countsMono <- GetAssayData(SObjSubset, slot = "counts")
phenoDataMono <- new("AnnotatedDataFrame", data = SObjSubset@meta.data)
featureDataMono <- new(
    "AnnotatedDataFrame",
    data = data.frame(
        "gene" = rownames(SObjSubset),
        "gene_short_name" = rownames(SObjSubset),
        row.names = rownames(SObjSubset)
    )
)

monoObj <- newCellDataSet(
    countsMono,
    phenoData = phenoDataMono,
    featureData = featureDataMono,
    expressionFamily = negbinomial.size()
)

monoObj <- estimateSizeFactors(monoObj)
monoObj <- estimateDispersions(monoObj, modelFormulaStr = "~orig.ident+group")

## Construct trajectory based on differential genes between Stem and TA in reference

varFeatures <- rownames(
    diffExprRefStemTA %>% filter((p_val_adj <= pAdjThreshold) & (abs(avg_log2FC) >= lfcThreshold))
)

monoObj <- setOrderingFilter(monoObj, varFeatures)

monoObj <- reduceDimension(
    monoObj,
    max_components = 2,
    reduction_method = 'DDRTree',
    residualModelFormulaStr = "~orig.ident+group+G2M_phase_UCell+S_phase_UCell",
    verbose = TRUE
)
monoObj <- orderCells(monoObj, reverse = FALSE)

saveRDS(
    monoObj,
    proot$find_file(
        config$outputpath, "data",
        paste0("monocle_", analysisLabel, "_", today, ".rds")
    )
)