### Libraries

# Using rbioc3-14

library(rprojroot)
library(yaml)
library(Seurat)
library(SingleR)
library(stringr)
library(ggplot2)
library(dplyr)
library(WriteXLS)
library(future)

clusterCol <- "RNA_snn_res.0.8"
refName <- "Haber_2017"

### Paths

root_criterion <- ".projroot"
proot <- rprojroot::has_file(root_criterion)

source(
    proot$find_file("experiments/scRNA/analysis/custom/functions/custom_functions.R")
)

params <- list()
params$config <- proot$find_file("experiments/scRNA/analysis/custom/hep62703_2023-03-24_A_harmony/config.yaml")
params$analparams <- proot$find_file("experiments/scRNA/analysis/custom/hep62703_2023-03-24_A_harmony/analparams_preprocessing_filters_2023-03-24.yaml")
params$date <- Sys.Date()

set.seed(7)

config <- yaml.load_file(proot$find_file(params$config))
today <- params$date

print(today)

config <- yaml.load_file(proot$find_file(params$config))
analparams <- yaml.load_file(proot$find_file(params$analparams))

# Parallelization

# Increase max size, otherwise it is exceeded for these data
options(future.globals.maxSize = 1200 * 1024^2)

future.seed=TRUE
plan("multicore", workers = 6)
plan()

### Load data

SObj <- readRDS(proot$find_file(analparams$Seurat_clustered))
scRef <- read.table(
    proot$find_file("experiments/scRNA/reference_data/GSE92332_AtlasFullLength_TPM.txt"),
    sep = "\t"
)
scRefLog <- log1p(scRef)

# Extract labels from reference and simplify

 refLabels <- do.call(
    rbind,
    strsplit(colnames(scRefLog), "_(?!.*_)", perl=TRUE)
 )[, 2]

refLabels[grep("^Tuft", refLabels)] <- "Tuft"
refLabels[grep("^Enterocyte\\.Progenitor", refLabels)] <- "Enterocyte.Progenitor"

print(unique(refLabels))

# Get counts from Seurat object

scCountsLog <- as.matrix(GetAssayData(SObj, slot = "data"))

# Run SingleR

labelsSingleR <- SingleR(
    test=scCountsLog,
    ref=scRefLog,
    labels=refLabels,
    de.method="wilcox"
)

labelsSingleRClusters <- SingleR(
    test=scCountsLog,
    ref=scRefLog,
    labels=refLabels,
    clusters=SObj@meta.data[[clusterCol]],
    de.method="wilcox"
)

# Add annotation to Seurat object

annotColClusters <- paste0("annot_SingleR_", clusterCol)

SObj[["annot_SingleR"]] <- labelsSingleR$labels
SObj[[annotColClusters]] <- NA

for (cluster in levels(SObj@meta.data[[clusterCol]])) {
    SObj@meta.data[SObj@meta.data[[clusterCol]] == cluster, annotColClusters] <- labelsSingleRClusters[cluster, "labels"]
}

# Combine Golbet and Paneth cells as Paneth cells do hardly occur

SObj$annot_SingleR_mod <- SObj$annot_SingleR
SObj$annot_SingleR_mod[SObj$annot_SingleR_mod %in% c("Goblet", "Paneth")] <- "Goblet.Paneth"

# Save data

saveRDS(
    labelsSingleR,
    proot$find_file(
        config$outputpath, "data",
        paste0("annot_SingleR_", refName, ".rds")
    )
)

saveRDS(
    labelsSingleRClusters,
    proot$find_file(
        config$outputpath, "data",
        paste0("annot_SingleR_", refName, "_", clusterCol, ".rds")
    )
)

saveRDS(SObj, file = proot$find_file(analparams[["Seurat_clustered"]]))