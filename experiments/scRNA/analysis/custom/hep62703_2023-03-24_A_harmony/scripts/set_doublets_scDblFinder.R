# THIS MUST BE RUN AFTER QC REPORT WITH FILTERS!

# Setup

#rbioc_3-14
library(rprojroot)
library(yaml)
library(knitr)

root_criterion <- ".projroot"
proot <- rprojroot::has_file(root_criterion)

params <- list()
params$config <- proot$find_file("experiments/scRNA/analysis/custom/hep62703_2023-03-24_A_harmony/config.yaml")
params$analparams <- proot$find_file("experiments/scRNA/analysis/custom/hep62703_2023-03-24_A_harmony/analparams_preprocessing_filters_2023-03-24.yaml")
params$date <- Sys.Date()

commonSeed <- 7

set.seed(commonSeed)

config <- yaml.load_file(proot$find_file(params$config))
analparams <- yaml.load_file(proot$find_file(params$analparams))

# libraries

libload <- c(
    "Seurat","ggplot2", "harmony", "Nebulosa", "future",
    "WriteXLS", "stringr", "dplyr", "BiocParallel", "SingleCellExperiment",
    "scDblFinder"
)

for(i in libload) {
    library(i, character.only = T)
}

# Parallelization

# Increase max size, otherwise it is exceeded for these data
options(future.globals.maxSize = 1200 * 1024^2)

future.seed=TRUE
plan("multicore", workers = 6)
plan()

# Load Seurat object and set scDblFinder results

SObj <- readRDS(proot$find_file(analparams$Seurat_rawall))

scDatasets <- readRDS(
    proot$find_file(config$outputpath, "data", paste0("scDblFinder_cluster-dbr-custom-no-neg_2023-03-29.rds"))
)

scDblFinderCols <- colnames(colData(scDatasets[[1]]))[grep("^scDblFinder", colnames(colData(scDatasets[[1]])))]

for (sample in names(scDatasets)) {
    SObj@meta.data[colnames(scDatasets[[sample]]), scDblFinderCols] <- as.data.frame(colData(scDatasets[[sample]]))[
       colnames(scDatasets[[sample]]), scDblFinderCols
    ]
}

SObj$keep_doublets_HTO <- SObj$keep_doublets
SObj$keep_doublets_comp <- ifelse(SObj$scDblFinder.class == "singlet", TRUE, FALSE)
SObj$keep_doublets_comp[is.na(SObj$keep_doublets_comp)] <- FALSE
SObj$keep_doublets <- SObj$keep_doublets_HTO & SObj$keep_doublets_comp
SObj$keep <- SObj$keep_doublets & SObj$keep_filters & SObj$keep_libs

# Save Seurat object

saveRDS(SObj, proot$find_file(analparams$Seurat_rawall))

analparams$doublet_data <- file.path(config$outputpath, "data", paste0("scDblFinder_cluster-dbr-custom-no-neg_2023-03-29.rds"))
analparams$last_modified_time <- format(Sys.time())
analparams$last_modified_entry <- "Seurat_rawall"

write_yaml(
    analparams,
    proot$find_file(params$analparams)
)