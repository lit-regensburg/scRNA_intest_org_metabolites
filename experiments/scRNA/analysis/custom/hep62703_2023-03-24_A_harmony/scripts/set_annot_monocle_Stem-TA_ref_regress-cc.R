### Setup

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

config <- yaml.load_file(proot$find_file(params$config))
analparams <- yaml.load_file(proot$find_file(params$analparams))
today <- params$date

# Set params

annotCol <- "annot_SingleR_mod"
annotColMonocle <- "annot_SingleR_mod_monocle_Stem.TA_ref_regress.cc"

pathMonoObj <- proot$find_file(
    config$outputpath, "data", "monocle_Stem-TA_ref_regress-cc_2024-03-20.rds"
)

# Random seed

set.seed(7)

# Libraries

libload <- c(
    "Seurat","ggplot2", "harmony", "Nebulosa", "future",
    "WriteXLS", "stringr", "dplyr", "monocle"
)

for(i in libload) {
    library(i, character.only = T)
}

## Parallelization

# Increase max size, otherwise it is exceeded for these data
options(future.globals.maxSize = 1200 * 1024^2)

future.seed=TRUE
plan("multicore", workers = 6)
plan()

### Load complete Seurat object and monocle object

SObj <- readRDS(proot$find_file(analparams$Seurat_clustered))
monoObj <- readRDS(pathMonoObj)

# Redo cell annotation with new clusters

annotationMap <- c(
    "1" = "Stem",
    "2" = "Stem.TA",
    "3" = "Stem.TA",
    "4" = "Stem.TA",
    "5" = "Stem.TA",
    "6" = "TA",
    "7" = "TA",
    "8" = "TA",
    "9" = "TA",
    "10" = "Stem.TA",
    "11" = "Stem.TA"
)

pData(monoObj)[["annot"]] <- sapply(
    pData(monoObj)[["State"]],
    function(x) annotationMap[as.character(x)]
)

# Add new annotation to complete Seurat object

SObj[[annotColMonocle]] <- SObj[[annotCol]]
SObj@meta.data[rownames(pData(monoObj)), annotColMonocle] <- pData(monoObj)[["annot"]]

## Save Seurat object

saveRDS(
    SObj,
    proot$find_file(analparams$Seurat_clustered)
)

analparams$last_modified_time <- format(Sys.time())
analparams$last_modified_entry <- "Seurat_clustered"

write_yaml(
    analparams,
    proot$find_file(params$analparams)
)