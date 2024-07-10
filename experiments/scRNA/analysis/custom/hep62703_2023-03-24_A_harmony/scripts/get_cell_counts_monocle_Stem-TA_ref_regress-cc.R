### Libraries

# Using rbioc3-14

library(rprojroot)
library(yaml)
library(Seurat)
library(SingleR)
library(stringr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(WriteXLS)
library(future)

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

# Params

identCol <- "annot_SingleR_mod_monocle_Stem.TA_ref_regress.cc"
relevantCelltypes <- c("Stem", "Stem.TA", "TA", "Enterocyte.Progenitor", "Enterocyte", "Endocrine", "Goblet.Paneth", "Tuft")

# Load data

SObj <- readRDS(proot$find_file(analparams$Seurat_clustered))

# Get celltype counts with metadata

metadata <- SObj@meta.data

countsLong <- metadata %>%
    dplyr::count(group, orig.ident, .data[[identCol]]) %>%
    dplyr::rename(celltype = .data[[identCol]]) %>%
    filter(celltype %in% relevantCelltypes) %>%
    mutate(sample = paste0(group, "_", orig.ident))

groups <- unique(countsLong$group)

countsLong <- countsLong %>%
    mutate(celltype = factor(celltype, levels = relevantCelltypes)) %>%
    mutate(group = factor(group, levels = c(groups[4:6], groups[1:3]))) %>%
    mutate(orig.ident = factor(orig.ident, levels = sort(unique(orig.ident)))) %>%
    arrange(group, orig.ident, celltype)

countsPivot <- pivot_wider(
    countsLong,
    id_cols = c("group", "orig.ident"),
    names_from = "celltype",
    values_from = "n"
)
countsPivot[is.na(countsPivot)] <- 0

# Save data

write.csv(
    countsPivot,
    proot$find_file(config$outputpath, "data", paste0("celltype_counts_monocle_Stem-TA_ref_regress-cc_", today, ".csv")),
    quote = FALSE,
    row.names = FALSE
)