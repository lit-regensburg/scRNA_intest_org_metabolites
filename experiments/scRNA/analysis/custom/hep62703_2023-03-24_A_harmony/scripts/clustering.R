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

set.seed(7)

config <- yaml.load_file(proot$find_file(params$config))
today <- params$date

print(today)

config <- yaml.load_file(proot$find_file(params$config))
analparams <- yaml.load_file(proot$find_file(params$analparams))

# libraries

libload <- c(
    "Seurat","ggplot2", "harmony", "Nebulosa", "future",
    "WriteXLS", "stringr", "dplyr"
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

# Load Seurat object and subset

SObj <- readRDS(proot$find_file(analparams$Seurat_rawall))

cellsSelected <- rownames(SObj@meta.data)[SObj@meta.data[["keep"]]]
SObj <- subset(SObj, cells = cellsSelected)

# Normalization, scaling, PCA

SObj <- NormalizeData(
  SObj,
  normalization.method = "LogNormalize",
  scale.factor = 1e6
)

SObj[["RNA"]]@data[1:10, 1:10]
feat<-rownames(SObj[["RNA"]]@data)
SObj <- SObj %>%
  FindVariableFeatures(nfeatures = analparams$scParam$VarFeat) %>%
  ScaleData(features=feat) %>% 
  RunPCA()

if(!is.null(analparams$scParam$PCdims)) {
    PCdimsmax <- analparams$scParam$PCdims
} else {
    PCdimsmax <- 30
    analparams$scParam$PCdims <- PCdimsmax
}

ElbowPlot(SObj,ndims=PCdimsmax)

# Harmony integration

if (is.null(analparams$scParam$integration)) {
    reduction <- "pca"
} else if (analparams$scParam$integration == "harmony") {
    SObj <- RunHarmony(SObj, group.by.vars = "orig.ident")
    reduction <- "harmony"
} else {
    stop("Integration parameter not supported")
}

# Clustering

clusteringResolutions <- c(0.3, 0.5, 0.8, 1.0, 1.2, 1.6)
names(clusteringResolutions) <- as.character(clusteringResolutions)

SObj <- FindNeighbors(SObj, reduction = reduction, dims = 1:analparams$scParam$PCdims)
SObj <- RunUMAP(SObj, reduction = reduction, dims = 1:analparams$scParam$PCdims)

for (res in names(clusteringResolutions)) {
    SObj <- FindClusters(
        SObj,
        resolution = clusteringResolutions[res]
    )
}
SObj$seurat_clusters <- SObj$RNA_snn_res.0.8
Idents(SObj) <- SObj$seurat_clusters

# Save Seurat object and config

SObjName <- paste0(
    "SeuratObj_clustered_", today, ".rds"
)

analparams[["Seurat_clustered"]] <- file.path(config$outputpath, "data", SObjName)
analparams$last_modified_time <- as.character(Sys.time())
analparams$last_modified_entry <- "Seurat_clustered"

saveRDS(SObj, file = proot$find_file(analparams[["Seurat_clustered"]]))

write_yaml(
    analparams,
    file=proot$find_file(params$analparams)
)