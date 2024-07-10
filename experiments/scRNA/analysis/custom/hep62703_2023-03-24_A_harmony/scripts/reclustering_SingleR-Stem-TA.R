### Setup

#rbioc_3-14
library(rprojroot)
library(yaml)
library(knitr)

root_criterion <- ".projroot"
proot <- rprojroot::has_file(root_criterion)

identCol <- "annot_SingleR_mod"
identSubset <- c("Stem", "TA") # Stem and TA from SingleR annotation

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

## libraries

libload <- c(
    "Seurat","ggplot2", "harmony", "Nebulosa", "future",
    "WriteXLS", "stringr", "dplyr"
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

### Load Seurat object and subset

SObjAll <- readRDS(proot$find_file(analparams$Seurat_clustered))

Idents(SObjAll) <- SObjAll[[identCol]]
SObj <- subset(SObjAll, idents = identSubset)

### Marker genes from differential expression between Stem and TA in reference (Haber et al., 2017)

scRef <- read.table(
    proot$find_file("experiments/scRNA/reference_data/GSE92332_AtlasFullLength_TPM.txt"),
    sep = "\t"
)

# Extract labels from reference and simplify

refLabels <- do.call(
    rbind,
    strsplit(colnames(scRef), "_(?!.*_)", perl=TRUE)
 )[, 2]

refLabels[grep("^Tuft", refLabels)] <- "Tuft"
refLabels[grep("^Enterocyte\\.Progenitor", refLabels)] <- "Enterocyte.Progenitor"

print(unique(refLabels))

# Create Seurat object from TPM data

scRefMetadata <- data.frame(
    "cell_id" = colnames(scRef),
    "celltype" = refLabels,
    row.names = "cell_id"
)

SObjRef <- CreateSeuratObject(
    counts = scRef,
    project = "sc_ref_Haber_2017",
    assay = "RNA",
    names.field = 1,
    names.delim = "_",
    meta.data = scRefMetadata
)
Idents(SObjRef) <- "celltype"

SObjRef <- NormalizeData(
  SObjRef,
  normalization.method = "LogNormalize",
  scale.factor = 1e6
)

diffExprRefStemTA <- FindMarkers(
    SObjRef,
    ident.1 = "Stem",
    ident.2 = "TA",
    test.use = "wilcox",
    only.pos = FALSE,
    min.pct = 0.1,
    logfc.threshold = -Inf,
    return.thresh = 1,
)
markersRefStemTA <- diffExprRefStemTA %>% filter((p_val_adj <= 0.05) & (abs(avg_log2FC) >= 0.25))

## Normalization, scaling, PCA

SObj <- NormalizeData(
  SObj,
  normalization.method = "LogNormalize",
  scale.factor = 1e6
)

feat <- rownames(SObj[["RNA"]]@data)

featPCA <- rownames(markersRefStemTA)[
    rownames(markersRefStemTA) %in% feat
]

SObj <- ScaleData(SObj, features = feat)
SObj <- RunPCA(SObj, features = featPCA)

if(!is.null(analparams$scParam$PCdims)) {
    PCdimsmax <- analparams$scParam$PCdims
} else {
    PCdimsmax <- 30
    analparams$scParam$PCdims <- PCdimsmax
}

ElbowPlot(SObj,ndims=PCdimsmax)

## Harmony integration

if (is.null(analparams$scParam$integration)) {
    reduction <- "pca"
} else if (analparams$scParam$integration == "harmony") {
    SObj <- RunHarmony(SObj, group.by.vars = "orig.ident")
    reduction <- "harmony"
} else {
    stop("Integration parameter not supported")
}

## Clustering

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

SObj[["seurat_clusters"]] <- SObj[["RNA_snn_res.0.5"]]

DimPlot(SObj, group.by="RNA_snn_res.0.5")
DimPlot(SObj, group.by=identCol)

## Save Seurat object subset

saveRDS(
    SObj,
    proot$find_file(
        config$outputpath,
        "data",
        paste0("SeuratObj_reclustered_SingleR-Stem-TA_", today, ".rds")
    )
)

# Save reference marker genes

write.csv(
    diffExprRefStemTA,
    proot$find_file(config$outputpath, "table", "diff_expr_Stem-TA_ref-Haber.csv"),
    row.names = TRUE,
    quote = FALSE
)

# plt <- DimPlot(SObj, group.by = "RNA_snn_res.0.3")

# ggsave(
#     proot$find_file(config$outputpath, "fig", "tmp", "reclustering_res-03.png"),
#     plt
# )