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

# Load Seurat object, modify and subset and create SingleCellExperiment

SObj <- readRDS(proot$find_file(analparams$Seurat_rawall))

SObj[["is_doublet"]] <- ifelse(SObj$HTO_classification.global == "Doublet", TRUE, FALSE)

# Select only non-negative cells

cellsSelected <- rownames(SObj@meta.data)[!(SObj$HTO_classification.global == "Negative")]
SObj <- subset(SObj, cells = cellsSelected)

SObj <- NormalizeData(
  SObj,
  normalization.method = "LogNormalize",
  scale.factor = 1e4
)

feat <- rownames(SObj[["RNA"]]@data)

SObjSplit <- SplitObject(SObj, split.by = "orig.ident")

# Cluster cells for each sample (without doublets) and compute probabilities of doublet formation
# from same sample and same cluster

SObjSplitClustered <- list()
fractionsGroup <- list()
fractionsCluster <- list()
probSameGroup <- c()
probSameCluster <- c()
fractionHtoDoublets <- c()
fractionDetectableDoublets <- c()

for (sample in names(SObjSplit)) {
    cellsSelected <- rownames(SObjSplit[[sample]]@meta.data)[
        !(SObjSplit[[sample]]$HTO_classification.global == "Doublet")
    ]
    SObjSplitClustered[[sample]] <- subset(SObjSplit[[sample]], cells = cellsSelected)

    SObjSplitClustered[[sample]] <- SObjSplitClustered[[sample]] %>%
        FindVariableFeatures(nfeatures = analparams$scParam$VarFeat) %>%
        ScaleData() %>%
        RunPCA()

    if(!is.null(analparams$scParam$PCdims)) {
        PCdimsmax <- analparams$scParam$PCdims
    } else {
        PCdimsmax <- 30
        analparams$scParam$PCdims <- PCdimsmax
    }

    SObjSplitClustered[[sample]] <- SObjSplitClustered[[sample]] %>%
        RunUMAP(reduction = "pca", dims = 1:analparams$scParam$PCdims) %>%
        FindNeighbors(reduction = "pca", dims = 1:analparams$scParam$PCdims) %>%
        FindClusters(resolution = 0.5)

    fractionsGroup[[sample]] <- SObjSplitClustered[[sample]]@meta.data %>%
        dplyr::count(group) %>%
        mutate(fraction = prop.table(n)) %>%
        mutate(fraction_sqrd = fraction^2)
    fractionsCluster[[sample]] <- SObjSplitClustered[[sample]]@meta.data %>%
        dplyr::count(seurat_clusters) %>%
        mutate(fraction = prop.table(n)) %>%
        mutate(fraction_sqrd = fraction^2)

    probSameGroup[sample] <- sum(fractionsGroup[[sample]] %>% pull(fraction_sqrd))
    probSameCluster[sample] <- sum(fractionsCluster[[sample]] %>% pull(fraction_sqrd))

    fractionHtoDoublets[sample] <- SObjSplit[[sample]]@meta.data %>%
        dplyr::count(HTO_classification.global) %>%
        mutate(fraction = prop.table(n)) %>%
        filter(HTO_classification.global == "Doublet") %>%
        pull(fraction)

    fractionDetectableDoublets[sample] <- (fractionHtoDoublets[sample]/(1 - probSameGroup[sample]))*(1 - probSameCluster[sample])
}

# Add clusters to Seurat object with Doublets

for (sample in names(SObjSplit)) {
    SObjSplit[[sample]][["seurat_clusters"]] <- max(
        as.numeric(as.character(SObjSplitClustered[[sample]]@meta.data[["seurat_clusters"]]))
    ) + 1

    SObjSplit[[sample]]@meta.data[
        rownames(SObjSplitClustered[[sample]]@meta.data), "seurat_clusters"
    ] <- as.numeric(as.character(SObjSplitClustered[[sample]]$seurat_clusters))

    SObjSplit[[sample]]$seurat_clusters <- factor(SObjSplit[[sample]]$seurat_clusters)

    print(all.equal(
        as.character(SObjSplit[[sample]]$seurat_clusters[rownames(SObjSplitClustered[[sample]]@meta.data)]),
        as.character(SObjSplitClustered[[sample]]$seurat_clusters)
    ))
}

# Create SingleCellExperiment for each sample

scDatasets <- list()

for (sample in names(SObjSplit)) {
    counts <- GetAssayData(SObjSplit[[sample]], slot = "counts")
    scDatasets[[sample]] <- SingleCellExperiment(
        list("counts" = counts),
        colData = SObjSplit[[sample]]@meta.data
    )
}

# Run scDblFinder

for (sample in names(SObjSplit)) {
    print(Sys.time())
    print(paste0("Sample: ", sample))
    scDatasets[[sample]] <- scDblFinder(
        scDatasets[[sample]],
        clusters = "seurat_clusters",
        dbr = fractionDetectableDoublets[sample],
        nfeatures = 3000,
        dims = 30,
        includePCs = 15,
        iter = 3,
        knownDoublets = "is_doublet",
        knownUse = "discard",
        verbose = TRUE
    )
}

# Add data to Seurat object

scDblFinderCols <- colnames(colData(scDatasets[[1]]))[grep("^scDblFinder", colnames(colData(scDatasets[[1]])))]

for (sample in names(SObjSplit)) {
    SObj@meta.data[colnames(scDatasets[[sample]]), scDblFinderCols] <- as.data.frame(colData(scDatasets[[sample]]))[
       colnames(scDatasets[[sample]]), scDblFinderCols
    ]
}

classificationTable <- SObj@meta.data %>%
    dplyr::group_by(orig.ident, HTO_classification.global, scDblFinder.class) %>%
    summarise(mean_scDbl_score = mean(scDblFinder.score), n = n(), mean_lib_size = mean(nCount_RNA))

classificationTable

SObj@meta.data %>%
    dplyr::count(orig.ident, keep_filters & keep_doublets, scDblFinder.class)

saveRDS(scDatasets, proot$find_file(config$outputpath, "data", paste0("scDblFinder_cluster-dbr-custom-no-neg_", today, ".rds")))