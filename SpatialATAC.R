rm(list=ls())
setwd("/Users/ruiqil/Documents/Spatial_ATAC-seq")

# load packages
library(ArchR)
library(Seurat)
library(grid)

library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)

library(ggplot2)
library(patchwork)
library(dplyr)

library(SeuratDisk)

threads = 8

addArchRThreads(threads = threads)
addArchRGenome("mm10")

# ME13_50um
inputFiles <- '/Users/ruiqil/Downloads/GSM5238386_ME13_50um.fragments.tsv.gz'
sampleNames <- 'ME13'
## Create ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0,
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000)
)
ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE
)
proj

## Select pixels in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)
data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"
image <- Read10X_Image(image.dir = file.path("/Users/ruiqil/Downloads/GSM5238386_ME13_50um_spatial"), 
                       filter.matrix = filter.matrix)

### rename
rownms = rownames(meta.data)
for (i in 1:length(rownms)) {
  rownms[i] = strsplit(strsplit(rownms[i], "#")[[1]][2], "-")[[1]][1]
}
rownames(meta.data) = rownms

meta.data.spatial <- meta.data[row.names(image@coordinates), ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]
proj_in_tissue

## Data normalization and dimensionality reduction 
proj_in_tissue <- addIterativeLSI(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

proj_in_tissue <- addClusters(
  input = proj_in_tissue,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)

proj_in_tissue <- addUMAP(
  ArchRProj = proj_in_tissue, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

plotEmbedding(ArchRProj = proj_in_tissue, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 1.5)

proj_in_tissue <- addImputeWeights(proj_in_tissue)


## Identify the marker genes for each cluster 
markersGS <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  testMethod = "wilcoxon"
)

markerList_pos <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")

markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}
markerGenes <- unlist(markerGenes)


## Call peaks
proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "Clusters")

pathToMacs2 <- "/Users/ruiqil/miniconda3/bin/Macs2"

proj_in_tissue <- addReproduciblePeakSet(
  ArchRProj = proj_in_tissue, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

proj_in_tissue <- addPeakMatrix(proj_in_tissue)

getAvailableMatrices(proj_in_tissue)
getPeakSet(proj_in_tissue)

if("Motif" %ni% names(proj_in_tissue@peakAnnotation)){
  proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)
}

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_in_tissue,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.1")


## ChromVAR Deviatons Enrichment
proj_in_tissue <- addBgdPeaks(proj_in_tissue, force = TRUE)

proj_in_tissue <- addDeviationsMatrix(
  ArchRProj = proj_in_tissue, 
  peakAnnotation = "Motif",
  force = TRUE
)

markersMotifs <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "MotifMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = 'z'
)

source('Data_visualization/getGeneScore_ArchR.R')
source('Data_visualization/SpatialPlot_new.R')
## Create Seurat object 
proj_in_tissue <- addImputeWeights(proj_in_tissue)
gene_score <- getGeneScore_ArchR(ArchRProj = proj_in_tissue, name = markerGenes, imputeWeights = getImputeWeights(proj_in_tissue))
# library(anndata)
# ad <- AnnData(X = t(gene_score))
# write_h5ad(
#   ad,
#   "genescore.h5ad"
# )
data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = gene_score, assay = assay, meta.data = meta.data)
image <- Read10X_Image(image.dir = file.path("/Users/ruiqil/Downloads/GSM5238386_ME13_50um_spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object

spatial.obj@meta.data[["in_tissue"]] <- object@images[["slice1"]]@coordinates[["tissue"]]
spatial.obj@meta.data[["array_row"]] <- object@images[["slice1"]]@coordinates[["row"]]
spatial.obj@meta.data[["array_col"]] <- object@images[["slice1"]]@coordinates[["col"]]
SaveH5Seurat(spatial.obj, filename = "data.h5Seurat", overwrite = TRUE)
Convert("data.h5Seurat", dest = "h5ad", overwrite = TRUE)

## Plot results
p1 <- SpatialPlot_new(spatial.obj, features = "nFrags",  pt.size.factor = 3, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1 + scale_x_reverse()

# n_clusters <- length(unique(projCUTA$Clusters))
# cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
# names(cols) <- paste0('C', seq_len(n_clusters))
# cols

# p <- SpatialPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'Clusters', pt.size.factor = 3, cols = cols, image.alpha = 0, stroke = 0)
# p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
# p

# feature <- 'Rarg'
# p <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 3, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
#   theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
# p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
# p + scale_x_reverse()

feature <- 'TSSEnrichment'
p2 <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 3, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p2$layers[[1]]$aes_params <- c(p2$layers[[1]]$aes_params, shape=22)
p2 + scale_x_reverse()

feature <- 'PromoterRatio'
p3 <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 3, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p3$layers[[1]]$aes_params <- c(p3$layers[[1]]$aes_params, shape=22)
p3 + scale_x_reverse()

feature <- 'NucleosomeRatio'
p4 <- SpatialPlot_new(spatial.obj, features = feature, pt.size.factor = 3, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q90") +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p4$layers[[1]]$aes_params <- c(p4$layers[[1]]$aes_params, shape=22)
p4 + scale_x_reverse()
library(gridExtra)
grid.arrange(p1 + scale_x_reverse(), p2 + scale_x_reverse(), p3 + scale_x_reverse(), p4 + scale_x_reverse(), nrow=2)



PlotClusters = data.frame(rownames = proj_in_tissue@cellColData@rownames, Clusters = proj_in_tissue@cellColData@listData[["Clusters"]])
ref = data.frame(order = 1:length(proj_in_tissue@cellColData@rownames), rownames = spatial.obj@meta.data[["cellID_archr"]])
PlotClusters = merge(ref, PlotClusters, by = 'rownames')
PlotClusters = PlotClusters[order(PlotClusters$order),]
spatial.obj@meta.data[["Clusters"]] = PlotClusters$Clusters

n_clusters <- length(unique(PlotClusters$Clusters))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
names(cols) <- paste0('C', seq_len(n_clusters))
cols
p <- SpatialPlot_new(spatial.obj, label.size = 5, group.by = "Clusters", pt.size.factor = 4, cols = cols, image.alpha = 0., stroke = 0)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p + scale_x_reverse()

p <- SpatialPlot_new(spatial.obj, label.size = 5, group.by = "Clusters", pt.size.factor = 2, cols = cols, image.alpha = 1, stroke = 0)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p

