# Top 20 markers per cluster.
# Also including dot plots of these genes.
# Michael Steinbaugh
# 2018-11-14

# Clusters 0,2,3 are probably CMs and the rest are non-CMs.

library(basejump)  # v0.8.1
library(bcbioSingleCell)  # v0.3.1
library(Seurat)  # v2.3.4
library(pointillism)  # v0.2.0
library(tidyverse)

options("basejump.save.dir" = file.path("data", Sys.Date()))
theme_set(theme_paperwhite())

# Load the seurat object.
loadDataAsName(seurat = wildtype, dir = "data/2018-10-22")

# Calculate the markers.
all_markers <- FindAllMarkers(seurat)
saveData(all_markers)

# Load the corresponding bcbioSingleCell object to get metadata.
loadDataAsName(bcb = bcb_pool_filtered, dir = "data/2017-12-01")
x <- seurat %>% .@data %>% colnames()
y <- bcb %>% .@assays %>% .[[1]] %>% colnames()
stopifnot(all(x %in% y))

# bcbioSingleCell v0.0.24
# Homo sapiens Ensembl 90.
bcb %>% .@metadata %>% .$version
bcb %>% .@metadata %>% .$organism
bcb %>% .@metadata %>% .$ensemblVersion
# We didn't use GenomicRanges in this version of the analysis, but we can
# pull the matching ones from Ensembl.
rowData <- bcb %>% .@elementMetadata

rowRanges <- makeGRangesFromEnsembl("Homo sapiens", release = 90L)
saveData(rowRanges)

# Check that the IDs match.
stopifnot(length(setdiff(rownames(rowData),names(rowRanges))) == 0)

# Need to remap the names to match the gene symbols. I may simplify this step
# in the `SeuratMarkers()` call in a future update.
ranges <- rowRanges
names(ranges) <- make.unique(as.character(ranges$geneName))

# Check that we can slot these back into the seurat object correctly.
stopifnot(all(rownames(seurat) %in% names(ranges)))
stopifnot(all(all_markers$gene %in% names(ranges)))
ranges <- ranges[rownames(seurat)]
seurat@misc$rowRanges <- ranges
assignAndSaveData(name = "wildtype_resave_with_rowRanges", object = seurat)

plotTSNE(seurat)

seurat_markers <- SeuratMarkersPerCluster(object = all_markers, ranges = ranges)
saveData(seurat_markers)

top_markers <- topMarkers(seurat_markers, n = 20L)
saveData(top_markers)

export(x = top_markers, file = "top_20_markers_per_cluster.csv")

# Loop across the clusters and generate dot plots.
invisible(lapply(
    X = levels(clusterID(seurat)),
    FUN = function(cluster) {
        genes <- top_markers %>% filter(cluster == !!cluster) %>% pull(name)
        plotDot(
            object = seurat,
            genes = genes,
            perSample = TRUE,
            color = ggplot2::scale_colour_gradient2(
                low = "orange",
                mid = "gray75",
                midpoint = 0L,
                high = "purple"
            )
        )
        ggsave(
            filename = paste0("wildtype_plot_dot_", cluster, ".pdf"),
            plot = last_plot(),
            width = 12L,
            height = 4L
        )
    }
))





# Everything else
