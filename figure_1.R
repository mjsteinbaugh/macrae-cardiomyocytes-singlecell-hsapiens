# Run this script directly on O2.
# R library: /home/mjs65/R/library/3.5-bioc-release-20181022
# Only Seurat is used for calculation below.

library(readr)  # v1.1.1
library(basejump)  # v0.7.3
library(pointillism)  # v0.1.3
library(Seurat)  # v2.3.4

theme_set(theme_paperwhite())

load("data/2017-12-12/seurat.rda")

data_dir <- file.path("data", Sys.Date())
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

results_dir <- file.path("results", Sys.Date())
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

width <- 10
height <- 10

# Subset and re-cluster wild-type samples ======================================
# Guide for cell selection
# https://satijalab.org/seurat/interaction_vignette.html
# `WhichCells()`, `SubsetData()`
colnames(seurat@meta.data)
table(seurat@meta.data$genotype)

# Subset the wild-type samples (M1, M7).
# n = 1784
wildtype <- seurat %>%
    SubsetData(
        subset.name = "genotype",
        accept.value = "wildtype",
        do.clean = TRUE
    ) %>%
    NormalizeData(
        normalization.method = "LogNormalize",
        scale.factor = 10000
    ) %>%
    FindVariableGenes(
        mean.function = ExpMean,
        dispersion.function = LogVMR,
        do.plot = FALSE
    ) %>%
    ScaleData(
        model.use = "linear",
        vars.to.regress = c("nUMI", "S.Score", "G2M.Score")
    ) %>%
    RunPCA(do.print = FALSE) %>%
    ProjectPCA(do.print = FALSE)

Seurat::PCElbowPlot(wildtype)
dims_use <- pointillism::plotPCElbow(wildtype)
# Using 15 PCs.

wildtype <- wildtype %>%
    FindClusters(
        dims.use = dims_use,
        resolution = c(0.4, 0.6, 0.8),
        save.SNN = TRUE,
        force.recalc = TRUE
    ) %>%
    RunTSNE(
        dims.use = dims_use,
        do.fast = TRUE
    )
# Set the cluster identities to the resolution 0.6 mappings.
wildtype <- SetAllIdent(wildtype, id = "res.0.6")
saveData(wildtype, dir = data_dir)

# t-SNE plots ==================================================================
# Default Seurat t-SNE.
tsne1 <- Seurat::DimPlot(wildtype, reduction.use = "tsne")
ggsave(
    filename = file.path(results_dir, "wildtype_seurat_tsne.pdf"),
    width = width,
    height = height
)

# Friendlier ggplot return that's easier to customize.
tsne2 <- pointillism::plotTSNE(wildtype)
ggsave(
    filename = file.path(results_dir, "wildtype_custom_tsne.pdf"),
    width = width,
    height = height
)
saveData(tsne1, tsne2, dir = data_dir)

# Marker plots =================================================================
marker_tsne_dir <- file.path(results_dir, "marker_tsne")
dir.create(marker_tsne_dir, recursive = TRUE, showWarnings = FALSE)

dot_plot_dir <- file.path(results_dir, "dot_plot")
dir.create(dot_plot_dir, recursive = TRUE, showWarnings = FALSE)

dot_colors <- c("orange", "purple")
tsne_colors <- c("gray75", "purple")

# MYH6, SHOX2 ------------------------------------------------------------------
genes <- c("MYH6", "SHOX2")
Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(
        dot_plot_dir,
        paste0(paste(genes, collapse = "_"), ".pdf")
    ),
    width = 4,
    height = 4
)

Seurat::FeaturePlot(
    object = wildtype,
    features.plot = genes,
    cols.use = tsne_colors
)
ggsave(
    filename = file.path(
        marker_tsne_dir,
        paste0(paste(genes, collapse = "_"), ".pdf")
    ),
    width = 10,
    height = 5
)

# Bulk RNA-seq DEG -------------------------------------------------------------
stem <- "bulk_rnaseq_deg"
input <- read_lines(
    file = file.path(
        "genes",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# [1] "PWAR6"   "PWAR5"   "MIR34-A" "CACNA1S" "GABRR1"  "KCNK5"  
# [7] "CLCN1"   "MMEL1"

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 10,
    height = 4
)

# Cardiomyocyte markers > ion channels -----------------------------------------
stem <- "ion_channels"
input <- read_lines(
    file = file.path(
        "genes",
        "cardiomyocyte_markers",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# [1] "KCNK9" "KCNA6" "RYR1"

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 8,
    height = 4
)

# Cardiomyocyte markers > protein kinases --------------------------------------
stem <- "protein_kinases"
input <- read_lines(
    file = file.path(
        "genes",
        "cardiomyocyte_markers",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# character(0)

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 10,
    height = 4
)

# Cardiomyocyte markers > structural proteins ----------------------------------
stem <- "structural_proteins"
input <- read_lines(
    file = file.path(
        "genes",
        "cardiomyocyte_markers",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# [1] "CTNT"

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 8,
    height = 4
)

# Cardiomyocyte markers > transcription factors --------------------------------
stem <- "transcription_factors"
input <- read_lines(
    file = file.path(
        "genes",
        "cardiomyocyte_markers",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# [1] "NKX2.5"  "PITX2C"  "NODAL"   "LEFTY-1" "LEFTY-2" "NKX3.2" 
# [7] "FGF-8"

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 6,
    height = 4
)

# Negative fraction > endocardial lineage --------------------------------------
stem <- "endocardial_lineage"
input <- read_lines(
    file = file.path(
        "genes",
        "negative_fraction",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# [1] "NFATC11" "TIE2"

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 4,
    height = 4
)

# Negative fraction > epicardial lineage ---------------------------------------
stem <- "epicardial_lineage"
input <- read_lines(
    file = file.path(
        "genes",
        "negative_fraction",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# [1] "WT-1"   "ZO1"    "SNAIL2"

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 4,
    height = 4
)

# Negative fraction > fibroblsts -----------------------------------------------
stem <- "fibroblasts"
input <- read_lines(
    file = file.path(
        "genes",
        "negative_fraction",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# [1] "WT-1"

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 4,
    height = 4
)

# Negative fraction > retinoic acid pathway ------------------------------------
stem <- "retinoic_acid_pathway"
input <- read_lines(
    file = file.path(
        "genes",
        "negative_fraction",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# [1] "ADH1A"   "CYP26C1" "HOXB1"   "WNT8A"

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 8,
    height = 4
)

# Pacemaker population > SA node cells -----------------------------------------
stem <- "SA_node_cells"
input <- read_lines(
    file = file.path(
        "genes",
        "pacemaker_population",
        paste0(stem, ".txt")
    )
)
genes <- input %>%
    .[. %in% rownames(wildtype@data)] %>%
    unique() %>%
    sort(decreasing = TRUE)

# Return missing symbols.
setdiff(input, genes)
# [1] "ADH1A"   "CYP26C1" "HOXB1"   "WNT8A"

Seurat::DotPlot(
    object = wildtype,
    genes.plot = genes,
    cols.use = dot_colors,
    plot.legend = TRUE
)
ggsave(
    filename = file.path(dot_plot_dir, paste0(stem, ".pdf")),
    width = 6,
    height = 4
)

# Session info =================================================================
sessionInfo <- sessioninfo::session_info()
saveData(sessionInfo, dir = data_dir)
