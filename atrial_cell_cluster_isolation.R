library(bcbioSingleCell)
library(tidyverse)

dataDir <- file.path("data", Sys.Date())

load("data/2017-12-01/bcb_pool_filtered.rda")
load("data/2017-12-12/seurat.rda")

# Replot the tSNE
plotTSNE(seurat, pointSize = 0.5, pointAlpha = 0.8)

# Check that the dimensions match ==============================================
identical(
    dim(bcb_pool_filtered),
    dim(seurat@data)
)
# TRUE

identical(
    colnames(bcb_pool_filtered),
    colnames(seurat@data)
)
# TRUE

# Note that Seurat uses the gene symbols, with the original Ensembl ids stashed
# as the names of the rownames
identical(
    rownames(bcb_pool_filtered),
    names(rownames(seurat@data))
)
# TRUE

# Import incoming markers (2018-01-30 list) ====================================
library(basejump)
library(readxl)
xlsx_file <- file.path("meta", "markers_2018-01-30.xlsx")
sheets <- excel_sheets(xlsx_file)
# [1] "cardiomyocyte"          "left_right_pattern"     "epicardial_lineage"    
# [4] "retinoic_acid_pathway"  "sheer_stress"           "atp_electron_transport"
# [7] "bulk_rnaseq_de"  
new_markers <- lapply(
    X = sheets,
    FUN = function(sheet) {
        read_xlsx(xlsx_file, sheet = sheet) %>%
            pull(symbol)
    }
)
names(new_markers) <- sheets
saveData(new_markers, dir = dataDir)

# Now we need to check that all of these marker symbols map to Ensembl gene
# identifiers, and are present in our dataset. For missing symbols, let's refer
# to our list of synonyms downloaded from Ensembl, saved in the basejump
# package.
# Here `annotable` contains the information for genes we detected, and
# `full_annotable` contains all genes from Ensembl
annotable <- annotable(bcb_pool_filtered)
full_annotable <- metadata(bcb_pool_filtered)$annotable
synonyms <- basejump::synonyms$homoSapiens
saveData(full_annotable, synonyms, dir = dataDir)

seurat_gene_ids <- rownames(seurat@data)
undetected_markers <- lapply(
    X = new_markers,
    FUN = function(x) {
        x[which(!x %in% seurat_gene_ids)]
    }
)
# 2018-01-31 unmapped symbols we need to update
# $cardiomyocyte
# [1] "CTNT"   "NKX2.5" "SCN5a"  "VIMP"  
# 
# $left_right_pattern
# [1] "Pitx2c"  "nodal"   "Lefty-1" "Lefty-2" "NKX3.2"  "FGF-8"  
# 
# $epicardial_lineage
# [1] "WT-1"   "ZO1"    "SNAIL2"
# 
# $retinoic_acid_pathway
# [1] "ADH1A"   "CYP26C1" "HOXB1"   "WNT8a"   "WNT5a"  
# 
# $sheer_stress
# [1] "IP3R"
# 
# $atp_electron_transport
# character(0)
# 
# $bulk_rnaseq_de
# [1] "PWAR6"   "PWAR5"   "mir34-a" "CACNA1S" "GABRR1"  "KCNK5"   "CLCN1"   "mmel1"
saveData(undetected_markers, dir = dataDir)

# Ensure all of the undetected markers symbols map correctly to the Ensembl
# gene identifiers, but that they're just not expressed.
lapply(
    X = undetected_markers,
    FUN = function(x) {
        all(x %in% full_annotable[, "symbol"])
    })

new_markers_detected <- lapply(
    new_markers,
    FUN = function(x) {
        # Keep the input order Zaniar used in the original spreadsheet
        x[which(x %in% seurat_gene_ids)]
    }
)
saveData(new_markers_detected, dir = dataDir)

# New marker plots =============================================================
plotsDir <- file.path("results", "marker_plots", Sys.Date())
dir.create(plotsDir, recursive = TRUE)
# Loop across the aggregate markers and export in PNG format
lapply(
    X = names(new_markers_detected),
    FUN = function(name) {
        png(file.path(plotsDir, paste0(name, ".png")), width = 1000, height = 1000)
        p <- plotMarkerTSNE(seurat, genes = new_markers_detected[[name]])
        # Have to use `print` inside an `lapply()` call
        print(p)
        dev.off()
    }
)

gene2symbol <- gene2symbol(bcb_pool_filtered)
markers <- readCellTypeMarkersFile(
    "meta/cellTypeMarkers.csv",
    gene2symbol = gene2symbol)

atrial_markers <- markers %>%
    filter(cell == "Atrial") %>%
    pull("symbol")
atrial_markers
# "HEPACAM" "MYH6"    "MYL4"    "NKX2-5"  "NPPA"    "SCN5A"   "TNNT2"

counts <- seurat@data

# Get the summaries for each marker gene
lapply(seq_along(atrial_markers), function(a) {
    gene <- atrial_markers[[a]]
    summary(counts[gene, , drop = TRUE])
}) %>%
    set_names(atrial_markers)

# First, let's see how these genes are marking
plotMarkerTSNE(seurat, genes = atrial_markers, title = "Atrial markers")

# Let's look at the normalized expression data from Seurat
data <- FetchData(seurat, vars.all = atrial_markers)


data <- bcbioSingleCell:::.fetchGeneData.seurat(seurat, genes = atrial_markers)
geomean_x <- Matrix::colMeans(t(data))
head(geomean)

geomean <- colMeans(t(data))
geomean_2 <- rowMeans(data)

# Now, let's fetch the expression data in a more useful format. This builds
# upon the `Seurat::FetchData()` function but summarizes the counts into
# an `expression` column, and also calculates the `geomean`
data2 <- fetchTSNEExpressionData(seurat, genes = atrial_markers)
head(data2$geomean)

assignAndSaveData(
    name = "atrial_marker_expression_data",
    object = data,
    dir = dataDir
)
