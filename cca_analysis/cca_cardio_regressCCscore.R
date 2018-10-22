library(dplyr)
library(tibble)
library(Seurat)


load("/n/data1/cores/bcbio/PIs/calum_macrae/macrae_scRNASeq_cardiomyocytes/data/2017-12-12/seurat.rda")

seurat@meta.data$KO <- as.factor(ifelse(as.character(seurat@meta.data$genotype) == "wildtype", "WT","KO"))

seurat_WT_metadata <- seurat@meta.data %>% tibble::rownames_to_column(var = "cell_id") %>% dplyr::filter(KO == "WT")
seurat_KO_metadata <- seurat@meta.data %>% tibble::rownames_to_column(var = "cell_id") %>% dplyr::filter(KO == "KO")


seurat_WT <- SubsetData(seurat, cells.use = seurat_WT_metadata$cell_id, subset.raw = T,do.clean = T)
seurat_KO <- SubsetData(seurat, cells.use = seurat_KO_metadata$cell_id, subset.raw = T,do.clean = T)

seurat_WT <- NormalizeData(object = seurat_WT, normalization.method = "LogNormalize", 
                           scale.factor = 10000)

seurat_WT <- ScaleData(
  object = seurat_WT,
  model.use = "linear")



seurat_KO <- NormalizeData(object = seurat_KO, normalization.method = "LogNormalize", 
                           scale.factor = 10000)

seurat_KO <- ScaleData(
  object = seurat_KO,
  model.use = "linear")



ccm <- read.table("./ccm.csv",sep = ",", header = T)
sGenes <- ccm %>%
  dplyr::filter(phase == "S") %>%
  pull("geneName")
g2mGenes <- ccm %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("geneName")


seurat_WT <- CellCycleScoring(
  seurat_WT,
  g2m.genes = g2mGenes,
  s.genes = sGenes)

seurat_KO <- CellCycleScoring(
  seurat_KO,
  g2m.genes = g2mGenes,
  s.genes = sGenes)



seurat_WT <- ScaleData(object = seurat_WT, vars.to.regress = c("nUMI", "S.Score", "G2M.Score"), 
                       display.progress = FALSE)



seurat_KO <- ScaleData(object = seurat_KO, vars.to.regress = c("nUMI", "S.Score", "G2M.Score"), 
                       display.progress = FALSE)


# Gene selection for input to CCA
seurat_WT <- FindVariableGenes(seurat_WT, do.plot = F)
seurat_KO <- FindVariableGenes(seurat_KO, do.plot = F)
g.1 <- head(rownames(seurat_WT@hvg.info), 1000)
g.2 <- head(rownames(seurat_KO@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(seurat_WT@scale.data))
genes.use <- intersect(genes.use, rownames(seurat_KO@scale.data))


cardio.combined <- RunCCA(seurat_WT, seurat_KO, genes.use = genes.use, num.cc = 30)

cardio.combined <- RunPCA(object = cardio.combined, pc.genes = genes.use, pcs.compute = 100, do.print = FALSE)



# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = cardio.combined, reduction.use = "cca", group.by = "KO", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = cardio.combined, features.plot = "CC1", group.by = "KO", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = cardio.combined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)



p3 <- MetageneBicorPlot(cardio.combined, grouping.var = "KO", dims.eval = 1:30, 
                        display.progress = FALSE)



DimHeatmap(object = cardio.combined, reduction.type = "cca", cells.use = 500, 
           dim.use = 1:9, do.balanced = TRUE)





cardio.combined <- AlignSubspace(cardio.combined, reduction.type = "cca", grouping.var = "KO", 
                                 dims.align = 1:25)



cardio.combined <- RunTSNE(cardio.combined, reduction.use = "cca.aligned", dims.use = 1:25, 
                           do.fast = T)
cardio.combined <- FindClusters(cardio.combined, reduction.type = "cca.aligned", 
                                resolution = 0.6, dims.use = 1:25)

save(cardio.combined,file = "./cardio_combined.RData")
