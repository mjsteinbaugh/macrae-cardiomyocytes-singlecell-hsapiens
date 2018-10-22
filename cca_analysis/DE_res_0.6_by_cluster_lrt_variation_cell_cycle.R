library(dplyr)
library(tibble)
library(DESeq2)
library(BiocParallel)

for (cluster_n in c(1,4,6,11,10)) {
    load("/n/data1/cores/bcbio/PIs/calum_macrae/macrae_scRNASeq_cardiomyocytes/cca_analysis/cardio_combined_res_0.6.RData")
    
    metadata = cardio.combined@meta.data
    counts = cardio.combined@raw.data
    counts = counts[Matrix::rowSums(counts >= 5) >= 5, ]
    
    se = SummarizedExperiment(assays = list(counts = as.matrix(counts)),
                              colData = metadata)
    
    
    design <- ~nUMI + mitoRatio + nGene + S.Score + G2M.Score + KO
    reduced <- ~nUMI + mitoRatio + nGene + KO
    
    ddsraw = DESeqDataSet(se, design = design)
    
    ddsraw_cardio <- ddsraw[,ddsraw$res.0.6 %in% cluster_n]
    
    dds_cardio_lrt = DESeq(ddsraw_cardio, test = "LRT", full = design, reduced = reduced,
                           sfType = "poscounts", minmu = 1e-6, minRep = Inf,
                           parallel = TRUE,BPPARAM = MulticoreParam(8))

        
    save(dds_cardio_lrt, file = paste0("/n/data1/cores/bcbio/PIs/calum_macrae/macrae_scRNASeq_cardiomyocytes/cca_analysis/DE_KO_WT_results/by_cluster/variations/deseqObject_cardio_lrt_cluster_",cluster_n,"_variation_cell_cycle.Rdata"))
    rm(list = ls())
}
