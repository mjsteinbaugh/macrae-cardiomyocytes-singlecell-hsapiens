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
    
    
    design <- ~nUMI + mitoRatio  + Phase + KO
    reduced <- ~nUMI + mitoRatio + Phase
    
    ddsraw = DESeqDataSet(se, design = design)
    
    ddsraw_cardio <- ddsraw[,ddsraw$res.0.6 %in% cluster_n]
    
    dds_cardio_lrt = DESeq(ddsraw_cardio, test = "LRT", full = design, reduced = reduced,
                           sfType = "poscounts", minmu = 1e-6, minRep = Inf,
                           parallel = TRUE,BPPARAM = MulticoreParam(8))

    contrast = c("KO","KO","WT")
    dds_cardio_lrt_results_unshrunken = results(dds_cardio_lrt, contrast=contrast, cooksCutoff = FALSE) 
dds_cardio_lrt_results_shrunken = lfcShrink(dds_cardio_lrt,contrast = contrast, res = dds_cardio_lrt_results_unshrunken) 
    save(dds_cardio_lrt, dds_cardio_lrt_results_unshrunken,dds_cardio_lrt_results_shrunken,file = paste0("/n/data1/cores/bcbio/PIs/calum_macrae/macrae_scRNASeq_cardiomyocytes/cca_analysis/DE_KO_WT_results/by_cluster/lfcShrink/deseqObject_cardio_lrt_cluster_",cluster_n,".Rdata"))
    rm(list = ls())
}
