library(dplyr)
library(tibble)
library(DESeq2)
library(BiocParallel)

load("/n/data1/cores/bcbio/PIs/calum_macrae/macrae_scRNASeq_cardiomyocytes/cca_analysis/cardio_combined_res_0.6.RData")





metadata = cardio.combined@meta.data
counts = cardio.combined@raw.data
counts = counts[Matrix::rowSums(counts >= 5) >= 5, ]

se = SummarizedExperiment(assays=list(counts=as.matrix(counts)),
                          colData=metadata)


design <- ~nUMI + mitoRatio + nGene + S.Score + G2M.Score + KO
reduced <- ~nUMI + mitoRatio + nGene + S.Score + G2M.Score

ddsraw = DESeqDataSet(se, design=design)

ddsraw_cardio <- ddsraw[,ddsraw$res.0.6 %in% c(1,4,6,11,10) ]


dds_cardio_full = DESeq(ddsraw_cardio, test="Wald",
                   sfType="poscounts", minmu=1e-6, minRep=Inf,
		parallel=TRUE,BPPARAM = MulticoreParam(8))



dds_cardio_lrt = DESeq(ddsraw_cardio, test="LRT", full=design, reduced=reduced,
                         sfType="poscounts", minmu=1e-6, minRep=Inf,
		parallel=TRUE,BPPARAM = MulticoreParam(8))


contrast = c("KO","KO","WT")

dds_cardio_full_results = results(dds_cardio_full, contrast=contrast, cooksCutoff=FALSE,parallel=TRUE,BPPARAM = MulticoreParam(8))
dds_cardio_lrt_results = results(dds_cardio_lrt, contrast=contrast, cooksCutoff=FALSE,parallel=TRUE,BPPARAM = MulticoreParam(8))

save(dds_cardio_full,dds_cardio_lrt, file = "/n/data1/cores/bcbio/PIs/calum_macrae/macrae_scRNASeq_cardiomyocytes/cca_analysis/deseqObjects_cardio.Rdata")
save(dds_cardio_lrt_results, file = "/n/data1/cores/bcbio/PIs/calum_macrae/macrae_scRNASeq_cardiomyocytes/cca_analysis/dds_cardio_results.Rdata")
