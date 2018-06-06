# v0.1.14
library(bcbioSingleCell)
bcb <- bcbioSingleCell(
    uploadDir = file.path(
        "bcbio_harvard_indrop_v3",
        "final"
    ),
    organism = "Homo sapiens",
    interestingGroups = "genotype",
    sampleMetadataFile = file.path(
        "meta",
        "multiplexed_sample_metadata.xlsx"
    ),
    gffFile = file.path(
        "meta",
        "Homo_sapiens.GRCh38.90.chr_patch_hapl_scaff.gtf"
    ),
    ensemblRelease = 90
)
# Back up all data inside bcbio object
flatFiles <- flatFiles(bcb)
saveData(bcb, flatFiles, dir = "data")
