library(bcbioSingleCell)  # v0.1.18
bcb <- bcbioSingleCell(
    uploadDir = "bcbio_harvard_indrop_v3/final",
    organism = "Homo sapiens",
    interestingGroups = "genotype",
    sampleMetadataFile = "meta/multiplexed_sample_metadata.xlsx",
    genomeBuild = "GRCh38",
    ensemblRelease = 90
)
# Back up all data inside bcbio object
flatFiles <- flatFiles(bcb)
saveData(bcb, flatFiles, dir = file.path("data", Sys.Date()))
