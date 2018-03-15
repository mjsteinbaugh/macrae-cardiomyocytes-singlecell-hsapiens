# Latest version of this script is available here:
# script <- system.file(
#     file.path("R_scripts", "loadSingleCell.R"),
#     package = "bcbioSingleCell")
# file.edit(script)

library(bcbioSingleCell)

# Genome build is mismatched between GRCh38 (hg38) and GRCh37 (hg19) in YAML.
# Manually pass in the GRCh38 annotable as a fix.
bcb <- loadSingleCell(
    uploadDir = file.path(
        "data",
        "bcbio_harvard_indrop_v3",
        "final"),
    interestingGroups = "genotype",
    sampleMetadataFile = file.path(
        "meta",
        "multiplexed_sample_metadata.xlsx"),
    gtfFile = file.path(
        "meta",
        "Homo_sapiens.GRCh38.90.chr_patch_hapl_scaff.gtf"),
    annotable = annotable("Homo sapiens", release = 90),
    ensemblVersion = 90)

# Back up all data inside bcbio object
flatFiles <- flatFiles(bcb)

saveData(bcb, flatFiles, dir = "data")
