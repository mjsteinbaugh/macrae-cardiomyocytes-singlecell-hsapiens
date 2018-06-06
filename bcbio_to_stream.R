library(bcbioSingleCell)
loadDataAsName(bcb = bcb_pool_filtered, dir = "data/2017-12-01")
validObject(bcb)
show(bcb)
filtered_counts <- bcb@assays[["assay"]]
saveData(filtered_counts, dir = file.path("data", Sys.Date()))
dim(filtered_counts)
# Convert to TSV file
filtered_counts %>%
    as.matrix() %>%
    write.table(
        file = "filtered_counts.tsv",
        sep = "\t",
        quote = FALSE,
        col.names = NA
    )
