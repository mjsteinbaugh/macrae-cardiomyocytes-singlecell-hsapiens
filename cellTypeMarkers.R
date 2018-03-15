library(bcbioSingleCell)
library(tidyverse)

annotable <- annotable("Homo sapiens", release = 90)
saveData(annotable, dir = "data")

myh_genes <- annotable %>%
    filter(str_detect(symbol, "^MYH")) %>%
    select(ensgene, symbol)
write_csv(myh_genes, "meta/myh_genes.csv")

synonyms <- basejump::synonyms$hsapiens
saveData(synonyms, dir = "data")
