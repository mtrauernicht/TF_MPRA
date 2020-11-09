## Create barcodes with length 13 using R
## Do this via bash because it takes forever

library(DNABarcodes)

barc <- create.dnabarcodes(n = 13, dist = 3, filter.triplets = T,metric = "seqlev",filter.gc = T, filter.self_complementary = T, cores = 24)

write.csv2(barc, file = "/DATA/usr/m.trauernicht/projects/SuRE-TF/data/library_design/output/barc-13.csv")
