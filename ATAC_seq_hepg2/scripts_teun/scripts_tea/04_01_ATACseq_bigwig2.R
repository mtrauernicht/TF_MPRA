# Libraries ---------------------------------------------------------------

library(here)
library(glue)
library(data.table)

library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(GenomeInfoDb)


# Directories and data ----------------------------------------------------

files_df <- readRDS("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/rds/bamfile_atacseq_metadata2.rds")
bw_dir <- "/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/bigwig/20210610"
size_factors <- readRDS("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/bigwig/ATAC_size_factors.rds") %>% rownames_to_column("sample")

# Split up pools ----------------------------------------------------------


files_df <- files_df[files_df$run == "6420",]

files <- split(files_df$tabix_file, interaction(files_df$celltype, files_df$treatment, files_df$replicate, sep = "_"))
files <- files[lapply(files,length)>0]

seqinfo <- SeqinfoForUCSCGenome("hg38")
seqinfo <- keepStandardChromosomes(seqinfo, "Homo_sapiens")
seqinfo <- seqinfo[names(seqinfo)[-length(seqinfo)]]
seqnames <- seqnames(seqinfo)
gsize <- sum(seqlengths(seqinfo))

for (group in names(files)) {
  f <- files[[group]]
  outfile <- paste0(bw_dir, group, "_3.bw")
  if (file.exists(outfile)) {
    next
  }
  
  data <- lapply(f, fread)
  data <- rbindlist(data)

  
  data <- with(data, GRanges(V1, IRanges(V2, V3)))
  
  factor <- size_factors$sizefactor[size_factors$sample==group]
  
  data <- coverage(data) / factor
  data <- data[names(data) %in% seqnames]
  export.bw(data, paste0("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/bigwig/", group, "_3.bw"))
}
