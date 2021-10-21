# Libraries ---------------------------------------------------------------

library(here)
library(glue)
library(data.table)

library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(GenomeInfoDb)
library(dplyr)
library(tibble)


# Directories and data ----------------------------------------------------

files_df <- readRDS("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/rds/bamfile_atacseq_metadata2.rds")
bw_dir <- "/DATA/usr/m.trauernicht/projects/SuRE-TF/ATAC_seq/bigwig/bamCoverage_rpm_peaks/"
size_factor <- readRDS("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/bigwig/ATAC_size_factors.rds") %>%
  rownames_to_column("sample")

# Split up pools ----------------------------------------------------------

files_df <- files_df[files_df$run == "6420",]

files <- split(files_df$bam_file, interaction(files_df$celltype, files_df$treatment, files_df$replicate, sep = "_"))
files <- files[lapply(files,length)>0]

seqinfo <- SeqinfoForUCSCGenome("hg38")
seqinfo <- keepStandardChromosomes(seqinfo, "Homo_sapiens")
seqinfo <- seqinfo[names(seqinfo)[-length(seqinfo)]]
seqnames <- seqnames(seqinfo)
gsize <- sum(seqlengths(seqinfo))

bamCoverage <- "/DATA/usr/t.filipovska/software/Miniconda3/envs/tf_activity/bin/bamCoverage"

bam_coverage <- function(file, factor, name) {
  cmd <- glue(
    "{bamCoverage} -b {file} --scaleFactor {factor} -bs 1 -p 20 --centerReads -o {name}", 
  )
  system(cmd)
}


for (group in names(files)) {
  f <- files[[group]]
  outfile <- paste0(bw_dir, group, ".bw")
  if (file.exists(outfile)) {
    next
  }
  factor <- 1 / (size_factor$rip[size_factor$sample == group] / 1e6)
  bam_coverage(f, factor, outfile)
}
