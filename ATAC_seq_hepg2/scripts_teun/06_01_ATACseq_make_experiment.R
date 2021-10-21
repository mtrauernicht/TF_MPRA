# Libraries ---------------------------------------------------------------

library(SummarizedExperiment)
library(GenomicRanges)
library(here)
library(rtracklayer)
library(data.table)

# Data --------------------------------------------------------------------

meta <- readRDS("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/rds/bamfile_atacseq_metadata.rds")

peaks <- import("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/bed_peaks/6420_peaks.narrowPeak")
peaks <- as(peaks, "GNCList")


# Metadata formatting -----------------------------------------------------

meta <- meta[meta$run == "6420", ]
meta <- meta[order(interaction(meta$celltype, meta$treatment, meta$replicate)),]
#rep <- tapply(meta$celltype, interaction(meta$celltype, meta$treatment, meta$replicate))
#meta$rep <- unlist(rep)
meta$alias <- paste(meta$celltype, meta$treatment, meta$replicate, sep = "_")
meta$celltype <- factor(meta$celltype, levels = c("HepG2", "K562"))
meta$treatment <- factor(meta$treatment, levels = c("DMSO", "T0901317", "Wy-14643", "CDCA", "Rifampicin"))
meta$replicate <- factor(meta$replicate, levels=c("R1", "R2", "R3"))

# Count -------------------------------------------------------------------

mat <- vapply(meta$tabix_file, function(file) {
  dat <- data.table::fread(file)
  dat <- with(dat, GRanges(V1, IRanges(V2, V3)))
  olaps <- findOverlaps(dat, peaks)
  olaps <- data.table(olap = to(olaps))
  olaps <- olaps[, .N, by = olap]
  olaps <- olaps[.(seq_along(peaks)), on = "olap"]
  return(olaps$N)
}, integer(NROW(peaks)))
mat[is.na(mat)] <- 0
colnames(mat) <- meta$alias


# Summarise Experiment ----------------------------------------------------

exp <- SummarizedExperiment(
  assays = mat,
  rowRanges = as(peaks, "GRanges"),
  colData = meta
)

saveRDS(exp, here("rds", "experiment_6420.rds"))
