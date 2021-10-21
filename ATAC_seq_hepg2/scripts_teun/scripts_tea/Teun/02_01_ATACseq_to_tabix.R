# Libraries ---------------------------------------------------------------

library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(S4Vectors)
library(data.table, exclude = c("first", "second", "shift", "last"))
library(here)
library(glue)

# Files and directories ---------------------------------------------------

files <- readRDS(here("rds", "bamfile_atacseq_metadata.rds"))
dir <- "/DATA/projects/OSN/teun/processed_data/"
dir_in  <- glue("{dir}bam/atacseq")
dir_out <- glue("{dir}tabix/atacseq")

bamfiles <- files$bam_file
bamfiles <- BamFileList(bamfiles, yieldSize = 1e6L, asMates = TRUE)

# Blacklist from https://www.encodeproject.org/files/ENCFF547MET/
blacklist <- fread("/DATA/projects/OSN/teun/references/ENCFF547MET.bed")
blacklist <- with(blacklist, GRanges(V1, IRanges(V2, V3)))
blacklist <- as(blacklist, "GNCList")

# Setup parameters --------------------------------------------------------

setDTthreads(1)

param <- ScanBamParam(
  flag = scanBamFlag(
    isPaired = TRUE, isProperPair = TRUE, isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE, hasUnmappedMate = FALSE
  ),
  mapqFilter = 10L
)

# Functions ---------------------------------------------------------------

shift_fragments <- function(frags, positive = 4, negative = -5) {
  read1 <- tn5shift(GenomicAlignments::first(frags),  positive, negative)
  read2 <- tn5shift(GenomicAlignments::second(frags), positive, negative)
  punion(read1, read2, fill.gap = TRUE, ignore.strand = TRUE)
}

tn5shift <- function(x, positive = 4, negative = -5) {
  x <- resize(granges(x), fix = "start", 1)
  offset <- strand(x) == "+"
  offset <- Rle(
    ifelse(runValue(offset), positive, negative),
    runLength(offset)
  )
  GenomicRanges::shift(x, decode(offset))
}

write_bed <- function(granges, file) {
  granges <- data.table(
    chrom = decode(seqnames(granges)),
    start = start(granges),
    end   = end(granges)
  )
  fwrite(granges, file, append = TRUE, quote = FALSE, sep = "\t",
         col.names = FALSE)
  file
}

finalise_tabix <- function(file) {
  sorter <- glue("sort -k1,1 -k2,2n -o {file} {file}")
  system(sorter)
  zipper <- glue("bgzip {file}")
  system(zipper)
  file <- paste0(file, ".gz")
  index <- indexTabix(file, seq = 1, start = 2, end = 3)
  TabixFile(file, index)
}

read_bam <- function(file) {
  readGAlignmentPairs(file, param = param)
}

# Convert -----------------------------------------------------------------

bedfiles <- path(bamfiles)
bedfiles <- paste0(dir_out, "/", basename(bedfiles))
bedfiles <- gsub(".bam$", ".bed", bedfiles)
tabixfiles <- paste0(bedfiles, ".gz")

for(i in seq_along(bamfiles)) {
  bamfile <- bamfiles[[i]]
  outfile <- bedfiles[[i]]
  
  if (file.exists(tabixfiles[[i]])) {
    next
  }
  
  open(bamfile)
  while(length({data <- read_bam(bamfile)})) {
    data    <- shift_fragments(data)
    data    <- data[width(data) > 19]
    data    <- data[!overlapsAny(data, blacklist)]
    outfile <- write_bed(data, outfile)
  }
  close(bamfile)
  
  tbx <- finalise_tabix(outfile)
}

files$tabix_file <- tabixfiles
saveRDS(files, here("rds", "bamfile_atacseq_metadata.rds"))