# README:
#
# This script makes BAM files from fastq files. There are some fastq files that
# need to be merged as they are resequenced samples. Meanwhile, it collects
# metadata about the BAM files and saves this as a table.

# Libraries ---------------------------------------------------------------

library(here)
library(glue)
library(S4Vectors)
library(IRanges)
library(withr)


# Files -------------------------------------------------------------------

bam_dir <- "/DATA/usr/m.trauernicht/projects/SuRE-TF/ATAC_seq_hepg2/public_data/analysis/bam"
genome <- "/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/data/hg38bwaidx"
picard <- "/DATA/usr/t.filipovska/software/picard.jar"
table_out <- here("rds", "bamfile_atacseq_metadata.rds")

files <- read.table("ATAC_seq_hepg2/public_data/analysis/atac_fastq_table.tsv", 
                    sep = "\t", stringsAsFactors = FALSE, header = TRUE)

# Declare merges ----------------------------------------------------------

# Merge function
 merge_this <- function(this, into, list = merges) {
   list[[into]] <- c(list[[into]], list[[this]])
   list[[this]] <- NULL
   return(list)
 }
 
# # Setup named indices
 merges <- as.list(setNames(
   seq_len(nrow(files)), 
   with(files, glue("{run}_{basename}"))
 ))

# Declare merges
# merge_table <- tibble::tribble(
#   ~reseq, ~orig,
#   "5859_17_NANOG_2H_GTAGAGGA-GCGATCTA_S17",  "5831_12_NANOG_2H_GTAGAGGA-GCGATCTA_S12",
#   "5859_18_NANOG_24H_AGGCAGAA-TCTACTCT_S18", "5831_15_NANOG_24H_AGGCAGAA-TCTACTCT_S15",
#   "5859_19_OCT4_0H_TAGGCATG-GCGATCTA_S19",   "5831_6_OCT4_0H_TAGGCATG-GCGATCTA_S8",
#   "5859_20_OCT4_2H_CTCTCTAC-GCGATCTA_S20",   "5831_7_OCT4_2H_CTCTCTAC-GCGATCTA_S9"
# )

# Make merge indices
#for (i in seq_len(nrow(merge_table))) {
  #merges <- merge_this(
   # this = merge_table$reseq[i],
    #into = merge_table$orig[i],
    #list = merges
#  )
#}

# Build BAM table ---------------------------------------------------------

bam_table <- files[unlist(heads(merges, 1)), ]
bam_table$read1 <- NULL
bam_table$read2 <- NULL
bam_table$bam_file <- as.character(with(
  bam_table, 
  glue("{bam_dir}/{run}_{basename}.bam")
))
bam_table <- as(bam_table, "DataFrame")

read1 <- CharacterList(lapply(merges, function(i) {
  files$read1[i]
}))
read2 <- CharacterList(lapply(merges, function(i) {
  files$read2[i]
}))
fastq <- DataFrame(
  read1 = read1,
  read2 = read2
)
bam_table$fastq_files <- fastq

# Make commands -----------------------------------------------------------

cmd_bwa <- function(read1, read2, output) {
              path <- "/DATA/usr/m.trauernicht/software/miniconda3/bin"
              glue("{path}/bwa mem -M -t 20 {genome} {read1} {read2}",
              "{path}/samtools view -h -b -q 10",
              "{path}/samtools sort -o {output}",
    .sep = " | "
  )
  
}


cmd_picard <- function(input, output) {
  logfile <- gsub(".bam$", "_log.txt", output)
  
  path <- "/DATA/usr/m.trauernicht/software/miniconda3/bin"
  glue("{path}/java -jar {picard} MarkDuplicates",
    "REMOVE_DUPLICATES=true",
    "ASSUME_SORTED=true",
    "INPUT={input}",
    "OUTPUT={output}",
    "METRICS_FILE={logfile}",
    .sep = " "
  )
}


# Do mapping --------------------------------------------------------------

# I can for loop here because the looping won't be the time consuming part ;)
for (i in seq_len(nrow(bam_table))) {
  
  # Declare files
  # Fastq files
  r1 <- bam_table$fastq_files$read1[[i]]
  r2 <- bam_table$fastq_files$read2[[i]]
  stopifnot(length(r1) == length(r2))
  
  # (temporary) bam files
  out_bam <- bam_table$bam_file[[i]]
  if (file.exists(out_bam)) {
    next
  }
  tmpfiles1 <- gsub("_1.+", "_temp.bam", r1)
  tmpfiles1 <- paste0(dirname(out_bam), "/", basename(tmpfiles1))

  # Perform mapping
  map <- cmd_bwa(r1, r2, tmpfiles1)
  for (j in seq_along(r1)) {
    if (file.exists(tmpfiles1[j])) {
      next
    }
    system(map[j])
  }
  
  # Merge what needs to be merged
  if (length(tmpfiles1) > 1) {
    tmpfiles2 <- as.character(glue("{dirname(out_bam)}/temp.bam"))
    cmd_merge <- paste0(tmpfiles1, collapse = " ")
    cmd_merge <- glue("{path}/samtools merge {tmpfiles2} {cmd_merge}")
    system(cmd_merge)
  } else {
    tmpfiles2 <- tmpfiles1
  }
  
  # De-duplicate
  dedup <- cmd_picard(tmpfiles2, out_bam)
  system(dedup)
  
  # Index
  path <- "/DATA/usr/m.trauernicht/software/miniconda3/bin"
  
  index <- glue("{path}/samtools index {out_bam}")
  system(index)
  
  # Clean up temporary files
  cleanup <- unique(c(tmpfiles1, tmpfiles2))
  for (file in cleanup) {
    unlink(file)
  }
}


# Collect log info --------------------------------------------------------

logfiles <- gsub(".bam$", "_log.txt", bam_table$bam_file)
logs <- lapply(logfiles, function(logfile) {
  txt <- readLines(logfile)
  empty <- which(txt == "")
  metrics <- which(startsWith(txt, "## METRICS CLASS"))[[1]] + 1
  empty <- empty[empty > metrics][1] - 1
  read.table(text = txt[metrics:empty], sep = "\t", header = TRUE)
})
logs <- do.call(rbind, logs)[-1]
colnames(logs) <- tolower(colnames(logs))
logs <- as(logs, "DataFrame")
bam_table$logs <- logs

# Save metadata -----------------------------------------------------------

saveRDS(bam_table, table_out)


