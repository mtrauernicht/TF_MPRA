# Library statements ------------------------------------------------------

library(here)
library(data.table)
library(dplyr, include.only = "case_when")
library(readr, include.only = "write_tsv")
library(stringr, include.only = "str_extract")
library(glue)

# Declare samples ---------------------------------------------------------

df <- tibble::tribble(
  ~run, ~basename, ~celltype, ~treatment, ~replicate,
  # Run 6420
  "6420", "10_HepG2_CDCA_R2_CGAGGCTG-ATCTACAC_S1", "HepG2", "CDCA", "R2",
  "6420", "11_HepG2_Rifampicin_R2_AAGAGGCA-ATCTACAC_S3", "HepG2", "Rifampicin", "R2",
  "6420", "12_K562_DMSO_R2_GTAGAGGA-ATCTACAC_S6", "K562", "DMSO", "R2",
  "6420", "13_HepG2_DMSO_R3_CAGAGAGG-CTCTCTAT_S13", "HepG2", "DMSO", "R3", 
  "6420", "14_HepG2_T0901317_R3_CAGAGAGG-TATCCTCT_S14", "HepG2","T0901317", "R3", 
  "6420", "15_HepG2_Wy-14643_R3_CAGAGAGG-AGAGTAGA_S15", "HepG2", "Wy-14643", "R3", 
  "6420", "16_HepG2_CDCA_R3_CAGAGAGG-GTAAGGAG_S16", "HepG2", "CDCA", "R3",
  "6420", "17_HepG2_Rifampicin_R3_CAGAGAGG-ACTGCATA_S17", "HepG2", "Rifampicin", "R3",
  "6420", "18_K562_DMSO_R3_CAGAGAGG-AAGGAGTA_S18", "K562", "DMSO", "R3",
  "6420", "1_HepG2_DMSO_R1_TAAGGCGA-ATCTACAC_S2", "HepG2","DMSO", "R1",
  "6420", "2_HepG2_T0901317_R1_CGTACTAG-ATCTACAC_S5", "HepG2","T0901317", "R1",
  "6420", "3_HepG2_Wy-14643_R1_AGGCAGAA-ATCTACAC_S4", "HepG2", "Wy-14643", "R1",
  "6420", "4_HepG2_CDCA_R1_TCCTGAGC-ATCTACAC_S8", "HepG2", "CDCA", "R1",
  "6420", "5_HepG2_Rifampicin_R1_GGACTCCT-ATCTACAC_S9", "HepG2", "Rifampicin", "R1",
  "6420", "6_K562_DMSO_R1_TAGGCATG-ATCTACAC_S7", "K562", "DMSO", "R1",
  "6420", "7_HepG2_DMSO_R2_CTCTCTAC-ATCTACAC_S10", "HepG2","DMSO", "R2",
  "6420", "8_HepG2_T0901317_R2_CAGAGAGG-ATCTACAC_S11", "HepG2","T0901317", "R2",
  "6420", "9_HepG2_Wy-14643_R2_GCTACGCT-ATCTACAC_S12", "HepG2", "Wy-14643", "R2"
)

# Files -------------------------------------------------------------------

# List all files from the runs indicated above
runs <- paste0(sort(unique(df$run)), collapse = "|")
gcffiles <- list.files("/shared/gcf", pattern = runs, recursive = TRUE, full.names = TRUE)
base_gcf <- basename(gcffiles)

# Find the matching fastq files
read1 <- gcffiles[pmatch(with(df, glue("{run}_{basename}_R1")), base_gcf)]
read1 <- read1[!is.na(read1)]
read2 <- gcffiles[pmatch(with(df, glue("{run}_{basename}_R2")), base_gcf)]
read2 <- read2[!is.na(read2)]

# Check if all samples have been found
check <- any(is.na(read1) | is.na(read2))
stopifnot(!check)

# Check filename validity -------------------------------------------------

names <- c(read1, read2)
names <- gsub("R1|R2", "", names)
names <- matrix(names, ncol = 2)

# If removing R1 and R2, all parallel filenames should be equal
check <- all(names[, 1] == names[, 2])
stopifnot(check)

# Symlinking files --------------------------------------------------------

dir <- "/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis"
new_read1 <- glue("{dir}/{basename(read1)}")
new_read2 <- glue("{dir}/{basename(read2)}")

file_is_linked <- function(x) {
  f <- list.files(dirname(x)[1], full.names = TRUE)
  x %in% f
}

check <- file_is_linked(new_read1) & file_is_linked(new_read2)

cmd_r1 <- glue("ln -s {read1} {new_read1}")
cmd_r2 <- glue("ln -s {read2} {new_read2}")

for (i in seq_along(check)[!check]) {
  system(cmd_r1[i])
  system(cmd_r2[i])
}

# Update data with fastq locations ----------------------------------------

df <- transform(
  df,
  read1 = new_read1,
  read2 = new_read2
)

# Export ------------------------------------------------------------------

write_tsv(
  df, "/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/atac_fastq_table.tsv"
  
)
