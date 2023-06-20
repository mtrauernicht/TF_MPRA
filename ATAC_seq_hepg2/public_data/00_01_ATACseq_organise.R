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
  "GSE169929", "SRR14103347", "MCF7", "DMSO", "R1"#,
  # "GSE169929", "SRR14103348", "MCF7", "DMSO", "R2",
  # "GSE187226", "SRR16809040", "HCT116", "DMSO", "R1",
  # "GSE187226", "SRR16809041", "HCT116", "DMSO", "R2", 
  # "GSE211724", "SRR21152045", "NPC","DMSO", "R1", 
  # "GSE108513", "SRR6418071", "HEK293", "DMSO", "R1", 
  # "GSE108513", "SRR6418073", "HEK293", "DMSO", "R2",
  # "GSE169955", "SRR14103420", "A549", "DMSO", "R1",
  # "GSE169955", "SRR14103424", "A549", "DMSO", "R2",
  # "GSE169955", "SRR14103428", "A549","DMSO", "R3",
  # "GSE121840", "SRR8171322", "U2OS","DMSO", "R1",
  # "GSE121840", "SRR8171326", "U2OS", "DMSO", "R2"
)

# Files -------------------------------------------------------------------

# List all files from the runs indicated above
runs <- paste0(sort(unique(df$run)), collapse = "|")
gcffiles <- list.files("/DATA/usr/m.trauernicht/data/ATAC", pattern = runs, recursive = F, full.names = TRUE)
base_gcf <- basename(gcffiles)

# Find the matching fastq files
read1 <- gcffiles[pmatch(with(df, glue("{run}_{basename}_1")), base_gcf)]
read1 <- read1[!is.na(read1)]
read2 <- gcffiles[pmatch(with(df, glue("{run}_{basename}_2")), base_gcf)]
read2 <- read2[!is.na(read2)]

# Check if all samples have been found
check <- any(is.na(read1) | is.na(read2))
stopifnot(!check)

# Check filename validity -------------------------------------------------

names <- c(read1, read2)
names <- gsub("1|2", "", names)
names <- matrix(names, ncol = 2)

# If removing R1 and R2, all parallel filenames should be equal
check <- all(names[, 1] == names[, 2])
stopifnot(check)

# Symlinking files --------------------------------------------------------

dir <- "/DATA/usr/m.trauernicht/data/ATAC"
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
  df, "/DATA/usr/m.trauernicht/projects/SuRE-TF/ATAC_seq_hepg2/public_data/analysis/atac_fastq_table.tsv"
  
)
