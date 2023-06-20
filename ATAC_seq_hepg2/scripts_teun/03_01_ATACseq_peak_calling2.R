# Libraries ---------------------------------------------------------------

library(here)
library(glue)

path_bgzip <- "/DATA/usr/t.filipovska/software/Miniconda3/pkgs/tabix-0.2.6-ha92aebf_0/bin"

# Directories and data ----------------------------------------------------

files_df <- readRDS("/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/rds/bamfile_atacseq_metadata2.rds")
macs2 <- "/DATA/usr/t.filipovska/software/Miniconda3/envs/tf_activity/bin/macs2"
peak_dir <- "/DATA/usr/m.trauernicht/projects/SuRE-TF/ATAC_seq_hepg2/bed_peaks"

exps <- split(files_df$tabix_file, files_df$run)
#exps[names(exps) != "technical"]

# Functions ---------------------------------------------------------------

merge_tabixes <- function(files, out_file = tempfile(fileext = ".bed.gz")) {
  files <- paste0(files, collapse = " ")
  cmd <- glue("cat {files} > {out_file}")
  system(cmd)
  cmd <- glue("{path_bgzip}/bgzip -d {out_file}")
  system(cmd)
  gsub(".gz$", "", out_file)
}

call_macs2 <- function(file, name) {
  cmd <- glue(
    "{macs2} callpeak -t {file} -f BEDPE -g hs -n {name}", # change mm to hs for human
    " --nomodel --outdir {peak_dir} --keep-dup all"
  )
  system(cmd)
}


for (expname in names(exps)) {
  tabixes <- exps[[expname]]
  peakfile <- paste0(peak_dir, expname, "_peaks_2.narrowPeak")
  if (file.exists(peakfile)) {
    next
  }
  temp <- merge_tabixes(tabixes)
  call_macs2(temp, expname)
  unlink(temp)
}

