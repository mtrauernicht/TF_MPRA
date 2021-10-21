# Libraries ---------------------------------------------------------------

library(here)
library(glue)

# Directories and data ----------------------------------------------------

files_df <- readRDS(here("rds", "bamfile_atacseq_metadata.rds"))
macs2 <- "/home/t.vd.brand/.local/bin/macs2"
peak_dir <- "/DATA/projects/OSN/teun/processed_data/bed_peaks/"

exps <- split(files_df$tabix_file, files_df$exp)
exps[names(exps) != "technical"]

# Functions ---------------------------------------------------------------

merge_tabixes <- function(files, out_file = tempfile(fileext = ".bed.gz")) {
  files <- paste0(files, collapse = " ")
  cmd <- glue("cat {files} > {out_file}")
  system(cmd)
  cmd <- glue("bgzip -d {out_file}")
  system(cmd)
  gsub(".gz$", "", out_file)
}

call_macs2 <- function(file, name) {
  cmd <- glue(
    "{macs2} callpeak -t {file} -f BEDPE -g mm -n {name}", # change mm to hg for human
    " --nomodel --outdir {peak_dir} --keep-dup all"
  )
  system(cmd)
}


for (expname in names(exps)) {
  tabixes <- exps[[expname]]
  peakfile <- paste0(peak_dir, expname, "_peaks.narrowPeak")
  if (file.exists(peakfile)) {
    next
  }
  temp <- merge_tabixes(tabixes)
  call_macs2(temp, expname)
  unlink(temp)
}
