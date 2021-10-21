# Activate conda environment to load installed packages (create this environemt and install packages first)
source activate tf-activity

# Change directory to working directory (directory with snakemake file and all other scripts)
cd /DATA/usr/m.trauernicht/projects/SuRE_deep_scan_trp53_gr/pDNA_insert_seq/analysis_gcf6382/

# Run snakemake script using 12 cores
snakemake --cores 12
