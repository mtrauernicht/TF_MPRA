from os.path import join

# Globals ---------------------------------------------------------------------

# Full path to working directory
W_DIR = '/DATA/usr/m.trauernicht/projects/SuRE_deep_scan_trp53_gr/'

# Expression and copy number (ECN)---------------------------------------------

# Full path to cDNA and pDNA raw data folder (fastq)
ECN_DIR = W_DIR + 'data/pDNA_insert_seq/fastq_files/'

# Extract cDNA and pDNA sample names
ECN, = glob_wildcards(join(ECN_DIR, '{ecn,[^/]+}_R1_001.fastq.gz'))
print(ECN)

# Pattern for SE read
S1 = '{ecn}_R1_001.fastq.gz'
