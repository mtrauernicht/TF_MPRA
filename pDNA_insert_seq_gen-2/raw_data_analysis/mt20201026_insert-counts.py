# This script will be used to extract inserts from raw fastq.gz files


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MT@NKI
# pDNA insert seq

# Extract inserts from fastq files
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load modules
import pysam
import regex

# define seq immediately downstream of barcode (truncated to 18nt)
# allow for 1 mismatch
downstream_seq = '(' + 'CATCGTCGCATCCAAGAG' + '){e<1}'

# open output file
tsv_out = open(snakemake.output["tsv"], "w")

# open input fastq stream
with pysam.FastxFile(snakemake.input["fq"]) as fq_in:

    # iterate over reads
    for read in fq_in:

        # extract read sequence
        seq = read.sequence

        # identify downstream seq position
        match = regex.search(downstream_seq, seq, regex.BESTMATCH)

        # if no match, skip to next
        if match is None:
            continue

        # extract insert
        end_bc = match.span()[0]
        insert = seq[0:end_bc]

        # if insert intact (>50 bp) and no N in insert, write to file
        if((len(insert) > 150) and ('N' not in insert)):
            # write to output file
            tsv_out.write(insert + '\n')

tsv_out.close()
