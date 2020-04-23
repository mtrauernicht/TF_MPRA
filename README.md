# SuRE_TF_1

**Introduction:**
I modified the SuRE plasmid that Joris van Arensbergen used to probe genome-wide autonomous promoter activity in human cells (doi: 10.1038/nbt.3754) to measure many TF activities in parallel. This repository contains oligo design and data analysis of the first library.

**TF reporter library design:**\
The designed first library contains:
- ~18,000 TF reporters, each with 4 identical TF binding sites, followed by a minP and a barcode in the transcription unit
- 29 TFs
- 10 or 5 bp spacing between the TF binding sites
- 10 or 21 bp distance from the TF binding sites to the minimal promoter
- 3 different minimal promoters
- 8 barcodes per TF reporter

All TF reporters were designed using FIMO. This way, the spacings were designed to be inactive, while the TF binding sites were ensured to be active.

**Experimental setting:**\
In a first experiment, the library was transfected into cells using 7 different conditions:
1. mESC: 2i+LIF
2. mESC: 2i-LIF
3. mESC: LIF+PD
4. mESC: LIF+CH
5. mESC: -vitA
6. mESC: N2B27 (neural differentiation medium)
7. NPC

**Sequencing data analysis:**
- Raw sequencing data were processed by counting the barcodes and clustering the barcodes using starcode