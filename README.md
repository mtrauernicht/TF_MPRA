# TF reporter MPRA

**Introduction:**\
Transcriptional reporters for 86 TFs were systematically designed and probed in 9 cell types and dozens of TF perturbation conditions with the aim to identify the most sensitive reporters per TF.

**TF reporter library design:**\
The designed first library contains:
- General reporter design: with 4 identical TF binding sites followed by a minimal promoter and a barcode in the transcription unit
- In total reporters for 86 different TFs
- 10 or 5 bp spacing between the TF binding sites
- 10 or 21 bp distance from the TF binding sites to the minimal promoter
- 3 different minimal promoters
- 5-8 barcodes per TF reporter
- In total ~36,000 uniquely barcoded reporters per TF

All TF reporters were designed using FIMO. This way, the spacings were designed to be inactive, while the TF binding sites were ensured to be intact.

**Experimental setup:**\
The library was transfected into:
- 9 different cell types
- 75 TF perturbation conditions
- TF knockdown
- TF overexpression
- Signaling pathway perturbation
- TF degradation

The cells were grown in these condition for 24h after transfection before RNA isolation. Barcoded transcripts were then reverse transcribed and amplified before sequencing.

**Sequencing data analysis:**
- Raw sequencing data were processed by counting the barcodes and clustering the barcodes using starcode.
- Barcode counts in the cDNA were normalized by the counts in the plasmid library.
