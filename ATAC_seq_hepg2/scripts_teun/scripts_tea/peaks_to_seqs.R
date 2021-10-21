

peaks <- rowRanges(exp)
library(BSgenome.Hsapiens.UCSC.hg38)

fasta <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, peaks)
names(fasta) <- peaks$name

ShortRead::writeFasta(fasta, file = "/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/peak_sequences.fa")


##Annotating open regions

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)

exp_annotate <- annotatePeak(peaks, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)

print(exp_annotate)

plotAnnoPie(exp_annotate)

pdf(file="/DATA/usr/t.filipovska/projects/ATAC-seq_TF-footprinting/analysis/Plots/test.pdf")
print(plotAnnoPie(exp_annotate))
dev.off()

## show overlap between annotation in cases where a peak falls under multiple gene regions.
## this plot shows the order of preference of ChIP seeker
upsetplot(exp_annotate)
