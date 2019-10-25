library(data.table)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(OrganismDbi)
library("org.Hs.eg.db") # remember to install it if you don't have it already

#txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- loadDb("/home/bioinf/Documents/ATAC-seq/atac_pipeline/GRCH38_ENSEMBL.db")
peaks <- readPeakFile("/home/bioinf/Documents/Courses/MedBioInfo_courses/Algorithms_bioinf/Algorithms_Bioinformatics/project_data/Activation_mpbs.bed", header=FALSE)
peakAnno <- annotatePeak(peaks, tssRegion = c(-3000,3000),assignGenomicAnnotation=TRUE, level="gene",TxDb=txdb,overlap="TSS")

annon.out <- data.frame(peakAnno)
annon.out <- annon.out[!(colnames(annon.out) %like% "X")]
ensemblIds <- as.character(annon.out$geneId)
symbols <- mapIds(org.Hs.eg.db, keys = ensemblIds, keytype = "ENSEMBL", column="SYMBOL", multiVals="first")
annon.out$geneId <- symbols
write.table(annon.out,file = "/home/bioinf/Documents/Courses/MedBioInfo_courses/Algorithms_bioinf/Algorithms_Bioinformatics/project_data/Activation_nearest_gene.bed",quote = FALSE, sep='\t',row.names = FALSE)
