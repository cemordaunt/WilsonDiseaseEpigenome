# WD Liver and Blood DMR Fibrosis Overlap ####
# From Hotta et al.
# Charles Mordaunt
# 2/22/18

# Packages ####
library(GenomicRanges)
library(ChIPpeakAnno)
library(rtracklayer)
library(R.utils)

# Data ####
fibrosis <- read.csv("Tables/2017 - Hotta - suppl table 1.csv", header=TRUE, stringsAsFactors = FALSE) # Assuming hg19
fibrosis$Start <- as.integer(gsub(",", "", fibrosis$Start))
fibrosis$End <- as.integer(gsub(",", "", fibrosis$End))
fibrosisAll <- fibrosis[,c("Chr", "Start", "End")]
fibrosisAll <- subset(fibrosisAll, !is.na(fibrosisAll$Start))
fibrosisHyper <- subset(fibrosis, Direction=="hypermethylated", select = c("Chr", "Start", "End"))
fibrosisHypo <- subset(fibrosis, Direction=="hypomethylated", select = c("Chr", "Start", "End"))
fibrosisBack <- read.delim("UCSC Tracks/HM450_hg38.bed", sep="\t", header=FALSE, stringsAsFactors = FALSE) #hg38 already

liverAll <- read.delim("WD_Specific_DMRs.bed", sep="\t", header=FALSE, stringsAsFactors = FALSE)
liverHyper <- read.delim("WD_Specific_hyper_DMRs.bed", sep="\t", header=FALSE, stringsAsFactors = FALSE)
liverHypo <- read.delim("WD_Specific_hypo_DMRs.bed", sep="\t", header=FALSE, stringsAsFactors = FALSE)
liverBack <- read.delim("Healthy_Ctrl_vs_WD_Liver_Subsetted_DMR_Background.bed", header=FALSE, stringsAsFactors = FALSE)

setwd("C:/Users/Booboo/Charles/Documents/Programming/Wilson's Disease Blood Batch 2")
blood <- read.delim("Tables/WD Specific DMR Stats with annotation 2.txt", sep="\t", header = TRUE, stringsAsFactors = FALSE)
blood <- blood[,c("chr", "start", "end", "direction")]
blood <- unique(blood)
bloodAll <- blood[,c("chr", "start", "end")]
bloodHyper <- subset(blood, direction=="hyper", select=c("chr", "start", "end"))
bloodHypo <- subset(blood, direction=="hypo", select=c("chr", "start", "end"))
bloodBack <- read.delim("UCSC Tracks/HCvWD Blood 2 Subsetted Background DMRs.bed", sep="\t", header=FALSE, stringsAsFactors = FALSE)
rm(blood, fibrosis)

# Make Genomic Ranges ####
setwd("C:/Users/Booboo/Charles/Documents/Programming/Wilson's Disease Liver")
GR_fibrosisAll <- GRanges(seqnames = fibrosisAll$Chr, ranges=IRanges(start=fibrosisAll$Start, end=fibrosisAll$End))
GR_fibrosisHyper <- GRanges(seqnames = fibrosisHyper$Chr, ranges=IRanges(start=fibrosisHyper$Start, end=fibrosisHyper$End))
GR_fibrosisHypo <- GRanges(seqnames = fibrosisHypo$Chr, ranges=IRanges(start=fibrosisHypo$Start, end=fibrosisHypo$End))
GR_fibrosisBack <- GRanges(seqnames = fibrosisBack$V1, ranges=IRanges(start=fibrosisBack$V2, end=fibrosisBack$V3))

GR_liverAll <- GRanges(seqnames = liverAll$V1, ranges=IRanges(start=liverAll$V2, end=liverAll$V3))
GR_liverHyper <- GRanges(seqnames = liverHyper$V1, ranges=IRanges(start=liverHyper$V2, end=liverHyper$V3))
GR_liverHypo <- GRanges(seqnames = liverHypo$V1, ranges=IRanges(start=liverHypo$V2, end=liverHypo$V3))
GR_liverBack <- GRanges(seqnames = liverBack$V1, ranges=IRanges(start=liverBack$V2, end=liverBack$V3))

GR_bloodAll <- GRanges(seqnames = bloodAll$chr, ranges=IRanges(start=bloodAll$start, end=bloodAll$end))
GR_bloodHyper <- GRanges(seqnames = bloodHyper$chr, ranges=IRanges(start=bloodHyper$start, end=bloodHyper$end))
GR_bloodHypo <- GRanges(seqnames = bloodHypo$chr, ranges=IRanges(start=bloodHypo$start, end=bloodHypo$end))
GR_bloodBack <- GRanges(seqnames = bloodBack$V1, ranges=IRanges(start=bloodBack$V2, end=bloodBack$V3))

# Liftover Fibrosis DMRs from hg19 to hg38 ####
# Get chain file
url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
download.file(url, basename(url))
gunzip(basename(url))
chain <- import.chain("hg19ToHg38.over.chain")

# liftOver from hg38 to hg19
seqlevelsStyle(GR_fibrosisAll) <- "UCSC"
seqlevelsStyle(GR_fibrosisHyper) <- "UCSC"
seqlevelsStyle(GR_fibrosisHypo) <- "UCSC"

GR_fibrosisAll_hg38 <- unlist(liftOver(GR_fibrosisAll, chain))
GR_fibrosisHyper_hg38 <- unlist(liftOver(GR_fibrosisHyper, chain))
GR_fibrosisHypo_hg38 <- unlist(liftOver(GR_fibrosisHypo, chain))

# If regions are overlapping, make them non-overlapping (disjoint)
if(!isDisjoint(GR_fibrosisAll_hg38)){GR_fibrosisAll_hg38 <- disjoin(GR_fibrosisAll_hg38)} 
if(!isDisjoint(GR_fibrosisHyper_hg38)){GR_fibrosisHyper_hg38 <- disjoin(GR_fibrosisHyper_hg38)} 
if(!isDisjoint(GR_fibrosisHypo_hg38)){GR_fibrosisHypo_hg38 <- disjoin(GR_fibrosisHypo_hg38)} 
if(!isDisjoint(GR_fibrosisBack)){GR_fibrosisBack <- disjoin(GR_fibrosisBack)} 

# Liver Venn Diagrams ####
# Background and 450K Probes
pdf(file="Figures/Fibrosis Liver Background and 450K Probes Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_fibrosisBack, GR_liverBack), NameOfPeaks = c("Probes", "WD_Liver_Background"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), cat.fontfamily="sans",fontfamily="sans")
dev.off()

GR_fibrosis_liver_back <- subsetByOverlaps(GR_liverBack, GR_fibrosisBack) # Background regions that overlap 450K probes
length(GR_fibrosis_liver_back) #73537

# All
pdf(file="Figures/Fibrosis Liver All DMR Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_fibrosisAll_hg38, GR_liverAll), NameOfPeaks = c("Fibrosis", "WD_Liver"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), totalTest = length(GR_fibrosis_liver_back), cat.fontfamily="sans",fontfamily="sans")
dev.off()

# Hypergeometric Test pvalue
as.data.frame(venn$p.value)$pval #2.465073e-22

# Hyper
pdf(file="Figures/Fibrosis Liver Hyper DMR Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_fibrosisHyper_hg38, GR_liverHyper), NameOfPeaks = c("Fibrosis", "WD_Liver"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), totalTest = length(GR_fibrosis_liver_back), cat.fontfamily="sans",fontfamily="sans")
dev.off()

# Hypergeometric Test pvalue
as.data.frame(venn$p.value)$pval #9.417433e-25

GR_liverHyperFibrosis <- subsetByOverlaps(GR_liverHyper, GR_fibrosisHyper_hg38)
liverHyperFibrosis <- as.data.frame(GR_liverHyperFibrosis)
write.table(liverHyperFibrosis, "WD Liver Hyper DMRs Overlapping with Fibrosis DMRs.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Hypo
pdf(file="Figures/Fibrosis Liver Hypo DMR Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_fibrosisHypo_hg38, GR_liverHypo), NameOfPeaks = c("Fibrosis", "WD_Liver"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), totalTest = length(GR_fibrosis_liver_back), cat.fontfamily="sans",fontfamily="sans")
dev.off()

# Hypergeometric Test pvalue
as.data.frame(venn$p.value)$pval #7.436024e-17

GR_liverHypoFibrosis <- subsetByOverlaps(GR_liverHypo, GR_fibrosisHypo_hg38)
liverHypoFibrosis <- as.data.frame(GR_liverHypoFibrosis)
write.table(liverHypoFibrosis, "WD Liver Hypo DMRs Overlapping with Fibrosis DMRs.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Blood Venn Diagrams ####
# Background and 450K Probes
pdf(file="Figures/Fibrosis Blood Background and 450K Probes Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_fibrosisBack, GR_bloodBack), NameOfPeaks = c("Probes", "WD_Blood_Background"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), cat.fontfamily="sans",fontfamily="sans")
dev.off()

GR_fibrosis_blood_back <- subsetByOverlaps(GR_bloodBack, GR_fibrosisBack) # Background regions that overlap 450K probes
length(GR_fibrosis_blood_back) #80360

# All
pdf(file="Figures/Fibrosis Blood All DMR Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_fibrosisAll_hg38, GR_bloodAll), NameOfPeaks = c("Fibrosis", "WD_Blood"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), totalTest = length(GR_fibrosis_blood_back), cat.fontfamily="sans",fontfamily="sans")
dev.off()

# Hypergeometric Test pvalue
as.data.frame(venn$p.value)$pval #0.5950856

# Hyper
pdf(file="Figures/Fibrosis Blood Hyper DMR Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_fibrosisHyper_hg38, GR_bloodHyper), NameOfPeaks = c("Fibrosis", "WD_Blood"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), totalTest = length(GR_fibrosis_blood_back), cat.fontfamily="sans",fontfamily="sans")
dev.off()

# Hypergeometric Test pvalue
as.data.frame(venn$p.value)$pval #1

GR_bloodHyperFibrosis <- subsetByOverlaps(GR_bloodHyper, GR_fibrosisHyper_hg38)
length(GR_bloodHyperFibrosis) #0

# Hypo
pdf(file="Figures/Fibrosis Blood Hypo DMR Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_fibrosisHypo_hg38, GR_bloodHypo), NameOfPeaks = c("Fibrosis", "WD_Blood"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), totalTest = length(GR_fibrosis_blood_back), cat.fontfamily="sans",fontfamily="sans")
dev.off()

# Hypergeometric Test pvalue
as.data.frame(venn$p.value)$pval #0.04936673

GR_bloodHypoFibrosis <- subsetByOverlaps(GR_bloodHypo, GR_fibrosisHypo_hg38)
bloodHypoFibrosis <- as.data.frame(GR_bloodHypoFibrosis)
write.table(bloodHypoFibrosis, "WD Blood Hypo DMRs Overlapping with Fibrosis DMRs.txt", sep="\t", quote=FALSE, row.names=FALSE)
