# WD Liver TF DMR Overlap Analysis ####
# Charles Mordaunt
# 1/16/18

# Packages ####
library(GenomicRanges)
library(VennDiagram)
library(ChIPpeakAnno)
library(ggplot2)
library(scales)
library(reshape2)

# Load Data ####
dmrs <- read.delim("WD_Specific_DMRs.bed", sep="\t", header=FALSE)
colnames(dmrs) <- c("chr", "start", "end")

dmrs_hyper <- read.delim("WD_Specific_hyper_DMRs.bed", sep="\t", header=FALSE)
colnames(dmrs_hyper) <- c("chr", "start", "end")

dmrs_hypo <- read.delim("WD_Specific_hypo_DMRs.bed", sep="\t", header=FALSE)
colnames(dmrs_hypo) <- c("chr", "start", "end")

hnf4a <- read.delim("C:/Users/Booboo/Charles/Documents/Programming/LOLA/ENCODE_Liver_TF_ChIPseq/ENCFF639EJK.bed", sep="\t", header=FALSE)
hnf4a <- hnf4a[,1:3]
colnames(hnf4a) <- c("chr", "start", "end")
write.table(hnf4a, "Tables/HNF4A Liver ChIPseq.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

rxra <- read.delim("C:/Users/Booboo/Charles/Documents/Programming/LOLA/ENCODE_Liver_TF_ChIPseq/ENCFF420GLO.bed", sep="\t", header=FALSE)
rxra <- rxra[,1:3]
colnames(rxra) <- c("chr", "start", "end")
write.table(rxra, "Tables/RXRA Liver ChIPseq.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

foxa1 <- read.delim("C:/Users/Booboo/Charles/Documents/Programming/LOLA/ENCODE_Liver_TF_ChIPseq/ENCFF741PQE.bed", sep="\t", header=FALSE)
foxa1 <- foxa1[,1:3]
colnames(foxa1) <- c("chr", "start", "end")
# Used male sample
write.table(foxa1, "Tables/FOXA1 Liver ChIPseq.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

foxa2 <- read.delim("C:/Users/Booboo/Charles/Documents/Programming/LOLA/ENCODE_Liver_TF_ChIPseq/ENCFF168JTA.bed", sep="\t", header=FALSE)
foxa2 <- foxa2[,1:3]
colnames(foxa2) <- c("chr", "start", "end")
# Used male sample
write.table(foxa2, "Tables/FOXA2 Liver ChIPseq.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

dmrInfo <- read.delim("WD Specific DMR Stats with annotation 2.txt", sep="\t", stringsAsFactors = FALSE)
dmrInfo <- dmrInfo[,1:13]
dmrInfo <- base::unique(dmrInfo)

# Make Genomic Ranges ####
GR_dmrs <- GRanges(seqnames = dmrs$chr, ranges=IRanges(start=dmrs$start, end=dmrs$end))
GR_dmrs_hyper <- GRanges(seqnames = dmrs_hyper$chr, ranges=IRanges(start=dmrs_hyper$start, end=dmrs_hyper$end))
GR_dmrs_hypo <- GRanges(seqnames = dmrs_hypo$chr, ranges=IRanges(start=dmrs_hypo$start, end=dmrs_hypo$end))
GR_hnf4a <- GRanges(seqnames = hnf4a$chr, ranges=IRanges(start=hnf4a$start, end=hnf4a$end))
GR_rxra <- GRanges(seqnames = rxra$chr, ranges=IRanges(start=rxra$start, end=rxra$end))
GR_foxa1 <- GRanges(seqnames = foxa1$chr, ranges=IRanges(start=foxa1$start, end=foxa1$end))
GR_foxa2 <- GRanges(seqnames = foxa2$chr, ranges=IRanges(start=foxa2$start, end=foxa2$end))

# Make Sure Regions are Disjoint ####
if(!isDisjoint(GR_dmrs)){GR_dmrs <- disjoin(GR_dmrs)} 
if(!isDisjoint(GR_dmrs_hyper)){GR_dmrs_hyper <- disjoin(GR_dmrs_hyper)} 
if(!isDisjoint(GR_dmrs_hypo)){GR_dmrs_hypo <- disjoin(GR_dmrs_hypo)} 
if(!isDisjoint(GR_hnf4a)){GR_hnf4a <- disjoin(GR_hnf4a)} 
if(!isDisjoint(GR_rxra)){GR_rxra <- disjoin(GR_rxra)} 
if(!isDisjoint(GR_foxa1)){GR_foxa1 <- disjoin(GR_foxa1)} 
if(!isDisjoint(GR_foxa2)){GR_foxa2 <- disjoin(GR_foxa2)} 

# Overlaps ####

# All TF Binding Sites Venn Diagram
#pdf(file="Figures/All TF Binding Site Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_hnf4a, GR_rxra, GR_foxa1, GR_foxa2), NameOfPeaks = c("HNF4A", "RXRA", "FOXA1", "FOXA2"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen", "orange"), cat.pos = c(0,0,0,0), 
                        cat.dist = c(0.05, 0.05, 0.1, 0.1))
#dev.off()
hyper_results <- as.data.frame(venn$p.value)


# All DMRs
GR_hnf4a_dmrs <- subsetByOverlaps(query=GR_dmrs, subject=GR_hnf4a) #463 DMRs
GR_rxra_dmrs <- subsetByOverlaps(query=GR_dmrs, subject=GR_rxra) #462 DMRs
GR_foxa1_dmrs <- subsetByOverlaps(query=GR_dmrs, subject=GR_foxa1) #138 DMRs
GR_foxa2_dmrs <- subsetByOverlaps(query=GR_dmrs, subject=GR_foxa2) #171 DMRs

pdf(file="Figures/All TF DMRs Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_hnf4a_dmrs, GR_rxra_dmrs, GR_foxa1_dmrs, GR_foxa2_dmrs), NameOfPeaks = c("HNF4A", "RXRA", "FOXA1", "FOXA2"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen", "orange"), cat.pos = c(0,0,0,0), 
                        cat.dist = c(0.05, 0.05, 0.1, 0.1), totalTest = nrow(dmrs))
dev.off()
hyper_results_dmrs <- as.data.frame(venn$p.value)
hyper_results_dmrs$TF1 <- factor(c("FOXA1", "RXRA", "RXRA", "HNF4A", "HNF4A", "HNF4A"), levels=c("HNF4A", "RXRA", "FOXA1", "FOXA2"), ordered=TRUE)
hyper_results_dmrs$TF2 <- factor(c("FOXA2", "FOXA2", "FOXA1", "FOXA2", "FOXA1", "RXRA"), levels=c("HNF4A", "RXRA", "FOXA1", "FOXA2"), ordered=TRUE)
hyper_results_dmrs$logPvalue <- -log10(hyper_results_dmrs$pval)
hyper_results_dmrs$qval <- p.adjust(hyper_results_dmrs$pval, method="fdr")
hyper_results_dmrs$logQvalue <- -log10(hyper_results_dmrs$qval)

# logqvalue Heatmap
gg <- ggplot(data = hyper_results_dmrs)
gg +
        geom_tile(aes(x = TF1, y = TF2, fill = logQvalue)) +
        scale_fill_gradientn("-log(q-value)\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits=c(0,max(hyper_results_dmrs$logQvalue))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.25, 0.77), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 15, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_blank(), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/TF DMR Overlap All DMR Hypergeometric Heatmap.png", dpi = 600, width = 6, height = 5, units = "in")

pdf(file="Figures/All TF DMRs Overlap Venn with DMRs.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_hnf4a_dmrs, GR_rxra_dmrs, GR_foxa1_dmrs, GR_foxa2_dmrs, GR_dmrs), 
                        NameOfPeaks = c("HNF4A", "RXRA", "FOXA1", "FOXA2", "All DMRs"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen", "orange", "purple"), 
                        cat.pos = c(0,0,0,0,0), cat.dist = c(0.05, 0.05, 0.1, 0.1, 0.05))
dev.off()

# Hyper DMRs
GR_hnf4a_dmrs_hyper <- subsetByOverlaps(query=GR_dmrs_hyper, subject=GR_hnf4a) #447 dmrs_hyper
GR_rxra_dmrs_hyper <- subsetByOverlaps(query=GR_dmrs_hyper, subject=GR_rxra) #432 dmrs_hyper
GR_foxa1_dmrs_hyper <- subsetByOverlaps(query=GR_dmrs_hyper, subject=GR_foxa1) #135 dmrs_hyper
GR_foxa2_dmrs_hyper <- subsetByOverlaps(query=GR_dmrs_hyper, subject=GR_foxa2) #165 dmrs_hyper
GR_allTF_dmrs_hyper <- subsetByOverlaps(subsetByOverlaps(subsetByOverlaps(GR_hnf4a_dmrs_hyper, GR_rxra_dmrs_hyper), 
                                                         GR_foxa1_dmrs_hyper), GR_foxa2_dmrs_hyper) # 82 dmrs with all TFs
allTF_dmrs_hyper <- as.data.frame(GR_allTF_dmrs_hyper)
allTF_dmrs_hyper <- merge(x=allTF_dmrs_hyper, y=dmrInfo, by.x=c("seqnames", "start", "end"), by.y=c("chr", "start", "end"), 
                          all.x=TRUE, all.y=FALSE)

pdf(file="Figures/All TF Hyper DMRs Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_hnf4a_dmrs_hyper, GR_rxra_dmrs_hyper, GR_foxa1_dmrs_hyper, GR_foxa2_dmrs_hyper), NameOfPeaks = c("HNF4A", "RXRA", "FOXA1", "FOXA2"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen", "orange"), cat.pos = c(0,0,0,0), 
                        cat.dist = c(0.05, 0.05, 0.1, 0.1), totalTest = nrow(dmrs_hyper), fontfamily="sans", cat.fontfamily="sans")
dev.off()
hyper_results_dmrs_hyper <- as.data.frame(venn$p.value)
hyper_results_dmrs_hyper$TF1 <- factor(c("FOXA1", "RXRA", "RXRA", "HNF4A", "HNF4A", "HNF4A"), levels=c("HNF4A", "RXRA", "FOXA1", "FOXA2"), ordered=TRUE)
hyper_results_dmrs_hyper$TF2 <- factor(c("FOXA2", "FOXA2", "FOXA1", "FOXA2", "FOXA1", "RXRA"), levels=c("HNF4A", "RXRA", "FOXA1", "FOXA2"), ordered=TRUE)
hyper_results_dmrs_hyper$logPvalue <- -log10(hyper_results_dmrs_hyper$pval)
hyper_results_dmrs_hyper$qval <- p.adjust(hyper_results_dmrs_hyper$pval, method="fdr")
hyper_results_dmrs_hyper$logQvalue <- -log10(hyper_results_dmrs_hyper$qval)

# logqvalue Heatmap
gg <- ggplot(data = hyper_results_dmrs_hyper)
gg +
        geom_tile(aes(x = TF1, y = TF2, fill = logQvalue)) +
        scale_fill_gradientn("-log(q-value)\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits=c(0,max(hyper_results_dmrs_hyper$logQvalue))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.25, 0.77), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 15, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_blank(), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/TF DMR Overlap Hyper DMR Hypergeometric Heatmap.png", dpi = 600, width = 6, height = 5, units = "in")

pdf(file="Figures/All TF Hyper DMRs Overlap Venn with DMRs.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_hnf4a_dmrs_hyper, GR_rxra_dmrs_hyper, GR_foxa1_dmrs_hyper, GR_foxa2_dmrs_hyper, GR_dmrs_hyper), 
                        NameOfPeaks = c("HNF4A", "RXRA", "FOXA1", "FOXA2", "Hyper DMRs"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen", "orange", "purple"), cat.pos = c(0,0,0,0,0), 
                        cat.dist = c(0.05, 0.05, 0.1, 0.1, 0.05))
dev.off()

# Hypo DMRs
GR_hnf4a_dmrs_hypo <- subsetByOverlaps(query=GR_dmrs_hypo, subject=GR_hnf4a) #16 dmrs_hypo
GR_rxra_dmrs_hypo <- subsetByOverlaps(query=GR_dmrs_hypo, subject=GR_rxra) #30 dmrs_hypo
GR_foxa1_dmrs_hypo <- subsetByOverlaps(query=GR_dmrs_hypo, subject=GR_foxa1) #3 dmrs_hypo
GR_foxa2_dmrs_hypo <- subsetByOverlaps(query=GR_dmrs_hypo, subject=GR_foxa2) #6 dmrs_hypo

pdf(file="Figures/All TF Hypo DMRs Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_hnf4a_dmrs_hypo, GR_rxra_dmrs_hypo, GR_foxa1_dmrs_hypo, GR_foxa2_dmrs_hypo), NameOfPeaks = c("HNF4A", "RXRA", "FOXA1", "FOXA2"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen", "orange"), cat.pos = c(0,0,0,0), 
                        cat.dist = c(0.05, 0.05, 0.1, 0.1), totalTest = nrow(dmrs_hypo))
dev.off()

hyper_results_dmrs_hypo <- as.data.frame(venn$p.value)
hyper_results_dmrs_hypo$TF1 <- factor(c("FOXA1", "RXRA", "RXRA", "HNF4A", "HNF4A", "HNF4A"), levels=c("HNF4A", "RXRA", "FOXA1", "FOXA2"), ordered=TRUE)
hyper_results_dmrs_hypo$TF2 <- factor(c("FOXA2", "FOXA2", "FOXA1", "FOXA2", "FOXA1", "RXRA"), levels=c("HNF4A", "RXRA", "FOXA1", "FOXA2"), ordered=TRUE)
hyper_results_dmrs_hypo$logPvalue <- -log10(hyper_results_dmrs_hypo$pval)
hyper_results_dmrs_hypo$qval <- p.adjust(hyper_results_dmrs_hypo$pval, method="fdr")
hyper_results_dmrs_hypo$logQvalue <- -log10(hyper_results_dmrs_hypo$qval)

# logqvalue Heatmap
gg <- ggplot(data = hyper_results_dmrs_hypo)
gg +
        geom_tile(aes(x = TF1, y = TF2, fill = logQvalue)) +
        scale_fill_gradientn("-log(q-value)\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", 
                             limits=c(0,max(hyper_results_dmrs_hypo$logQvalue))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.25, 0.77), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 15, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_blank(), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/TF DMR Overlap Hypo DMR Hypergeometric Heatmap.png", dpi = 600, width = 6, height = 5, units = "in")

pdf(file="Figures/All TF Hypo DMRs Overlap Venn with DMRs.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_hnf4a_dmrs_hypo, GR_rxra_dmrs_hypo, GR_foxa1_dmrs_hypo, GR_foxa2_dmrs_hypo, GR_dmrs_hypo), 
                        NameOfPeaks = c("HNF4A", "RXRA", "FOXA1", "FOXA2", "Hypo DMRs"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen", "orange", "purple"), 
                        cat.pos = c(0,0,0,0,0), cat.dist = c(0.05, 0.05, 0.1, 0.1, 0.05))
dev.off()