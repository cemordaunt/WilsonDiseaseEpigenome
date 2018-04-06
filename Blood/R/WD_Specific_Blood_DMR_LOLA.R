# WD Blood LOLA Enrichment Analysis ####
# WD Specific DMRs
# Charles Mordaunt
# 3/21/18

# Load Packages ############
library(ggplot2)
library(RColorBrewer)
library(scales)
library(reshape2)

# Load Data ####
index <- read.delim(file = "C:/Users/Booboo/Charles/Documents/Programming/MARBLES Cord Blood/Tables/LOLA/LOLA Roadmap ChromHMM index.txt", sep = "\t", header = TRUE)
stats <- read.delim(file = "Tables/WD_Blood_2_Specific_DMRs_LOLA.tsv", sep = "\t", header = TRUE)
stats$userSet[stats$userSet == 1] <- "allDMRs"
stats$userSet[stats$userSet == 2] <- "hyperDMRs"
stats$userSet[stats$userSet == 3] <- "hypoDMRs"

# All DMRs ####
chromhmm_all <- subset(stats, collection == "Roadmap_ChromHMM" & userSet == "allDMRs")
roadmap_all <- subset(stats, collection == "roadmap_epigenomics" & userSet == "allDMRs")

# Make Stats Table for ChromHMM ####
# Merge Tables
match <- match(chromhmm_all$filename, index$filename)
chromhmm_all$type <- index$type[match]
chromhmm_all$order <- index$order[match]
chromhmm_all$color <- index$color[match]
chromhmm_all <- chromhmm_all[,c("dbSet", "antibody", "cellType", "tissue", "type", "pValueLog", "qValue", "logOddsRatio", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", 
                                "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]
colnames(chromhmm_all)[2] <- "chromState"

# Adjust chromhmm_all for region sets with < 5 overlaps 
table(chromhmm_all$support < 5) # 483 region sets 
chromhmm_all$logOddsRatio[chromhmm_all$support < 5] <- 0
chromhmm_all$pValueLog[chromhmm_all$support < 5] <- 0
chromhmm_all$qValue[chromhmm_all$support < 5] <- 1

# Replace inf logP with 1.5*max
table(is.na(chromhmm_all$pValueLog)) # All false
chromhmm_all$pValueLog[is.infinite(chromhmm_all$pValueLog)] <- NA
max_logp <- max(chromhmm_all$pValueLog, na.rm = TRUE)
chromhmm_all$pValueLog[is.na(chromhmm_all$pValueLog)] <- 1.5*max_logp

# Adjust pValue
chromhmm_all$pValue <- 10^(-chromhmm_all$pValueLog)
chromhmm_all$qValueLog <- -log10(chromhmm_all$qValue)

# Replace qValueLog with 1.5*max
table(is.na(chromhmm_all$qValueLog)) # All false
chromhmm_all$qValueLog[is.infinite(chromhmm_all$qValueLog)] <- NA
max_logq <- max(chromhmm_all$qValueLog, na.rm = TRUE)
chromhmm_all$qValueLog[is.na(chromhmm_all$qValueLog)] <- 1.5*max_logq

rm(match, max_logp, max_logq)

# Add percent DMRs
numDMRs <- chromhmm_all$support[1] + chromhmm_all$c[1]
chromhmm_all$pct_DMRs <- chromhmm_all$support*100/numDMRs

# Rearrange Columns
chromhmm_all <- chromhmm_all[,c("dbSet", "chromState", "cellType", "tissue", "type", "pValue", "qValue", "pValueLog", "qValueLog", "logOddsRatio", "support", 
                                "pct_DMRs", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]
chromhmm_all$order <- factor(chromhmm_all$order, levels = sort(unique(chromhmm_all$order), decreasing = TRUE), ordered = TRUE)
chromhmm_all$chromState <- sapply(strsplit(as.character(chromhmm_all$chromState), split='_', fixed=TRUE), function(x) (x[2]))
chromhmm_all$chromState <- factor(chromhmm_all$chromState, levels = c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZnfRpts", "Het", "TssBiv", "BivFlnk",
                                                                      "EnhBiv", "ReprPC", "ReprPCwk", "Quies"), ordered = TRUE)
chromhmm_all <- chromhmm_all[order(chromhmm_all$order, chromhmm_all$chromState),]
tissues <- rev(as.character(unique(chromhmm_all$tissue)))
chromhmm_all$tissue <- factor(chromhmm_all$tissue, levels = tissues, ordered = TRUE)
tissue_colors <- rev(as.character(unique(chromhmm_all$color)))

write.table(x = chromhmm_all, file = "Tables/WD_Blood_2_Specific_DMRs_All_ChromHMM_Enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Make Stats table for Roadmap ####
# Merge Tables
EID <- strsplit(as.character(roadmap_all$filename), "-")
EID <- sapply(EID, function(x){x[1]})
roadmap_all$EID <- EID
match <- match(roadmap_all$EID, index$EID)
roadmap_all$type <- index$type[match]
roadmap_all$order <- index$order[match]
roadmap_all$color <- index$color[match]
roadmap_all$cellType <- index$cellType[match]
roadmap_all$tissue <- index$tissue[match]

roadmap_all$antibody <- as.character(roadmap_all$antibody)
roadmap_all <- subset(roadmap_all, antibody %in% c("H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"))
roadmap_all <- roadmap_all[,c("dbSet", "antibody", "cellType", "tissue", "type", "pValueLog", "qValue", "logOddsRatio", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", 
                              "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]

# Adjust roadmap_all for region sets with < 5 overlaps 
table(roadmap_all$support < 5) # 62 region sets 
roadmap_all$logOddsRatio[roadmap_all$support < 5] <- 0
roadmap_all$pValueLog[roadmap_all$support < 5] <- 0
roadmap_all$qValue[roadmap_all$support < 5] <- 1

# Replace inf logP with 1.5*max
table(is.na(roadmap_all$pValueLog)) # All false
roadmap_all$pValueLog[is.infinite(roadmap_all$pValueLog)] <- NA
max_logp <- max(roadmap_all$pValueLog, na.rm = TRUE)
roadmap_all$pValueLog[is.na(roadmap_all$pValueLog)] <- 1.5*max_logp

# Adjust pValue
roadmap_all$pValue <- 10^(-roadmap_all$pValueLog)
roadmap_all$qValueLog <- -log10(roadmap_all$qValue)

# Replace qValueLog with 1.5*max
table(is.na(roadmap_all$qValueLog)) # All false
roadmap_all$qValueLog[is.infinite(roadmap_all$qValueLog)] <- NA
max_logq <- max(roadmap_all$qValueLog, na.rm = TRUE)
roadmap_all$qValueLog[is.na(roadmap_all$qValueLog)] <- 1.5*max_logq

rm(match, max_logp, max_logq)

# Add percent DMRs
numDMRs <- roadmap_all$support[1] + roadmap_all$c[1]
roadmap_all$pct_DMRs <- roadmap_all$support*100/numDMRs

# Rearrange Columns
roadmap_all <- roadmap_all[,c("dbSet", "antibody", "cellType", "tissue", "type", "pValue", "qValue", "pValueLog", "qValueLog", "logOddsRatio", "support", 
                              "pct_DMRs", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]
roadmap_all$order <- factor(roadmap_all$order, levels = sort(unique(roadmap_all$order), decreasing = TRUE), ordered = TRUE)
roadmap_all$antibody <- factor(roadmap_all$antibody, levels = c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"), ordered = TRUE)
roadmap_all <- roadmap_all[order(roadmap_all$order, roadmap_all$antibody),]
tissues <- rev(as.character(unique(roadmap_all$tissue)))
roadmap_all$tissue <- factor(roadmap_all$tissue, levels = tissues, ordered = TRUE)
tissue_colors <- rev(as.character(unique(roadmap_all$color)))

write.table(x = roadmap_all, file = "Tables/WD_Blood_2_Specific_DMRs_all_Roadmap_Enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Make Stats Table for Liver TFs ####
liver_TFs_all <- subset(stats, collection == "encode_liver_chipseq" & userSet == "allDMRs")

# Adjust liver_TFs_all for region sets with < 5 overlaps 
table(liver_TFs_all$support < 5) # 5 region sets 
liver_TFs_all$logOddsRatio[liver_TFs_all$support < 5] <- 0
liver_TFs_all$pValueLog[liver_TFs_all$support < 5] <- 0
liver_TFs_all$qValue[liver_TFs_all$support < 5] <- 1

# Replace inf logP with 1.5*max
table(is.na(liver_TFs_all$pValueLog)) # All false
table(is.infinite(liver_TFs_all$pValueLog)) # All false

# Adjust pValue
liver_TFs_all$pValue <- 10^(-liver_TFs_all$pValueLog)
liver_TFs_all$qValueLog <- -log10(liver_TFs_all$qValue)

# Replace qValueLog with 1.5*max
table(is.na(liver_TFs_all$qValueLog)) # All false
table(is.infinite(liver_TFs_all$qValueLog)) # All false

# Add percent DMRs
numDMRs <- liver_TFs_all$support[1] + liver_TFs_all$c[1]
liver_TFs_all$pct_DMRs <- liver_TFs_all$support*100/numDMRs
backDMRs <- sum(liver_TFs_all$support[1], liver_TFs_all$b[1], liver_TFs_all$c[1], liver_TFs_all$d[1]) 
liver_TFs_all$pct_BackDMRs <- (liver_TFs_all$support + liver_TFs_all$b)*100/backDMRs

# Reorder columns
liver_TFs_all <- liver_TFs_all[,c("dbSet", "antibody", "cellType", "pValue", "qValue", "pValueLog", "qValueLog", "logOddsRatio", "support", 
                                  "pct_DMRs", "pct_BackDMRs", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "size", "filename")]

# Fix Antibody Names
liver_TFs_all$antibody <- factor(liver_TFs_all$antibody, levels=sort(unique(liver_TFs_all$antibody)), ordered=TRUE)
write.table(x = liver_TFs_all, file = "Tables/WD_Blood_2_Specific_DMRs_all_LOLA_Liver_TF_Enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Chromhmm Figures ####################

# logOR Heatmap
gg <- ggplot(data = chromhmm_all)
gg +
        geom_tile(aes(x = chromState, y = order, fill = logOddsRatio)) +
        scale_fill_gradientn("-logOR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,20)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR ChromHMM logOR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# logFDR Heatmap
gg <- ggplot(data = chromhmm_all)
gg +
        geom_tile(aes(x = chromState, y = order, fill = qValueLog)) +
        scale_fill_gradientn("-logFDR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,58)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR ChromHMM logFDR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# Heatmap Legend
gg <- ggplot(data = chromhmm_all)
gg +
        geom_tile(aes(x = 16, y = order, fill = tissue)) +
        scale_fill_manual(name = "Tissue", values = tissue_colors) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(4.7, 0.21), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,22,6.3,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_blank(),
              axis.title = element_blank(), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR ChromHMM legend by tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# Stacked Barchart of percent each chromstate
percent <- aggregate(x = chromhmm_all$pct_DMRs, by = list(chromhmm_all$chromState), FUN = "mean")
colnames(percent) <- c("chromState", "Mean_Percent")
chromColors <- c(rgb(255,0,0, maxColorValue = 255), rgb(255,69,0, maxColorValue = 255), rgb(50,205,50, maxColorValue = 255), 
                 rgb(0,128,0, maxColorValue = 255), rgb(0,100,0, maxColorValue = 255), rgb(194,225,5, maxColorValue = 255), 
                 rgb(255,255,0, maxColorValue = 255), rgb(102,205,170, maxColorValue = 255), rgb(138,145,208, maxColorValue = 255), 
                 rgb(205,92,92, maxColorValue = 255), rgb(233,150,122, maxColorValue = 255), rgb(189,183,107, maxColorValue = 255),
                 rgb(128,128,128, maxColorValue = 255), rgb(192,192,192, maxColorValue = 255), rgb(255,255,255, maxColorValue = 255))

backDMRs <- sum(chromhmm_all$support[1], chromhmm_all$b[1], chromhmm_all$c[1], chromhmm_all$d[1]) 
percentBackDMRs <- (chromhmm_all$support + chromhmm_all$b)*100/backDMRs
meanBackDMRs <- aggregate(x = percentBackDMRs, by = list(chromhmm_all$chromState), FUN = "mean")
colnames(meanBackDMRs) <- c("chromState", "Background")
meanBackDMRs$DMRs <- percent$Mean_Percent
meanBackDMRs$DMRs <- meanBackDMRs$DMRs*100/sum(meanBackDMRs$DMRs)
meanBackDMRs$Background <- meanBackDMRs$Background*100/sum(meanBackDMRs$Background)

percentDMRs_m <- melt(meanBackDMRs, id.vars="chromState")
percentDMRs_m$variable <- factor(percentDMRs_m$variable, levels=c("Background", "DMRs"), ordered=TRUE)

gg <- ggplot(data = percentDMRs_m)
gg +
        geom_col(aes(x = variable, y = value, fill = chromState), color = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.45, 0.64), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,12,1,0.5), "lines"), 
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 19, color = "black"),
              axis.ticks = element_line(size = 0.6, color = "black"),
              axis.text.x = element_text(angle=45, hjust=0.9),
              axis.title.x = element_blank(), 
              legend.title = element_text(size = 18)) +
        ylab("Mean % DMRs Across Tissues") +
        scale_fill_manual(name = "Chromatin State", values=chromColors) +
        scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(expand = c(0,0))
ggsave("Figures/LOLA WD Blood 2 Specific DMR and Background ChromHMM Mean Percent chromState.png", dpi = 600, width = 6, height = 7, units = "in")

# Stacked Barchart of percent each chromstate in liver

liver <- chromhmm_all[chromhmm_all$cellType == "Liver", c("chromState", "support", "b", "c", "d")]
liver$DMRs <- liver$support*100/numDMRs
liver$DMRs <- liver$DMRs*100/sum(liver$DMRs)
liver$Background <- (liver$support + liver$b)*100/backDMRs
liver$Background <- liver$Background*100/sum(liver$Background)

liver_m <- melt(liver, id.vars=c("chromState", "support", "b", "c", "d"))
liver_m$variable <- factor(liver_m$variable, levels=c("Background", "DMRs"), ordered=TRUE)

gg <- ggplot(data = liver_m)
gg +
        geom_col(aes(x = variable, y = value, fill = chromState), color = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.45, 0.64), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,12,1,0.5), "lines"), 
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 19, color = "black"),
              axis.ticks = element_line(size = 0.6, color = "black"),
              axis.text.x = element_text(angle=45, hjust=0.9),
              axis.title.x = element_blank(), 
              legend.title = element_text(size = 18)) +
        ylab("% DMRs in Liver") +
        scale_fill_manual(name = "Chromatin State", values=chromColors) +
        scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(expand = c(0,0))
ggsave("Figures/LOLA WD Blood 2 Specific DMR and Background ChromHMM Liver Percent chromState.png", dpi = 600, width = 6, height = 7, units = "in")

# Stacked Barchart of percent each chromstate in HSC & Bcell

HSCB <- chromhmm_all[chromhmm_all$tissue == "HSC & B-cell", c("chromState", "support", "b", "c", "d")]
HSCB$DMRs <- HSCB$support*100/numDMRs
HSCB$Background <- (HSCB$support + HSCB$b)*100/backDMRs

percent <- aggregate(x = HSCB$DMRs, by = list(HSCB$chromState), FUN = "mean")
colnames(percent) <- c("chromState", "Mean_Percent")

meanBackDMRs <- aggregate(x = HSCB$Background, by = list(HSCB$chromState), FUN = "mean")
colnames(meanBackDMRs) <- c("chromState", "Background")

meanBackDMRs$DMRs <- percent$Mean_Percent
meanBackDMRs$DMRs <- meanBackDMRs$DMRs*100/sum(meanBackDMRs$DMRs)
meanBackDMRs$Background <- meanBackDMRs$Background*100/sum(meanBackDMRs$Background)

percentDMRs_m <- melt(meanBackDMRs, id.vars="chromState")
percentDMRs_m$variable <- factor(percentDMRs_m$variable, levels=c("Background", "DMRs"), ordered=TRUE)

gg <- ggplot(data = percentDMRs_m)
gg +
        geom_col(aes(x = variable, y = value, fill = chromState), color = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.45, 0.64), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,12,1,0.5), "lines"), 
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 19, color = "black"),
              axis.ticks = element_line(size = 0.6, color = "black"),
              axis.text.x = element_text(angle=45, hjust=0.9),
              axis.title.x = element_blank(), 
              legend.title = element_text(size = 18)) +
        ylab("% DMRs in HSC & B-cells") +
        scale_fill_manual(name = "Chromatin State", values=chromColors) +
        scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(expand = c(0,0))
ggsave("Figures/LOLA WD Blood 2 Specific DMR and Background ChromHMM HSC B-cells Percent chromState.png", dpi = 600, width = 6, height = 7, units = "in")

# roadmap_all Figures ####################

# logOR Heatmap
gg <- ggplot(data = roadmap_all)
gg +
        geom_tile(aes(x = antibody, y = order, fill = logOddsRatio)) +
        scale_fill_gradientn("-logOR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,16)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR roadmap_all logOR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# logFDR Heatmap
gg <- ggplot(data = roadmap_all)
gg +
        geom_tile(aes(x = antibody, y = order, fill = qValueLog)) +
        scale_fill_gradientn("-logFDR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,75)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR roadmap_all logFDR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# Liver_TFs_all Figures ####
# Volcano Plot
top <- data.frame(antibody=liver_TFs_all$antibody, logOddsRatio=liver_TFs_all$logOddsRatio, qValueLog=liver_TFs_all$qValueLog)
top <- subset(top, rank(-top$qValueLog) <= 15)
top$antibody <- as.character(top$antibody)

gg <- ggplot(data = liver_TFs_all)
gg +
        geom_point(aes(x = logOddsRatio, y = qValueLog), color="#3366CC", size=2) +
        geom_text(data=top, aes(x=logOddsRatio, y=qValueLog, label=antibody), 
                  size=4, check_overlap = TRUE, hjust=0, nudge_x=0.05, nudge_y=0) +
        theme_bw(base_size = 24) +
        labs(x="-log(Odds Ratio)", y="-log(q-value)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title = element_text(size=15, color="Black"), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) +
        scale_x_continuous(breaks=pretty_breaks(n=5)) +
        coord_cartesian(xlim=c(0,3.8))
ggsave("Figures/LOLA WD Blood 2 Specific DMR Liver TF Volcano Plot.png", dpi = 600, width = 7, height = 5, units = "in")

# Hyper DMRs ####
chromhmm_hyper <- subset(stats, collection == "Roadmap_ChromHMM" & userSet == "hyperDMRs")
roadmap_hyper <- subset(stats, collection == "roadmap_epigenomics" & userSet == "hyperDMRs")

# Make Stats Table for ChromHMM ####
# Merge Tables
match <- match(chromhmm_hyper$filename, index$filename)
chromhmm_hyper$type <- index$type[match]
chromhmm_hyper$order <- index$order[match]
chromhmm_hyper$color <- index$color[match]
chromhmm_hyper <- chromhmm_hyper[,c("dbSet", "antibody", "cellType", "tissue", "type", "pValueLog", "qValue", "logOddsRatio", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", 
                                    "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]
colnames(chromhmm_hyper)[2] <- "chromState"

# Adjust chromhmm_hyper for region sets with < 5 overlaps 
table(chromhmm_hyper$support < 5) # 596 region sets 
chromhmm_hyper$logOddsRatio[chromhmm_hyper$support < 5] <- 0
chromhmm_hyper$pValueLog[chromhmm_hyper$support < 5] <- 0
chromhmm_hyper$qValue[chromhmm_hyper$support < 5] <- 1

# Replace inf logP with 1.5*max
table(is.na(chromhmm_hyper$pValueLog)) # All false
chromhmm_hyper$pValueLog[is.infinite(chromhmm_hyper$pValueLog)] <- NA
max_logp <- max(chromhmm_hyper$pValueLog, na.rm = TRUE)
chromhmm_hyper$pValueLog[is.na(chromhmm_hyper$pValueLog)] <- 1.5*max_logp

# Adjust pValue
chromhmm_hyper$pValue <- 10^(-chromhmm_hyper$pValueLog)
chromhmm_hyper$qValueLog <- -log10(chromhmm_hyper$qValue)

# Replace qValueLog with 1.5*max
table(is.na(chromhmm_hyper$qValueLog)) # All false
chromhmm_hyper$qValueLog[is.infinite(chromhmm_hyper$qValueLog)] <- NA
max_logq <- max(chromhmm_hyper$qValueLog, na.rm = TRUE)
chromhmm_hyper$qValueLog[is.na(chromhmm_hyper$qValueLog)] <- 1.5*max_logq

rm(match, max_logp, max_logq)

# Add percent DMRs
numDMRs <- chromhmm_hyper$support[1] + chromhmm_hyper$c[1]
chromhmm_hyper$pct_DMRs <- chromhmm_hyper$support*100/numDMRs

# Rearrange Columns
chromhmm_hyper <- chromhmm_hyper[,c("dbSet", "chromState", "cellType", "tissue", "type", "pValue", "qValue", "pValueLog", "qValueLog", "logOddsRatio", "support", 
                                    "pct_DMRs", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]
chromhmm_hyper$order <- factor(chromhmm_hyper$order, levels = sort(unique(chromhmm_hyper$order), decreasing = TRUE), ordered = TRUE)
chromhmm_hyper$chromState <- sapply(strsplit(as.character(chromhmm_hyper$chromState), split='_', fixed=TRUE), function(x) (x[2]))
chromhmm_hyper$chromState <- factor(chromhmm_hyper$chromState, levels = c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZnfRpts", "Het", "TssBiv", "BivFlnk",
                                                                          "EnhBiv", "ReprPC", "ReprPCwk", "Quies"), ordered = TRUE)
chromhmm_hyper <- chromhmm_hyper[order(chromhmm_hyper$order, chromhmm_hyper$chromState),]
tissues <- rev(as.character(unique(chromhmm_hyper$tissue)))
chromhmm_hyper$tissue <- factor(chromhmm_hyper$tissue, levels = tissues, ordered = TRUE)
tissue_colors <- rev(as.character(unique(chromhmm_hyper$color)))

write.table(x = chromhmm_hyper, file = "Tables/WD_Blood_2_Specific_DMRs_Hyper_ChromHMM_Enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Make Stats table for Roadmap ####
# Merge Tables
EID <- strsplit(as.character(roadmap_hyper$filename), "-")
EID <- sapply(EID, function(x){x[1]})
roadmap_hyper$EID <- EID
match <- match(roadmap_hyper$EID, index$EID)
roadmap_hyper$type <- index$type[match]
roadmap_hyper$order <- index$order[match]
roadmap_hyper$color <- index$color[match]
roadmap_hyper$cellType <- index$cellType[match]
roadmap_hyper$tissue <- index$tissue[match]

roadmap_hyper$antibody <- as.character(roadmap_hyper$antibody)
roadmap_hyper <- subset(roadmap_hyper, antibody %in% c("H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"))
roadmap_hyper <- roadmap_hyper[,c("dbSet", "antibody", "cellType", "tissue", "type", "pValueLog", "qValue", "logOddsRatio", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", 
                                  "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]

# Adjust roadmap_hyper for region sets with < 5 overlaps 
table(roadmap_hyper$support < 5) # 77 region sets 
roadmap_hyper$logOddsRatio[roadmap_hyper$support < 5] <- 0
roadmap_hyper$pValueLog[roadmap_hyper$support < 5] <- 0
roadmap_hyper$qValue[roadmap_hyper$support < 5] <- 1

# Replace inf logP with 1.5*max
table(is.na(roadmap_hyper$pValueLog)) # All false
roadmap_hyper$pValueLog[is.infinite(roadmap_hyper$pValueLog)] <- NA
max_logp <- max(roadmap_hyper$pValueLog, na.rm = TRUE)
roadmap_hyper$pValueLog[is.na(roadmap_hyper$pValueLog)] <- 1.5*max_logp

# Adjust pValue
roadmap_hyper$pValue <- 10^(-roadmap_hyper$pValueLog)
roadmap_hyper$qValueLog <- -log10(roadmap_hyper$qValue)

# Replace qValueLog with 1.5*max
table(is.na(roadmap_hyper$qValueLog)) # All false
roadmap_hyper$qValueLog[is.infinite(roadmap_hyper$qValueLog)] <- NA
max_logq <- max(roadmap_hyper$qValueLog, na.rm = TRUE)
roadmap_hyper$qValueLog[is.na(roadmap_hyper$qValueLog)] <- 1.5*max_logq

rm(match, max_logp, max_logq)

# Add percent DMRs
numDMRs <- roadmap_hyper$support[1] + roadmap_hyper$c[1]
roadmap_hyper$pct_DMRs <- roadmap_hyper$support*100/numDMRs

# Rearrange Columns
roadmap_hyper <- roadmap_hyper[,c("dbSet", "antibody", "cellType", "tissue", "type", "pValue", "qValue", "pValueLog", "qValueLog", "logOddsRatio", "support", 
                                  "pct_DMRs", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]
roadmap_hyper$order <- factor(roadmap_hyper$order, levels = sort(unique(roadmap_hyper$order), decreasing = TRUE), ordered = TRUE)
roadmap_hyper$antibody <- factor(roadmap_hyper$antibody, levels = c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"), ordered = TRUE)
roadmap_hyper <- roadmap_hyper[order(roadmap_hyper$order, roadmap_hyper$antibody),]
tissues <- rev(as.character(unique(roadmap_hyper$tissue)))
roadmap_hyper$tissue <- factor(roadmap_hyper$tissue, levels = tissues, ordered = TRUE)
tissue_colors <- rev(as.character(unique(roadmap_hyper$color)))

write.table(x = roadmap_hyper, file = "Tables/WD_Blood_2_Specific_DMRs_hyper_Roadmap_Enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Make Stats Table for Liver TFs ####
liver_TFs_hyper <- subset(stats, collection == "encode_liver_chipseq" & userSet == "hyperDMRs")

# Adjust liver_TFs_hyper for region sets with < 5 overlaps 
table(liver_TFs_hyper$support < 5) # 5 region sets 
liver_TFs_hyper$logOddsRatio[liver_TFs_hyper$support < 5] <- 0
liver_TFs_hyper$pValueLog[liver_TFs_hyper$support < 5] <- 0
liver_TFs_hyper$qValue[liver_TFs_hyper$support < 5] <- 1

# Replace inf logP with 1.5*max
table(is.na(liver_TFs_hyper$pValueLog)) # All false
table(is.infinite(liver_TFs_hyper$pValueLog)) # All false ##
liver_TFs_hyper$pValueLog[is.infinite(liver_TFs_hyper$pValueLog)] <- NA
max_logp <- max(liver_TFs_hyper$pValueLog, na.rm = TRUE)
liver_TFs_hyper$pValueLog[is.na(liver_TFs_hyper$pValueLog)] <- 1.5*max_logp

# Adjust pValue
liver_TFs_hyper$pValue <- 10^(-liver_TFs_hyper$pValueLog)
liver_TFs_hyper$qValueLog <- -log10(liver_TFs_hyper$qValue)

# Replace qValueLog with 1.5*max
table(is.na(liver_TFs_hyper$qValueLog)) # All false
table(is.infinite(liver_TFs_hyper$qValueLog)) # All false
liver_TFs_hyper$qValueLog[is.infinite(liver_TFs_hyper$qValueLog)] <- NA
max_logq <- max(liver_TFs_hyper$qValueLog, na.rm = TRUE)
liver_TFs_hyper$qValueLog[is.na(liver_TFs_hyper$qValueLog)] <- 1.5*max_logq

# Add percent DMRs
numDMRs <- liver_TFs_hyper$support[1] + liver_TFs_hyper$c[1]
liver_TFs_hyper$pct_DMRs <- liver_TFs_hyper$support*100/numDMRs
backDMRs <- sum(liver_TFs_hyper$support[1], liver_TFs_hyper$b[1], liver_TFs_hyper$c[1], liver_TFs_hyper$d[1]) 
liver_TFs_hyper$pct_BackDMRs <- (liver_TFs_hyper$support + liver_TFs_hyper$b)*100/backDMRs

# Reorder columns
liver_TFs_hyper <- liver_TFs_hyper[,c("dbSet", "antibody", "cellType", "pValue", "qValue", "pValueLog", "qValueLog", "logOddsRatio", "support", 
                                      "pct_DMRs", "pct_BackDMRs", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "size", "filename")]

# Fix Antibody Names
liver_TFs_hyper$antibody <- factor(liver_TFs_hyper$antibody, levels=sort(unique(liver_TFs_hyper$antibody)), ordered=TRUE)
write.table(x = liver_TFs_hyper, file = "Tables/WD_Blood_2_Specific_DMRs_hyper_LOLA_Liver_TF_Enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Chromhmm Figures ####################

# logOR Heatmap
gg <- ggplot(data = chromhmm_hyper)
gg +
        geom_tile(aes(x = chromState, y = order, fill = logOddsRatio)) +
        scale_fill_gradientn("-logOR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,20)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific Hyper DMR ChromHMM logOR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# logFDR Heatmap
gg <- ggplot(data = chromhmm_hyper)
gg +
        geom_tile(aes(x = chromState, y = order, fill = qValueLog)) +
        scale_fill_gradientn("-logFDR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,58)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific Hyper DMR ChromHMM logFDR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# Stacked Barchart of percent each chromstate
percent <- aggregate(x = chromhmm_hyper$pct_DMRs, by = list(chromhmm_hyper$chromState), FUN = "mean")
colnames(percent) <- c("chromState", "Mean_Percent")
chromColors <- c(rgb(255,0,0, maxColorValue = 255), rgb(255,69,0, maxColorValue = 255), rgb(50,205,50, maxColorValue = 255), 
                 rgb(0,128,0, maxColorValue = 255), rgb(0,100,0, maxColorValue = 255), rgb(194,225,5, maxColorValue = 255), 
                 rgb(255,255,0, maxColorValue = 255), rgb(102,205,170, maxColorValue = 255), rgb(138,145,208, maxColorValue = 255), 
                 rgb(205,92,92, maxColorValue = 255), rgb(233,150,122, maxColorValue = 255), rgb(189,183,107, maxColorValue = 255),
                 rgb(128,128,128, maxColorValue = 255), rgb(192,192,192, maxColorValue = 255), rgb(255,255,255, maxColorValue = 255))

backDMRs <- sum(chromhmm_hyper$support[1], chromhmm_hyper$b[1], chromhmm_hyper$c[1], chromhmm_hyper$d[1]) 
percentBackDMRs <- (chromhmm_hyper$support + chromhmm_hyper$b)*100/backDMRs
meanBackDMRs <- aggregate(x = percentBackDMRs, by = list(chromhmm_hyper$chromState), FUN = "mean")
colnames(meanBackDMRs) <- c("chromState", "Background")
meanBackDMRs$DMRs <- percent$Mean_Percent
meanBackDMRs$DMRs <- meanBackDMRs$DMRs*100/sum(meanBackDMRs$DMRs)
meanBackDMRs$Background <- meanBackDMRs$Background*100/sum(meanBackDMRs$Background)

percentDMRs_m <- melt(meanBackDMRs, id.vars="chromState")
percentDMRs_m$variable <- factor(percentDMRs_m$variable, levels=c("Background", "DMRs"), ordered=TRUE)

gg <- ggplot(data = percentDMRs_m)
gg +
        geom_col(aes(x = variable, y = value, fill = chromState), color = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.45, 0.64), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,12,1,0.5), "lines"), 
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 19, color = "black"),
              axis.ticks = element_line(size = 0.6, color = "black"),
              axis.text.x = element_text(angle=45, hjust=0.9),
              axis.title.x = element_blank(), 
              legend.title = element_text(size = 18)) +
        ylab("Mean % DMRs Across Tissues") +
        scale_fill_manual(name = "Chromatin State", values=chromColors) +
        scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(expand = c(0,0))
ggsave("Figures/LOLA WD Blood 2 Specific Hyper DMR and Background ChromHMM Mean Percent chromState.png", dpi = 600, width = 6, height = 7, units = "in")

# Stacked Barchart of percent each chromstate in liver

liver <- chromhmm_hyper[chromhmm_hyper$cellType == "Liver", c("chromState", "support", "b", "c", "d")]
liver$DMRs <- liver$support*100/numDMRs
liver$DMRs <- liver$DMRs*100/sum(liver$DMRs)
liver$Background <- (liver$support + liver$b)*100/backDMRs
liver$Background <- liver$Background*100/sum(liver$Background)

liver_m <- melt(liver, id.vars=c("chromState", "support", "b", "c", "d"))
liver_m$variable <- factor(liver_m$variable, levels=c("Background", "DMRs"), ordered=TRUE)

gg <- ggplot(data = liver_m)
gg +
        geom_col(aes(x = variable, y = value, fill = chromState), color = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.45, 0.64), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,12,1,0.5), "lines"), 
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 19, color = "black"),
              axis.ticks = element_line(size = 0.6, color = "black"),
              axis.text.x = element_text(angle=45, hjust=0.9),
              axis.title.x = element_blank(), 
              legend.title = element_text(size = 18)) +
        ylab("% DMRs in Liver") +
        scale_fill_manual(name = "Chromatin State", values=chromColors) +
        scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(expand = c(0,0))
ggsave("Figures/LOLA WD Blood 2 Specific Hyper DMR and Background ChromHMM Liver Percent chromState.png", dpi = 600, width = 6, height = 7, units = "in")

# Stacked Barchart of percent each chromstate in HSC & Bcell

HSCB <- chromhmm_hyper[chromhmm_hyper$tissue == "HSC & B-cell", c("chromState", "support", "b", "c", "d")]
HSCB$DMRs <- HSCB$support*100/numDMRs
HSCB$Background <- (HSCB$support + HSCB$b)*100/backDMRs

percent <- aggregate(x = HSCB$DMRs, by = list(HSCB$chromState), FUN = "mean")
colnames(percent) <- c("chromState", "Mean_Percent")

meanBackDMRs <- aggregate(x = HSCB$Background, by = list(HSCB$chromState), FUN = "mean")
colnames(meanBackDMRs) <- c("chromState", "Background")

meanBackDMRs$DMRs <- percent$Mean_Percent
meanBackDMRs$DMRs <- meanBackDMRs$DMRs*100/sum(meanBackDMRs$DMRs)
meanBackDMRs$Background <- meanBackDMRs$Background*100/sum(meanBackDMRs$Background)

percentDMRs_m <- melt(meanBackDMRs, id.vars="chromState")
percentDMRs_m$variable <- factor(percentDMRs_m$variable, levels=c("Background", "DMRs"), ordered=TRUE)

gg <- ggplot(data = percentDMRs_m)
gg +
        geom_col(aes(x = variable, y = value, fill = chromState), color = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.45, 0.64), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,12,1,0.5), "lines"), 
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 19, color = "black"),
              axis.ticks = element_line(size = 0.6, color = "black"),
              axis.text.x = element_text(angle=45, hjust=0.9),
              axis.title.x = element_blank(), 
              legend.title = element_text(size = 18)) +
        ylab("% DMRs in HSC & B-cells") +
        scale_fill_manual(name = "Chromatin State", values=chromColors) +
        scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(expand = c(0,0))
ggsave("Figures/LOLA WD Blood 2 Specific Hyper DMR and Background ChromHMM HSC B-cells Percent chromState.png", dpi = 600, width = 6, height = 7, units = "in")

# roadmap_hyper Figures ####################

# logOR Heatmap
gg <- ggplot(data = roadmap_hyper)
gg +
        geom_tile(aes(x = antibody, y = order, fill = logOddsRatio)) +
        scale_fill_gradientn("-logOR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,16)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR roadmap_hyper logOR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# logFDR Heatmap
gg <- ggplot(data = roadmap_hyper)
gg +
        geom_tile(aes(x = antibody, y = order, fill = qValueLog)) +
        scale_fill_gradientn("-logFDR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,75)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR roadmap_hyper logFDR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# Liver_TFs_hyper Figures ####
# Volcano Plot
top <- data.frame(antibody=liver_TFs_hyper$antibody, logOddsRatio=liver_TFs_hyper$logOddsRatio, qValueLog=liver_TFs_hyper$qValueLog)
top <- subset(top, rank(-top$qValueLog) <= 15)
top$antibody <- as.character(top$antibody)

gg <- ggplot(data = liver_TFs_hyper)
gg +
        geom_point(aes(x = logOddsRatio, y = qValueLog), color="#3366CC", size=2) +
        geom_text(data=top, aes(x=logOddsRatio, y=qValueLog, label=antibody), 
                  size=4, check_overlap = TRUE, hjust=0, nudge_x=0.05, nudge_y=0) +
        #annotate("text", label=c("RXRA", "YY1"), x=c(19.3,6.7), y=c(352,110), size=4) +
        theme_bw(base_size = 24) +
        labs(x="-log(Odds Ratio)", y="-log(q-value)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title = element_text(size=15, color="Black"), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) +
        scale_x_continuous(breaks=pretty_breaks(n=5)) +
        coord_cartesian(xlim=c(0,4.5))
ggsave("Figures/LOLA WD Blood 2 Specific Hyper DMR Liver TF Volcano Plot.png", dpi = 600, width = 7, height = 5, units = "in")

# Hypo DMRs ####
chromhmm_hypo <- subset(stats, collection == "Roadmap_ChromHMM" & userSet == "hypoDMRs")
roadmap_hypo <- subset(stats, collection == "roadmap_epigenomics" & userSet == "hypoDMRs")

# Make Stats Table for ChromHMM ####
# Merge Tables
match <- match(chromhmm_hypo$filename, index$filename)
chromhmm_hypo$type <- index$type[match]
chromhmm_hypo$order <- index$order[match]
chromhmm_hypo$color <- index$color[match]
chromhmm_hypo <- chromhmm_hypo[,c("dbSet", "antibody", "cellType", "tissue", "type", "pValueLog", "qValue", "logOddsRatio", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", 
                                  "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]
colnames(chromhmm_hypo)[2] <- "chromState"

# Adjust chromhmm_hypo for region sets with < 5 overlaps 
table(chromhmm_hypo$support < 5) # 997 region sets 
chromhmm_hypo$logOddsRatio[chromhmm_hypo$support < 5] <- 0
chromhmm_hypo$pValueLog[chromhmm_hypo$support < 5] <- 0
chromhmm_hypo$qValue[chromhmm_hypo$support < 5] <- 1

# Replace inf logP with 1.5*max
table(is.na(chromhmm_hypo$pValueLog)) # All false
chromhmm_hypo$pValueLog[is.infinite(chromhmm_hypo$pValueLog)] <- NA
max_logp <- max(chromhmm_hypo$pValueLog, na.rm = TRUE)
chromhmm_hypo$pValueLog[is.na(chromhmm_hypo$pValueLog)] <- 1.5*max_logp

# Adjust pValue
chromhmm_hypo$pValue <- 10^(-chromhmm_hypo$pValueLog)
chromhmm_hypo$qValueLog <- -log10(chromhmm_hypo$qValue)

# Replace qValueLog with 1.5*max
table(is.na(chromhmm_hypo$qValueLog)) # All false
chromhmm_hypo$qValueLog[is.infinite(chromhmm_hypo$qValueLog)] <- NA
max_logq <- max(chromhmm_hypo$qValueLog, na.rm = TRUE)
chromhmm_hypo$qValueLog[is.na(chromhmm_hypo$qValueLog)] <- 1.5*max_logq

rm(match, max_logp, max_logq)

# Add percent DMRs
numDMRs <- chromhmm_hypo$support[1] + chromhmm_hypo$c[1]
chromhmm_hypo$pct_DMRs <- chromhmm_hypo$support*100/numDMRs

# Rearrange Columns
chromhmm_hypo <- chromhmm_hypo[,c("dbSet", "chromState", "cellType", "tissue", "type", "pValue", "qValue", "pValueLog", "qValueLog", "logOddsRatio", "support", 
                                  "pct_DMRs", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]
chromhmm_hypo$order <- factor(chromhmm_hypo$order, levels = sort(unique(chromhmm_hypo$order), decreasing = TRUE), ordered = TRUE)
chromhmm_hypo$chromState <- sapply(strsplit(as.character(chromhmm_hypo$chromState), split='_', fixed=TRUE), function(x) (x[2]))
chromhmm_hypo$chromState <- factor(chromhmm_hypo$chromState, levels = c("TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZnfRpts", "Het", "TssBiv", "BivFlnk",
                                                                        "EnhBiv", "ReprPC", "ReprPCwk", "Quies"), ordered = TRUE)
chromhmm_hypo <- chromhmm_hypo[order(chromhmm_hypo$order, chromhmm_hypo$chromState),]
tissues <- rev(as.character(unique(chromhmm_hypo$tissue)))
chromhmm_hypo$tissue <- factor(chromhmm_hypo$tissue, levels = tissues, ordered = TRUE)
tissue_colors <- rev(as.character(unique(chromhmm_hypo$color)))

write.table(x = chromhmm_hypo, file = "Tables/WD_Blood_2_Specific_DMRs_Hypo_ChromHMM_Enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Make Stats table for Roadmap ####
# Merge Tables
EID <- strsplit(as.character(roadmap_hypo$filename), "-")
EID <- sapply(EID, function(x){x[1]})
roadmap_hypo$EID <- EID
match <- match(roadmap_hypo$EID, index$EID)
roadmap_hypo$type <- index$type[match]
roadmap_hypo$order <- index$order[match]
roadmap_hypo$color <- index$color[match]
roadmap_hypo$cellType <- index$cellType[match]
roadmap_hypo$tissue <- index$tissue[match]

roadmap_hypo$antibody <- as.character(roadmap_hypo$antibody)
roadmap_hypo <- subset(roadmap_hypo, antibody %in% c("H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"))
roadmap_hypo <- roadmap_hypo[,c("dbSet", "antibody", "cellType", "tissue", "type", "pValueLog", "qValue", "logOddsRatio", "support", "rnkPV", "rnkLO", "rnkSup", "maxRnk", 
                                "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]

# Adjust roadmap_hypo for region sets with < 5 overlaps 
table(roadmap_hypo$support < 5) # 151 region sets 
roadmap_hypo$logOddsRatio[roadmap_hypo$support < 5] <- 0
roadmap_hypo$pValueLog[roadmap_hypo$support < 5] <- 0
roadmap_hypo$qValue[roadmap_hypo$support < 5] <- 1

# Replace inf logP with 1.5*max
table(is.na(roadmap_hypo$pValueLog)) # All false
roadmap_hypo$pValueLog[is.infinite(roadmap_hypo$pValueLog)] <- NA
max_logp <- max(roadmap_hypo$pValueLog, na.rm = TRUE)
roadmap_hypo$pValueLog[is.na(roadmap_hypo$pValueLog)] <- 1.5*max_logp

# Adjust pValue
roadmap_hypo$pValue <- 10^(-roadmap_hypo$pValueLog)
roadmap_hypo$qValueLog <- -log10(roadmap_hypo$qValue)

# Replace qValueLog with 1.5*max
table(is.na(roadmap_hypo$qValueLog)) # All false
roadmap_hypo$qValueLog[is.infinite(roadmap_hypo$qValueLog)] <- NA
max_logq <- max(roadmap_hypo$qValueLog, na.rm = TRUE)
roadmap_hypo$qValueLog[is.na(roadmap_hypo$qValueLog)] <- 1.5*max_logq

rm(match, max_logp, max_logq)

# Add percent DMRs
numDMRs <- roadmap_hypo$support[1] + roadmap_hypo$c[1]
roadmap_hypo$pct_DMRs <- roadmap_hypo$support*100/numDMRs

# Rearrange Columns
roadmap_hypo <- roadmap_hypo[,c("dbSet", "antibody", "cellType", "tissue", "type", "pValue", "qValue", "pValueLog", "qValueLog", "logOddsRatio", "support", 
                                "pct_DMRs", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "size", "filename", "order", "color")]
roadmap_hypo$order <- factor(roadmap_hypo$order, levels = sort(unique(roadmap_hypo$order), decreasing = TRUE), ordered = TRUE)
roadmap_hypo$antibody <- factor(roadmap_hypo$antibody, levels = c("H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"), ordered = TRUE)
roadmap_hypo <- roadmap_hypo[order(roadmap_hypo$order, roadmap_hypo$antibody),]
tissues <- rev(as.character(unique(roadmap_hypo$tissue)))
roadmap_hypo$tissue <- factor(roadmap_hypo$tissue, levels = tissues, ordered = TRUE)
tissue_colors <- rev(as.character(unique(roadmap_hypo$color)))

write.table(x = roadmap_hypo, file = "Tables/WD_Blood_2_Specific_DMRs_hypo_Roadmap_Enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Make Stats Table for Liver TFs ####
liver_TFs_hypo <- subset(stats, collection == "encode_liver_chipseq" & userSet == "hypoDMRs")

# Adjust liver_TFs_hypo for region sets with < 5 overlaps 
table(liver_TFs_hypo$support < 5) # 24 region sets 
liver_TFs_hypo$logOddsRatio[liver_TFs_hypo$support < 5] <- 0
liver_TFs_hypo$pValueLog[liver_TFs_hypo$support < 5] <- 0
liver_TFs_hypo$qValue[liver_TFs_hypo$support < 5] <- 1

# Replace inf logP with 1.5*max
table(is.na(liver_TFs_hypo$pValueLog)) # All false
table(is.infinite(liver_TFs_hypo$pValueLog)) # All false
liver_TFs_hypo$pValueLog[is.infinite(liver_TFs_hypo$pValueLog)] <- NA
max_logp <- max(liver_TFs_hypo$pValueLog, na.rm = TRUE)
liver_TFs_hypo$pValueLog[is.na(liver_TFs_hypo$pValueLog)] <- 1.5*max_logp

# Adjust pValue
liver_TFs_hypo$pValue <- 10^(-liver_TFs_hypo$pValueLog)
liver_TFs_hypo$qValueLog <- -log10(liver_TFs_hypo$qValue)

# Replace qValueLog with 1.5*max
table(is.na(liver_TFs_hypo$qValueLog)) # All false
table(is.infinite(liver_TFs_hypo$qValueLog)) # All false
liver_TFs_hypo$qValueLog[is.infinite(liver_TFs_hypo$qValueLog)] <- NA
max_logq <- max(liver_TFs_hypo$qValueLog, na.rm = TRUE)
liver_TFs_hypo$qValueLog[is.na(liver_TFs_hypo$qValueLog)] <- 1.5*max_logq

# Add percent DMRs
numDMRs <- liver_TFs_hypo$support[1] + liver_TFs_hypo$c[1]
liver_TFs_hypo$pct_DMRs <- liver_TFs_hypo$support*100/numDMRs
backDMRs <- sum(liver_TFs_hypo$support[1], liver_TFs_hypo$b[1], liver_TFs_hypo$c[1], liver_TFs_hypo$d[1]) 
liver_TFs_hypo$pct_BackDMRs <- (liver_TFs_hypo$support + liver_TFs_hypo$b)*100/backDMRs

# Reorder columns
liver_TFs_hypo <- liver_TFs_hypo[,c("dbSet", "antibody", "cellType", "pValue", "qValue", "pValueLog", "qValueLog", "logOddsRatio", "support", 
                                    "pct_DMRs", "pct_BackDMRs", "rnkPV", "rnkLO", "rnkSup", "maxRnk", "meanRnk", "b", "c", "d", "size", "filename")]

# Fix Antibody Names
liver_TFs_hypo$antibody <- factor(liver_TFs_hypo$antibody, levels=sort(unique(liver_TFs_hypo$antibody)), ordered=TRUE)
write.table(x = liver_TFs_hypo, file = "Tables/WD_Blood_2_Specific_DMRs_hypo_LOLA_Liver_TF_Enrichment.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Chromhmm Figures ####################

# logOR Heatmap
gg <- ggplot(data = chromhmm_hypo)
gg +
        geom_tile(aes(x = chromState, y = order, fill = logOddsRatio)) +
        scale_fill_gradientn("-logOR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,20)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific Hypo DMR ChromHMM logOR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# logFDR Heatmap
gg <- ggplot(data = chromhmm_hypo)
gg +
        geom_tile(aes(x = chromState, y = order, fill = qValueLog)) +
        scale_fill_gradientn("-logFDR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,58)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific Hypo DMR ChromHMM logFDR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# Stacked Barchart of percent each chromstate
percent <- aggregate(x = chromhmm_hypo$pct_DMRs, by = list(chromhmm_hypo$chromState), FUN = "mean")
colnames(percent) <- c("chromState", "Mean_Percent")
chromColors <- c(rgb(255,0,0, maxColorValue = 255), rgb(255,69,0, maxColorValue = 255), rgb(50,205,50, maxColorValue = 255), 
                 rgb(0,128,0, maxColorValue = 255), rgb(0,100,0, maxColorValue = 255), rgb(194,225,5, maxColorValue = 255), 
                 rgb(255,255,0, maxColorValue = 255), rgb(102,205,170, maxColorValue = 255), rgb(138,145,208, maxColorValue = 255), 
                 rgb(205,92,92, maxColorValue = 255), rgb(233,150,122, maxColorValue = 255), rgb(189,183,107, maxColorValue = 255),
                 rgb(128,128,128, maxColorValue = 255), rgb(192,192,192, maxColorValue = 255), rgb(255,255,255, maxColorValue = 255))

backDMRs <- sum(chromhmm_hypo$support[1], chromhmm_hypo$b[1], chromhmm_hypo$c[1], chromhmm_hypo$d[1]) 
percentBackDMRs <- (chromhmm_hypo$support + chromhmm_hypo$b)*100/backDMRs
meanBackDMRs <- aggregate(x = percentBackDMRs, by = list(chromhmm_hypo$chromState), FUN = "mean")
colnames(meanBackDMRs) <- c("chromState", "Background")
meanBackDMRs$DMRs <- percent$Mean_Percent
meanBackDMRs$DMRs <- meanBackDMRs$DMRs*100/sum(meanBackDMRs$DMRs)
meanBackDMRs$Background <- meanBackDMRs$Background*100/sum(meanBackDMRs$Background)

percentDMRs_m <- melt(meanBackDMRs, id.vars="chromState")
percentDMRs_m$variable <- factor(percentDMRs_m$variable, levels=c("Background", "DMRs"), ordered=TRUE)

gg <- ggplot(data = percentDMRs_m)
gg +
        geom_col(aes(x = variable, y = value, fill = chromState), color = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.45, 0.64), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,12,1,0.5), "lines"), 
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 19, color = "black"),
              axis.ticks = element_line(size = 0.6, color = "black"),
              axis.text.x = element_text(angle=45, hjust=0.9),
              axis.title.x = element_blank(), 
              legend.title = element_text(size = 18)) +
        ylab("Mean % DMRs Across Tissues") +
        scale_fill_manual(name = "Chromatin State", values=chromColors) +
        scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(expand = c(0,0))
ggsave("Figures/LOLA WD Blood 2 Specific Hypo DMR and Background ChromHMM Mean Percent chromState.png", dpi = 600, width = 6, height = 7, units = "in")

# Stacked Barchart of percent each chromstate in liver

liver <- chromhmm_hypo[chromhmm_hypo$cellType == "Liver", c("chromState", "support", "b", "c", "d")]
liver$DMRs <- liver$support*100/numDMRs
liver$DMRs <- liver$DMRs*100/sum(liver$DMRs)
liver$Background <- (liver$support + liver$b)*100/backDMRs
liver$Background <- liver$Background*100/sum(liver$Background)

liver_m <- melt(liver, id.vars=c("chromState", "support", "b", "c", "d"))
liver_m$variable <- factor(liver_m$variable, levels=c("Background", "DMRs"), ordered=TRUE)

gg <- ggplot(data = liver_m)
gg +
        geom_col(aes(x = variable, y = value, fill = chromState), color = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.45, 0.64), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,12,1,0.5), "lines"), 
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 19, color = "black"),
              axis.ticks = element_line(size = 0.6, color = "black"),
              axis.text.x = element_text(angle=45, hjust=0.9),
              axis.title.x = element_blank(), 
              legend.title = element_text(size = 18)) +
        ylab("% DMRs in Liver") +
        scale_fill_manual(name = "Chromatin State", values=chromColors) +
        scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(expand = c(0,0))
ggsave("Figures/LOLA WD Blood 2 Specific Hypo DMR and Background ChromHMM Liver Percent chromState.png", dpi = 600, width = 6, height = 7, units = "in")

# Stacked Barchart of percent each chromstate in HSC & Bcell

HSCB <- chromhmm_hypo[chromhmm_hypo$tissue == "HSC & B-cell", c("chromState", "support", "b", "c", "d")]
HSCB$DMRs <- HSCB$support*100/numDMRs
HSCB$Background <- (HSCB$support + HSCB$b)*100/backDMRs

percent <- aggregate(x = HSCB$DMRs, by = list(HSCB$chromState), FUN = "mean")
colnames(percent) <- c("chromState", "Mean_Percent")

meanBackDMRs <- aggregate(x = HSCB$Background, by = list(HSCB$chromState), FUN = "mean")
colnames(meanBackDMRs) <- c("chromState", "Background")

meanBackDMRs$DMRs <- percent$Mean_Percent
meanBackDMRs$DMRs <- meanBackDMRs$DMRs*100/sum(meanBackDMRs$DMRs)
meanBackDMRs$Background <- meanBackDMRs$Background*100/sum(meanBackDMRs$Background)

percentDMRs_m <- melt(meanBackDMRs, id.vars="chromState")
percentDMRs_m$variable <- factor(percentDMRs_m$variable, levels=c("Background", "DMRs"), ordered=TRUE)

gg <- ggplot(data = percentDMRs_m)
gg +
        geom_col(aes(x = variable, y = value, fill = chromState), color = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_blank(),
              legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.45, 0.64), 
              legend.background = element_blank(), legend.text = element_text(size = 15, color = "Black"),
              plot.margin = unit(c(1,12,1,0.5), "lines"), 
              axis.text = element_text(size = 15, color = "black"),
              axis.title.y = element_text(size = 19, color = "black"),
              axis.ticks = element_line(size = 0.6, color = "black"),
              axis.text.x = element_text(angle=45, hjust=0.9),
              axis.title.x = element_blank(), 
              legend.title = element_text(size = 18)) +
        ylab("% DMRs in HSC & B-cells") +
        scale_fill_manual(name = "Chromatin State", values=chromColors) +
        scale_y_continuous(expand=c(0,0), breaks = c(0, 25, 50, 75, 100)) +
        scale_x_discrete(expand = c(0,0))
ggsave("Figures/LOLA WD Blood 2 Specific hypo DMR and Background ChromHMM HSC B-cells Percent chromState.png", dpi = 600, width = 6, height = 7, units = "in")

# roadmap_hypo Figures ####################

# logOR Heatmap
gg <- ggplot(data = roadmap_hypo)
gg +
        geom_tile(aes(x = antibody, y = order, fill = logOddsRatio)) +
        scale_fill_gradientn("-logOR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,16)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR roadmap_hypo logOR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# logFDR Heatmap
gg <- ggplot(data = roadmap_hypo)
gg +
        geom_tile(aes(x = antibody, y = order, fill = qValueLog)) +
        scale_fill_gradientn("-logFDR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,75)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR roadmap_hypo logFDR by Tissue.png", dpi = 600, width = 5, height = 7, units = "in")

# Liver_TFs_hypo Figures ####
# Volcano Plot
top <- data.frame(antibody=liver_TFs_hypo$antibody, logOddsRatio=liver_TFs_hypo$logOddsRatio, qValueLog=liver_TFs_hypo$qValueLog)
#top <- subset(top, qValueLog > -log(0.05))
top$antibody <- as.character(top$antibody)

gg <- ggplot(data = liver_TFs_hypo)
gg +
        geom_point(aes(x = logOddsRatio, y = qValueLog), color="#3366CC", size=2) +
        geom_text(data=top, aes(x=logOddsRatio, y=qValueLog, label=antibody), 
                  size=4, check_overlap = TRUE, hjust=0, nudge_x=0.1, nudge_y=0) +
        #annotate("text", label=c("CTCF"), x=c(1.9), y=c(7.51), size=4) +
        theme_bw(base_size = 24) +
        labs(x="-log(Odds Ratio)", y="-log(q-value)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title = element_text(size=15, color="Black"), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        coord_cartesian(xlim=c(0,2))
ggsave("Figures/LOLA WD Blood 2 Specific Hypo DMR Liver TF Volcano Plot.png", dpi = 600, width = 7, height = 5, units = "in")

# Max Values for Scales
# logOR
max(roadmap_all$logOddsRatio, roadmap_hyper$logOddsRatio, roadmap_hypo$logOddsRatio) #16
max(chromhmm_all$logOddsRatio, chromhmm_hyper$logOddsRatio, chromhmm_hypo$logOddsRatio) #20

# FDR
max(roadmap_all$qValueLog, roadmap_hyper$qValueLog, roadmap_hypo$qValueLog) #75
max(chromhmm_all$qValueLog, chromhmm_hyper$qValueLog, chromhmm_hypo$qValueLog) #58

# liver TFs Hyper vs Hypo ####
# logOddsRatio Volcano Plot
liver_TFs_hyperhypo_OR <- data.frame(antibody=liver_TFs_hyper$antibody, hyper_logOddsRatio=liver_TFs_hyper$logOddsRatio,
                                     hypo_logOddsRatio=liver_TFs_hypo$logOddsRatio[match(liver_TFs_hyper$dbSet, liver_TFs_hypo$dbSet)])
top <- subset(liver_TFs_hyperhypo_OR, rank(-liver_TFs_hyperhypo_OR$hyper_logOddsRatio) <= 10)
top$antibody <- as.character(top$antibody)

gg <- ggplot(data = liver_TFs_hyperhypo_OR)
gg +
        geom_point(aes(x = hyper_logOddsRatio, y = hypo_logOddsRatio), color="#3366CC", size=2) +
        geom_text(data=top, aes(x=hyper_logOddsRatio, y=hypo_logOddsRatio, label=antibody), 
                  size=3.5, check_overlap = TRUE, hjust=0, nudge_x=0.1, nudge_y=0) +
        #annotate("text", label=c("HNF4G", "FOXA2", "FOXA1", "JUND"), x=c(20.9, 14.4, 10.9, 11), y=c(0.2, 0.8, 0, 2.6), size=3.5, hjust=0) +
        theme_bw(base_size = 24) +
        labs(x="Hyper DMR -log(Odds Ratio)", y="Hypo DMR -log(Odds Ratio)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title = element_text(size=15, color="Black"), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) +
        scale_x_continuous(breaks=pretty_breaks(n=5)) +
        scale_y_continuous(breaks=pretty_breaks(n=5)) +
        coord_cartesian(xlim=c(0,5), ylim=c(0,2))
ggsave("Figures/LOLA WD Blood 2 Specific Hyper vs Hypo DMR liver TF logOR Volcano Plot.png", dpi = 600, width = 7, height = 5, units = "in")

# qValueLog Volcano Plot

liver_TFs_hyperhypo_q <- data.frame(antibody=liver_TFs_hyper$antibody, hyper_qValueLog=liver_TFs_hyper$qValueLog,
                                    hypo_qValueLog=liver_TFs_hypo$qValueLog[match(liver_TFs_hyper$dbSet, liver_TFs_hypo$dbSet)])
top <- subset(liver_TFs_hyperhypo_q, rank(-liver_TFs_hyperhypo_q$hyper_qValueLog) <= 10)
top$antibody <- as.character(top$antibody)

gg <- ggplot(data = liver_TFs_hyperhypo_q)
gg +
        geom_point(aes(x = hyper_qValueLog, y = hypo_qValueLog), color="#3366CC", size=2) +
        geom_text(data=top, aes(x=hyper_qValueLog, y=hypo_qValueLog, label=antibody), 
                  size=3.5, check_overlap = TRUE, hjust=0, nudge_x=0.1, nudge_y=0) +
        #annotate("text", label=c("RXRA", "HNF4G", "REST", "FOXA2", "YY1"), x=c(352, 188, 172, 130, 107), y=c(30, 3, 19, 18, 13), size=3.5, hjust=0) +
        theme_bw(base_size = 24) +
        labs(x="Hyper DMR -log(q-value)", y="Hypo DMR -log(q-value)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title = element_text(size=15, color="Black"), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) +
        scale_x_continuous(breaks=pretty_breaks(n=5)) +
        scale_y_continuous(breaks=pretty_breaks(n=5)) +
        coord_cartesian(xlim=c(0,10), ylim=c(0,0.5))
ggsave("Figures/LOLA WD Blood 2 Specific Hyper vs Hypo DMR liver TF logFDR Volcano Plot.png", dpi = 600, width = 7, height = 5, units = "in")

# logOddsRatio Heatmap
liver_TFs_heat_OR <- data.frame(antibody=liver_TFs_all$antibody, All=liver_TFs_all$logOddsRatio, 
                                Hyper=liver_TFs_hyper$logOddsRatio[match(liver_TFs_all$dbSet, liver_TFs_hyper$dbSet)], 
                                Hypo=liver_TFs_hypo$logOddsRatio[match(liver_TFs_all$dbSet, liver_TFs_hypo$dbSet)])
liver_TFs_heat_OR <- melt(liver_TFs_heat_OR, id.vars="antibody")
liver_TFs_heat_OR <- aggregate(liver_TFs_heat_OR$value, by=list(antibody=liver_TFs_heat_OR$antibody, variable=liver_TFs_heat_OR$variable), 
                               FUN=mean) # Aggregate duplicates by averaging them
TF_order <- subset(liver_TFs_heat_OR, variable=="Hyper")
TF_order <- TF_order[order(TF_order$x),] # Order TFs by enrichment in hypermethylated DMRs
TF_order <- as.character(TF_order$antibody)
liver_TFs_heat_OR$antibody <- factor(liver_TFs_heat_OR$antibody, levels=TF_order, ordered=TRUE)

gg <- ggplot(data = liver_TFs_heat_OR)
gg +
        geom_tile(aes(x = variable, y = antibody, fill = x)) +
        scale_fill_gradientn("-logOR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_text(size = 15, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR Liver TF OR Heatmap.png", dpi = 600, width = 5, height = 7, units = "in")

# FDR Heatmap
liver_TFs_heat_FDR <- data.frame(antibody=liver_TFs_all$antibody, All=liver_TFs_all$qValueLog, 
                                 Hyper=liver_TFs_hyper$qValueLog[match(liver_TFs_all$dbSet, liver_TFs_hyper$dbSet)], 
                                 Hypo=liver_TFs_hypo$qValueLog[match(liver_TFs_all$dbSet, liver_TFs_hypo$dbSet)])
liver_TFs_heat_FDR <- melt(liver_TFs_heat_FDR, id.vars="antibody")
liver_TFs_heat_FDR <- aggregate(liver_TFs_heat_FDR$value, by=list(antibody=liver_TFs_heat_FDR$antibody, variable=liver_TFs_heat_FDR$variable), 
                                FUN=mean) # Aggregate duplicates by averaging them
TF_order <- subset(liver_TFs_heat_FDR, variable=="Hyper")
TF_order <- TF_order[order(TF_order$x),] # Order TFs by enrichment in hypermethylated DMRs
TF_order <- as.character(TF_order$antibody)
liver_TFs_heat_FDR$antibody <- factor(liver_TFs_heat_FDR$antibody, levels=TF_order, ordered=TRUE)

gg <- ggplot(data = liver_TFs_heat_FDR)
gg +
        geom_tile(aes(x = variable, y = antibody, fill = x)) +
        scale_fill_gradientn("-logFDR\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.25, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,7,1,1), "lines"), 
              axis.text.y = element_text(size = 15, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18))
ggsave("Figures/LOLA WD Blood 2 Specific DMR Liver TF FDR Heatmap.png", dpi = 600, width = 5, height = 7, units = "in")
