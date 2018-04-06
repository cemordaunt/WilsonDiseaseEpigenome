# WD Blood and Liver Gene-Drug Enrichment ####
# Charles Mordaunt
# 2/24/18

# Packages ####
library(GeneOverlap)
library(reshape)
library(ggplot2)
library(scales)

# Data ####
# Gene Drug Interactions
liver_drug <- read.delim("Tables/Liver Background Gene Drug Interactions.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)
blood_drug <- read.delim("Tables/Blood Background Gene Drug Interactions.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)
blood_liver_drug <- read.delim("Tables/Blood and Liver Background Gene Drug Interactions.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE)

# Blood DMR Genes
blood_hyper <- as.character(unlist(read.delim("Tables/WD Blood Hyper DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
blood_hypo <- as.character(unlist(read.delim("Tables/WD Blood Hypo DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
blood_back <- as.character(unlist(read.delim("Tables/GREAT/WD Blood Background GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
blood_all <- unique(c(blood_hyper, blood_hypo))

# Liver DMR Genes
liver_hyper <- as.character(unlist(read.delim("Tables/WD Liver Hyper DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
liver_hypo <- as.character(unlist(read.delim("Tables/WD Liver Hypo DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
liver_back <- as.character(unlist(read.delim("Tables/GREAT/WD Liver Distal Background Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
liver_all <- unique(c(liver_hyper, liver_hypo))

# Blood and Liver DMR Genes
blood_liver_hyper <- intersect(blood_hyper, liver_hyper)
blood_liver_hypo <- intersect(blood_hypo, liver_hypo)
blood_liver_back <- intersect(blood_back, liver_back)
blood_liver_all <- unique(c(blood_liver_hyper, blood_liver_hypo))
#write.table(sort(unique(blood_liver_back)), "Tables/GREAT/WD Blood and Liver Distal Background Genes.txt", quote=FALSE, sep="\t", 
#            row.names=FALSE, col.names=FALSE)

# Make Gene Lists ####
# DMR Genes
blood_genes <- list("Blood_All"=blood_all, "Blood_Hyper"=blood_hyper, "Blood_Hypo"=blood_hypo)
liver_genes <- list("Liver_All"=liver_all, "Liver_Hyper"=liver_hyper, "Liver_Hypo"=liver_hypo)
blood_liver_genes <- list("Blood_Liver_All"=blood_liver_all, "Blood_Liver_Hyper"=blood_liver_hyper, "Blood_Liver_Hypo"=blood_liver_hypo)

# Blood Drug Genes
blood_drug_names <- unique(unlist(blood_drug$drug))
blood_drug_genes <- lapply(blood_drug_names, function(x){unique(unlist(blood_drug$search_term[blood_drug$drug==x]))})
names(blood_drug_genes) <- blood_drug_names

# Liver Drug Genes
liver_drug_names <- unique(unlist(liver_drug$drug))
liver_drug_genes <- lapply(liver_drug_names, function(x){unique(unlist(liver_drug$search_term[liver_drug$drug==x]))})
names(liver_drug_genes) <- liver_drug_names

# Blood and Liver Drug Genes
blood_liver_drug_names <- unique(unlist(blood_liver_drug$drug))
blood_liver_drug_genes <- lapply(blood_liver_drug_names, function(x){unique(unlist(blood_liver_drug$search_term[blood_liver_drug$drug==x]))})
names(blood_liver_drug_genes) <- blood_liver_drug_names

# Blood Gene Overlap Tests ####
blood_gom <- newGOM(gsetA=blood_drug_genes, gsetB=blood_genes, genome.size = length(blood_back))
blood_int <- getMatrix(blood_gom, "intersection")
blood_OR <- getMatrix(blood_gom, "odds.ratio")
blood_p <- getMatrix(blood_gom, "pval") 
blood_jac <- getMatrix(blood_gom, "Jaccard")
blood_sum <- as.data.frame(cbind(blood_int, blood_OR, blood_p, blood_jac))
blood_cols <- paste(colnames(blood_sum), rep(c("Intersection", "oddsRatio", "pValue", "Jaccard"), each=3), sep="_")
colnames(blood_sum) <- blood_cols
blood_sum$Blood_All_qValue <- p.adjust(blood_sum$Blood_All_pValue, "fdr")
blood_sum$Blood_Hyper_qValue <- p.adjust(blood_sum$Blood_Hyper_pValue, "fdr")
blood_sum$Blood_Hypo_qValue <- p.adjust(blood_sum$Blood_Hypo_pValue, "fdr")
# No Enrichment with q < 1

blood_sum_top <- subset(blood_sum, Blood_All_pValue < 0.005 | Blood_Hyper_pValue < 0.005 | Blood_Hypo_pValue < 0.005)
blood_sum_top <- subset(blood_sum_top, Blood_All_Intersection > 1 | Blood_Hyper_Intersection > 1 | Blood_Hypo_Intersection > 1)
blood_logp <- data.frame("Drug"=rownames(blood_sum_top), "All"=-log10(blood_sum_top$Blood_All_pValue), 
                         "Hyper"=-log10(blood_sum_top$Blood_Hyper_pValue), "Hypo"=-log10(blood_sum_top$Blood_Hypo_pValue), 
                         stringsAsFactors = FALSE)
blood_logp$Drug[blood_logp$Drug == "CHEMBL1213492"] <- "GIVINOSTAT"
blood_logp$Drug[blood_logp$Drug == "CHEMBL1230989"] <- "GAMMA-IMINO-ATP"
blood_logp$Drug[blood_logp$Drug == "CHEMBL406845"] <- "PP2"
blood_logp$Drug[blood_logp$Drug == "CHEMBL100014"] <- "PIVANEX"
blood_logp$Drug[blood_logp$Drug == "CHEMBL3110004"] <- "TMP269"
blood_logp$Drug <- factor(blood_logp$Drug, levels=blood_logp$Drug[order(blood_logp$All)], ordered=TRUE)
blood_logp_m <- melt(blood_logp, id.vars="Drug")
colnames(blood_logp_m) <- c("Drug", "DMRs", "logpval")

# logpvalue Heatmap
gg <- ggplot(data = blood_logp_m)
gg +
        geom_tile(aes(x = DMRs, y = Drug, fill = logpval)) +
        scale_fill_gradientn("-log(p-value)\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000",
                             limits=c(0,3.6)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.4, 0.7), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 12, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_blank(), 
              legend.title = element_text(size = 17), legend.text = element_text(size = 17),
              plot.title = element_text(size = 18))
ggsave("Figures/Blood DMR Gene Drug Enrichment logp Heatmap.png", dpi = 600, width = 6.8, height = 3.8, units = "in")

# Liver Gene Overlap Tests ####
liver_gom <- newGOM(gsetA=liver_drug_genes, gsetB=liver_genes, genome.size = length(liver_back))
liver_int <- getMatrix(liver_gom, "intersection")
liver_OR <- getMatrix(liver_gom, "odds.ratio")
liver_p <- getMatrix(liver_gom, "pval") 
liver_jac <- getMatrix(liver_gom, "Jaccard")
liver_sum <- as.data.frame(cbind(liver_int, liver_OR, liver_p, liver_jac))
liver_cols <- paste(colnames(liver_sum), rep(c("Intersection", "oddsRatio", "pValue", "Jaccard"), each=3), sep="_")
colnames(liver_sum) <- liver_cols
liver_sum$Liver_All_qValue <- p.adjust(liver_sum$Liver_All_pValue, "fdr")
liver_sum$Liver_Hyper_qValue <- p.adjust(liver_sum$Liver_Hyper_pValue, "fdr")
liver_sum$Liver_Hypo_qValue <- p.adjust(liver_sum$Liver_Hypo_pValue, "fdr")
# No Enrichment with q < 1

liver_sum_top <- subset(liver_sum, Liver_All_pValue < 0.005 | Liver_Hyper_pValue < 0.005 | Liver_Hypo_pValue < 0.005)
liver_sum_top <- subset(liver_sum_top, Liver_All_Intersection > 1 | Liver_Hyper_Intersection > 1 | Liver_Hypo_Intersection > 1)
liver_logp <- data.frame("Drug"=rownames(liver_sum_top), "All"=-log10(liver_sum_top$Liver_All_pValue), 
                         "Hyper"=-log10(liver_sum_top$Liver_Hyper_pValue), "Hypo"=-log10(liver_sum_top$Liver_Hypo_pValue), 
                         stringsAsFactors = FALSE)
liver_logp$Drug[liver_logp$Drug == "CHEMBL3277946"] <- "ADRIAMYCINOL"
liver_logp$Drug[liver_logp$Drug == "CHEMBL1851943"] <- "PRACINOSTAT"
liver_logp$Drug[liver_logp$Drug == "CHEMBL3110004"] <- "TMP269"

liver_logp$Drug <- factor(liver_logp$Drug, levels=liver_logp$Drug[order(liver_logp$All)], ordered=TRUE)
liver_logp_m <- melt(liver_logp, id.vars="Drug")
colnames(liver_logp_m) <- c("Drug", "DMRs", "logpval")

# logpvalue Heatmap
gg <- ggplot(data = liver_logp_m)
gg +
        geom_tile(aes(x = DMRs, y = Drug, fill = logpval)) +
        scale_fill_gradientn("-log(p-value)\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000",
                             limits=c(0,3.6)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.4, 0.87), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 12, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_blank(), 
              legend.title = element_text(size = 17), legend.text = element_text(size = 17),
              plot.title = element_text(size = 18))
ggsave("Figures/Liver DMR Gene Drug Enrichment logp Heatmap.png", dpi = 600, width = 6.5, height = 7.5, units = "in")

# Blood_Liver Gene Overlap Tests ####
blood_liver_gom <- newGOM(gsetA=blood_liver_drug_genes, gsetB=blood_liver_genes, genome.size = length(blood_liver_back))
blood_liver_int <- getMatrix(blood_liver_gom, "intersection")
blood_liver_OR <- getMatrix(blood_liver_gom, "odds.ratio")
blood_liver_p <- getMatrix(blood_liver_gom, "pval") 
blood_liver_jac <- getMatrix(blood_liver_gom, "Jaccard")
blood_liver_sum <- as.data.frame(cbind(blood_liver_int, blood_liver_OR, blood_liver_p, blood_liver_jac))
blood_liver_cols <- paste(colnames(blood_liver_sum), rep(c("Intersection", "oddsRatio", "pValue", "Jaccard"), each=3), sep="_")
colnames(blood_liver_sum) <- blood_liver_cols
blood_liver_sum$Blood_Liver_All_qValue <- p.adjust(blood_liver_sum$Blood_Liver_All_pValue, "fdr")
blood_liver_sum$Blood_Liver_Hyper_qValue <- p.adjust(blood_liver_sum$Blood_Liver_Hyper_pValue, "fdr")
blood_liver_sum$Blood_Liver_Hypo_qValue <- p.adjust(blood_liver_sum$Blood_Liver_Hypo_pValue, "fdr")
# No Enrichment with q < 1

blood_liver_sum_top <- subset(blood_liver_sum, Blood_Liver_All_pValue < 0.005 | Blood_Liver_Hyper_pValue < 0.005 | Blood_Liver_Hypo_pValue < 0.01)
blood_liver_sum_top <- subset(blood_liver_sum_top, Blood_Liver_All_Intersection > 1 | Blood_Liver_Hyper_Intersection > 1 | Blood_Liver_Hypo_Intersection > 1)
blood_liver_logp <- data.frame("Drug"=rownames(blood_liver_sum_top), "All"=-log10(blood_liver_sum_top$Blood_Liver_All_pValue), 
                         "Hyper"=-log10(blood_liver_sum_top$Blood_Liver_Hyper_pValue), "Hypo"=-log10(blood_liver_sum_top$Blood_Liver_Hypo_pValue), 
                         stringsAsFactors = FALSE)
blood_liver_logp$Drug[blood_liver_logp$Drug == "CHEMBL3110004"] <- "TMP269"
blood_liver_logp$Drug[blood_liver_logp$Drug == "CHEMBL100014"] <- "PIVANEX"
blood_liver_logp$Drug[blood_liver_logp$Drug == "CHEMBL1851943"] <- "PRACINOSTAT"
blood_liver_logp$Drug[blood_liver_logp$Drug == "CHEMBL191482"] <- "AR-42"

blood_liver_logp$Drug <- factor(blood_liver_logp$Drug, levels=blood_liver_logp$Drug[order(blood_liver_logp$All)], ordered=TRUE)
blood_liver_logp_m <- melt(blood_liver_logp, id.vars="Drug")
colnames(blood_liver_logp_m) <- c("Drug", "DMRs", "logpval")

# logpvalue Heatmap
gg <- ggplot(data = blood_liver_logp_m)
gg +
        geom_tile(aes(x = DMRs, y = Drug, fill = logpval)) +
        scale_fill_gradientn("-log(p-value)\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000",
                             limits=c(0,3.6)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.5, 0.77), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 15, color = "Black", angle = 90, hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 12, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_blank(), 
              legend.title = element_text(size = 17), legend.text = element_text(size = 17),
              plot.title = element_text(size = 18))
ggsave("Figures/Blood_Liver DMR Gene Drug Enrichment logp Heatmap.png", dpi = 600, width = 6, height = 5, units = "in")

# HDAC Drugs ####
# Liver
liver_drug_top <- as.character(rownames(liver_sum_top))
liver_drug_genes_top <- subset(liver_drug_genes, names(liver_drug_genes) %in% liver_drug_top)
table(names(liver_drug_genes_top) == liver_drug_top) #All True
liver_sum_top$HDAC5 <- sapply(liver_drug_genes_top, function(x){"HDAC5" %in% x})
liver_sum_top$AnyHDAC <- sapply(liver_drug_genes_top, function(x){TRUE %in% grepl("HDAC", x=x)})
liver_sum_top_HDAC5 <- rownames(liver_sum_top)[liver_sum_top$HDAC5]
liver_drug_hdac5 <- subset(liver_drug, drug %in% liver_sum_top_HDAC5 & search_term == "HDAC5")

# Blood
blood_drug_top <- as.character(rownames(blood_sum_top))
blood_drug_genes_top <- subset(blood_drug_genes, names(blood_drug_genes) %in% blood_drug_top)
table(names(blood_drug_genes_top) == blood_drug_top) #All True
blood_sum_top$HDAC5 <- sapply(blood_drug_genes_top, function(x){"HDAC5" %in% x})
blood_sum_top$AnyHDAC <- sapply(blood_drug_genes_top, function(x){TRUE %in% grepl("HDAC", x=x)})
blood_sum_top_HDAC5 <- rownames(blood_sum_top)[blood_sum_top$HDAC5]
blood_drug_hdac5 <- subset(blood_drug, drug %in% blood_sum_top_HDAC5 & search_term == "HDAC5")

# Blood and Liver
blood_liver_drug_top <- as.character(rownames(blood_liver_sum_top))
blood_liver_drug_genes_top <- subset(blood_liver_drug_genes, names(blood_liver_drug_genes) %in% blood_liver_drug_top)
table(names(blood_liver_drug_genes_top) == blood_liver_drug_top) #All True
blood_liver_sum_top$HDAC5 <- sapply(blood_liver_drug_genes_top, function(x){"HDAC5" %in% x})
blood_liver_sum_top$AnyHDAC <- sapply(blood_liver_drug_genes_top, function(x){TRUE %in% grepl("HDAC", x=x)})
blood_liver_sum_top_HDAC5 <- rownames(blood_liver_sum_top)[blood_liver_sum_top$HDAC5]
blood_liver_drug_hdac5 <- subset(blood_liver_drug, drug %in% blood_liver_sum_top_HDAC5 & search_term == "HDAC5")
