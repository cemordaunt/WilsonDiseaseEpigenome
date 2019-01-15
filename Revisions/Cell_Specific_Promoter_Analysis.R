# WD DMR Cell Type-Specific Promoter Analysis ####
# Data from Taghdouini et al 2015
# 12/26/18

# Packages ####
library(dplyr)
library(biomaRt)
library(reshape2)
library(ggplot2)
library(scales)
library(GenomicRanges)

# Functions ####
methByDiagnosis <- function(data, samples){
        # takes input region percent methylation data and sample info and runs ANOVA by diagnosis followed by Tukey HSD test
        meth <- data[,samples$SequencingID]
        if(length(levels(samples$Diagnosis)) == 2){
                group1 <- levels(samples$Diagnosis)[1]
                group2 <- levels(samples$Diagnosis)[2]
                stats <- apply(meth, 1, function(x){
                        tukeyStats <- aov(x ~ samples$Diagnosis) %>% TukeyHSD %>% .[[1]]
                        temp <- c(mean(x[samples$Diagnosis == group1]), mean(x[samples$Diagnosis == group2]), 
                                  sd(x[samples$Diagnosis == group1]), sd(x[samples$Diagnosis == group2]), 
                                  tukeyStats[,"diff"], tukeyStats[,"p adj"])
                        return(temp)
                })
                stats <- stats %>% t %>% as.data.frame
                colnames(stats) <- c(paste(rep(c("mean", "sd"), each = 2), c(group1, group2), sep = "_"),
                                     paste(rep(c("diff", "pvalue"), each = 1), c(paste(group2, "vs", group1, sep = "")), sep = "_"))
                sig <- stats[,paste("pvalue_", group2, "vs", group1, sep = "")] < 0.05
                stats <- cbind(data[,c("entrezID", "chr", "start", "end", "strand", "geneName", "cellType")], stats, sig)
                colnames(stats)[length(colnames(stats))] <- paste("Sig", c(paste(group2, "vs", group1, sep = "")), sep = "_")
                return(stats)
        } else {
                if(length(levels(samples$Diagnosis)) == 3){
                        group1 <- levels(samples$Diagnosis)[1]
                        group2 <- levels(samples$Diagnosis)[2]
                        group3 <- levels(samples$Diagnosis)[3]
                        stats <- apply(meth, 1, function(x){
                                tukeyStats <- aov(x ~ samples$Diagnosis) %>% TukeyHSD %>% .[[1]]
                                temp <- c(mean(x[samples$Diagnosis == group1]), mean(x[samples$Diagnosis == group2]), mean(x[samples$Diagnosis == group3]),
                                        sd(x[samples$Diagnosis == group1]), sd(x[samples$Diagnosis == group2]), sd(x[samples$Diagnosis == group3]),
                                        tukeyStats[,"diff"], tukeyStats[,"p adj"])
                                return(temp)
                        })
                        stats <- stats %>% t %>% as.data.frame
                        colnames(stats) <- c(paste(rep(c("mean", "sd"), each = 3), c(group1, group2, group3), sep = "_"),
                                        paste(rep(c("diff", "pvalue"), each = 3), 
                                                c(paste(group2, "vs", group1, sep = ""), paste(group3, "vs", group1, sep = ""),
                                                paste(group3, "vs", group2, sep = "")), sep = "_"))
                        sig <- cbind(stats[,paste("pvalue_", group2, "vs", group1, sep = "")] < 0.05,
                                stats[,paste("pvalue_", group3, "vs", group1, sep = "")] < 0.05,
                                stats[,paste("pvalue_", group3, "vs", group2, sep = "")] < 0.05)
                        stats <- cbind(data[,c("entrezID", "chr", "start", "end", "strand", "geneName", "cellType")], stats, sig)
                        colnames(stats)[(length(colnames(stats))-2):length(colnames(stats))] <- paste("Sig", c(paste(group2, "vs", group1, sep = ""), 
                                                                                                        paste(group3, "vs", group1, sep = ""),
                                                                                                        paste(group3, "vs", group2, sep = "")), sep = "_")
                        return(stats)
                }
        }
}

# Data ####
hep <- read.delim(file = "Tables/HEP specific entrez IDs.txt", sep = "\n", header = FALSE, stringsAsFactors = FALSE) %>% 
        unlist %>% as.integer %>% unique %>% sort
hsc <- read.delim(file = "Tables/HSC specific entrez IDs.txt", sep = "\n", header = FALSE, stringsAsFactors = FALSE) %>% 
        unlist %>% as.integer %>% unique %>% sort
lsec <- read.delim(file = "Tables/LSEC specific entrez IDs.txt", sep = "\n", header = FALSE, stringsAsFactors = FALSE) %>% 
        unlist %>% as.integer %>% unique %>% sort
entrezIDs <- list("hep" = hep, "hsc" = hsc, "lsec" = lsec)

# Get hg38 gene annotations ####
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
ensemblAtt <- listAttributes(ensembl)
geneInfo <- lapply(entrezIDs, function(x){
        getBM(attributes = c("entrezgene", "external_gene_name", "chromosome_name", "start_position", 
                             "end_position", "strand"),
              filters = "entrezgene", values = x, mart = ensembl)
})
geneInfo <- lapply(geneInfo, function(x){subset(x, chromosome_name %in% c(1:22, "X", "Y"))}) # Remove misc chroms

# Add promoter locations ####
for(i in 1:length(geneInfo)){
        geneInfo[[i]]$tss <- ifelse(geneInfo[[i]]$strand == 1, geneInfo[[i]]$start_position, geneInfo[[i]]$end_position)
        geneInfo[[i]]$promoter_start <- ifelse(geneInfo[[i]]$strand == 1, geneInfo[[i]]$tss - 3000, geneInfo[[i]]$tss - 1000)
        geneInfo[[i]]$promoter_end <- ifelse(geneInfo[[i]]$strand == 1, geneInfo[[i]]$tss + 1000, geneInfo[[i]]$tss + 3000)
}

# Make bed files ####
# HEP
hepBed <- geneInfo$hep[,c("chromosome_name", "promoter_start", "promoter_end", "entrezgene")]
hepBed$chromosome_name <- paste("chr", hepBed$chromosome_name, sep = "")
hepBed <- hepBed[order(hepBed$chromosome_name, hepBed$promoter_start),]

# HSC
hscBed <- geneInfo$hsc[,c("chromosome_name", "promoter_start", "promoter_end", "entrezgene")]
hscBed$chromosome_name <- paste("chr", hscBed$chromosome_name, sep = "")
hscBed <- hscBed[order(hscBed$chromosome_name, hscBed$promoter_start),]

# LSEC
lsecBed <- geneInfo$lsec[,c("chromosome_name", "promoter_start", "promoter_end", "entrezgene")]
lsecBed$chromosome_name <- paste("chr", lsecBed$chromosome_name, sep = "")
lsecBed <- lsecBed[order(lsecBed$chromosome_name, lsecBed$promoter_start),]

# Combined
combined <- rbind(hepBed, hscBed, lsecBed)
combined <- combined[order(combined$chromosome_name, combined$promoter_start),]

# Error message from AvgMeth.2col.pl ####
# Warning: duplicate found, skipping: chr16       15392754        15396754        255027
# Warning: duplicate found, skipping: chr7        116950238       116954238       7982

# Fix Duplicates
combined <- combined %>% unique

# 2 entrezIDs annotated to the same gene and location
# entrezID 7982 is correctly annotated
# entrezID 93655 is not correctly annotated, remove this one
combined <- combined[!combined$entrezgene == 93655,]
write.table(combined, file = "Tables/Combined_specific_promoters.bed", sep = "\t", quote = FALSE, row.names = FALSE)
rm(ensembl, ensemblAtt, hepBed, hscBed, lsecBed, hep, hsc, i, lsec)
# AvgMeth.2col.pl now runs with no error message

# Methylation Analysis ####
# Data ####
samples <- read.csv(file = "Tables/Supplemental_Table_S1_for_stats.csv", header = TRUE, stringsAsFactors = FALSE)
countData <- read.delim(file = "Tables/specific_promoter_methylation.txt.2col", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
meth <- countData[,grepl("meth", colnames(countData))]
cov <- countData[,grepl("total", colnames(countData))]
table(is.na(cov)) # 760
nacount <- rowSums(is.na(cov))
table(nacount) # Out of 2530 regions in 21 samples
#    0    1    2    3    4    5    6    8    9   10   12   13   14   15   17   18   19   20 
# 2465    3    4    6    5    1    5    2    3    1    2    1    1    3    3    8    8    9
permeth <- meth * 100 / cov
permeth <- cbind(countData[,c("Chromosome", "Start", "End", "Name")], permeth)
permeth <- permeth[nacount == 0,]
table(is.na(permeth)) # All FALSE, 2465 regions remaining
colnames(permeth) <- gsub("_meth", "", colnames(permeth), fixed = TRUE)
colnames(permeth) <- gsub("VMDK", "VMDK00", colnames(permeth), fixed = TRUE)
rm(cov, meth, nacount)

# Annotate
infoCombined <- rbind(geneInfo$hep, geneInfo$hsc, geneInfo$lsec)
infoCombined$CellType <- c(rep("HEP", nrow(geneInfo$hep)), rep("HSC", nrow(geneInfo$hsc)), rep("LSEC", nrow(geneInfo$lsec)))
infoCombined <- subset(infoCombined, entrezgene %in% permeth$Name)

# Fix Duplicates ####
table(duplicated(permeth$Name)) #6 duplicated
table(duplicated(infoCombined$entrezgene)) # 7 duplicated

# Duplicate in permeth:         220074, 3963, 9997, 6314, 54984, 10597
# Duplicate in infoCombined:    220074, 3963, 9997, 6314, 54984, 10597, 727866

# Fix permeth
permeth <- subset(permeth, !(Name == 220074 & Start == 72102924))
permeth <- subset(permeth, !(Name == 3963 & Start == 38786211))
permeth <- subset(permeth, !(Name == 9997 & Start == 50524606))
permeth <- subset(permeth, !(Name == 6314 & Start == 63895399))
permeth <- subset(permeth, !(Name == 54984 & Start == 10838847))
permeth <- subset(permeth, !(Name == 10597 & Chromosome == "chrX"))
permeth <- subset(permeth, !Name == 727866)

# Fix infoCombined
infoCombined <- subset(infoCombined, !(entrezgene == 220074 & promoter_start == 72102924))
infoCombined <- subset(infoCombined, !(entrezgene == 3963 & promoter_start == 38786211))
infoCombined <- subset(infoCombined, !(entrezgene == 9997 & promoter_start == 50524606))
infoCombined <- subset(infoCombined, !(entrezgene == 6314 & promoter_start == 63895399))
infoCombined <- subset(infoCombined, !(entrezgene == 54984 & promoter_start == 10838847))
infoCombined <- subset(infoCombined, !(entrezgene == 10597 & chromosome_name == "X"))
infoCombined <- subset(infoCombined, !entrezgene == 727866)

# Merge Tables
permeth <- permeth[order(permeth$Name),]
infoCombined <- infoCombined[order(infoCombined$entrezgene),]
table(permeth$Name == infoCombined$entrezgene) # All TRUE, 2458 promoters
table(duplicated(permeth$Name)) # All FALSE
table(duplicated(infoCombined$entrezgene)) # All FALSE
permeth <- merge(x = permeth, y = infoCombined[,c("entrezgene", "external_gene_name", "strand", "CellType")], by.x = "Name", 
                 by.y = "entrezgene", all = FALSE, sort = FALSE)
sampleNames <- colnames(permeth)[grepl("VMDK", colnames(permeth), fixed = TRUE)]
permeth <- permeth[,c("Name", "Chromosome", "Start", "End", "strand", "external_gene_name", "CellType", sampleNames)]
colnames(permeth) <- c("entrezID", "chr", "start", "end", "strand", "geneName", "cellType", sampleNames)
write.table(permeth, file = "Tables/Specific Promoter Methylation with Gene Info.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(combined, countData, entrezIDs, geneInfo, infoCombined)

# Compare Methylation by Diagnosis ####
samples <- samples[match(sampleNames, samples$SequencingID),]
table(samples$SequencingID == sampleNames) # All TRUE
samples$Diagnosis <- factor(samples$Diagnosis, levels = c("Healthy", "NAFLD", "WD"))
statsMethByDiagnosis <- methByDiagnosis(data = permeth, samples = samples)
write.table(statsMethByDiagnosis, file = "Tables/Specific Promoter Methylation Stats by Diagnosis.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Summary Significant
sigMethByDiagnosis <- rbind(table(statsMethByDiagnosis$diff_NAFLDvsHealthy > 0 & statsMethByDiagnosis$Sig_NAFLDvsHealthy,
                                  statsMethByDiagnosis$diff_NAFLDvsHealthy < 0 & statsMethByDiagnosis$Sig_NAFLDvsHealthy,
                                  statsMethByDiagnosis$cellType) %>% as.data.frame,
                            table(statsMethByDiagnosis$diff_WDvsHealthy > 0 & statsMethByDiagnosis$Sig_WDvsHealthy,
                                  statsMethByDiagnosis$diff_WDvsHealthy < 0 & statsMethByDiagnosis$Sig_WDvsHealthy,
                                  statsMethByDiagnosis$cellType) %>% as.data.frame,
                            table(statsMethByDiagnosis$diff_WDvsNAFLD > 0 & statsMethByDiagnosis$Sig_WDvsNAFLD,
                                  statsMethByDiagnosis$diff_WDvsNAFLD < 0 & statsMethByDiagnosis$Sig_WDvsNAFLD,
                                  statsMethByDiagnosis$cellType) %>% as.data.frame)
colnames(sigMethByDiagnosis) <- c("HyperSig", "HypoSig", "cellType", "Freq")
sigMethByDiagnosis <- subset(sigMethByDiagnosis, !(HyperSig =="TRUE" & HypoSig == "TRUE"))
sigMethByDiagnosis$Comparison <- rep(c("NAFLDvsHealthy", "WDvsHealthy", "WDvsNAFLD"), each = 9)
sigMethByDiagnosis$Total <- rep(rep(table(statsMethByDiagnosis$cellType) %>% as.numeric, each = 3), 3)
sigMethByDiagnosis$Direction <- rep(c("NotSig", "HyperSig", "HypoSig"), 9)
sigMethByDiagnosis$Percent <- sigMethByDiagnosis$Freq * 100 / sigMethByDiagnosis$Total
sigMethByDiagnosis <- sigMethByDiagnosis[c("Comparison", "cellType", "Direction", "Freq", "Total", "Percent")]
sigMethByDiagnosis$Comparison <- factor(sigMethByDiagnosis$Comparison, levels = c("NAFLDvsHealthy", "WDvsHealthy", "WDvsNAFLD"))
sigMethByDiagnosis$Direction <- factor(sigMethByDiagnosis$Direction, levels = c("NotSig", "HypoSig", "HyperSig"))
write.table(sigMethByDiagnosis, "Tables/Specific Promoter Differential Significance Summary by Diagnosis.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Compare methylation by Early Diagnosis ####
sampleNamesEarly <- c("VMDK001C", "VMDK001D", "VMDK004C", "VMDK005C", "VMDK006B", "VMDK006C", "VMDK004B","VMDK005B", 
                      "VMDK002B", "VMDK006A", "VMDK002D", "VMDK003C", "VMDK003D")
samplesEarly <- subset(samples, SequencingID %in% sampleNamesEarly)
samplesEarly <- samplesEarly[match(sampleNamesEarly, samplesEarly$SequencingID),]
table(samplesEarly$SequencingID == sampleNamesEarly) # All TRUE
samplesEarly$Diagnosis <- factor(samplesEarly$Diagnosis, levels = c("Healthy", "NAFLD", "WD"))
statsMethByEarlyDiagnosis <- methByDiagnosis(data = permeth, samples = samplesEarly)
write.table(statsMethByEarlyDiagnosis, file = "Tables/Specific Promoter Methylation Stats by Early Diagnosis.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Summary Significant
sigMethByEarlyDiagnosis <- rbind(table(statsMethByEarlyDiagnosis$diff_NAFLDvsHealthy > 0 & statsMethByEarlyDiagnosis$Sig_NAFLDvsHealthy,
                                  statsMethByEarlyDiagnosis$diff_NAFLDvsHealthy < 0 & statsMethByEarlyDiagnosis$Sig_NAFLDvsHealthy,
                                  statsMethByEarlyDiagnosis$cellType) %>% as.data.frame,
                            table(statsMethByEarlyDiagnosis$diff_WDvsHealthy > 0 & statsMethByEarlyDiagnosis$Sig_WDvsHealthy,
                                  statsMethByEarlyDiagnosis$diff_WDvsHealthy < 0 & statsMethByEarlyDiagnosis$Sig_WDvsHealthy,
                                  statsMethByEarlyDiagnosis$cellType) %>% as.data.frame,
                            table(statsMethByEarlyDiagnosis$diff_WDvsNAFLD > 0 & statsMethByEarlyDiagnosis$Sig_WDvsNAFLD,
                                  statsMethByEarlyDiagnosis$diff_WDvsNAFLD < 0 & statsMethByEarlyDiagnosis$Sig_WDvsNAFLD,
                                  statsMethByEarlyDiagnosis$cellType) %>% as.data.frame)
colnames(sigMethByEarlyDiagnosis) <- c("HyperSig", "HypoSig", "cellType", "Freq")
sigMethByEarlyDiagnosis <- subset(sigMethByEarlyDiagnosis, !(HyperSig =="TRUE" & HypoSig == "TRUE"))
sigMethByEarlyDiagnosis$Comparison <- rep(c("NAFLDvsHealthy", "WDvsHealthy", "WDvsNAFLD"), each = 9)
sigMethByEarlyDiagnosis$Total <- rep(rep(table(statsMethByEarlyDiagnosis$cellType) %>% as.numeric, each = 3), 3)
sigMethByEarlyDiagnosis$Direction <- rep(c("NotSig", "HyperSig", "HypoSig"), 9)
sigMethByEarlyDiagnosis$Percent <- sigMethByEarlyDiagnosis$Freq * 100 / sigMethByEarlyDiagnosis$Total
sigMethByEarlyDiagnosis <- sigMethByEarlyDiagnosis[c("Comparison", "cellType", "Direction", "Freq", "Total", "Percent")]
sigMethByEarlyDiagnosis$Comparison <- factor(sigMethByEarlyDiagnosis$Comparison, levels = c("NAFLDvsHealthy", "WDvsHealthy", "WDvsNAFLD"))
sigMethByEarlyDiagnosis$Direction <- factor(sigMethByEarlyDiagnosis$Direction, levels = c("NotSig", "HypoSig", "HyperSig"))
write.table(sigMethByEarlyDiagnosis, "Tables/Specific Promoter Differential Significance Summary by Early Diagnosis.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Compare methylation by WD Severity ####
samplesSeverity <- subset(samples, Diagnosis == "WD")
samplesSeverity$Diagnosis <- samplesSeverity$Diagnosis %>% as.character
samplesSeverity$Diagnosis[samplesSeverity$SequencingID %in% c("VMDK004B","VMDK005B", "VMDK002B")] <- "EarlyWD"
samplesSeverity$Diagnosis[samplesSeverity$SequencingID %in% c("VMDK003B", "VMDK002A", "VMDK004A", "VMDK001A", "VMDK001B",
                                                              "VMDK005A", "VMDK003A")] <- "AdvancedWD"
samplesSeverity$Diagnosis <- factor(samplesSeverity$Diagnosis, levels = c("EarlyWD", "AdvancedWD"))
statsMethBySeverity <- methByDiagnosis(data = permeth, samples = samplesSeverity)
write.table(statsMethBySeverity, file = "Tables/Specific Promoter Methylation Stats by Severity.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Summary Significant
sigMethBySeverity <- table(statsMethBySeverity$diff_AdvancedWDvsEarlyWD > 0 & statsMethBySeverity$Sig_AdvancedWDvsEarlyWD,
                                statsMethBySeverity$diff_AdvancedWDvsEarlyWD < 0 & statsMethBySeverity$Sig_AdvancedWDvsEarlyWD,
                                statsMethBySeverity$cellType) %>% as.data.frame
colnames(sigMethBySeverity) <- c("HyperSig", "HypoSig", "cellType", "Freq")
sigMethBySeverity <- subset(sigMethBySeverity, !(HyperSig =="TRUE" & HypoSig == "TRUE"))
sigMethBySeverity$Comparison <- "AdvancedWDvsEarlyWD"
sigMethBySeverity$Total <- rep(table(statsMethBySeverity$cellType) %>% as.numeric, each = 3)
sigMethBySeverity$Direction <- rep(c("NotSig", "HyperSig", "HypoSig"), 3)
sigMethBySeverity$Percent <- sigMethBySeverity$Freq * 100 / sigMethBySeverity$Total
sigMethBySeverity <- sigMethBySeverity[c("Comparison", "cellType", "Direction", "Freq", "Total", "Percent")]
sigMethBySeverity$Direction <- factor(sigMethBySeverity$Direction, levels = c("NotSig", "HypoSig", "HyperSig"))
write.table(sigMethBySeverity, "Tables/Specific Promoter Differential Significance Summary by Severity.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Stacked barplots ####
# Methylation Significance and Direction by Diagnosis
g <- ggplot(sigMethByDiagnosis)
g +
        geom_col(aes(x = cellType, y = Percent, fill = Direction), size = 1, position = position_stack()) +
        facet_grid(cols = vars(Comparison)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.7, 1.13), panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(size = 0.9), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.title=element_blank(), axis.text = element_text(color = "black", size = 18), 
              axis.text.x = element_text(), strip.background = element_blank(),
              legend.background = element_blank(), plot.margin = unit(c(2,1,1,1), "lines")) +
        scale_fill_manual(breaks = c("HyperSig", "HypoSig", "NotSig"), 
                          values = c("HyperSig" = "#FF3366", "HypoSig" = "#3366CC", "NotSig" = "#009933")) +
        xlab("Cell Type") +
        ylab("Promoters (%)") +
        scale_y_continuous(breaks = pretty_breaks(n = 7), expand = c(0.004,0))
ggsave("Figures/Specific Promoter Differential Significance Summary by Diagnosis.png", dpi = 600, width = 8, height = 7, units = "in")

# Methylation Significance and Direction by Early Diagnosis
g <- ggplot(sigMethByEarlyDiagnosis)
g +
        geom_col(aes(x = cellType, y = Percent, fill = Direction), size = 1, position = position_stack()) +
        facet_grid(cols = vars(Comparison)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.7, 1.13), panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(size = 0.9), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.title=element_blank(), axis.text = element_text(color = "black", size = 18), 
              axis.text.x = element_text(), strip.background = element_blank(),
              legend.background = element_blank(), plot.margin = unit(c(2,1,1,1), "lines")) +
        scale_fill_manual(breaks = c("HyperSig", "HypoSig", "NotSig"), 
                          values = c("HyperSig" = "#FF3366", "HypoSig" = "#3366CC", "NotSig" = "#009933")) +
        xlab("Cell Type") +
        ylab("Promoters (%)") +
        scale_y_continuous(breaks = pretty_breaks(n = 7), expand = c(0.004,0))
ggsave("Figures/Specific Promoter Differential Significance Summary by Early Diagnosis.png", dpi = 600, width = 8, height = 7, units = "in")

# Methylation Significance and Direction by Severity
g <- ggplot(sigMethBySeverity)
g +
        geom_col(aes(x = cellType, y = Percent, fill = Direction), size = 1, position = position_stack()) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.3, 1.07), panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(size = 0.9), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.title=element_blank(), axis.text = element_text(color = "black", size = 18), 
              axis.text.x = element_text(), strip.background = element_blank(),
              legend.background = element_blank(), plot.margin = unit(c(2.5,1,1,1), "lines")) +
        scale_fill_manual(breaks = c("HyperSig", "HypoSig", "NotSig"), 
                          values = c("HyperSig" = "#FF3366", "HypoSig" = "#3366CC", "NotSig" = "#009933")) +
        xlab("Cell Type") +
        ylab("Promoters (%)") +
        scale_y_continuous(breaks = pretty_breaks(n = 7), expand = c(0.004,0))
ggsave("Figures/Specific Promoter Differential Significance Summary by Severity.png", dpi = 600, width = 4, height = 6, units = "in")

rm(g, permeth_rel, permeth_rel_agg, permeth_rel_agg_hep, permeth_rel_agg_hsc, 
   permeth_rel_agg_lsec, permeth_rel_agg_sub, sigMethByDiagnosis, sumMethByDiagnosis,
   sampleNames, x, permeth, permeth_rel_sub)

# Overlap Promoters with DMRs ####
# Load Regions
WD_DMRs <- read.delim(file = "Tables/WD Specific Liver DMRs.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
        .[,c("chr", "start", "end", "direction")]
WDE_DMRs <- read.delim(file = "Tables/WD Early Specific DMRs.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
        .[,c("chr", "start", "end", "direction")]
WDA_DMRs <- read.delim(file = "Tables/WD Adv vs Early Liver DMRs.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
        .[,c("chr", "start", "end", "direction")]
promoters <- statsMethByDiagnosis[,c("chr", "start", "end", "entrezID", "cellType")]

# Make GRanges
makeGRanges <- function(data){GRanges(seqnames = data$chr, ranges = IRanges(start = data$start, end = data$end))}

GR_WD_DMRs <- makeGRanges(WD_DMRs)
GR_WD_DMRs_Hyper <- makeGRanges(subset(WD_DMRs, direction == "hyper"))
GR_WD_DMRs_Hypo <- makeGRanges(subset(WD_DMRs, direction == "hypo"))

GR_WDE_DMRs <- makeGRanges(WDE_DMRs)
GR_WDE_DMRs_Hyper <- makeGRanges(subset(WDE_DMRs, direction == "hyper"))
GR_WDE_DMRs_Hypo <- makeGRanges(subset(WDE_DMRs, direction == "hypo"))

GR_WDA_DMRs <- makeGRanges(WDA_DMRs)
GR_WDA_DMRs_Hyper <- makeGRanges(subset(WDA_DMRs, direction == "hyper"))
GR_WDA_DMRs_Hypo <- makeGRanges(subset(WDA_DMRs, direction == "hypo"))

GR_promoters_hep <- makeGRanges(subset(promoters, cellType == "HEP"))
GR_promoters_hsc <- makeGRanges(subset(promoters, cellType == "HSC"))
GR_promoters_lsec <- makeGRanges(subset(promoters, cellType == "LSEC"))

# Intersect Regions
GR_DMRs <- list("WD_Hyper" = GR_WD_DMRs_Hyper, "WD_Hypo" = GR_WD_DMRs_Hypo, "WDE_Hyper" = GR_WDE_DMRs_Hyper, 
             "WDE_Hypo" = GR_WDE_DMRs_Hypo, "WDA_Hyper" = GR_WDA_DMRs_Hyper, "WDA_Hypo" = GR_WDA_DMRs_Hypo)
GR_promoters <- list("HEP" = GR_promoters_hep, "HSC" = GR_promoters_hsc, "LSEC" = GR_promoters_lsec)

overlapSum <- matrix(nrow = length(GR_DMRs), ncol = length(GR_promoters))
for(i in 1:length(GR_DMRs)){
        for(j in 1:length(GR_promoters)){
                overlapSum[i,j] <- length(intersect(GR_DMRs[[i]], GR_promoters[[j]]))
        }
}
overlapSum <- as.data.frame(overlapSum)
colnames(overlapSum) <- names(GR_promoters)
overlapSum$DMRs <- names(GR_DMRs)
overlapSum$Total <- sapply(GR_DMRs, length)
overlapSum <- overlapSum[,c("DMRs", "HEP", "HSC", "LSEC", "Total")]
write.table(overlapSum, file = "Tables/DMR Promoter Overlap Summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Stacked Barplot of Overlapping Regions
overlapTest <- outer(X = seq_along(GR_DMRs), Y = seq_along(GR_promoters), 
                     FUN = Vectorize(function(i,j){c(names(GR_DMRs)[i], names(GR_promoters)[j], 
                                                     overlapsAny(GR_DMRs[[i]],GR_promoters[[j]]) %>% table %>% as.character)}))
overlapTest <- lapply(overlapTest, function(x){if(length(x) == 3){x <- c(x, '0')} else {x <- x}}) %>% 
        as.data.frame %>% t %>% as.data.frame(stringsAsFactors = FALSE)
rownames(overlapTest) <- 1:nrow(overlapTest)
colnames(overlapTest) <- c("DMRs", "Promoters", "OverlapFalse", "OverlapTrue")
overlapTest$OverlapFalse <- overlapTest$OverlapFalse %>% as.integer
overlapTest$OverlapTrue <- overlapTest$OverlapTrue %>% as.integer
DMRnames <- names(GR_DMRs)
NoOverlapSum <- NULL
for(i in 1:length(DMRnames)){
        NoOverlap <- length(GR_DMRs[[i]]) - sum(overlapTest$OverlapTrue[overlapTest$DMRs == DMRnames[i]])
        temp <- c(DMRnames[i], "None", 0, NoOverlap)
        NoOverlapSum <- rbind(NoOverlapSum, temp)
}
NoOverlapSum <- NoOverlapSum %>% as.data.frame(stringsAsFactors = FALSE)
rownames(NoOverlapSum) <- 1:nrow(NoOverlapSum)
colnames(NoOverlapSum) <- colnames(overlapTest)
NoOverlapSum$OverlapFalse <- NoOverlapSum$OverlapFalse %>% as.integer
NoOverlapSum$OverlapTrue <- NoOverlapSum$OverlapTrue %>% as.integer
overlapTest <- rbind(overlapTest, NoOverlapSum)
overlapTest$Comparison <- rep(rep(c("WD", "WDE", "WDA"), each = 2), 4)
overlapTest$Direction <- rep(c("Hyper", "Hypo"), 12)
overlapTest$DMRs <- factor(overlapTest$DMRs, levels = names(GR_DMRs))
overlapTest$Promoters <- factor(overlapTest$Promoters, levels = c("None", "HEP", "HSC", "LSEC"))
overlapTest$Comparison <- factor(overlapTest$Comparison, levels = c("WD", "WDE", "WDA"))
overlapTest$Direction <- factor(overlapTest$Direction, levels = c("Hyper", "Hypo"))

g <- ggplot(overlapTest)
g +
        geom_col(aes(x = Direction, y = OverlapTrue, fill = Promoters), size = 1, position = position_fill()) +
        facet_grid(cols = vars(Comparison)) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.58, 1.13), panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(size = 0.9), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), panel.grid.minor = element_blank(), 
              legend.title=element_text(size = 20), axis.text = element_text(color = "black", size = 18), 
              axis.text.x = element_text(), strip.background = element_blank(), legend.text = element_text(size = 18),
              legend.background = element_blank(), plot.margin = unit(c(2,1,1,1), "lines")) +
        scale_fill_manual(name = "Promoters", breaks = c("HEP", "HSC", "LSEC", "None"), 
                          values = c("HEP" = "#FF3366", "HSC" = "#3366CC", "LSEC" = "#9933CC", "None" = "#009933")) +
        xlab("Direction") +
        ylab("Overlapping DMRs (%)") +
        scale_y_continuous(breaks = pretty_breaks(n = 7), expand = c(0.004,0), 
                           trans = trans_new(name = "temp", transform = function(x) x / 100, inverse = function(x) x * 100))
ggsave("Figures/Specific Promoter DMR Overlap.png", dpi = 600, width = 8, height = 7, units = "in")
