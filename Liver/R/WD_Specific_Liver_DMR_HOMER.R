# WD Specific Liver DMR HOMER TF Analysis ####
# Charles Mordaunt
# 12/14/17

# Packages ####
library(ggplot2)
library(reshape2)
library(reshape)
library(scales)

# Known Motif Enrichment Results

# All DMRs
homerAll <- read.delim("Homer/All DMRs/knownResults.txt", sep="\t", stringsAsFactors=FALSE)
colnames(homerAll) <- c("MotifName", "Consensus", "Pvalue", "LogPvalue", "qValue", "NumTargetSeq", "PerTargetSeq", 
                        "NumBackgroundSeq", "PerBackgroundSeq")
MotifName <- strsplit(homerAll$MotifName, "/")
MotifName1 <- sapply(MotifName, function(x) x[1])
MotifName1 <- strsplit(MotifName1, "(", fixed=TRUE)
TF <- sapply(MotifName1, function(x) x[1])
TFtype <- sapply(MotifName1, function(x) x[2])
TFtype <- as.character(sapply(TFtype, function(x) gsub(")", "", x, fixed=TRUE)))
TFtype <- as.character(sapply(TFtype, function(x) gsub("?", "", x, fixed=TRUE)))
TFtype[TFtype == ",Zf"] <- "Zf"
MotifName2 <- sapply(MotifName, function(x) x[2])
MotifName2 <- strsplit(MotifName2, "-")
CellType <- sapply(MotifName2, function(x) x[1])
Antibody <- sapply(MotifName2, function(x) x[2])
Reference <- sapply(MotifName2, function(x) x[4])
Reference <- as.character(sapply(Reference, function(x) gsub("Seq(", "", x, fixed=TRUE)))
Reference <- as.character(sapply(Reference, function(x) gsub("Chip(", "", x, fixed=TRUE)))
Reference <- as.character(sapply(Reference, function(x) gsub(")", "", x, fixed=TRUE)))
Reference <- as.character(sapply(Reference, function(x) gsub("Seq", "", x, fixed=TRUE)))
homerAll <- cbind(TF, TFtype, CellType, Antibody, Reference, homerAll)
homerAll$PerTargetSeq <- as.numeric(sapply(homerAll$PerTargetSeq, function(x) gsub("%", "", x, fixed=TRUE)))
homerAll$PerBackgroundSeq <- as.numeric(sapply(homerAll$PerBackgroundSeq, function(x) gsub("%", "", x, fixed=TRUE)))
homerAll$LogPvalue <- -homerAll$LogPvalue
homerAll$Pvalue <- 10^-homerAll$LogPvalue
homerAll$qValue <- p.adjust(homerAll$Pvalue, "fdr")
homerAll$LogQvalue <- -log10(homerAll$qValue)
homerAll$FoldEnrichment <- homerAll$PerTargetSeq/homerAll$PerBackgroundSeq
homerAll$Enriched <- sapply(homerAll$qValue, function(x) ifelse(x < 0.05, TRUE, FALSE))
homerAll <- homerAll[order(homerAll$MotifName),]
homerAll$ID <- 1:dim(homerAll)[1]
table(homerAll$Enriched) #258 Enriched Factors at FDR < 0.05

# Hyper DMRs
homerHyper <- read.delim("Homer/Hyper DMRs/knownResults.txt", sep="\t", stringsAsFactors=FALSE)
colnames(homerHyper) <- c("MotifName", "Consensus", "Pvalue", "LogPvalue", "qValue", "NumTargetSeq", "PerTargetSeq", 
                          "NumBackgroundSeq", "PerBackgroundSeq")
MotifName <- strsplit(homerHyper$MotifName, "/")
MotifName1 <- sapply(MotifName, function(x) x[1])
MotifName1 <- strsplit(MotifName1, "(", fixed=TRUE)
TF <- sapply(MotifName1, function(x) x[1])
TFtype <- sapply(MotifName1, function(x) x[2])
TFtype <- as.character(sapply(TFtype, function(x) gsub(")", "", x, fixed=TRUE)))
TFtype <- as.character(sapply(TFtype, function(x) gsub("?", "", x, fixed=TRUE)))
TFtype[TFtype == ",Zf"] <- "Zf"
MotifName2 <- sapply(MotifName, function(x) x[2])
MotifName2 <- strsplit(MotifName2, "-")
CellType <- sapply(MotifName2, function(x) x[1])
Antibody <- sapply(MotifName2, function(x) x[2])
Reference <- sapply(MotifName2, function(x) x[4])
Reference <- as.character(sapply(Reference, function(x) gsub("Seq(", "", x, fixed=TRUE)))
Reference <- as.character(sapply(Reference, function(x) gsub("Chip(", "", x, fixed=TRUE)))
Reference <- as.character(sapply(Reference, function(x) gsub(")", "", x, fixed=TRUE)))
Reference <- as.character(sapply(Reference, function(x) gsub("Seq", "", x, fixed=TRUE)))
homerHyper <- cbind(TF, TFtype, CellType, Antibody, Reference, homerHyper)
homerHyper$PerTargetSeq <- as.numeric(sapply(homerHyper$PerTargetSeq, function(x) gsub("%", "", x, fixed=TRUE)))
homerHyper$PerBackgroundSeq <- as.numeric(sapply(homerHyper$PerBackgroundSeq, function(x) gsub("%", "", x, fixed=TRUE)))
homerHyper$LogPvalue <- -homerHyper$LogPvalue
homerHyper$Pvalue <- 10^-homerHyper$LogPvalue
homerHyper$qValue <- p.adjust(homerHyper$Pvalue, "fdr")
homerHyper$LogQvalue <- -log10(homerHyper$qValue)
homerHyper$FoldEnrichment <- homerHyper$PerTargetSeq/homerHyper$PerBackgroundSeq
homerHyper$Enriched <- sapply(homerHyper$qValue, function(x) ifelse(x < 0.05, TRUE, FALSE))
homerHyper <- homerHyper[order(homerHyper$MotifName),]
homerHyper$ID <- 1:dim(homerHyper)[1]
table(homerHyper$Enriched) #205 Enriched Factors at FDR < 0.05

# Hypo DMRs
homerHypo <- read.delim("Homer/Hypo DMRs/knownResults.txt", sep="\t", stringsAsFactors=FALSE)
colnames(homerHypo) <- c("MotifName", "Consensus", "Pvalue", "LogPvalue", "qValue", "NumTargetSeq", "PerTargetSeq", 
                         "NumBackgroundSeq", "PerBackgroundSeq")
MotifName <- strsplit(homerHypo$MotifName, "/")
MotifName1 <- sapply(MotifName, function(x) x[1])
MotifName1 <- strsplit(MotifName1, "(", fixed=TRUE)
TF <- sapply(MotifName1, function(x) x[1])
TFtype <- sapply(MotifName1, function(x) x[2])
TFtype <- as.character(sapply(TFtype, function(x) gsub(")", "", x, fixed=TRUE)))
TFtype <- as.character(sapply(TFtype, function(x) gsub("?", "", x, fixed=TRUE)))
TFtype[TFtype == ",Zf"] <- "Zf"
MotifName2 <- sapply(MotifName, function(x) x[2])
MotifName2 <- strsplit(MotifName2, "-")
CellType <- sapply(MotifName2, function(x) x[1])
Antibody <- sapply(MotifName2, function(x) x[2])
Reference <- sapply(MotifName2, function(x) x[4])
Reference <- as.character(sapply(Reference, function(x) gsub("Seq(", "", x, fixed=TRUE)))
Reference <- as.character(sapply(Reference, function(x) gsub("Chip(", "", x, fixed=TRUE)))
Reference <- as.character(sapply(Reference, function(x) gsub(")", "", x, fixed=TRUE)))
Reference <- as.character(sapply(Reference, function(x) gsub("Seq", "", x, fixed=TRUE)))
homerHypo <- cbind(TF, TFtype, CellType, Antibody, Reference, homerHypo)
homerHypo$PerTargetSeq <- as.numeric(sapply(homerHypo$PerTargetSeq, function(x) gsub("%", "", x, fixed=TRUE)))
homerHypo$PerBackgroundSeq <- as.numeric(sapply(homerHypo$PerBackgroundSeq, function(x) gsub("%", "", x, fixed=TRUE)))
homerHypo$LogPvalue <- -homerHypo$LogPvalue
homerHypo$Pvalue <- 10^-homerHypo$LogPvalue
homerHypo$qValue <- p.adjust(homerHypo$Pvalue, "fdr")
homerHypo$LogQvalue <- -log10(homerHypo$qValue)
homerHypo$FoldEnrichment <- homerHypo$PerTargetSeq/homerHypo$PerBackgroundSeq
homerHypo$Enriched <- sapply(homerHypo$qValue, function(x) ifelse(x < 0.05, TRUE, FALSE))
homerHypo <- homerHypo[order(homerHypo$MotifName),]
homerHypo$ID <- 1:dim(homerHypo)[1]
table(homerHypo$Enriched) #169 Enriched Factors at FDR < 0.05

# Volcano Plot All DMRs
top <- data.frame(TF=homerAll$TF, FoldEnrichment=homerAll$FoldEnrichment, LogQvalue=homerAll$LogQvalue)
top <- subset(top, rank(-top$LogQvalue) <= 20)
top$TF <- as.character(top$TF)
top$TF[top$TF == "Erra"] <- "ERRa"
top$TF[top$TF == "Etv2"] <- "ETV2"
top$TF[top$TF == "Rfx6"] <- "RFX6"
top$TF[top$TF == "Foxa2"] <- "FOXA2"
top$TF[top$TF == "Fli1"] <- "FLI1"
top$TF[top$TF == "AP-2gamma"] <- "AP-2g"
top$TF
#[1] "NF1"          "THRb"         "ERRa"         "ERG"          "ZNF711"       "GABPA"        "NF1-halfsite" "ETS1"        
#[9] "ETV1"         "ZFX"          "FOXM1"        "COUP-TFII"    "ETV2"         "PPARE"        "HNF4a"        "RFX6"        
#[17] "FOXA2"        "FLI1"         "STAT4"        "AP-2g"       

gg <- ggplot(data = homerAll)
gg +
        geom_point(aes(x = FoldEnrichment, y = LogQvalue), color="#3366CC", size=2) +
        geom_text(data=top, aes(x=FoldEnrichment, y=LogQvalue, label=TF), 
                  size=3.5, check_overlap = TRUE, hjust=0, nudge_x=0.05, nudge_y=0) +
        theme_bw(base_size = 24) +
        labs(x="Fold Enrichment", y="-log(q-value)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title = element_text(size=15, color="Black"), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) +
        scale_x_continuous(breaks=pretty_breaks(n=4), limits=c(0,3))
ggsave("Figures/Homer WD Specific DMR TF Motif Enrichment.png", dpi = 600, width = 7, height = 5, units = "in")

# Volcano Plot Hyper DMRs
topHyper <- data.frame(TF=homerHyper$TF, FoldEnrichment=homerHyper$FoldEnrichment, LogQvalue=homerHyper$LogQvalue)
topHyper <- subset(topHyper, rank(-topHyper$LogQvalue) <= 20)
topHyper$TF <- as.character(topHyper$TF)
topHyper$TF[topHyper$TF == "Erra"] <- "ERRa"
topHyper$TF[topHyper$TF == "Rfx6"] <- "RFX6"
topHyper$TF[topHyper$TF == "Foxa2"] <- "FOXA2"
topHyper$TF[topHyper$TF == "Foxa3"] <- "FOXA3"
topHyper$TF[topHyper$TF == "Nr5a2"] <- "LRH1"
topHyper$TF[topHyper$TF == "Tlx?"] <- "TLX"
topHyper$TF[topHyper$TF == "Esrrb"] <- "ESRRb"
topHyper$TF
#[1] "HNF4a"        "NF1-halfsite" "ERRa"         "THRb"         "NF1"          "FOXA2"        "COUP-TFII"    "FOXM1"       
#[9] "ZNF711"       "FOXA1"        "RXR"          "PPARE"        "ZNF322"       "ZFX"          "FOXA3"        "FOXA1"       
#[17] "LRH1"         "TLX"          "RFX6"         "ESRRb"       

gg <- ggplot(data = homerHyper)
gg +
        geom_point(aes(x = FoldEnrichment, y = LogQvalue), color="#3366CC", size=2) +
        geom_text(data=topHyper, aes(x=FoldEnrichment, y=LogQvalue, label=TF), 
                  size=3.5, check_overlap = TRUE, hjust=0, nudge_x=0.05, nudge_y=0) +
        theme_bw(base_size = 24) +
        labs(x="Fold Enrichment", y="-log(q-value)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title = element_text(size=15, color="Black"), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) +
        scale_x_continuous(breaks=pretty_breaks(n=4), limits=c(0,3))
ggsave("Figures/Homer WD Specific Hyper DMR TF Motif Enrichment.png", dpi = 600, width = 7, height = 5, units = "in")

# Volcano Plot Hypo DMRs
topHypo <- data.frame(TF=homerHypo$TF, FoldEnrichment=homerHypo$FoldEnrichment, LogQvalue=homerHypo$LogQvalue)
topHypo <- subset(topHypo, rank(-topHypo$LogQvalue) <= 20)
topHypo$TF <- as.character(topHypo$TF)
topHypo$TF[topHypo$TF == "Atf3"] <- "ATF3"
topHypo$TF[topHypo$TF == "Ets1-distal"] <- "ETS1-distal"
topHypo$TF[topHypo$TF == "Etv2"] <- "ETV2"
topHypo$TF[topHypo$TF == "Fli1"] <- "FLI1"
topHypo$TF[topHypo$TF == "Fosl2"] <- "FOSL2"
topHypo$TF[topHypo$TF == "Fra1"] <- "FRA1"
topHypo$TF[topHypo$TF == "Fra2"] <- "FRA2"
topHypo$TF[topHypo$TF == "Jun-AP1"] <- "JUN-AP1"
topHypo$TF[topHypo$TF == "JunB"] <- "JUNB"
topHypo$TF
# [1] "ETS1"            "ERG"             "ETV2"            "FLI1"            "ETV1"            "EHF"             "GABPA"          
# [8] "EWS:ERG-fusion"  "ELF3"            "EWS:FLI1-fusion" "ETS1-distal"     "PU.1"            "FRA2"            "BATF"           
# [15] "FRA1"            "AP-1"            "JUNB"            "FOSL2"           "ATF3"            "JUN-AP1"

table(topHypo$TF %in% top$TF) #6
table(topHyper$TF %in% top$TF) #12
table(topHypo$TF %in% topHyper$TF) #0

enrichedAll <- subset(homerAll, qValue < 1e-10) #50 / 364
enrichedHyper <- subset(homerHyper, qValue < 1e-10) #28 / 364
enrichedHypo <- subset(homerHypo, qValue < 1e-10) #25 / 364

table(enrichedAll$MotifName %in% enrichedHyper$MotifName) #20
table(enrichedAll$MotifName %in% enrichedHypo$MotifName) #22
table(enrichedHyper$MotifName %in% enrichedHypo$MotifName) #0 None Overlap between hyper and hypo


gg <- ggplot(data = homerHypo)
gg +
        geom_point(aes(x = FoldEnrichment, y = LogQvalue), color="#3366CC", size=2) +
        geom_text(data=topHypo, aes(x=FoldEnrichment, y=LogQvalue, label=TF), 
                  size=3.5, check_overlap = TRUE, hjust=0, nudge_x=0.05, nudge_y=0) +
        theme_bw(base_size = 24) +
        labs(x="Fold Enrichment", y="-log(q-value)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title = element_text(size=15, color="Black"), 
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) #+
#scale_x_continuous(breaks=pretty_breaks(n=4), limits=c(0,3))
ggsave("Figures/Homer WD Specific Hypo DMR TF Motif Enrichment.png", dpi = 600, width = 7, height = 5, units = "in")

# Scatter Plot Hyper vs Hypo DMR log q-value
LogQ <- data.frame(TF=homerAll$TF, LogQvalueHyper=homerHyper$LogQvalue, LogQvalueHypo=homerHypo$LogQvalue)
topLogQ <- subset(LogQ, rank(-LogQ$LogQvalueHyper) <= 10 | rank(-LogQ$LogQvalueHypo) <= 10)
topLogQ$TF <- as.character(topLogQ$TF)
topLogQ$TF[topLogQ$TF == "Erra"] <- "ERRa"
topLogQ$TF[topLogQ$TF == "Foxa2"] <- "FOXA2"
topLogQ$TF[topLogQ$TF == "Etv2"] <- "ETV2"
topLogQ$TF[topLogQ$TF == "Fli1"] <- "FLI1"

topLogQ$TF
# [1] "COUP-TFII"       "EHF"             "ELF3"            "ERG"             "ERRa"            "ETS1"            "ETV1"           
# [8] "ETV2"            "EWS:ERG-fusion"  "EWS:FLI1-fusion" "FLI1"            "FOXA1"           "FOXA2"           "FOXM1"          
# [15] "GABPA"           "HNF4a"           "NF1-halfsite"    "NF1"             "THRb"            "ZNF711"  

gg <- ggplot(LogQ)
gg +
        geom_point(aes(x = LogQvalueHyper, y = LogQvalueHypo), color="#3366CC", size=2) +
        geom_text(data=topLogQ, aes(x=LogQvalueHyper, y=LogQvalueHypo, label=TF), 
                  size=3.5, check_overlap = TRUE, hjust=0, nudge_x=0.5, nudge_y=0) +
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
        coord_cartesian(xlim=c(0,46), ylim=c(0,46))
ggsave("Figures/Homer WD Specific Hyper vs Hypo DMR TF Motif logqvalue.png", dpi = 600, width = 7, height = 5, units = "in")

# LOLA Enriched TFs Hyper DMRs
lolaTFs <- c("HNF4a", "RXR", "Foxa2", "FOXA1", "Max", "Atf3", "JunD", "REST-NRSF", "YY1", "COUP-TFII", "GABPA", "ZBTB33", "CTCF")
homerHyper_lola <- subset(homerHyper, TF %in% lolaTFs)
# Atf3      COUP-TFII CTCF      FOXA1     FOXA1     Foxa2     GABPA     HNF4a     JunD      Max       REST-NRSF RXR       YY1      
# ZBTB33
homerHyper_lola <- subset(homerHyper_lola, !(TF == "FOXA1" & CellType == "MCF7")) #Remove duplicate FOXA1
homerHyper_lola$TF <- as.character(homerHyper_lola$TF)
homerHyper_lola$TF[homerHyper_lola$TF == "HNF4a"] <- "HNF4A"
homerHyper_lola$TF[homerHyper_lola$TF == "Foxa2"] <- "FOXA2"
homerHyper_lola$TF[homerHyper_lola$TF == "COUP-TFII"] <- "NR2F2"
homerHyper_lola$TF[homerHyper_lola$TF =="RXR"] <- "RXRA"
homerHyper_lola$TF[homerHyper_lola$TF =="Max"] <- "MAX"
homerHyper_lola$TF[homerHyper_lola$TF =="JunD"] <- "JUND"
homerHyper_lola$TF[homerHyper_lola$TF =="Atf3"] <- "ATF3"
homerHyper_lola$TF[homerHyper_lola$TF =="REST-NRSF"] <- "REST"

homerHyper_lola <- homerHyper_lola[order(homerHyper_lola$LogQvalue),]
homerHyper_lola$TF <- factor(homerHyper_lola$TF, levels=unique(homerHyper_lola$TF), ordered=TRUE)
gg <- ggplot(data = homerHyper_lola)
gg +
        geom_col(aes(x = TF, y = LogQvalue), color="#3366CC", fill="white", size=1.2, position="dodge") +
        coord_flip(ylim=c(0,45)) +
        theme_bw(base_size = 24) +
        labs(y="-log(q-value)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(size = 1.25), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title.x = element_text(size=15, color="Black"), 
              axis.title.y = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) +
        scale_y_continuous(breaks=pretty_breaks(n=4))
ggsave("Figures/Homer WD Specific Hyper DMR LOLA Liver TF Motif Enrichment.png", dpi = 600, width = 5, height = 7, units = "in")

# LOLA Enriched TFs Hypo DMRs
lolaTFs <- c("HNF4a", "RXR", "Foxa2", "FOXA1", "Max", "Atf3", "JunD", "REST-NRSF", "YY1", "COUP-TFII", "GABPA", "ZBTB33", "CTCF")
homerHypo_lola <- subset(homerHypo, TF %in% lolaTFs)
# Atf3      COUP-TFII CTCF      FOXA1     FOXA1     Foxa2     GABPA     HNF4a     JunD      Max       REST-NRSF RXR       YY1      
# ZBTB33
homerHypo_lola <- subset(homerHypo_lola, !(TF == "FOXA1" & CellType == "MCF7")) #Remove duplicate FOXA1
homerHypo_lola$TF <- as.character(homerHypo_lola$TF)
homerHypo_lola$TF[homerHypo_lola$TF == "HNF4a"] <- "HNF4A"
homerHypo_lola$TF[homerHypo_lola$TF == "Foxa2"] <- "FOXA2"
homerHypo_lola$TF[homerHypo_lola$TF == "COUP-TFII"] <- "NR2F2"
homerHypo_lola$TF[homerHypo_lola$TF =="RXR"] <- "RXRA"
homerHypo_lola$TF[homerHypo_lola$TF =="Max"] <- "MAX"
homerHypo_lola$TF[homerHypo_lola$TF =="JunD"] <- "JUND"
homerHypo_lola$TF[homerHypo_lola$TF =="Atf3"] <- "ATF3"
homerHypo_lola$TF[homerHypo_lola$TF =="REST-NRSF"] <- "REST"

homerHypo_lola <- homerHypo_lola[match(homerHyper_lola$TF, homerHypo_lola$TF),]
homerHypo_lola$TF <- factor(homerHypo_lola$TF, levels=unique(homerHypo_lola$TF), ordered=TRUE)
gg <- ggplot(data = homerHypo_lola)
gg +
        geom_col(aes(x = TF, y = LogQvalue), color="#3366CC", fill="white", size=1.2, position="dodge") +
        coord_flip(ylim=c(0,45)) +
        theme_bw(base_size = 24) +
        labs(y="-log(q-value)") +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks.x = element_line(size = 1.25), legend.key = element_blank(), axis.ticks.y=element_blank(),
              panel.grid.minor = element_blank(), legend.position = c(1.2, 0.84), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(size = 15, color = "Black"),
              axis.title.x = element_text(size=15, color="Black"), 
              axis.title.y = element_blank(),
              legend.title = element_text(size = 18),
              plot.title = element_text(size = 18)) +
        scale_y_continuous(breaks=pretty_breaks(n=4))
ggsave("Figures/Homer WD Specific Hypo DMR LOLA Liver TF Motif Enrichment.png", dpi = 600, width = 5, height = 7, units = "in")
