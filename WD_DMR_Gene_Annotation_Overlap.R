# WD DMR Gene Annotation and Overlap ####
# Charles Mordaunt
# 1/6/18

# Packages ####
library(ggplot2)
library(reshape2)
library(scales)
library(ggbiplot)
library(GenomicRanges)
library(dplyr)
library(VennDiagram)
library(GeneOverlap)
library(biomaRt)
library(grid)

# Functions ####
getDMRgenes <- function(DMRstats, regDomains){
        cat("[getDMRgenes] Making GRanges.\n")
        GR_regDomains <- GRanges(seqnames = regDomains$gene_chr, ranges=IRanges(start=regDomains$distal_start, end=regDomains$distal_end))
        GR_DMRstats <- GRanges(seqnames = DMRstats$chr, ranges=IRanges(start=DMRstats$start, end=DMRstats$end))
        cat("[getDMRgenes] Adding overlapping genes to DMRs.\n")
        overlaps <- as.data.frame(findOverlaps(GR_DMRstats, GR_regDomains))
        rm(GR_regDomains, GR_DMRstats)
        DMRstats_genes <- cbind("DMRid"=DMRstats$DMRid[overlaps$queryHits], regDomains[overlaps$subjectHits,], row.names=NULL)
        rm(overlaps)
        DMRstats_genes <- merge(DMRstats, DMRstats_genes, by="DMRid", all=TRUE, sort=FALSE)
        cat("[getDMRgenes] Getting DMR positions relative to genes.\n")
        cat("[getDMRgenes] Positions added:\t")
        DMRstats_genes_pos <- NULL
        for(i in 1:nrow(DMRstats_genes)){
                if(i %% 500 == 0){cat(i, "\t")}
                temp <- NULL
                temp <- DMRstats_genes[i,]
                if(is.na(temp$gene_strand)){
                        temp$distanceToTSS <- NA
                        temp$position <- NA
                } else {
                        if(temp$gene_strand == "+"){
                                if(temp$start < temp$gene_start & temp$end < temp$gene_start){temp$distanceToTSS <- temp$end - temp$gene_start} #upstream of TSS (-)
                                if(temp$start <= temp$gene_start & temp$end >= temp$gene_start){temp$distanceToTSS <- 0} #overlapping TSS
                                if(temp$start > temp$gene_start & temp$end > temp$gene_start){temp$distanceToTSS <- temp$start - temp$gene_start} #downstream of TSS (+)
                                if(temp$end < temp$gene_start){temp$position <- "upstream"}
                                if(temp$end > temp$gene_start & temp$start < temp$gene_end){temp$position <- "gene_body"}
                                if(temp$start > temp$gene_end){temp$position <- "downstream"}
                        }
                        if(temp$gene_strand == "-"){
                                if(temp$start > temp$gene_end & temp$end > temp$gene_end){temp$distanceToTSS <- temp$gene_end - temp$start} #upstream of TSS (-)
                                if(temp$start <= temp$gene_end & temp$end >= temp$gene_end){temp$distanceToTSS <- 0} #overlapping TSS
                                if(temp$start < temp$gene_end & temp$end < temp$gene_end){temp$distanceToTSS <- temp$gene_end - temp$end} #downstream of TSS (+)
                                if(temp$start > temp$gene_end){temp$position <- "upstream"}
                                if(temp$start < temp$gene_end & temp$end > temp$gene_start){temp$position <- "gene_body"}
                                if(temp$end < temp$gene_start){temp$position <- "downstream"}
                        }
                }
                DMRstats_genes_pos <- rbind(DMRstats_genes_pos, temp)
        }
        cat("\n[getDMRgenes] Complete!\n")
        return(DMRstats_genes_pos)
}

getDMRgeneList <- function(DMRstats, regDomains){
        GR_regDomains <- GRanges(seqnames = regDomains$gene_chr, ranges=IRanges(start=regDomains$distal_start, end=regDomains$distal_end))
        GR_DMRstats <- GRanges(seqnames = DMRstats$chr, ranges=IRanges(start=DMRstats$start, end=DMRstats$end))
        overlaps <- as.data.frame(findOverlaps(GR_DMRstats, GR_regDomains))
        DMRstats_genes <- regDomains[overlaps$subjectHits,]
        geneList <- sort(unique(DMRstats_genes$gene_name))
        return(geneList)
}

# Annotate DMRs with Genes and Get Distance ####
# Table 1
DMRs <- read.delim(file = "Tables/liver gold DMRs.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
DMRs$DMRid <- paste("DMR", 1:nrow(DMRs), sep = "_")
regDomains <- read.delim("Tables/Regulatory domains hg38.txt", header=TRUE, stringsAsFactors = FALSE, sep="\t")
DMRs_genes <- getDMRgenes(DMRstats = DMRs, regDomains = regDomains)
write.table(x = DMRs_genes, file = "Tables/liver gold DMRs with Genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)
rm(DMRs, DMRs_genes)

# Table S5
progDMRs <- read.delim(file = "Tables/liver progression gold DMRs.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
progDMRs$DMRid <- paste("DMR", 1:nrow(progDMRs), sep = "_")
progDMRs_genes <- getDMRgenes(DMRstats = progDMRs, regDomains = regDomains)
write.table(x = progDMRs_genes, file = "Tables/liver progression gold DMRs with Genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)
rm(progDMRs, progDMRs_genes)

# All DMRs
allDMRs <- read.delim(file = "Tables/WD DMRs All.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
allDMRs$DMRid <- paste("DMR", 1:nrow(allDMRs), sep = "_")
allDMRs_genes <- getDMRgenes(DMRstats = allDMRs, regDomains = regDomains)
allDMRs_genes$position[allDMRs_genes$distanceToTSS == 0] <- "TSS"
allDMRs_genes$DMRid <- factor(allDMRs_genes$DMRid, levels = unique(allDMRs_genes$DMRid))
allDMRs_collapse <- aggregate(formula = cbind(gene_name, distanceToTSS, position) ~ DMRid, data = allDMRs_genes, 
                              FUN = function(x) paste(x, collapse = ", "), simplify = TRUE)
allDMRs_table <- merge(x = allDMRs, y = allDMRs_collapse, by = "DMRid", all.x = TRUE, all.y = FALSE, sort = FALSE)
write.table(allDMRs_table, file = "Tables/WD DMRs All with Genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rm(allDMRs_collapse, allDMRs_table, wd, info, allDMRs)

# Adaboost DMRs
adaDMRs <- read.delim("Tables/Adaboost DMRs.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
adaDMRs$DMRid <- paste("DMR", 1:nrow(adaDMRs), sep = "_")
adaDMRs_genes <- getDMRgenes(DMRstats = adaDMRs, regDomains = regDomains)
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
geneNames <- getBM(attributes = c("external_gene_name", "description"), filters = "external_gene_name", 
                   values = adaDMRs_genes$gene_name, mart = ensembl)
geneNames$description <- strsplit(x = geneNames$description, split = " [", fixed = TRUE)
geneNames$description <- sapply(geneNames$description, function(x) x[1])
adaDMRs_genes <- merge(x = adaDMRs_genes, y = geneNames, by.x = "gene_name", by.y = "external_gene_name", all.x = TRUE, 
                       all.y = TRUE, sort = FALSE)
adaDMRs_genes <- adaDMRs_genes[,c("DMRid", "chr", "start", "end", "meanDiff..H...N.", "Direction..H...N.", 
                                  "Adaboost.Importance", "gene_name", "description", "distanceToTSS", "position")]
adaDMRs_genes$DMRid <- factor(adaDMRs_genes$DMRid, levels = sort(unique(adaDMRs_genes$DMRid)))
adaDMRs_genes$description[is.na(adaDMRs_genes$description)] <- adaDMRs_genes$gene_name[is.na(adaDMRs_genes$description)]
adaDMRs_collapse <- aggregate(formula = cbind(gene_name, description, distanceToTSS, position) ~ DMRid, 
                              data = adaDMRs_genes, FUN = function(x) paste(x, collapse = ", "), simplify = TRUE)
adaDMRs_table <- merge(x = adaDMRs, y = adaDMRs_collapse, by = "DMRid", all.x = TRUE, all.y = FALSE, sort = FALSE)
write.table(adaDMRs_table, file = "Tables/Adaboost DMRs with Genes.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DMR Liver and Blood Gene Overlaps ####
# Get DMR Genes
liverDMRs <- read.delim("Tables/WD_Specific_Liver_DMR_Info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
liverHyper <- liverDMRs$gene_name[liverDMRs$meanDiff > 0]
liverHyper <- strsplit(x = liverHyper, split = ", ", fixed = FALSE) %>% unlist %>% unique %>% sort
liverHypo <- liverDMRs$gene_name[liverDMRs$meanDiff < 0]
liverHypo <- strsplit(x = liverHypo, split = ", ", fixed = FALSE) %>% unlist %>% unique %>% sort

bloodDMRs <- read.delim("Tables/WD_Specific_Blood_DMR_Info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
bloodHyper <- bloodDMRs$gene_name[bloodDMRs$meanDiff > 0]
bloodHyper <- strsplit(x = bloodHyper, split = ", ", fixed = FALSE) %>% unlist %>% unique %>% sort
bloodHypo <- bloodDMRs$gene_name[bloodDMRs$meanDiff < 0]
bloodHypo <- strsplit(x = bloodHypo, split = ", ", fixed = FALSE) %>% unlist %>% unique %>% sort

liverBloodHyper <- intersect(liverHyper, bloodHyper) %>% unique %>% sort %>% as.data.frame
write.table(liverBloodHyper, "Tables/Liver and Blood Hypermethylated DMR Genes.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

liverBloodHypo <- intersect(liverHypo, bloodHypo) %>% unique %>% sort %>% as.data.frame
write.table(liverBloodHypo, "Tables/Liver and Blood Hypomethylated DMR Genes.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

# Venn Diagrams
venn.diagram(list("Blood"=bloodHyper, "Liver"=liverHyper), "Figures/WD Blood Liver Hyper DMR GREAT Gene Overlap.png", 
             height=8, width=10, imagetype="png", units="in", fontfamily="sans", cat.fontfamily="sans", 
             fill=c("lightblue", "lightpink"), cex=3, lwd=4, cat.cex=3, cat.pos=c(180,180), 
             main="Hypermethylated DMR Genes", main.fontfamily="sans", main.cex=3)

venn.diagram(list("Blood"=bloodHypo, "Liver"=liverHypo), "Figures/WD Blood Liver Hypo DMR GREAT Gene Overlap.png", 
             height=8, width=10, imagetype="png", units="in", fontfamily="sans", cat.fontfamily="sans", 
             fill=c("lightblue", "lightpink"), cex=3, lwd=4, cat.cex=3, cat.pos=c(180,180), 
             main="Hypomethylated DMR Genes", main.fontfamily="sans", main.cex=3, ext.text = FALSE)

# Background
liverBack <- read.delim(file = "Tables/Healthy_Ctrl_vs_WD_Liver_Subsetted_DMR_Background.bed", sep = "\t", 
                        header = FALSE, stringsAsFactors = FALSE)
colnames(liverBack) <- c("chr", "start", "end")
liverBackGenes <- getDMRgeneList(DMRstats = liverBack, regDomains = regDomains)

bloodBack <- read.delim(file = "Tables/HCvWD Blood 2 Subsetted Background DMRs.bed", sep = "\t", 
                        header = FALSE, stringsAsFactors = FALSE)
colnames(bloodBack) <- c("chr", "start", "end")
bloodBackGenes <- getDMRgeneList(DMRstats = bloodBack, regDomains = regDomains)

liverBloodBackGenes <- intersect(liverBackGenes, bloodBackGenes) %>% unique %>% sort # 24125 genes
write.table(as.data.frame(bloodBackGenes), "Tables/Liver and Blood Background Genes.txt", sep = "\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)

# Stats
hyper_gom <- newGOM(list("blood" = bloodHyper, "liver" = liverHyper), genome.size = length(liverBloodBackGenes))
getMatrix(hyper_gom, "odds.ratio") # 5.019832
getMatrix(hyper_gom, "pval") # 1.145865e-25 

hypo_gom <- newGOM(list("blood" = bloodHypo, "liver" = liverHypo), genome.size = length(liverBloodBackGenes))
getMatrix(hypo_gom, "odds.ratio") # 3.587828
getMatrix(hypo_gom, "pval") # 1.985941e-06 

# Human and Mouse DMR Gene Overlap ####
# Human List
earlyDMRs <- read.delim("Tables/WD_Specific_Early_Liver_DMR_Info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
earlyHyper <- earlyDMRs$gene_name[earlyDMRs$meanDiff > 0]
earlyHyper <- strsplit(x = earlyHyper, split = ", ", fixed = FALSE) %>% unlist %>% unique %>% sort
earlyHypo <- earlyDMRs$gene_name[earlyDMRs$meanDiff < 0]
earlyHypo <- strsplit(x = earlyHypo, split = ", ", fixed = FALSE) %>% unlist %>% unique %>% sort
human_genes <- list("Liver_All" = sort(unique(c(liverHyper, liverHypo))), "Liver_Hyper" = liverHyper, "Liver_Hypo" = liverHypo, 
                    "Liver_Early_All" = sort(unique(c(earlyHyper, earlyHypo))), "Liver_Early_Hyper" = earlyHyper, 
                    "Liver_Early_Hypo" = earlyHypo, "Blood_All" = sort(unique(c(bloodHyper, bloodHypo))), 
                    "Blood_Hyper" = bloodHyper, "Blood_Hypo" = bloodHypo)

earlyBack <- read.delim(file = "Tables/HCvWDE Liver Subsetted Background DMRs.bed", sep = "\t", 
                        header = FALSE, stringsAsFactors = FALSE)
colnames(earlyBack) <- c("chr", "start", "end")
earlyBackGenes <- getDMRgeneList(DMRstats = earlyBack, regDomains = regDomains)
humanBack <- intersect(liverBloodBackGenes, earlyBackGenes) %>% unique %>% sort # 23856 genes

# Mouse List
mouseDMRs <- read.delim("Tables/Mouse_Liver_DMR_Info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mouseDMRs <- subset(mouseDMRs, Comparison == "WTctrl vs txJctrl")
mouseHyper <- mouseDMRs$Genes.from.GREAT[mouseDMRs$meanDiff > 0]
mouseHyper <- strsplit(mouseHyper, split = " ", fixed = TRUE) %>% unlist
mouseHyper <- mouseHyper[!grepl(pattern = "(", x = mouseHyper, fixed = TRUE)] %>% unique %>% sort
mouseHypo <- mouseDMRs$Genes.from.GREAT[mouseDMRs$meanDiff < 0]
mouseHypo <- strsplit(mouseHypo, split = " ", fixed = TRUE) %>% unlist
mouseHypo <- mouseHypo[!grepl(pattern = "(", x = mouseHypo, fixed = TRUE)] %>% unique %>% sort

# Convert Mouse List to Human homologs
human <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
mouse <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

mouseHyperToHuman <- getLDS(attributes = "external_gene_name", filters = "external_gene_name", values = mouseHyper, 
                            mart = mouse, attributesL = "external_gene_name", martL = human, uniqueRows = TRUE)
mouseHyperToHuman <- mouseHyperToHuman$Gene.name.1 %>% unique %>% sort
mouseHypoToHuman <- getLDS(attributes = "external_gene_name", filters = "external_gene_name", values = mouseHypo, 
                            mart = mouse, attributesL = "external_gene_name", martL = human, uniqueRows = TRUE)
mouseHypoToHuman <- mouseHypoToHuman$Gene.name.1 %>% unique %>% sort
mouse_genes <- list("Mouse_All" = sort(unique(c(mouseHyperToHuman, mouseHypoToHuman))), "Mouse_Hyper" = mouseHyperToHuman, 
                    "Mouse_Hypo" = mouseHypoToHuman)

# Overlap
humanMouseGOM <- newGOM(gsetA = human_genes, gsetB = mouse_genes, genome.size = length(humanBack))
oddsRatio <- getMatrix(humanMouseGOM, "odds.ratio") %>% melt
colnames(oddsRatio) <- c("HumanDMR", "MouseDMR", "oddsRatio")

intersects <- getMatrix(humanMouseGOM, "intersection") %>% melt
oddsRatio$intersects <- intersects$value

pvalue <- getMatrix(humanMouseGOM, "pval") %>% melt
oddsRatio$pvalue <- pvalue$value
oddsRatio$qvalue <- p.adjust(p = oddsRatio$pvalue, method = "fdr")
oddsRatio$Sig <- oddsRatio$qvalue < 0.05

oddsRatio$HumanDMR <- gsub(pattern = "_", replacement = " ", x = as.character(oddsRatio$HumanDMR), fixed = TRUE) 
oddsRatio$HumanDMR <- factor(oddsRatio$HumanDMR, levels = rev(c("Liver All", "Liver Hyper", "Liver Hypo", "Liver Early All",
                                                                "Liver Early Hyper", "Liver Early Hypo", "Blood All",
                                                                "Blood Hyper", "Blood Hypo")))
oddsRatio$MouseDMR <- gsub(pattern = "Mouse_", replacement = " ", x = as.character(oddsRatio$MouseDMR), fixed = TRUE) %>% as.factor
oddsRatio$Sig <- factor(oddsRatio$Sig, levels = c("TRUE", "FALSE"))

intersects_genes <- getNestedList(humanMouseGOM, "intersection")
oddsRatio$Genes <- c(sapply(intersects_genes[[1]], function(x) paste(x, collapse = ", ")),
                    sapply(intersects_genes[[2]], function(x) paste(x, collapse = ", ")),
                    sapply(intersects_genes[[3]], function(x) paste(x, collapse = ", ")))
write.table(oddsRatio, file = "Tables/Human Mouse DMR Gene Overlap Analysis.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Odds Ratio Heatmap
gg <- ggplot(data = oddsRatio)
gg +
        geom_tile(aes(x = MouseDMR, y = HumanDMR, fill = oddsRatio)) +
        geom_text(aes(x = MouseDMR, y = HumanDMR, label = intersects), color="white", size=6) +
        geom_text(aes(x = MouseDMR, y = HumanDMR, alpha = Sig), label = "*", color = "white", size = 10, nudge_x = 0.28, 
                  nudge_y = -0.15) +
        scale_fill_gradientn("Odds Ratio", colors = c("Black",  "#FF0000", "#FF0000"), values = c(0,0.5,1), na.value = "black", 
                             limits = c(1, 12.83), breaks = pretty_breaks(n = 6)) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1), legend.key = element_blank(), legend.text=element_text(size=14),
              panel.grid.minor = element_blank(), legend.position = c(1.32, 0.85), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,8,1,1), "lines"), 
              axis.text.x = element_text(size = 14, color = "Black", angle = 0, hjust = 0.5, vjust = 1),
              axis.text.y = element_text(size = 14, color = "Black", angle = 0, hjust = 1, vjust = 0.5),
              axis.title = element_text(size = 16), legend.title = element_text(size = 16)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        scale_alpha_manual(breaks = c("TRUE", "FALSE"), values = c("TRUE" = 1, "FALSE" = 0), guide = FALSE) +
        xlab("Mouse Liver") +
        ylab("Human")
ggsave("Figures/Human Mouse WD DMR Gene Overlap Odds Ratio Heatmap.png", dpi = 600, width = 6, height = 6, units = "in")

# Overlapping Genes
intersects_genes_df <- NULL
for(i in 1:length(intersects_genes)){
        for(j in 1:length(intersects_genes[[i]])){
                temp <- NULL
                temp <- data.frame("MouseDMRs"=rep(names(intersects_genes)[i], length(intersects_genes[[i]][[j]])),
                                   "HumanDMRs"=rep(names(intersects_genes[[i]])[j], length(intersects_genes[[i]][[j]])),
                                   "Gene"=intersects_genes[[i]][[j]])
                intersects_genes_df <- rbind(intersects_genes_df, temp)
        }
}
write.table(intersects_genes_df, "Tables/Mouse Human DMRs Overlapping Genes.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# DMR Drug Target Enrichment ####
# Get WD Progression DMR and Background Genes
progDMRs <- read.delim("Tables/WD_Progression_Liver_DMR_Info.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
progHyper <- progDMRs$gene_name[progDMRs$meanDiff > 0]
progHyper <- strsplit(x = progHyper, split = ", ", fixed = FALSE) %>% unlist %>% unique %>% sort
progHypo <- progDMRs$gene_name[progDMRs$meanDiff < 0]
progHypo <- strsplit(x = progHypo, split = ", ", fixed = FALSE) %>% unlist %>% unique %>% sort
DMRgenes <- append(human_genes, values = list("Liver_Prog_Hyper" = progHyper, "Liver_Prog_Hypo" = progHypo))
DMRgenes <- list("Liver_All" = sort(unique(c(DMRgenes$Liver_Hyper, DMRgenes$Liver_Hypo))), 
                    "Liver_Hyper" = DMRgenes$Liver_Hyper, "Liver_Hypo" = DMRgenes$Liver_Hypo,
                    "Liver_Early_All" = sort(unique(c(DMRgenes$Liver_Early_Hyper, DMRgenes$Liver_Early_Hypo))), 
                    "Liver_Early_Hyper" = DMRgenes$Liver_Early_Hyper, "Liver_Early_Hypo" = DMRgenes$Liver_Early_Hypo,
                    "Liver_Prog_All" = sort(unique(c(DMRgenes$Liver_Prog_Hyper, DMRgenes$Liver_Prog_Hypo))), 
                    "Liver_Prog_Hyper" = DMRgenes$Liver_Prog_Hyper, "Liver_Prog_Hypo" = DMRgenes$Liver_Prog_Hypo,
                    "Blood_All" = sort(unique(c(DMRgenes$Blood_Hyper, DMRgenes$Blood_Hypo))), 
                    "Blood_Hyper" = DMRgenes$Blood_Hyper, "Blood_Hypo" = DMRgenes$Blood_Hypo)

progBack <- read.delim(file = "Tables/Liver WDEvWDA Subsetted Background.bed", sep = "\t", 
                        header = FALSE, stringsAsFactors = FALSE)
colnames(progBack) <- c("chr", "start", "end")
progBackGenes <- getDMRgeneList(DMRstats = progBack, regDomains = regDomains)
DMRback <- intersect(humanBack, progBackGenes) %>% unique %>% sort
rm(human_genes, progDMRs, regDomains, progHyper, progHypo, getDMRgeneList, getDMRgenes, bloodBackGenes, earlyBackGenes,
   humanBack, liverBackGenes, liverBloodBackGenes, progBackGenes, progBack)

# Get Drug gene lists for background
drug <- read.delim(paste("~/Documents/Programming/Wilson's Disease Blood Batch 2/Tables/Blood and Liver Background Gene Drug Interactions.tsv", sep = ""), sep="\t", header=TRUE, stringsAsFactors = FALSE)
table(drug$search_term %in% DMRback)
# FALSE  TRUE 
#  1197 25564 
drug <- subset(drug, search_term %in% DMRback)
drug_names <- unique(unlist(drug$drug))
drug_genes <- lapply(drug_names, function(x){unique(unlist(drug$search_term[drug$drug==x]))})
names(drug_genes) <- drug_names

# Overlap
DMRdrugGOM <- newGOM(gsetA = DMRgenes, gsetB = drug_genes, genome.size = length(DMRback))
DMRdrugSum <- rbind(getMatrix(DMRdrugGOM, "intersection"), getMatrix(DMRdrugGOM, "odds.ratio"), 
                    getMatrix(DMRdrugGOM, "pval")) %>% t %>% as.data.frame
colnames(DMRdrugSum) <- paste(rep(c("intersection", "oddsRatio", "pvalue"), each = length(DMRgenes)), 
                              colnames(DMRdrugSum), sep = "_")
fdr <- sapply(DMRdrugSum[,grepl("pvalue", colnames(DMRdrugSum), fixed = TRUE)], p.adjust, method = "fdr")
colnames(fdr) <- paste("qvalue", names(DMRgenes), sep = "_")
DMRdrugSum <- cbind(DMRdrugSum, fdr)
withIntersect <- rowSums(DMRdrugSum[,grepl("intersection", colnames(DMRdrugSum), fixed = TRUE)]) > 0 # 3153
withPvalue <- apply(DMRdrugSum[,grepl("pvalue", colnames(DMRdrugSum), fixed = TRUE)], 1, min) < 0.005 # 66
DMRdrugSum <- subset(DMRdrugSum, withIntersect & withPvalue) # 66
DMRdrugSum$Drug <- rownames(DMRdrugSum)
DMRdrugSum <- DMRdrugSum[,c(ncol(DMRdrugSum), 1:(ncol(DMRdrugSum) - 1))]
write.table(DMRdrugSum, file = "Tables/DMR Drug Enrichment Analysis.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Get Hits for each DMR set
testIntersect <- DMRdrugSum[,grepl("intersection", colnames(DMRdrugSum), fixed = TRUE)] > 0
testPvalue <- DMRdrugSum[,grepl("pvalue", colnames(DMRdrugSum), fixed = TRUE)] < 0.005
testHit <- as.data.frame(testIntersect & testPvalue)
colnames(testHit) <- names(DMRgenes)
DMRdrugHits <- sapply(testHit, function(x){DMRdrugSum$Drug[x]})

DMRdrugCollapse <- NULL
for(i in 1:length(DMRdrugHits)){
        temp <- DMRdrugSum[DMRdrugHits[[i]], grepl(pattern = names(DMRdrugHits)[i], x = colnames(DMRdrugSum), fixed = TRUE)]
        colnames(temp) <- c("intersection", "oddsRatio", "pvalue", "qvalue")
        if(length(temp$intersection) > 0){
                temp$DMRset <- names(DMRdrugHits)[i]
                temp$Drug <- rownames(temp)
                rownames(temp) <- 1:nrow(temp)
        }
        DMRdrugCollapse <- rbind(DMRdrugCollapse, temp)
}
DMRdrugCollapse <- DMRdrugCollapse[,c("DMRset", "Drug", "intersection", "oddsRatio", "pvalue", "qvalue")]

# Add Genes
DMRdrugGenes <- getNestedList(DMRdrugGOM, "intersection")
DMRdrugGeneList <- list()
for(i in 1:nrow(DMRdrugCollapse)){
        DMRdrugGeneList[[i]] <- DMRdrugGenes[[DMRdrugCollapse$Drug[i]]][[DMRdrugCollapse$DMRset[i]]]
}
table(DMRdrugCollapse$intersection == sapply(DMRdrugGeneList, length)) # All TRUE
DMRdrugGeneList <- sapply(DMRdrugGeneList, paste, collapse = ", ")
DMRdrugCollapse$Genes <- DMRdrugGeneList
write.table(DMRdrugCollapse, file = "Tables/DMR Drug Enrichment Analysis Summary.txt", sep = "\t", quote = FALSE, row.names = FALSE)

