# WD Blood 2 DMR Analysis ####
# Charles Mordaunt
# 1/16/18

# Packages ####
library(GenomicRanges)
library(VennDiagram)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggdendro)
library(ggplot2)
library(reshape2)
library(grid)
library(scales)
library(ggbiplot)
library(R.utils)
library(rtracklayer)
library(LOLA)
library(reshape)
library(devtools)
library(plyr)
library(sm)
library(GeneOverlap)

# Functions ####
combineFiles <- function(chroms, prefix, suffix){
        DMRs <- NULL
        for(i in 1:length(chroms)){
                temp <- NULL
                if(file.exists(paste(prefix,chroms[i],suffix, sep=""))){
                        temp <- read.delim(paste(prefix,chroms[i],suffix, sep=""), header = TRUE, sep = "\t", stringsAsFactors=FALSE)
                        DMRs <- rbind(DMRs, temp)
                }
        }
        DMRs
}

covCutoff <- function(numCtrl, numExp, minDiff){
        precisionCtrl <- (1/numCtrl)*2
        precisionExp <- (1/numExp)*2
        precisionSum <- precisionCtrl + precisionExp
        coverage <- ceiling(precisionSum/minDiff)
        coverage
}

DMRttest <- function(meth_info, numCtrl, numExp, reads){
        meth <- meth_info[,16:ncol(meth_info)]
        DMRs <- 1:nrow(meth)
        ttest <- matrix(nrow = length(DMRs), ncol = 5)
        for(i in DMRs){
                association <- NULL
                association <- t.test(meth[i,(numCtrl+1):(numCtrl+numExp)], meth[i,1:numCtrl])
                ttest[i,1] <- association$estimate[1] - association$estimate[2]
                ttest[i,2] <- association$conf.int[1]
                ttest[i,3] <- association$conf.int[2]
                ttest[i,4] <- association$statistic
                ttest[i,5] <- association$p.value
        }
        colnames(ttest) <- c("meanDiff", "confIntL", "confIntR", "tstat", "pValue")
        ttest <- as.data.frame(ttest)
        ttest$pValueFDR <- p.adjust(ttest$pValue, "fdr")
        ttest$pValueBonf <- p.adjust(ttest$pValue, "bonf")
        
        minreads <- apply(reads, 1, min)
        ttest$minReads <- minreads
        meanreads <- apply(reads, 1, mean)
        ttest$meanReads <- meanreads
        
        ttest <- cbind(meth_info[,1:15], ttest[,2:ncol(ttest)])
        ttest
}

makeGRange <- function(ttest, direction=c("all", "hyper", "hypo")){
        if(direction == "all"){
                GR <- GRanges(seqnames = ttest$chr, ranges=IRanges(start=ttest$start, end=ttest$end))
        }
        if(direction == "hyper"){
                ttest <- subset(ttest, meanDiff > 0)
                GR <- GRanges(seqnames = ttest$chr, ranges=IRanges(start=ttest$start, end=ttest$end))
        }
        if(direction == "hypo"){
                ttest <- subset(ttest, meanDiff < 0)
                GR <- GRanges(seqnames = ttest$chr, ranges=IRanges(start=ttest$start, end=ttest$end))
        }
        GR
}

GRangetoBED <- function(GR, writeFile=TRUE, fileName){
        GR <- as.data.frame(GR)
        BED <- GR[,c("seqnames", "start", "end")]
        BED <- BED[order(BED$seqnames, BED$start),]
        if(writeFile){
                write.table(BED, fileName , sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        }
        GR
}

prepGREAT <- function(DMRs, Background, writeFile=TRUE, fileName, writeBack=FALSE, backName){
        seqlevelsStyle(DMRs) <- "UCSC"
        seqlevelsStyle(Background) <- "UCSC"
        chain <- import.chain("hg38ToHg19.over.chain")
        DMRs_hg19 <- unlist(liftOver(DMRs, chain))
        Background_hg19 <- unlist(liftOver(Background, chain))
        if(!isDisjoint(DMRs_hg19)){DMRs_hg19 <- disjoin(DMRs_hg19)} 
        if(!isDisjoint(Background_hg19)){Background_hg19 <- disjoin(Background_hg19)} 
        if(length(Background_hg19) > 1000000){cat("\nWarning: Need to reduce background to < 1M regions\n")}
        DMRs_hg19 <- redefineUserSets(DMRs_hg19, Background_hg19)
        
        # DMRs
        chroms <- c(paste("chr",1:22,sep=""), "chrM") # Excluded chrX,Y
        DMRs_hg19 <- unlist(GRangesList(DMRs_hg19))
        DMRs_hg19 <- as.data.frame(DMRs_hg19)[,c("seqnames", "start", "end")]
        colnames(DMRs_hg19) <- c("chr", "start", "end")
        DMRs_hg19$chr <- as.character(DMRs_hg19$chr)
        DMRs_hg19 <- DMRs_hg19[order(DMRs_hg19$chr, DMRs_hg19$start),]
        DMRs_hg19 <- unique(subset(DMRs_hg19, chr %in% chroms))
        
        # Background
        Background_hg19 <- as.data.frame(Background_hg19)[,c("seqnames", "start", "end")]
        colnames(Background_hg19) <- c("chr", "start", "end")
        Background_hg19$chr <- as.character(Background_hg19$chr)
        Background_hg19 <- Background_hg19[order(Background_hg19$chr, Background_hg19$start),]
        Background_hg19 <- unique(subset(Background_hg19, chr %in% chroms))
        
        cat(table(DMRs_hg19$chr %in% Background_hg19$chr & DMRs_hg19$start %in% Background_hg19$start & DMRs_hg19$end %in% Background_hg19$end), "DMRs in Background ")
        cat("out of", nrow(DMRs_hg19), "total DMRs.\n")
        if(writeFile){
                write.table(DMRs_hg19, fileName, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
                cat("DMRs written to", fileName, "\n")
                
        }
        if(writeBack){
                write.table(Background_hg19, backName, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
                cat("Background written to", backName)
        }
}

mydplot <- function(ddata, row=!col, col=!row, labels=col) {
        ## plot a dendrogram
        yrange <- range(ddata$segments$y)
        yd <- yrange[2] - yrange[1]
        nc <- max(nchar(as.character(ddata$labels$label)))
        tangle <- if(row) { 0 } else { 90 }
        tshow <- col
        p <- ggplot() +
                geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend), lwd = 0.45) +
                labs(x = NULL, y = NULL) + theme_dendro()
        if(row) {
                p <- p +
                        scale_x_continuous(expand=c(0.5/length(ddata$labels$x),0)) +
                        coord_flip()
        } else {
                p <- p +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, color = "black"))
        }
        return(p)
}

g_legend<-function(a.gplot){
        ## from
        ## http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

ggheatmap.show <- function(L, col.width=0.2, row.width=0.2) {
        grid.newpage()
        top.layout <- grid.layout(nrow = 4, ncol = 4,
                                  widths = unit(c(0.02,0.96-row.width,row.width,0.02), "null"),
                                  heights = unit(c(0.02,col.width,0.96-col.width,0.02), "null"))
        pushViewport(viewport(layout=top.layout))
        if(col.width>0)
                print(L$col, vp=viewport(layout.pos.col=2, layout.pos.row=2))
        if(row.width>0)
                print(L$row, vp=viewport(layout.pos.col=3, layout.pos.row=3))
        ## print centre without legend
        print(L$centre +
                      theme(axis.line=element_blank(),
                            axis.text.x=element_blank(),axis.text.y=element_blank(),
                            axis.ticks=element_blank(),
                            axis.title.x=element_blank(),axis.title.y=element_blank(),
                            legend.position="none",
                            panel.background=element_blank(),
                            panel.border=element_blank(),panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),plot.background=element_blank()),
              vp=viewport(layout.pos.col=2, layout.pos.row=3))
        ## add legend
        legend <- g_legend(L$centre +
                                   theme(legend.title = element_blank(), 
                                         legend.text = element_text(size = 15)))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=2))
        grid.draw(legend)
        upViewport(0)
}

ggheatmap2 <- function(x, custom.label,
                       hm.colours=my.colours, my.values) {
        if(is.null(colnames(x)))
                colnames(x) <- sprintf("col%s",1:ncol(x))
        if(is.null(rownames(x)))
                rownames(x) <- sprintf("row%s",1:nrow(x))
        ## plot a heatmap
        ## x is an expression matrix
        row.hc <- hclust(dist(x), "ward.D")
        col.hc <- hclust(dist(t(x)), "ward.D")
        row.dendro <- dendro_data(as.dendrogram(row.hc),type="rectangle")
        col.dendro <- dendro_data(as.dendrogram(col.hc),type="rectangle")
        
        ## dendro plots
        col.plot <- mydplot(col.dendro, col=TRUE, labels=FALSE) +
                scale_x_continuous(breaks = 1:ncol(x), labels = custom.label[col.hc$order], expand=c(0,0)) +
                theme(plot.margin = unit(c(0,0.35,0.5,0.1), "lines"))
        row.plot <- mydplot(row.dendro, row=TRUE, labels=FALSE) +
                theme(plot.margin = unit(rep(0, 4), "lines"))
        
        ## order of the dendros
        col.ord <- match(col.dendro$labels$label, colnames(x))
        row.ord <- match(row.dendro$labels$label, rownames(x))
        xx <- x[row.ord,col.ord]
        dimnames(xx) <- NULL
        xx <- melt(xx)
        
        centre.plot <- ggplot(xx, aes(X2,X1)) + geom_tile(aes(fill=value, colour=value)) +
                scale_fill_gradientn(colours = hm.colours, values = my.values) +
                scale_colour_gradientn(colours = hm.colours, values = my.values) +
                labs(x = NULL, y = NULL) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0),breaks = NULL) +
                theme(plot.margin = unit(rep(0, 4), "lines"))
        ret <- list(col=col.plot,row=row.plot,centre=centre.plot)
        invisible(ret)
}
# Combine Files ####
chroms <- c(paste("chr",1:22,sep=""), "chrM")

# DMRs
suffix <- "_silver_DMR_info.txt"
DMRs_HCvWD <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_HCvWD/HCvWD_", suffix = suffix)
DMRs_HCvDC <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_HCvDC/HCvDC_", suffix = suffix)
DMRs_DCvWD <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_DCvWD/DCvWD_", suffix = suffix)

# Methylation
suffix <- "_silver_DMR_methylation.txt"
meth_HCvWD <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_HCvWD/HCvWD_", suffix = suffix)
meth_HCvDC <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_HCvDC/HCvDC_", suffix = suffix)
meth_DCvWD <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_DCvWD/DCvWD_", suffix = suffix)

# Coverage
suffix <- "_silver_DMR_coverage.txt"
cov_HCvWD <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_HCvWD/HCvWD_", suffix = suffix)
cov_HCvDC <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_HCvDC/HCvDC_", suffix = suffix)
cov_DCvWD <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_DCvWD/DCvWD_", suffix = suffix)

# Background
backCov_HCvWD <- combineFiles(chroms=chroms, prefix="DMRs/DMRs_HCvWD/HCvWD_", suffix="_background_DMR_coverage.txt")

# Subset HCvWD DMRs ####
# Parameters
numCtrl <- 12
numExp <- 40
numSamples <- numCtrl + numExp
minDiff <- 0.05

# Merge Tables
meth_HCvWD_info <- merge(meth_HCvWD, DMRs_HCvWD, by = c("chr", "start", "end"))
meth_HCvWD_info <- meth_HCvWD_info[,c(1:3,(4+numSamples):ncol(meth_HCvWD_info),4:(3+numSamples))]
cov_HCvWD_info <- merge(cov_HCvWD, DMRs_HCvWD, by = c("chr", "start", "end"))
cov_HCvWD_info <- cov_HCvWD_info[,c(1:3,(4+numSamples):ncol(cov_HCvWD_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_HCvWD_info[,16:ncol(cov_HCvWD_info)])
table(is.na(reads)) # No NA values
if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #5 reads per sample

# T-test and subset
ttest_HCvWD <- DMRttest(meth_info=meth_HCvWD_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_HCvWD <- subset(ttest_HCvWD, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #11064 -> 2608
write.table(ttest_HCvWD, "Tables/HCvWD Blood 2 Subsetted DMRs Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Subset HCvDC DMRs ####
# Parameters
numCtrl <- 12
numExp <- 20
numSamples <- numCtrl + numExp
minDiff <- 0.05

# Merge Tables
meth_HCvDC_info <- merge(meth_HCvDC, DMRs_HCvDC, by = c("chr", "start", "end"))
meth_HCvDC_info <- meth_HCvDC_info[,c(1:3,(4+numSamples):ncol(meth_HCvDC_info),4:(3+numSamples))]
cov_HCvDC_info <- merge(cov_HCvDC, DMRs_HCvDC, by = c("chr", "start", "end"))
cov_HCvDC_info <- cov_HCvDC_info[,c(1:3,(4+numSamples):ncol(cov_HCvDC_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_HCvDC_info[,16:ncol(cov_HCvDC_info)])
table(is.na(reads)) # No NA values
if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #6 reads per sample

# T-test and subset
ttest_HCvDC <- DMRttest(meth_info=meth_HCvDC_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_HCvDC <- subset(ttest_HCvDC, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #14247 -> 3411
write.table(ttest_HCvDC, "Tables/HCvDC Blood 2 Subsetted DMRs Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Subset DCvWD DMRs ####
# Parameters
numCtrl <- 20
numExp <- 40
numSamples <- numCtrl + numExp
minDiff <- 0.05

# Merge Tables
meth_DCvWD_info <- merge(meth_DCvWD, DMRs_DCvWD, by = c("chr", "start", "end"))
meth_DCvWD_info <- meth_DCvWD_info[,c(1:3,(4+numSamples):ncol(meth_DCvWD_info),4:(3+numSamples))]
cov_DCvWD_info <- merge(cov_DCvWD, DMRs_DCvWD, by = c("chr", "start", "end"))
cov_DCvWD_info <- cov_DCvWD_info[,c(1:3,(4+numSamples):ncol(cov_DCvWD_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_DCvWD_info[,16:ncol(cov_DCvWD_info)])
table(is.na(reads)) # No NA values
if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #4 reads per sample

# T-test and subset
ttest_DCvWD <- DMRttest(meth_info=meth_DCvWD_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_DCvWD <- subset(ttest_DCvWD, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #8343 -> 3044
write.table(ttest_DCvWD, "Tables/DCvWD Blood 2 Subsetted DMRs Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)


# Subset HCvWD Background ####
# Parameters
numCtrl <- 12
numExp <- 40
numSamples <- numCtrl + numExp
minDiff <- 0.05

# Coverage
reads <- as.matrix(backCov_HCvWD[,4:ncol(backCov_HCvWD)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #5 reads per sample
backCov_HCvWD$minReads <- apply(reads, 1, min)

# Subset
backCov_HCvWD <- subset(backCov_HCvWD, minReads >= coverage) #621382 -> 277797
background_HCvWD <- backCov_HCvWD[,c("chr", "start", "end")]
background_HCvWD <- background_HCvWD[order(background_HCvWD$chr, background_HCvWD$start),]
write.table(background_HCvWD, "UCSC Tracks/HCvWD Blood 2 Subsetted Background DMRs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Overlap HCvWD, HCvDC, and DCvWD DMRs ####
# Make GRanges
GR_HCvWD <- makeGRange(ttest=ttest_HCvWD, direction="all")
GR_HCvWD_hyper <- makeGRange(ttest=ttest_HCvWD, direction="hyper")
GR_HCvWD_hypo <- makeGRange(ttest=ttest_HCvWD, direction="hypo")

GR_HCvDC <- makeGRange(ttest=ttest_HCvDC, direction="all")
GR_HCvDC_hyper <- makeGRange(ttest=ttest_HCvDC, direction="hyper")
GR_HCvDC_hypo <- makeGRange(ttest=ttest_HCvDC, direction="hypo")

GR_DCvWD <- makeGRange(ttest=ttest_DCvWD, direction="all")
GR_DCvWD_hyper <- makeGRange(ttest=ttest_DCvWD, direction="hyper")
GR_DCvWD_hypo <- makeGRange(ttest=ttest_DCvWD, direction="hypo")

# Make Venn Diagrams
pdf(file="Figures/Subsetted DMR Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_HCvWD, GR_HCvDC, GR_DCvWD), NameOfPeaks = c("HC_vs_WD", "HC_vs_DC", "DC_vs_WD"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen"), cat.pos = c(0,0,0), 
                        cat.dist = c(0.05, 0.05, 0.05))
dev.off()

pdf(file="Figures/Subsetted Hyper DMR Overlap Venn.pdf", width=10, height=9, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_HCvDC_hyper, GR_HCvWD_hyper, GR_DCvWD_hyper), NameOfPeaks = c("HC_vs_DC_Hyper", "HC_vs_WD_Hyper", "DC_vs_WD_Hyper"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 3, fill = c("lightpink", "lightblue", "lightgreen"), cat.pos = c(0,0,0), 
                        cat.dist = c(0.3, 0.2, 0.05), cat.col=c("white", "white", "white"), fontfamily="sans")
dev.off()

pdf(file="Figures/Subsetted Hypo DMR Overlap Venn.pdf", width=10, height=9, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_HCvWD_hypo, GR_HCvDC_hypo, GR_DCvWD_hypo), NameOfPeaks = c("HC_vs_WD_Hypo", "HC_vs_DC_Hypo", "DC_vs_WD_Hypo"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen"), cat.pos = c(0,0,0), 
                        cat.dist = c(0.05, 0.05, 0.05), fontfamily="sans")#, cat.col=c("white", "white", "white"))
dev.off()

# Get Specific DMRs
HCvWD_overlaps_DCvWD_hyper <- overlapsAny(query = GR_HCvWD_hyper, subject = GR_DCvWD_hyper, maxgap = 0, minoverlap = 1, type = "any") #189
HCvWD_overlaps_HCvDC_hyper <- overlapsAny(query = GR_HCvWD_hyper, subject = GR_HCvDC_hyper, maxgap = 0, minoverlap = 1, type = "any") #484
HCvWD_overlaps_DCvWD_hypo <- overlapsAny(query = GR_HCvWD_hypo, subject = GR_DCvWD_hypo, maxgap = 0, minoverlap = 1, type = "any") #75
HCvWD_overlaps_HCvDC_hypo <- overlapsAny(query = GR_HCvWD_hypo, subject = GR_HCvDC_hypo, maxgap = 0, minoverlap = 1, type = "any") #171

GR_WD_Spec_hyper <- GR_HCvWD_hyper[HCvWD_overlaps_DCvWD_hyper & !HCvWD_overlaps_HCvDC_hyper,] #187 DMRs
GR_WD_Spec_hypo <- GR_HCvWD_hypo[HCvWD_overlaps_DCvWD_hypo & !HCvWD_overlaps_HCvDC_hypo,] #75 DMRs
GR_WD_Spec <- c(GR_WD_Spec_hyper, GR_WD_Spec_hypo) #262 DMRs
isDisjoint(GR_WD_Spec) #TRUE (not overlapping)

# Make bed files
DMRs_WD_Spec <- GRangetoBED(GR=GR_WD_Spec, writeFile=TRUE, fileName="UCSC Tracks/WD_Blood_2_Specific_DMRs.bed")
DMRs_WD_Spec_hyper <- GRangetoBED(GR=GR_WD_Spec_hyper, writeFile=TRUE, fileName="UCSC Tracks/WD_Blood_2_Specific_Hyper_DMRs.bed")
DMRs_WD_Spec_hypo <- GRangetoBED(GR=GR_WD_Spec_hypo, writeFile=TRUE, fileName="UCSC Tracks/WD_Blood_2_Specific_Hypo_DMRs.bed")

# Annotate Specific DMRs
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene)
seqlevelsStyle(GR_WD_Spec) <- seqlevelsStyle(annoData)
anno <- annotatePeakInBatch(GR_WD_Spec, featureType = "TSS", AnnotationData=annoData, output="overlapping", FeatureLocForDistance="TSS",
                            bindingRegion=c(-5000, 1000), multiple = TRUE, PeakLocForDistance = "middle", select = "all")
anno$symbol <- xget(anno$feature, org.Hs.egSYMBOL)
anno_WD_Spec <- as.data.frame(cbind(data.frame(chr=seqnames(anno), start=start(anno), end=end(anno)), anno@elementMetadata))
anno_WD_Spec <- anno_WD_Spec[,c("chr", "start", "end", "feature.strand", "distance", "symbol")]

ttest_WD_Spec <- merge(ttest_HCvWD, DMRs_WD_Spec, by.x=c("chr", "start", "end"), by.y=c("seqnames", "start", "end"), all.x=FALSE)
ttest_WD_Specific_anno <- merge(ttest_WD_Spec, anno_WD_Spec, by = c("chr", "start", "end"), all.x = TRUE)
write.table(ttest_WD_Specific_anno, "Tables/WD Specific DMR Stats with annotation.txt", sep = "\t", quote = FALSE, row.names = FALSE)

ttest_WD_Specific_anno_only <- ttest_WD_Specific_anno[,c("chr", "start", "end", "CpGs", "width.x", "areaStat", "Ctrl_mean", "Exp_mean", 
                                                         "meanDiff", "confIntL", "confIntR", "direction", "Rel_FWER", 
                                                         "feature.strand", "distance", "symbol")]
ttest_WD_Specific_anno_only <- ttest_WD_Specific_anno_only[order(ttest_WD_Specific_anno_only$Rel_FWER),]
write.table(ttest_WD_Specific_anno_only, "Tables/WD Specific DMR Stats with annotation 2.txt", sep = "\t", quote = FALSE, row.names = FALSE)

hyperGenes <- sort(unique(ttest_WD_Specific_anno_only$symbol[ttest_WD_Specific_anno_only$direction == "hyper"]))
write.table(hyperGenes, "Tables/WD Specific DMR Hyper Genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

hypoGenes <- sort(unique(ttest_WD_Specific_anno_only$symbol[ttest_WD_Specific_anno_only$direction == "hypo"]))
write.table(hypoGenes, "Tables/WD Specific DMR Hypo Genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# GREAT WD Specific DMRs ####
# Get background file
GR_Background <- GRanges(seqnames = background_HCvWD$chr, ranges = IRanges(start = background_HCvWD$start, end = background_HCvWD$end))

# Get chain file
url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
download.file(url, basename(url))
gunzip(basename(url))

# Make Files for GREAT
prepGREAT(DMRs=GR_WD_Spec, Background=GR_Background, fileName="UCSC Tracks/WD_Blood_2_WD_Specific_DMRs_hg19_for_GREAT.bed", writeBack=TRUE,
          backName="UCSC Tracks/WD_Blood_2_HCvWD_Background_hg19_for_GREAT.bed")
prepGREAT(DMRs=GR_WD_Spec_hyper, Background=GR_Background, fileName="UCSC Tracks/WD_Blood_2_WD_Specific_Hyper_DMRs_hg19_for_GREAT.bed")
prepGREAT(DMRs=GR_WD_Spec_hypo, Background=GR_Background, fileName="UCSC Tracks/WD_Blood_2_WD_Specific_Hypo_DMRs_hg19_for_GREAT.bed")

# Comparison with WD Liver Specific DMRs ####
liverDMRs <- read.delim("C:/Users/Booboo/Charles/Documents/Programming/Wilson's Disease Liver/WD_Specific_DMRs.bed", sep="\t", header=FALSE)
colnames(liverDMRs) <- c("chr", "start", "end")
GR_liverDMRs <- GRanges(seqnames = liverDMRs$chr, ranges=IRanges(start=liverDMRs$start, end=liverDMRs$end))

liverDMRs_hyper <- read.delim("C:/Users/Booboo/Charles/Documents/Programming/Wilson's Disease Liver/WD_Specific_hyper_DMRs.bed", sep="\t", header=FALSE)
colnames(liverDMRs_hyper) <- c("chr", "start", "end")
GR_liverDMRs_hyper <- GRanges(seqnames = liverDMRs_hyper$chr, ranges=IRanges(start=liverDMRs_hyper$start, end=liverDMRs_hyper$end))

liverDMRs_hypo <- read.delim("C:/Users/Booboo/Charles/Documents/Programming/Wilson's Disease Liver/WD_Specific_hypo_DMRs.bed", sep="\t", header=FALSE)
colnames(liverDMRs_hypo) <- c("chr", "start", "end")
GR_liverDMRs_hypo <- GRanges(seqnames = liverDMRs_hypo$chr, ranges=IRanges(start=liverDMRs_hypo$start, end=liverDMRs_hypo$end))

liverBackground <- read.delim("C:/Users/Booboo/Charles/Documents/Programming/Wilson's Disease Liver/Healthy_Ctrl_vs_WD_Liver_Subsetted_DMR_Background.bed", sep="\t", header=FALSE)
colnames(liverBackground) <- c("chr", "start", "end")
GR_liverBackground <- GRanges(seqnames = liverBackground$chr, ranges=IRanges(start=liverBackground$start, end=liverBackground$end))

hyperInBoth <- as.data.frame(GR_WD_Spec_hyper[overlapsAny(query = GR_WD_Spec_hyper, subject = GR_liverDMRs_hyper, maxgap = 0, minoverlap = 1, type = "any"),])
hyperInBoth <- merge(hyperInBoth, ttest_WD_Specific_anno, by.x=c("seqnames", "start", "end"), by.y=c("chr", "start", "end"), all.x=TRUE)

hypoInBoth <- as.data.frame(GR_WD_Spec_hypo[overlapsAny(query = GR_WD_Spec_hypo, subject = GR_liverDMRs_hypo, maxgap = 0, minoverlap = 1, type = "any"),])
hypoInBoth <- merge(hypoInBoth, ttest_WD_Specific_anno, by.x=c("seqnames", "start", "end"), by.y=c("chr", "start", "end"), all.x=TRUE)

# Backgrounds
pdf(file="Figures/Blood Liver Background Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_liverBackground, GR_Background), NameOfPeaks = c("Liver", "Blood"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05))
dev.off()

# All DMRs
pdf(file="Figures/Blood Liver DMRs Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_liverDMRs, GR_WD_Spec), NameOfPeaks = c("Liver", "Blood"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05))
dev.off()

# Hyper DMRs
pdf(file="Figures/Blood Liver Hyper DMRs Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_liverDMRs_hyper, GR_WD_Spec_hyper), NameOfPeaks = c("Liver", "Blood"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05))
dev.off()

# Hypo DMRs
pdf(file="Figures/Blood Liver Hypo DMRs Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_liverDMRs_hypo, GR_WD_Spec_hypo), NameOfPeaks = c("Liver", "Blood"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05))
dev.off()

# Comparison of Blood and Liver GREAT DMR Genes ####
blood_hyper <- as.character(unlist(read.delim("Tables/WD Blood Hyper DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
blood_hypo <- as.character(unlist(read.delim("Tables/WD Blood Hypo DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
blood_back <- as.character(unlist(read.delim("Tables/GREAT/WD Blood Background GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
liver_hyper <- as.character(unlist(read.delim("Tables/WD Liver Hyper DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
liver_hypo <- as.character(unlist(read.delim("Tables/WD Liver Hypo DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))
liver_back <- as.character(unlist(read.delim("Tables/GREAT/WD Liver Distal Background Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))

blood_liver_hyper <- intersect(blood_hyper, liver_hyper)
blood_liver_hypo <- intersect(blood_hypo, liver_hypo)
blood_liver_back <- intersect(blood_back, liver_back)
write.table(blood_liver_hyper, "Tables/Blood and Liver Hyper DMR Overlap Genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(blood_liver_hypo, "Tables/Blood and Liver Hypo DMR Overlap Genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(blood_liver_back, "Tables/Blood and Liver Background Overlap Genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

venn.diagram(list("Blood"=blood_hyper, "Liver"=liver_hyper), "Figures/WD Blood Liver Hyper DMR GREAT Distal Gene Overlap.png", height=8, 
             width=10, imagetype="png", units="in", fontfamily="sans", cat.fontfamily="sans", fill=c("lightblue", "lightpink"), cex=3, 
             lwd=4, cat.cex=3, cat.pos=c(180,180), main="Hypermethylated DMR Genes", main.fontfamily="sans", main.cex=3)

venn.diagram(list("Blood"=blood_hypo, "Liver"=liver_hypo), "Figures/WD Blood Liver Hypo DMR GREAT Distal Gene Overlap.png", height=8, 
             width=10, imagetype="png", units="in", fontfamily="sans", cat.fontfamily="sans", fill=c("lightblue", "lightpink"), cex=3, 
             lwd=4, cat.cex=3, cat.pos=c(180,180), main="Hypomethylated DMR Genes", main.fontfamily="sans", main.cex=3, ext.text=FALSE)

hyperlist <- list("blood"=blood_hyper, "liver"=liver_hyper)
hyperlist_gom <- newGOM(hyperlist, genome.size = length(blood_liver_back))
getMatrix(hyperlist_gom, "odds.ratio") # 6.119585
getMatrix(hyperlist_gom, "pval") #6.772944e-40 

hypolist <- list("blood"=blood_hypo, "liver"=liver_hypo)
hypolist_gom <- newGOM(hypolist, genome.size = length(blood_liver_back))
getMatrix(hypolist_gom, "odds.ratio") # 3.656413
getMatrix(hypolist_gom, "pval") #4.10701e-08 

