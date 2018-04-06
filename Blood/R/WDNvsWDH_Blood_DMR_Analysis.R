# WD Blood WDN vs WDH DMR Analysis ####
# Charles Mordaunt
# 3/21/18

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
                theme(plot.margin = unit(c(0,2,0.5,1), "lines"))
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

# Heatmap with Pheno Data Functions
mydplot_pheno <- function(ddata, row=!col, col=!row, labels=col) {
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
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = "black"))
        }
        return(p)
}

g_legend_pheno <-function(a.gplot){
        ## from
        ## http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

ggheatmap.show_pheno <- function(L, widths=c(0.02,0.85,0.11,0.02), heights=c(0.02,0.2,0.06,0.7,0.02)){
        grid.newpage()
        top.layout <- grid.layout(nrow = 5, ncol = 4,
                                  widths = unit(widths, "null"),
                                  heights = unit(heights, "null"))
        pushViewport(viewport(layout=top.layout))
        print(L$col, vp=viewport(layout.pos.col=2, layout.pos.row=2))
        print(L$row, vp=viewport(layout.pos.col=3, layout.pos.row=4))
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
              vp=viewport(layout.pos.col=2, layout.pos.row=4))
        
        print(L$phenoData +
                      theme_bw(base_size = 24) +
                      theme(panel.grid.major = element_blank(), panel.border = element_blank(),
                            legend.key = element_blank(), legend.key.size = unit(1, "lines"),
                            panel.grid.minor = element_blank(), legend.position = "none", 
                            legend.background = element_blank(), legend.text = element_text(size = 12, color = "Black"),
                            plot.margin = unit(c(0,0,0,-0.45), "lines"), 
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.title = element_blank(), 
                            legend.title = element_blank(),
                            plot.title = element_blank()), vp=viewport(layout.pos.col=2, layout.pos.row=3))
        
        ## add heatmap legend
        legend <- g_legend_pheno(L$centre +
                                         theme(legend.title = element_blank(), 
                                               legend.text = element_text(size = 15),
                                               legend.background = element_blank(),
                                               legend.position = c(0.57, -0.48)))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=2))
        grid.draw(legend)
        upViewport(0)
        
        ## add pheno legend
        phenoLegend <- g_legend_pheno(L$phenoData +
                                              theme(legend.title = element_blank(), 
                                                    legend.text = element_text(size = 15),
                                                    legend.direction = "vertical",
                                                    legend.position = c(0.94, 0.89),
                                                    legend.background = element_blank()))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=3))
        grid.draw(phenoLegend)
        upViewport(0)
}

ggheatmap2_pheno <- function(x, phenoData, hm.colours=my.colours, my.values, low, high) {
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
        col.plot <- mydplot_pheno(col.dendro, col=TRUE, labels=FALSE) +
                theme(plot.margin = unit(c(0,-2.1,0,-2.2), "lines"),
                      axis.text.x = element_blank())
        row.plot <- mydplot_pheno(row.dendro, row=TRUE, labels=FALSE) +
                theme(plot.margin = unit(c(0,2,0,0), "lines"))
        
        ## order of the dendros
        col.ord <- match(col.dendro$labels$label, colnames(x))
        row.ord <- match(row.dendro$labels$label, rownames(x))
        xx <- x[row.ord,col.ord]
        dimnames(xx) <- NULL
        xx <- melt(xx)
        
        # Heatmap
        centre.plot <- ggplot(xx, aes(X2,X1)) + 
                geom_tile(aes(fill=value, colour=value)) +
                scale_fill_gradientn(colours = hm.colours, values = my.values, limits = c(low, high), na.value = "black") +
                scale_colour_gradientn(colours = hm.colours, values = my.values, limits = c(low, high), na.value = "black") +
                labs(x = NULL, y = NULL) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0),breaks = NULL) +
                theme(plot.margin = unit(rep(0, 4), "lines"))
        
        # phenoData
        sample.ord <- match(col.dendro$labels$label, as.character(phenoData$Sample))
        phenoData$Sample <- factor(as.character(phenoData$Sample), levels = as.character(phenoData$Sample)[sample.ord], ordered = TRUE)
        phenoData_m <- melt(phenoData, id.vars = "Sample")
        phenoData_m$variable <- factor(phenoData_m$variable, levels = rev(unique(phenoData_m$variable)), ordered = TRUE)
        phenoData.plot <- ggplot(phenoData_m, aes(Sample, variable)) +
                geom_tile(aes(fill=value, color=value)) +
                scale_x_discrete(expand=c(0,0)) +
                scale_y_discrete(expand=c(0,0)) +
                scale_color_manual(breaks = c("WDN", "WDH", "Male", "Female"), 
                                   values = c("WDN"="#3366CC", "WDH"="#FF3366", "Male"="#FFFF33", "Female"="#FF6633")) +
                scale_fill_manual(breaks = c("WDN", "WDH", "Male", "Female"), 
                                  values = c("WDN"="#3366CC", "WDH"="#FF3366", "Male"="#FFFF33", "Female"="#FF6633"))
        ret <- list(col=col.plot,row=row.plot,centre=centre.plot, phenoData=phenoData.plot)
        invisible(ret)
}

# Heatmap with Sorted DMRs Functions
mydplot_sort <- function(ddata, row=!col, col=!row, labels=col) {
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

g_legend_sort <-function(a.gplot){
        ## from
        ## http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

ggheatmap.show_sort <- function(L, col.width=0.2, row.width=0.2) {
        grid.newpage()
        top.layout <- grid.layout(nrow = 4, ncol = 4,
                                  widths = unit(c(0.02,0.96-row.width,row.width,0.02), "null"),
                                  heights = unit(c(0.02,col.width,0.96-col.width,0.02), "null"))
        pushViewport(viewport(layout=top.layout))
        if(col.width>0)
                print(L$col, vp=viewport(layout.pos.col=2, layout.pos.row=2))
        #if(row.width>0)
        #        print(L$row, vp=viewport(layout.pos.col=3, layout.pos.row=3))
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
        legend <- g_legend_sort(L$centre +
                                        theme(legend.title = element_blank(), 
                                              legend.text = element_text(size = 15)))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=2))
        grid.draw(legend)
        upViewport(0)
}

ggheatmap2_sort <- function(x, custom.label,
                            hm.colours=my.colours, my.values) {
        if(is.null(colnames(x)))
                colnames(x) <- sprintf("col%s",1:ncol(x))
        if(is.null(rownames(x)))
                rownames(x) <- sprintf("row%s",1:nrow(x))
        ## plot a heatmap
        ## x is an expression matrix
        #row.hc <- hclust(dist(x), "ward.D")
        col.hc <- hclust(dist(t(x)), "ward.D")
        #row.dendro <- dendro_data(as.dendrogram(row.hc),type="rectangle")
        col.dendro <- dendro_data(as.dendrogram(col.hc),type="rectangle")
        
        ## dendro plots
        col.plot <- mydplot_sort(col.dendro, col=TRUE, labels=FALSE) +
                scale_x_continuous(breaks = 1:ncol(x), labels = custom.label[col.hc$order], expand=c(0,0)) +
                theme(plot.margin = unit(c(0,2,0.5,1.5), "lines"))
        #row.plot <- mydplot_sort(row.dendro, row=TRUE, labels=FALSE) +
        #        theme(plot.margin = unit(rep(0, 4), "lines"))
        
        ## order of the dendros
        col.ord <- match(col.dendro$labels$label, colnames(x))
        row.ord <- 1:nrow(x)       #match(row.dendro$labels$label, rownames(x))
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
        ret <- list(col=col.plot,centre=centre.plot)
        invisible(ret)
}

# Combine Files ####
chroms <- c(paste("chr",1:22,sep=""), "chrM")
prefix <- "DMRs/DMRs_WDNvWDH/WDNvWDH_"
DMRs <- combineFiles(chroms=chroms, prefix=prefix, suffix ="_silver_DMR_info.txt")
meth <- combineFiles(chroms=chroms, prefix=prefix, suffix = "_silver_DMR_methylation.txt")
cov <- combineFiles(chroms=chroms, prefix=prefix, suffix = "_silver_DMR_coverage.txt")
backCov <- combineFiles(chroms=chroms, prefix=prefix, suffix="_background_DMR_coverage.txt")

# Subset DMRs ####
# Parameters
numCtrl <- 20
numExp <- 20
numSamples <- numCtrl + numExp
minDiff <- 0.05

# Merge Tables
meth_info <- merge(meth, DMRs, by = c("chr", "start", "end"))
meth_info <- meth_info[,c(1:3,(4+numSamples):ncol(meth_info),4:(3+numSamples))]
cov_info <- merge(cov, DMRs, by = c("chr", "start", "end"))
cov_info <- cov_info[,c(1:3,(4+numSamples):ncol(cov_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_info[,16:ncol(cov_info)])
table(is.na(reads)) # No NA values
if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #4 reads per sample

# T-test and subset
ttest <- DMRttest(meth_info=meth_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest <- subset(ttest, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #9712 -> 2860
write.table(ttest, "Tables/WD Blood 2 Subsetted WDNvWDH DMRs Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Subset Background ####
# Parameters
numCtrl <- 20
numExp <- 20
numSamples <- numCtrl + numExp
minDiff <- 0.05

# Coverage
reads <- as.matrix(backCov[,4:ncol(backCov)])
#table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #4 reads per sample
backCov$minReads <- apply(reads, 1, min)

# Subset
backCov <- subset(backCov, minReads >= coverage) #625588 -> 341892
background <- backCov[,c("chr", "start", "end")]
background <- background[order(background$chr, background$start),]
write.table(background, "UCSC Tracks/WD Blood 2 WDNvWDH Subsetted Background DMRs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Write Files ####
# Make GRanges
GR_DMRs_all <- makeGRange(ttest=ttest, direction="all")
GR_DMRs_hyper <- makeGRange(ttest=ttest, direction="hyper")
GR_DMRs_hypo <- makeGRange(ttest=ttest, direction="hypo")

# Make bed files
DMRs_all <- GRangetoBED(GR=GR_DMRs_all, writeFile=TRUE, fileName="UCSC Tracks/WD_Blood_2_WDNvWDH_DMRs.bed")
DMRs_hyper <- GRangetoBED(GR=GR_DMRs_hyper, writeFile=TRUE, fileName="UCSC Tracks/WD_Blood_2_WDNvWDH_Hyper_DMRs.bed")
DMRs_hypo <- GRangetoBED(GR=GR_DMRs_hypo, writeFile=TRUE, fileName="UCSC Tracks/WD_Blood_2_WDNvWDH_Hypo_DMRs.bed")

# GREAT File Prep ####
# Get background file
GR_Background <- GRanges(seqnames = background$chr, ranges = IRanges(start = background$start, end = background$end))

# Get chain file
#url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
#download.file(url, basename(url))
#gunzip(basename(url))

# Make Files for GREAT
prepGREAT(DMRs=GR_DMRs_all, Background=GR_Background, fileName="UCSC Tracks/WD_Blood_2_WDNvWDH_DMRs_hg19_for_GREAT.bed", writeBack=TRUE,
          backName="UCSC Tracks/WD_Blood_2_WDNvWDH_Background_hg19_for_GREAT.bed")
prepGREAT(DMRs=GR_DMRs_hyper, Background=GR_Background, fileName="UCSC Tracks/WD_Blood_2_WDNvWDH_Hyper_DMRs_hg19_for_GREAT.bed")
prepGREAT(DMRs=GR_DMRs_hypo, Background=GR_Background, fileName="UCSC Tracks/WD_Blood_2_WDNvWDH_Hypo_DMRs_hg19_for_GREAT.bed")

# Heatmap ####
meth_heat <- merge(x=DMRs_all, y=meth, by.x=c("seqnames", "start", "end"), by.y=c("chr", "start", "end"), all.x=TRUE)
meth_heat <- meth_heat[,6:ncol(meth_heat)]
methavg <- rowMeans(meth_heat, na.rm = TRUE)
methdiff <- meth_heat - methavg
methdiff <- as.matrix(methdiff)
colnames(methdiff) <- c(paste("WDN", 1:20,sep=""), paste("WDH", 1:20, sep=""))
methplot <- ggheatmap2(x = methdiff, custom.label = colnames(methdiff),
                       hm.colours = c("#0000FF", "#0000FF", "Black", "#FF0000", "#FF0000"), my.values = c(0,0,0.5,1,1))
pdf(file="Figures/WDNvWDH DMRs Heatmap.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show(methplot, col.width = 0.25, row.width = 0.12)
dev.off()

# Heatmap with pheno data
samples <- read.csv("Sample Info/WD Blood 2 Sample Info.csv", header=TRUE, stringsAsFactors = FALSE)
phenoData <- samples[,c("Sequencing.ID", "Diagnosis", "Sex", "Neurologic.form", "Hepatic.form")]
colnames(phenoData)[1] <- "Sample"
phenoData <- phenoData[match(colnames(meth_heat), phenoData$Sample),]
table(colnames(meth_heat) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels=unique(phenoData$Sample), ordered=TRUE)
phenoData$Diagnosis[phenoData$Diagnosis == "Healthy Control"] <- "HC"
phenoData$Diagnosis[phenoData$Neurologic.form==1] <- "WDN"
phenoData$Diagnosis[phenoData$Hepatic.form==1] <- "WDH"
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels=c("WDN", "WDH"), ordered=TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels=c("Male", "Female"), ordered=TRUE)
phenoData <- phenoData[,c("Sample", "Diagnosis", "Sex")]

methavg <- rowMeans(meth_heat, na.rm = TRUE)
methdiff <- meth_heat - methavg
methdiff <- as.matrix(methdiff) 
methplot <- ggheatmap2_pheno(x = methdiff, phenoData=phenoData, hm.colours = c("#0000FF", "#0000FF", "Black", "#FF0000", "#FF0000"), 
                             my.values = c(0,0,0.5,1,1), low = -0.5, high = 0.5)
pdf(file="Figures/WDNvWDH Blood DMRs Heatmap phenoData.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show_pheno(methplot)
dev.off()

# PCA Plot ####
data <- t(as.matrix(meth_heat))
diagnosis <- c(rep("WDN", 20), rep("WDH", 20))
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
plot(data.pca, type = "l")
summary(data.pca)
# PC1 15.9%
# PC2 4.8%

g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = diagnosis, ellipse = TRUE, circle = FALSE, 
              var.axes = FALSE, varname.abbrev = FALSE, choices = 1:2,ellipse.prob = 0.95)
g + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.82, 0.95), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank()) +
        coord_cartesian(xlim = c(-35, 25), ylim = c(-30,30)) +
        xlab("PC1 (16% of Variance)") +
        ylab("PC2 (5% of Variance)") +
        scale_color_manual(breaks = c("WDN", "WDH"), values = c("WDN"="#3366CC", "WDH"="#FF3366")) +
        geom_point(aes(color = diagnosis), size=3)
ggsave("Figures/WDNvWDH Blood 2 DMRs PCA plot.png", dpi = 600, width = 8, height = 8, units = "in")

# Compare covariates to DMR methylation ####
samples$Group <- factor(samples$Group, levels=c("Healthy Control", "Disease Control", "Wilson's Disease"), ordered=TRUE)
samples$Diagnosis <- factor(samples$Diagnosis, levels=c("Healthy Control", "NAFLD", "PSC", "Wilson's Disease"), ordered=TRUE)
samples$Sex <- factor(samples$Sex, levels=c("Male", "Female"), ordered=TRUE)
samples <- samples[,c("Sequencing.ID", "Group", "Diagnosis", "Sex", "Age", "Age.at.onset", "Height..cm.", "Weight..kg.", "BMI", 
                      "Neurologic.form", "Hepatic.form")]
colnames(samples) <- c("Sample", "Group", "Diagnosis", "Sex", "Age", "OnsetAge", "Height", "Weight", "BMI", "NeuroWD", "HepaticWD")
covAssoc <- matrix(nrow = dim(meth_heat)[1], ncol = 38)
covAssoc[,1:3] <- as.matrix(DMRs_all[,c("seqnames", "start", "end")])
samples <- samples[match(colnames(meth_heat), samples$Sample),]
table(colnames(meth_heat) == as.character(samples$Sample)) #All True

for(i in 1:dim(meth_heat)[1]){
        # FvsM
        ttest <- NULL
        ttest <- t.test(meth_heat[i,samples$Sample[samples$Sex == "Female"]], 
                        meth_heat[i,samples$Sample[samples$Sex == "Male"]])
        covAssoc[i,4] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,5] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,6] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,7] <- ttest$statistic #tstat
        covAssoc[i,8] <- ttest$p.value #pvalue
        
        # Age
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_heat[i,]) ~ samples$Age)
        sum <- summary(lm)
        covAssoc[i,9] <- sum$coefficients['samples$Age','Estimate'] #Effect
        covAssoc[i,10] <- sum$coefficients['samples$Age','Std. Error'] #Std Error
        covAssoc[i,11] <- sum$r.squared #Rsquared
        covAssoc[i,12] <- sum$coefficients['samples$Age','t value'] #tstat
        covAssoc[i,13] <- sum$coefficients['samples$Age', 'Pr(>|t|)'] #pvalue
        
        # OnsetAge
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_heat[i,]) ~ samples$OnsetAge)
        sum <- summary(lm)
        covAssoc[i,14] <- sum$coefficients['samples$OnsetAge','Estimate'] #Effect
        covAssoc[i,15] <- sum$coefficients['samples$OnsetAge','Std. Error'] #Std Error
        covAssoc[i,16] <- sum$r.squared #Rsquared
        covAssoc[i,17] <- sum$coefficients['samples$OnsetAge','t value'] #tstat
        covAssoc[i,18] <- sum$coefficients['samples$OnsetAge', 'Pr(>|t|)'] #pvalue
        
        # Height
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_heat[i,]) ~ samples$Height)
        sum <- summary(lm)
        covAssoc[i,19] <- sum$coefficients['samples$Height','Estimate'] #Effect
        covAssoc[i,20] <- sum$coefficients['samples$Height','Std. Error'] #Std Error
        covAssoc[i,21] <- sum$r.squared #Rsquared
        covAssoc[i,22] <- sum$coefficients['samples$Height','t value'] #tstat
        covAssoc[i,23] <- sum$coefficients['samples$Height', 'Pr(>|t|)'] #pvalue
        
        # Weight
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_heat[i,]) ~ samples$Weight)
        sum <- summary(lm)
        covAssoc[i,24] <- sum$coefficients['samples$Weight','Estimate'] #Effect
        covAssoc[i,25] <- sum$coefficients['samples$Weight','Std. Error'] #Std Error
        covAssoc[i,26] <- sum$r.squared #Rsquared
        covAssoc[i,27] <- sum$coefficients['samples$Weight','t value'] #tstat
        covAssoc[i,28] <- sum$coefficients['samples$Weight', 'Pr(>|t|)'] #pvalue
        
        # BMI
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_heat[i,]) ~ samples$BMI)
        sum <- summary(lm)
        covAssoc[i,29] <- sum$coefficients['samples$BMI','Estimate'] #Effect
        covAssoc[i,30] <- sum$coefficients['samples$BMI','Std. Error'] #Std Error
        covAssoc[i,31] <- sum$r.squared #Rsquared
        covAssoc[i,32] <- sum$coefficients['samples$BMI','t value'] #tstat
        covAssoc[i,33] <- sum$coefficients['samples$BMI', 'Pr(>|t|)'] #pvalue
        
        # Neuro vs Hepatic
        ttest <- NULL
        ttest <- t.test(meth_heat[i,samples$Sample[samples$NeuroWD == 1 & !is.na(samples$NeuroWD)]], 
                        meth_heat[i,samples$Sample[samples$HepaticWD == 1 & !is.na(samples$HepaticWD)]])
        covAssoc[i,34] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,35] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,36] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,37] <- ttest$statistic #tstat
        covAssoc[i,38] <- ttest$p.value #pvalue
}

covAssoc <- as.data.frame(covAssoc, stringsAsFactors = FALSE)
newCols <- c("chr", "start", "end", paste(rep(c("FvsM"), each=5), c("meanDiff", "confIntL", "confIntR", "tstat", "pvalue"), sep="_"),
             paste(rep(c("Age", "OnsetAge", "Height", "Weight", "BMI"), each=5), c("Effect", "StdErr", "Rsquared", "tstat", "pvalue"), sep="_"),
             paste(rep(c("WDNvsWDH"), each=5),c("meanDiff", "confIntL", "confIntR", "tstat", "pvalue"), sep="_"))
colnames(covAssoc) <- newCols
covAssoc$start <- as.integer(covAssoc$start)
covAssoc$end <- as.integer(covAssoc$end)
covAssoc[,4:dim(covAssoc)[2]] <- apply(covAssoc[,4:dim(covAssoc)[2]], 2, function(x) as.numeric(x))

# Adj p-values
covAssoc$FvsM_pAdj <- p.adjust(covAssoc$FvsM_pvalue, "bonf")
covAssoc$Age_pAdj <- p.adjust(covAssoc$Age_pvalue, "bonf")
covAssoc$OnsetAge_pAdj <- p.adjust(covAssoc$OnsetAge_pvalue, "bonf")
covAssoc$Height_pAdj <- p.adjust(covAssoc$Height_pvalue, "bonf")
covAssoc$Weight_pAdj <- p.adjust(covAssoc$Weight_pvalue, "bonf")
covAssoc$BMI_pAdj <- p.adjust(covAssoc$BMI_pvalue, "bonf")
covAssoc$WDNvsWDH_pAdj <- p.adjust(covAssoc$WDNvsWDH_pvalue, "bonf")
write.table(covAssoc, file="Tables/WDNvWDH Blood DMRs Covariate Associations.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Bonf Significant DMRs
table(covAssoc$FvsM_pAdj < 0.05)        # 0
table(covAssoc$Age_pAdj < 0.05)         # 0
table(covAssoc$OnsetAge_pAdj < 0.05)    # 0
table(covAssoc$Height_pAdj < 0.05)      # 0
table(covAssoc$Weight_pAdj < 0.05)      # 0
table(covAssoc$BMI_pAdj < 0.05)         # 0
table(covAssoc$WDNvsWDH_pAdj < 0.05)    # 3

# Heatmap Sorted by WDNvWDH logp-value ####
pvals <- covAssoc[,c("FvsM_pvalue", "Age_pvalue", "OnsetAge_pvalue", "Height_pvalue", "Weight_pvalue", "BMI_pvalue", "WDNvsWDH_pvalue")]
logPvals <- -log10(pvals)
logPvals <- as.matrix(logPvals)
newCols <- c("FvsM", "Age", "OnsetAge", "Height", "Weight", "BMI", "WDNvsWDH")
colnames(logPvals) <- newCols
logPvals_sort <- logPvals[order(logPvals[,"WDNvsWDH"]),]
plot <- ggheatmap2_sort(x = logPvals_sort, custom.label = newCols, hm.colours = c("Black", "#FF0000"), my.values = c(0,1))
pdf(file="Figures/WDNvWDH Blood DMRs Covariate Association Heatmap sorted.pdf", width=6, height=7, onefile = FALSE)
ggheatmap.show_sort(plot, col.width = 0.25, row.width = 0.1)
dev.off()
