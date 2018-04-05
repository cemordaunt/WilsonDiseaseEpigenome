# Liver WD Early Specific DMRs ####
# Charles Mordaunt
# 3/26/18

# Packages ####
library(GenomicRanges)
library(VennDiagram)
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
library(ChIPpeakAnno)

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

# Heatmap Functions
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
                                                    legend.position = c(0.94, 0.9),
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
                theme(plot.margin = unit(c(0,-0.5,0,-0.6), "lines"),
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
                scale_color_manual(breaks = c("HC", "Early WD", "Early DC", "Male", "Female"), 
                                   values = c("HC"="#3366CC", "Early WD"="#FF3366", "Early DC"="#009933","Male"="#FFFF33", "Female"="#FF6633")) +
                scale_fill_manual(breaks = c("HC", "Early WD", "Early DC", "Male", "Female"), 
                                  values = c("HC"="#3366CC", "Early WD"="#FF3366", "Early DC"="#009933","Male"="#FFFF33", "Female"="#FF6633"))
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

# DMRs
suffix <- "_silver_DMR_info.txt"
DMRs_HCvWDE <- combineFiles(chroms=chroms, prefix="DMRs_HCvWDE/HCvWDE_", suffix = suffix)
DMRs_HCvDCE <- combineFiles(chroms=chroms, prefix="DMRs_HCvDCE/HCvDCE_", suffix = suffix)
DMRs_DCEvWDE <- combineFiles(chroms=chroms, prefix="DMRs_DCEvWDE/DCEvWDE_", suffix = suffix)

# Methylation
suffix <- "_silver_DMR_methylation.txt"
meth_HCvWDE <- combineFiles(chroms=chroms, prefix="DMRs_HCvWDE/HCvWDE_", suffix = suffix)
meth_HCvDCE <- combineFiles(chroms=chroms, prefix="DMRs_HCvDCE/HCvDCE_", suffix = suffix)
meth_DCEvWDE <- combineFiles(chroms=chroms, prefix="DMRs_DCEvWDE/DCEvWDE_", suffix = suffix)

# Coverage
suffix <- "_silver_DMR_coverage.txt"
cov_HCvWDE <- combineFiles(chroms=chroms, prefix="DMRs_HCvWDE/HCvWDE_", suffix = suffix)
cov_HCvDCE <- combineFiles(chroms=chroms, prefix="DMRs_HCvDCE/HCvDCE_", suffix = suffix)
cov_DCEvWDE <- combineFiles(chroms=chroms, prefix="DMRs_DCEvWDE/DCEvWDE_", suffix = suffix)

# Background
backCov_HCvWDE <- combineFiles(chroms=chroms, prefix="DMRs_HCvWDE/HCvWDE_", suffix="_background_DMR_coverage.txt")

# Subset HCvWDE DMRs ####
# Parameters
numCtrl <- 6
numExp <- 3
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Merge Tables
meth_HCvWDE_info <- merge(meth_HCvWDE, DMRs_HCvWDE, by = c("chr", "start", "end"))
meth_HCvWDE_info <- meth_HCvWDE_info[,c(1:3,(4+numSamples):ncol(meth_HCvWDE_info),4:(3+numSamples))]
cov_HCvWDE_info <- merge(cov_HCvWDE, DMRs_HCvWDE, by = c("chr", "start", "end"))
cov_HCvWDE_info <- cov_HCvWDE_info[,c(1:3,(4+numSamples):ncol(cov_HCvWDE_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_HCvWDE_info[,16:ncol(cov_HCvWDE_info)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample

# T-test and subset
ttest_HCvWDE <- DMRttest(meth_info=meth_HCvWDE_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_HCvWDE <- subset(ttest_HCvWDE, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #6686 -> 926
write.table(ttest_HCvWDE, "Tables/HCvWDE Blood 2 Subsetted DMRs Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Subset HCvDCE DMRs ####
# Parameters
numCtrl <- 6
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Merge Tables
meth_HCvDCE_info <- merge(meth_HCvDCE, DMRs_HCvDCE, by = c("chr", "start", "end"))
meth_HCvDCE_info <- meth_HCvDCE_info[,c(1:3,(4+numSamples):ncol(meth_HCvDCE_info),4:(3+numSamples))]
cov_HCvDCE_info <- merge(cov_HCvDCE, DMRs_HCvDCE, by = c("chr", "start", "end"))
cov_HCvDCE_info <- cov_HCvDCE_info[,c(1:3,(4+numSamples):ncol(cov_HCvDCE_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_HCvDCE_info[,16:ncol(cov_HCvDCE_info)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #9 reads per sample

# T-test and subset
ttest_HCvDCE <- DMRttest(meth_info=meth_HCvDCE_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_HCvDCE <- subset(ttest_HCvDCE, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #6193 -> 1108
write.table(ttest_HCvDCE, "Tables/HCvDCE Blood 2 Subsetted DMRs Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Subset DCEvWDE DMRs ####
# Parameters
numCtrl <- 4
numExp <- 3
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Merge Tables
meth_DCEvWDE_info <- merge(meth_DCEvWDE, DMRs_DCEvWDE, by = c("chr", "start", "end"))
meth_DCEvWDE_info <- meth_DCEvWDE_info[,c(1:3,(4+numSamples):ncol(meth_DCEvWDE_info),4:(3+numSamples))]
cov_DCEvWDE_info <- merge(cov_DCEvWDE, DMRs_DCEvWDE, by = c("chr", "start", "end"))
cov_DCEvWDE_info <- cov_DCEvWDE_info[,c(1:3,(4+numSamples):ncol(cov_DCEvWDE_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_DCEvWDE_info[,16:ncol(cov_DCEvWDE_info)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #12 reads per sample

# T-test and subset
ttest_DCEvWDE <- DMRttest(meth_info=meth_DCEvWDE_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_DCEvWDE <- subset(ttest_DCEvWDE, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #6962 -> 809
write.table(ttest_DCEvWDE, "Tables/DCEvWDE Blood 2 Subsetted DMRs Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Subset HCvWDE Background ####
# Parameters
numCtrl <- 6
numExp <- 3
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Coverage
reads <- as.matrix(backCov_HCvWDE[,4:ncol(backCov_HCvWDE)])
#table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample
backCov_HCvWDE$minReads <- apply(reads, 1, min)

# Subset
backCov_HCvWDE <- subset(backCov_HCvWDE, minReads >= coverage) #611063 -> 233900
background_HCvWDE <- backCov_HCvWDE[,c("chr", "start", "end")]
background_HCvWDE <- background_HCvWDE[order(background_HCvWDE$chr, background_HCvWDE$start),]
write.table(background_HCvWDE, "UCSC Tracks/HCvWDE Blood 2 Subsetted Background DMRs.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Overlap HCvWD, HCvDC, and DCvWD DMRs ####
# Make GRanges
GR_HCvWDE <- makeGRange(ttest=ttest_HCvWDE, direction="all")
GR_HCvWDE_hyper <- makeGRange(ttest=ttest_HCvWDE, direction="hyper")
GR_HCvWDE_hypo <- makeGRange(ttest=ttest_HCvWDE, direction="hypo")

GR_HCvDCE <- makeGRange(ttest=ttest_HCvDCE, direction="all")
GR_HCvDCE_hyper <- makeGRange(ttest=ttest_HCvDCE, direction="hyper")
GR_HCvDCE_hypo <- makeGRange(ttest=ttest_HCvDCE, direction="hypo")

GR_DCEvWDE <- makeGRange(ttest=ttest_DCEvWDE, direction="all")
GR_DCEvWDE_hyper <- makeGRange(ttest=ttest_DCEvWDE, direction="hyper")
GR_DCEvWDE_hypo <- makeGRange(ttest=ttest_DCEvWDE, direction="hypo")

# Make Venn Diagrams
pdf(file="Figures/Subsetted WD Early DMR Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_HCvWDE, GR_HCvDCE, GR_DCEvWDE), NameOfPeaks = c("HC_vs_WDE", "HC_vs_DCE", "DCE_vs_WDE"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen"), cat.pos = c(0,0,0), 
                        cat.dist = c(0.05, 0.05, 0.05))
dev.off()

pdf(file="Figures/Subsetted WD Early Hyper DMR Overlap Venn.pdf", width=10, height=9, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_HCvDCE_hyper, GR_HCvWDE_hyper, GR_DCEvWDE_hyper), NameOfPeaks = c("HC_vs_DCE_Hyper", "HC_vs_WDE_Hyper", "DCE_vs_WDE_Hyper"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 3, fill = c("lightpink", "lightblue", "lightgreen"), cat.pos = c(0,0,0), 
                        cat.dist = c(0.3, 0.2, 0.05), cat.col=c("black", "black", "black"), fontfamily="sans")
dev.off()

pdf(file="Figures/Subsetted WD Early Hypo DMR Overlap Venn.pdf", width=10, height=9, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_HCvWDE_hypo, GR_HCvDCE_hypo, GR_DCEvWDE_hypo), NameOfPeaks = c("HC_vs_WDE_Hypo", "HC_vs_DCE_Hypo", "DCE_vs_WDE_Hypo"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 180, margin = 0.02,
                        cat.cex = 2, cex = 3, fill = c("lightblue", "lightpink", "lightgreen"), cat.pos = c(0,0,0), 
                        cat.dist = c(0.05, 0.05, 0.05), fontfamily="sans", cat.col=c("black", "black", "black"))
dev.off()

# Get Specific DMRs
HCvWDE_overlaps_DCEvWDE_hyper <- overlapsAny(query = GR_HCvWDE_hyper, subject = GR_DCEvWDE_hyper, maxgap = 0, minoverlap = 1, type = "any") #189
HCvWDE_overlaps_HCvDCE_hyper <- overlapsAny(query = GR_HCvWDE_hyper, subject = GR_HCvDCE_hyper, maxgap = 0, minoverlap = 1, type = "any") #484
HCvWDE_overlaps_DCEvWDE_hypo <- overlapsAny(query = GR_HCvWDE_hypo, subject = GR_DCEvWDE_hypo, maxgap = 0, minoverlap = 1, type = "any") #75
HCvWDE_overlaps_HCvDCE_hypo <- overlapsAny(query = GR_HCvWDE_hypo, subject = GR_HCvDCE_hypo, maxgap = 0, minoverlap = 1, type = "any") #171

GR_WDE_Spec_hyper <- GR_HCvWDE_hyper[HCvWDE_overlaps_DCEvWDE_hyper & !HCvWDE_overlaps_HCvDCE_hyper,] #124 DMRs
GR_WDE_Spec_hypo <- GR_HCvWDE_hypo[HCvWDE_overlaps_DCEvWDE_hypo & !HCvWDE_overlaps_HCvDCE_hypo,] #70 DMRs
GR_WDE_Spec <- c(GR_WDE_Spec_hyper, GR_WDE_Spec_hypo) #194 DMRs
isDisjoint(GR_WDE_Spec) #TRUE (not overlapping)

# Make bed files
DMRs_WDE_Spec <- GRangetoBED(GR=GR_WDE_Spec, writeFile=TRUE, fileName="UCSC Tracks/WD_Early_Blood_2_Specific_DMRs.bed")
DMRs_WDE_Spec_hyper <- GRangetoBED(GR=GR_WDE_Spec_hyper, writeFile=TRUE, fileName="UCSC Tracks/WD_Early_Blood_2_Specific_Hyper_DMRs.bed")
DMRs_WDE_Spec_hypo <- GRangetoBED(GR=GR_WDE_Spec_hypo, writeFile=TRUE, fileName="UCSC Tracks/WD_Early_Blood_2_Specific_Hypo_DMRs.bed")

# Specific DMR Stats
ttest_WDE_Spec <- merge(ttest_HCvWDE, DMRs_WDE_Spec, by.x=c("chr", "start", "end"), by.y=c("seqnames", "start", "end"), all.x=FALSE)
write.table(ttest_WDE_Spec, "Tables/WD Early Specific DMR Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# GREAT WD Early Specific DMRs ####
# Get background file
GR_Background <- GRanges(seqnames = background_HCvWDE$chr, ranges = IRanges(start = background_HCvWDE$start, end = background_HCvWDE$end))

# Get chain file
#url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
#download.file(url, basename(url))
#gunzip(basename(url))

# Make Files for GREAT
prepGREAT(DMRs=GR_WDE_Spec, Background=GR_Background, fileName="UCSC Tracks/WD_Early_Blood_2_WD_Specific_DMRs_hg19_for_GREAT.bed", writeBack=TRUE,
          backName="UCSC Tracks/WD_Early_Blood_2_HCvWD_Background_hg19_for_GREAT.bed")
prepGREAT(DMRs=GR_WDE_Spec_hyper, Background=GR_Background, fileName="UCSC Tracks/WD_Early_Blood_2_WD_Specific_Hyper_DMRs_hg19_for_GREAT.bed")
prepGREAT(DMRs=GR_WDE_Spec_hypo, Background=GR_Background, fileName="UCSC Tracks/WD_Early_Blood_2_WD_Specific_Hypo_DMRs_hg19_for_GREAT.bed")

# Heatmap ####
# Redo with methylation included for DC samples
meth <- merge(x=DMRs_WDE_Spec, y=meth_HCvWDE, by.x=c("seqnames", "start", "end"), by.y=c("chr", "start", "end"), all.x=TRUE)
meth <- meth[,6:ncol(meth)]
methavg <- rowMeans(meth, na.rm = TRUE)
methdiff <- meth - methavg
methdiff <- as.matrix(methdiff)
colnames(methdiff) <- c(paste("HC", 1:6,sep=""), paste("WDE", 1:3, sep=""))
methplot <- ggheatmap2(x = methdiff, custom.label = colnames(methdiff),
                       hm.colours = c("#0000FF", "#0000FF", "Black", "#FF0000", "#FF0000"), my.values = c(0,0,0.5,1,1))
pdf(file="Figures/HC vs WDE Specific DMRs Heatmap.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show(methplot, col.width = 0.25, row.width = 0.12)
dev.off()

# PCA Plot ####
# Redo with methylation included for DC Samples
data <- t(as.matrix(meth))
diagnosis <- c(rep("HC", 6), rep("WD Early", 3))
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
plot(data.pca, type = "l")
summary(data.pca)
# PC1 60.4%
# PC2 8.8%

g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = diagnosis, ellipse = TRUE, circle = FALSE, 
              var.axes = FALSE, varname.abbrev = FALSE, choices = 1:2,ellipse.prob = 0.95)
g + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.8, 0.95), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), legend.text=element_text(size=22)) +
        coord_cartesian(xlim = c(-10, 16), ylim = c(-13,13)) +
        xlab("PC1 (60% of Variance)") +
        ylab("PC2 (9% of Variance)") +
        scale_color_manual(breaks = c("HC", "WD Early"), values = c("HC"="#3366CC", "WD Early"="#FF3366")) +
        geom_point(aes(color = diagnosis), size=3)
ggsave("Figures/HC vs WDE Blood 2 Specific DMRs PCA plot.png", dpi = 600, width = 8, height = 8, units = "in")

# Comparison with WD Liver Specific DMRs ####
liverDMRs <- read.delim("WD_Specific_DMRs.bed", sep="\t", header=FALSE)
colnames(liverDMRs) <- c("chr", "start", "end")
GR_liverDMRs <- GRanges(seqnames = liverDMRs$chr, ranges=IRanges(start=liverDMRs$start, end=liverDMRs$end))

liverDMRs_hyper <- read.delim("WD_Specific_hyper_DMRs.bed", sep="\t", header=FALSE)
colnames(liverDMRs_hyper) <- c("chr", "start", "end")
GR_liverDMRs_hyper <- GRanges(seqnames = liverDMRs_hyper$chr, ranges=IRanges(start=liverDMRs_hyper$start, end=liverDMRs_hyper$end))

liverDMRs_hypo <- read.delim("WD_Specific_hypo_DMRs.bed", sep="\t", header=FALSE)
colnames(liverDMRs_hypo) <- c("chr", "start", "end")
GR_liverDMRs_hypo <- GRanges(seqnames = liverDMRs_hypo$chr, ranges=IRanges(start=liverDMRs_hypo$start, end=liverDMRs_hypo$end))

liverBackground <- read.delim("Healthy_Ctrl_vs_WD_Liver_Subsetted_DMR_Background.bed", sep="\t", header=FALSE)
colnames(liverBackground) <- c("chr", "start", "end")
GR_liverBackground <- GRanges(seqnames = liverBackground$chr, ranges=IRanges(start=liverBackground$start, end=liverBackground$end))

hyperInBoth <- as.data.frame(GR_WDE_Spec_hyper[overlapsAny(query = GR_WDE_Spec_hyper, subject = GR_liverDMRs_hyper, maxgap = 0, minoverlap = 1, type = "any"),])
hyperInBoth <- merge(hyperInBoth, ttest_WDE_Spec, by.x=c("seqnames", "start", "end"), by.y=c("chr", "start", "end"), all.x=TRUE)

hypoInBoth <- as.data.frame(GR_WDE_Spec_hypo[overlapsAny(query = GR_WDE_Spec_hypo, subject = GR_liverDMRs_hypo, maxgap = 0, minoverlap = 1, type = "any"),])
hypoInBoth <- merge(hypoInBoth, ttest_WDE_Spec, by.x=c("seqnames", "start", "end"), by.y=c("chr", "start", "end"), all.x=TRUE)

backInBoth <- as.data.frame(GR_Background[overlapsAny(query = GR_Background, subject = GR_liverBackground, maxgap = 0, minoverlap = 1, type = "any"),])


# Backgrounds
pdf(file="Figures/WD and WDE Liver Background Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_liverBackground, GR_Background), NameOfPeaks = c("WD", "WDE"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightpink", "lightblue"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05))
dev.off()

# All DMRs
pdf(file="Figures/WD and WDE Liver DMRs Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_liverDMRs, GR_WDE_Spec), NameOfPeaks = c("WD", "WDE"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightpink", "lightblue"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05))
dev.off()

# Hyper DMRs
pdf(file="Figures/WD and WDE Liver Hyper DMRs Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_liverDMRs_hyper, GR_WDE_Spec_hyper), NameOfPeaks = c("WD", "WDE"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 3, cex = 3, fill = c("lightpink", "lightblue"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), fontfamily="sans", cat.fontfamily="sans", main="Hypermethylated DMRs", 
                        main.fontfamily="sans", main.cex=3, ext.text=FALSE, cat.col=c("white", "white"), totalTest=nrow(backInBoth))
dev.off()
#p 1.650889e-93

# Hypo DMRs
pdf(file="Figures/WD and WDE Liver Hypo DMRs Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_liverDMRs_hypo, GR_WDE_Spec_hypo), NameOfPeaks = c("WD", "WDE"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 0, margin = 0.02,
                        cat.cex = 3, cex = 3, fill = c("lightpink", "lightblue"), cat.pos = c(0,0), 
                        cat.dist = c(0.05, 0.05), fontfamily="sans", cat.fontfamily="sans", main="Hypomethylated DMRs", 
                        main.fontfamily="sans", main.cex=3, ext.text=FALSE, cat.col=c("white", "white"), totalTest=nrow(backInBoth))
dev.off()
# p  4.181405e-33

# Comparison of Liver WDE and WD Specific GREAT DMR Genes ####
# WDE Specific
WDE_hyper <- sort(unique(as.character(unlist(read.delim("GREAT/Liver WD Early Specific Hyper DMRs GREAT Distal Gene List.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WDE_hypo <- sort(unique(as.character(unlist(read.delim("GREAT/Liver WD Early Specific Hypo DMRs GREAT Distal Gene List.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WDE_back <- read.delim("GREAT/Liver WD Early Specific Background GREAT Distal Gene List.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
WDE_back <- sort(unique(as.character(unlist(WDE_back[,1]))))

# WD Specific
WD_hyper <- sort(unique(as.character(unlist(read.delim("GREAT/WD Liver Hyper DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WD_hypo <- sort(unique(as.character(unlist(read.delim("GREAT/WD Liver Hypo DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WD_back <- sort(unique(as.character(unlist(read.delim("GREAT/WD Liver Distal Background Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))

WDE_WD_hyper <- intersect(WDE_hyper, WD_hyper)
WDE_WD_hypo <- intersect(WDE_hypo, WD_hypo)
WDE_WD_back <- intersect(WDE_back, WD_back)
write.table(WDE_WD_hyper, "Tables/WDE and WD Hyper DMR Overlap Genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(WDE_WD_hypo, "Tables/WDE and WD Hypo DMR Overlap Genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(WDE_WD_back, "Tables/WDE and WD Background Overlap Genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

venn.diagram(list("WDE"=WDE_hyper, "WD"=WD_hyper), "Figures/WD WDE WD Hyper DMR GREAT Distal Gene Overlap.png", height=8, 
             width=10, imagetype="png", units="in", fontfamily="sans", cat.fontfamily="sans", fill=c("lightblue", "lightpink"), cex=3, 
             lwd=4, cat.cex=3, cat.pos=c(190,180), main="Hypermethylated DMR Genes", main.fontfamily="sans", main.cex=3, cat.col=c("white", "white"))

venn.diagram(list("WDE"=WDE_hypo, "WD"=WD_hypo), "Figures/WD WDE WD Hypo DMR GREAT Distal Gene Overlap.png", height=8, 
             width=10, imagetype="png", units="in", fontfamily="sans", cat.fontfamily="sans", fill=c("lightblue", "lightpink"), cex=3, 
             lwd=4, cat.cex=3, cat.pos=c(190,180), main="Hypomethylated DMR Genes", main.fontfamily="sans", main.cex=3, ext.text=FALSE, cat.col=c("white", "white"))

hyperlist <- list("WDE"=WDE_hyper, "WD"=WD_hyper)
hyperlist_gom <- newGOM(hyperlist, genome.size = length(WDE_WD_back))
getMatrix(hyperlist_gom, "odds.ratio") # 14.89887
getMatrix(hyperlist_gom, "pval") #5.926662e-72 

hypolist <- list("WDE"=WDE_hypo, "WD"=WD_hypo)
hypolist_gom <- newGOM(hypolist, genome.size = length(WDE_WD_back))
getMatrix(hypolist_gom, "odds.ratio") # 8.92636
getMatrix(hypolist_gom, "pval") #3.347435e-26 

# Methylation Analysis ####

# Load Data
meth_WDE_Spec <- combineFiles(chroms=paste("chr",1:22, sep=""), prefix="WD_Early_Specific_Liver_DMR_Methylation/WD_Early_Specific_Blood_DMR_Liver_Methyl_",
                     suffix="_DMR_methylation.txt")
meth_WDE_Spec <- meth_WDE_Spec[order(meth_WDE_Spec$chr, meth_WDE_Spec$start),]
DMRs_WDE_Spec <- DMRs_WDE_Spec[order(DMRs_WDE_Spec$seqnames, DMRs_WDE_Spec$start),]
ttest_WDE_Spec <- ttest_WDE_Spec[order(ttest_WDE_Spec$chr, ttest_WDE_Spec$start),]
table(DMRs_WDE_Spec$start == meth_WDE_Spec$start) # All True
table(ttest_WDE_Spec$start == meth_WDE_Spec$start) # All True

# Heatmap
meth_WDE_Spec <- meth_WDE_Spec[,4:ncol(meth_WDE_Spec)]
methavg <- rowMeans(meth_WDE_Spec, na.rm = TRUE)
methdiff <- meth_WDE_Spec - methavg
methdiff <- as.matrix(methdiff)
colnames(methdiff) <- c(paste("HC", 1:6,sep=""), paste("WDE", 1:3, sep=""), paste("DCE", 1:4, sep=""))
methplot <- ggheatmap2(x = methdiff, custom.label = colnames(methdiff),
                       hm.colours = c("#0000FF", "#0000FF", "Black", "#FF0000", "#FF0000"), my.values = c(0,0,0.5,1,1))
pdf(file="Figures/WD Early Specific Liver DMRs All Samples Heatmap.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show(methplot, col.width = 0.25, row.width = 0.12)
dev.off()

# Heatmap with pheno data
info <- read.csv("Tables/WGBS Liver Sample Info.csv", header=TRUE)
info$EarlyGroup[info$Sequencing.ID %in% c("VMDK1C", "VMDK1D", "VMDK4C", "VMDK5C", "VMDK6B", "VMDK6C")] <- "HC"
info$EarlyGroup[info$Sequencing.ID %in% c("VMDK4B","VMDK5B", "VMDK2B")] <- "Early WD"
info$EarlyGroup[info$Sequencing.ID %in% c("VMDK6A", "VMDK2D", "VMDK3C", "VMDK3D")] <- "Early DC"
info$EarlyGroup[is.na(info$EarlyGroup)] <- "None"
phenoData <- info[,c("Sequencing.ID", "EarlyGroup", "Sex")]
colnames(phenoData)[1] <- "Sample"
phenoData <- subset(phenoData, EarlyGroup %in% c("HC", "Early WD", "Early DC"))
phenoData <- phenoData[match(colnames(meth_WDE_Spec), phenoData$Sample),]
table(colnames(meth_WDE_Spec) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels=unique(phenoData$Sample), ordered=TRUE)
phenoData$EarlyGroup <- factor(phenoData$EarlyGroup, levels=c("HC", "Early WD", "Early DC"), ordered=TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels=c("Male", "Female"), ordered=TRUE)

methavg <- rowMeans(meth_WDE_Spec, na.rm = TRUE)
methdiff <- meth_WDE_Spec - methavg
methdiff <- as.matrix(methdiff) 
methplot <- ggheatmap2_pheno(x = methdiff, phenoData=phenoData, hm.colours = c("#0000FF", "#0000FF", "Black", "#FF0000", "#FF0000"), 
                             my.values = c(0,0,0.5,1,1), low = min(methdiff), high = max(methdiff))
pdf(file="Figures/Liver WDE Specific DMRs Heatmap phenoData.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show_pheno(methplot)
dev.off()

# PCA Plot
data <- t(as.matrix(meth_WDE_Spec))
diagnosis <- c(rep("HC", 6), rep("Early WD", 3), rep("Early DC", 4))
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
plot(data.pca, type = "l")
summary(data.pca)
# PC1 56.49%
# PC2 7.3%

g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = diagnosis, ellipse = TRUE, circle = FALSE, 
              var.axes = FALSE, varname.abbrev = FALSE, choices = 1:2,ellipse.prob = 0.95)
g + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(0.87, 0.94), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), legend.text=element_text(size=22)) +
        coord_cartesian(xlim = c(-10, 20), ylim = c(-15,15)) +
        xlab("PC1 (56% of Variance)") +
        ylab("PC2 (7% of Variance)") +
        scale_color_manual(breaks = c("HC", "Early WD", "Early DC"), 
                           values = c("HC"="#3366CC", "Early WD"="#FF3366", "Early DC"="#009933")) +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        geom_point(aes(color = diagnosis), size=3)
ggsave("Figures/WD Early Specific DMRs Liver All Samples PCA plot.png", dpi = 600, width = 8, height = 8, units = "in")

# Compare covariates to DMR methylation ####
samples <- subset(info, EarlyGroup %in% c("HC", "Early WD", "Early DC"))
samples <- samples[,c("Sequencing.ID", "Sex", "Age", "BMI", "FibrosisStage", "Steatosis", "InflammationGrade", "EarlyGroup")]
colnames(samples)[1] <- "Sample"
samples$Sample <- as.character(samples$Sample)
samples$Sex <- factor(as.character(samples$Sex), levels=c("Male", "Female"), ordered=TRUE)
samples$EarlyGroup <- factor(samples$EarlyGroup, levels=c("HC", "Early WD", "Early DC"), ordered=TRUE)

covAssoc <- matrix(nrow = dim(meth_WDE_Spec)[1], ncol = 48)
covAssoc[,1:3] <- as.matrix(DMRs_WDE_Spec[,c("seqnames", "start", "end")])
samples <- samples[match(colnames(meth_WDE_Spec), samples$Sample),]
table(colnames(meth_WDE_Spec) == as.character(samples$Sample)) #All True

for(i in 1:dim(meth_WDE_Spec)[1]){
        # FvsM
        ttest <- NULL
        ttest <- t.test(meth_WDE_Spec[i,samples$Sample[samples$Sex == "Female"]], 
                        meth_WDE_Spec[i,samples$Sample[samples$Sex == "Male"]])
        covAssoc[i,4] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,5] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,6] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,7] <- ttest$statistic #tstat
        covAssoc[i,8] <- ttest$p.value #pvalue
        
        # Age
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_WDE_Spec[i,]) ~ samples$Age)
        sum <- summary(lm)
        covAssoc[i,9] <- sum$coefficients['samples$Age','Estimate'] #Effect
        covAssoc[i,10] <- sum$coefficients['samples$Age','Std. Error'] #Std Error
        covAssoc[i,11] <- sum$r.squared #Rsquared
        covAssoc[i,12] <- sum$coefficients['samples$Age','t value'] #tstat
        covAssoc[i,13] <- sum$coefficients['samples$Age', 'Pr(>|t|)'] #pvalue
        
        # BMI
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_WDE_Spec[i,]) ~ samples$BMI)
        sum <- summary(lm)
        covAssoc[i,14] <- sum$coefficients['samples$BMI','Estimate'] #Effect
        covAssoc[i,15] <- sum$coefficients['samples$BMI','Std. Error'] #Std Error
        covAssoc[i,16] <- sum$r.squared #Rsquared
        covAssoc[i,17] <- sum$coefficients['samples$BMI','t value'] #tstat
        covAssoc[i,18] <- sum$coefficients['samples$BMI', 'Pr(>|t|)'] #pvalue
        
        # FibrosisStage
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_WDE_Spec[i,]) ~ samples$FibrosisStage)
        sum <- summary(lm)
        covAssoc[i,19] <- sum$coefficients['samples$FibrosisStage','Estimate'] #Effect
        covAssoc[i,20] <- sum$coefficients['samples$FibrosisStage','Std. Error'] #Std Error
        covAssoc[i,21] <- sum$r.squared #Rsquared
        covAssoc[i,22] <- sum$coefficients['samples$FibrosisStage','t value'] #tstat
        covAssoc[i,23] <- sum$coefficients['samples$FibrosisStage', 'Pr(>|t|)'] #pvalue
        
        # Steatosis
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_WDE_Spec[i,]) ~ samples$Steatosis)
        sum <- summary(lm)
        covAssoc[i,24] <- sum$coefficients['samples$Steatosis','Estimate'] #Effect
        covAssoc[i,25] <- sum$coefficients['samples$Steatosis','Std. Error'] #Std Error
        covAssoc[i,26] <- sum$r.squared #Rsquared
        covAssoc[i,27] <- sum$coefficients['samples$Steatosis','t value'] #tstat
        covAssoc[i,28] <- sum$coefficients['samples$Steatosis', 'Pr(>|t|)'] #pvalue
        
        # InflammationGrade
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(meth_WDE_Spec[i,]) ~ samples$InflammationGrade)
        sum <- summary(lm)
        covAssoc[i,29] <- sum$coefficients['samples$InflammationGrade','Estimate'] #Effect
        covAssoc[i,30] <- sum$coefficients['samples$InflammationGrade','Std. Error'] #Std Error
        covAssoc[i,31] <- sum$r.squared #Rsquared
        covAssoc[i,32] <- sum$coefficients['samples$InflammationGrade','t value'] #tstat
        covAssoc[i,33] <- sum$coefficients['samples$InflammationGrade', 'Pr(>|t|)'] #pvalue
        
        # Early WD vs HC
        ttest <- NULL
        ttest <- t.test(meth_WDE_Spec[i,samples$Sample[samples$EarlyGroup == "Early WD"]], 
                        meth_WDE_Spec[i,samples$Sample[samples$EarlyGroup == "HC"]])
        covAssoc[i,34] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,35] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,36] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,37] <- ttest$statistic #tstat
        covAssoc[i,38] <- ttest$p.value #pvalue
        
        # Early WD vs Early DC
        ttest <- NULL
        ttest <- t.test(meth_WDE_Spec[i,samples$Sample[samples$EarlyGroup == "Early WD"]], 
                        meth_WDE_Spec[i,samples$Sample[samples$EarlyGroup == "Early DC"]])
        covAssoc[i,39] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,40] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,41] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,42] <- ttest$statistic #tstat
        covAssoc[i,43] <- ttest$p.value #pvalue
        
        # Early DC vs HC
        ttest <- NULL
        ttest <- t.test(meth_WDE_Spec[i,samples$Sample[samples$EarlyGroup == "Early DC"]], 
                        meth_WDE_Spec[i,samples$Sample[samples$EarlyGroup == "HC"]])
        covAssoc[i,44] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,45] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,46] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,47] <- ttest$statistic #tstat
        covAssoc[i,48] <- ttest$p.value #pvalue
        
}

covAssoc <- as.data.frame(covAssoc, stringsAsFactors = FALSE)
newCols <- c("chr", "start", "end", paste(rep(c("FvsM"), each=5), c("meanDiff", "confIntL", "confIntR", "tstat", "pvalue"), sep="_"),
             paste(rep(c("Age", "BMI", "FibrosisStage", "Steatosis", "InflammationGrade"), each=5), c("Effect", "StdErr", "Rsquared", "tstat", "pvalue"), sep="_"),
             paste(rep(c("WDEvsHC", "WDEvsDCE", "DCEvsHC"), each=5),c("meanDiff", "confIntL", "confIntR", "tstat", "pvalue"), sep="_"))
colnames(covAssoc) <- newCols
covAssoc$start <- as.integer(covAssoc$start)
covAssoc$end <- as.integer(covAssoc$end)
covAssoc[,4:dim(covAssoc)[2]] <- apply(covAssoc[,4:dim(covAssoc)[2]], 2, function(x) as.numeric(x))

# Adj p-values
covAssoc$FvsM_pAdj <- p.adjust(covAssoc$FvsM_pvalue, "bonf")
covAssoc$Age_pAdj <- p.adjust(covAssoc$Age_pvalue, "bonf")
covAssoc$BMI_pAdj <- p.adjust(covAssoc$BMI_pvalue, "bonf")
covAssoc$FibrosisStage_pAdj <- p.adjust(covAssoc$FibrosisStage_pvalue, "bonf")
covAssoc$Steatosis_pAdj <- p.adjust(covAssoc$Steatosis_pvalue, "bonf")
covAssoc$InflammationGrade_pAdj <- p.adjust(covAssoc$InflammationGrade_pvalue, "bonf")
covAssoc$WDEvsHC_pAdj <- p.adjust(covAssoc$WDEvsHC_pvalue, "bonf")
covAssoc$WDEvsDCE_pAdj <- p.adjust(covAssoc$WDEvsDCE_pvalue, "bonf")
covAssoc$DCEvsHC_pAdj <- p.adjust(covAssoc$DCEvsHC_pvalue, "bonf")

write.table(covAssoc, file="Tables/Liver WDE Specific DMRs Covariate Associations.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Bonf Significant DMRs
table(covAssoc$FvsM_pAdj < 0.05)              # 0
table(covAssoc$Age_pAdj < 0.05)               # 0
table(covAssoc$BMI_pAdj < 0.05)               # 0
table(covAssoc$FibrosisStage_pAdj < 0.05)     # 0
table(covAssoc$Steatosis_pAdj < 0.05)         # 0
table(covAssoc$InflammationGrade_pAdj < 0.05) # 0
table(covAssoc$WDEvsHC_pAdj < 0.05)           # 7 *
table(covAssoc$WDEvsDCE_pAdj < 0.05)          # 1 *
table(covAssoc$DCEvsHC_pAdj < 0.05)           # 0

# Heatmap Sorted by WDEvsHC logp-value
pvals <- covAssoc[,paste(c("FvsM", "Age", "BMI", "FibrosisStage", "Steatosis", "InflammationGrade", "WDEvsHC", "WDEvsDCE", "DCEvsHC"), "pvalue", sep="_")]
logPvals <- -log10(pvals)
logPvals <- as.matrix(logPvals)
newCols <- c("FvsM", "Age", "BMI", "Fibrosis\nStage", "Steatosis", "Inflammation\nGrade", "WDEvsHC", "WDEvsDCE", "DCEvsHC")
colnames(logPvals) <- newCols
logPvals_sort <- logPvals[order(logPvals[,"WDEvsHC"]),]
plot <- ggheatmap2_sort(x = logPvals_sort, custom.label = newCols, hm.colours = c("Black", "#FF0000"), my.values = c(0,1))
pdf(file="Figures/WD Early Specific Liver DMRs Covariate Association Heatmap sorted.pdf", width=6, height=7, onefile = FALSE)
ggheatmap.show_sort(plot, col.width = 0.25, row.width = 0.1)
dev.off()
