# WD Specific Liver DMR Analysis ####
# Charles Mordaunt

# Packages ########################
library(ggdendro)
library(ggplot2)
library(reshape2)
library(reshape)
library(grid)
library(scales)
library(devtools)
library(plyr)
library(ggbiplot)
library(sm)
library(GenomicRanges)
library(VennDiagram)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(WebGestaltR)

# Functions ####################
write.bed = function(df,file,header="",name="name"){
        #Generalized write bed function
        # Example Output:
        #  chr1    10469   10470   0.00-1  0       +       0       0       0,0,0
        #Area stat as name
        #color for direction (hyper = red, hypo = blue)
        colnum = 3
        if(name %in% colnames(df)) {namedat = round(df[,name], 3);colnum = 4} else{namedat = rep(".",length(df[,1]))}
        if("score" %in% colnames(df)) {score = df$score;colnum = 5} else{score = rep("0",length(df[,1]))}
        if("strand" %in% colnames(df)) {strand = df$strand;colnum = 6} else{strand = rep("+",length(df[,1]))}
        if("thickStart" %in% colnames(df)) {thickStart = df$thickStart;colnum = 7} else{thickStart = rep("0",length(df[,1]))}
        if("thickEnd" %in% colnames(df)) {thickEnd = df$thickEnd;colnum = 8} else{thickEnd = rep("0",length(df[,1]))}
        if("itemRgb" %in% colnames(df)) {itemRgb = df$itemRgb;colnum = 9} else{itemRgb = rep("0,0,0",length(df[,1]))}
        if("blockCount" %in% colnames(df)) {blockCount = df$blockCount;colnum = 10} else{blockCount = rep("0",length(df[,1]))}
        if("blockSizes" %in% colnames(df)) {blockSizes = df$blockSizes;colnum = 11} else{blockSizes = rep("0",length(df[,1]))}
        if("blockStarts" %in% colnames(df)) {blockStarts = df$blockStarts;colnum = 12} else{blockStarts = rep("0",length(df[,1]))}
        
        subdf = cbind(df$chr,df$start,df$end,namedat,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,blockStarts)
        
        if(header == ""){
                write.table(subdf[,1:colnum], file=file, quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE)
        }else{
                write(header,file=file)
                write.table(subdf[,1:colnum], file=file, quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE,append=TRUE)
        }
        cat(paste("Finished writing bedfile",file,"with",colnum,"columns\n"))
}
write.dmrs_bed = function(df,file,trackname,genome){
        #Area stat as name
        #color for direction (hyper = red, hypo = blue)
        df$itemRgb = ifelse(df$direction == "hypo","0,0,255","255,0,0")
        headerline = paste("track name=",trackname," description=",trackname," useScore=0 itemRgb=On genome=",genome,sep="")
        write.bed(df,file,header=headerline,name="areaStat")
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
                scale_x_continuous(breaks = 1:ncol(x), labels = custom.label[col.hc$order]) +
                theme(plot.margin = unit(c(0,-2,0.5,-2), "lines"))
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

# Combining DMR Files ######################################

# Combining DMRs into one table (Healthy Ctrl vs WD, Reclassification, Chr X, Y excluded)
chroms_HCvWD = c(paste("chr",1:22,sep=""), "chrM")
DMRs_HCvWD <- NULL
for(i in 1:length(chroms_HCvWD)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_healthy_v_WD_DMRs/Reclass_healthy_v_WD_",chroms_HCvWD[i],
                                 "_silver_DMR_info.txt", sep=""), header = TRUE, sep = "\t")
        DMRs_HCvWD <- rbind(DMRs_HCvWD, temp)
}
DMRs_HCvWD$chr <- as.character(DMRs_HCvWD$chr)

# Combining DMRs into one table (Healthy Ctrl vs Disease Ctrl, Reclassification, Chr X, Y excluded)
chroms_HCvDC = c(paste("chr",1:22,sep=""), "chrM")
DMRs_HCvDC <- NULL
for(i in 1:length(chroms_HCvDC)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_healthyCon_v_diseaseCon_DMRs/Reclass_healthyCon_v_diseaseCon_",chroms_HCvDC[i],
                                 "_silver_DMR_info.txt", sep=""), header = TRUE, sep = "\t")
        DMRs_HCvDC <- rbind(DMRs_HCvDC, temp)
}
DMRs_HCvDC$chr <- as.character(DMRs_HCvDC$chr)

# Combining DMRs into one table (Disease Ctrl vs WD, Reclassification, Chr X, Y excluded)
chroms_DCvWD = c(paste("chr",1:22,sep=""), "chrM")
DMRs_DCvWD <- NULL
for(i in 1:length(chroms_DCvWD)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_diseaseCon_vs_WD_DMRs/disCon_v_WD_",chroms_DCvWD[i],
                                 "_silver_DMR_info.txt", sep=""), header = TRUE, sep = "\t")
        DMRs_DCvWD <- rbind(DMRs_DCvWD, temp)
}
DMRs_DCvWD$chr <- as.character(DMRs_DCvWD$chr)

# Combining DMR Smoothed Methylation into one table (Healthy Ctrl vs WD, Reclassification, Chr X, Y excluded)
chroms_HCvWD = c(paste("chr",1:22,sep=""), "chrM")
meth_HCvWD <- NULL
for(i in 1:length(chroms_HCvWD)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_healthy_v_WD_DMRs/Reclass_healthy_v_WD_",chroms_HCvWD[i],
                                 "_silver_DMR_methylation.txt", sep=""), header = TRUE, sep = "\t")
        meth_HCvWD <- rbind(meth_HCvWD, temp)
}
meth_HCvWD$chr <- as.character(meth_HCvWD$chr)

# Combining DMR Smoothed Methylation into one table (Healthy Ctrl vs Disease Ctrl, Reclassification, Chr X, Y excluded)
chroms_HCvDC = c(paste("chr",1:22,sep=""), "chrM")
meth_HCvDC <- NULL
for(i in 1:length(chroms_HCvDC)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_healthyCon_v_diseaseCon_DMRs/Reclass_healthyCon_v_diseaseCon_",chroms_HCvDC[i],
                                 "_silver_DMR_methylation.txt", sep=""), header = TRUE, sep = "\t")
        meth_HCvDC <- rbind(meth_HCvDC, temp)
}
meth_HCvDC$chr <- as.character(meth_HCvDC$chr)

# Combining DMR Smoothed Methylation into one table (Disease Ctrl vs WD, Reclassification, Chr X, Y excluded)
chroms_DCvWD = c(paste("chr",1:22,sep=""), "chrM")
meth_DCvWD <- NULL
for(i in 1:length(chroms_DCvWD)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_diseaseCon_vs_WD_DMRs/disCon_v_WD_",chroms_DCvWD[i],
                                 "_silver_DMR_methylation.txt", sep=""), header = TRUE, sep = "\t")
        meth_DCvWD <- rbind(meth_DCvWD, temp)
}
meth_DCvWD$chr <- as.character(meth_DCvWD$chr)

# Combining DMR Coverage into one table (Healthy Ctrl vs WD, Reclassification, Chr X, Y excluded)
chroms_HCvWD = c(paste("chr",1:22,sep=""), "chrM")
cov_HCvWD <- NULL
for(i in 1:length(chroms_HCvWD)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_healthy_v_WD_DMRs/Reclass_healthy_v_WD_",chroms_HCvWD[i],
                                 "_silver_DMR_coverage.txt", sep=""), header = TRUE, sep = "\t")
        cov_HCvWD <- rbind(cov_HCvWD, temp)
}
cov_HCvWD$chr <- as.character(cov_HCvWD$chr)

# Combining DMR Coverage into one table (Healthy Ctrl vs Disease Ctrl, Reclassification, Chr X, Y excluded)
chroms_HCvDC = c(paste("chr",1:22,sep=""), "chrM")
cov_HCvDC <- NULL
for(i in 1:length(chroms_HCvDC)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_healthyCon_v_diseaseCon_DMRs/Reclass_healthyCon_v_diseaseCon_",chroms_HCvDC[i],
                                 "_silver_DMR_coverage.txt", sep=""), header = TRUE, sep = "\t")
        cov_HCvDC <- rbind(cov_HCvDC, temp)
}
cov_HCvDC$chr <- as.character(cov_HCvDC$chr)

# Combining DMR Coverage into one table (Disease Ctrl vs WD, Reclassification, Chr X, Y excluded)
chroms_DCvWD = c(paste("chr",1:22,sep=""), "chrM")
cov_DCvWD <- NULL
for(i in 1:length(chroms_DCvWD)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_diseaseCon_vs_WD_DMRs/disCon_v_WD_",chroms_DCvWD[i],
                                 "_silver_DMR_coverage.txt", sep=""), header = TRUE, sep = "\t")
        cov_DCvWD <- rbind(cov_DCvWD, temp)
}
cov_DCvWD$chr <- as.character(cov_DCvWD$chr)

# Combining Background DMRs into one table (Healthy Ctrl vs WD, Reclassification, Chr X, Y excluded)
chroms_HCvWD = c(paste("chr",1:22,sep=""), "chrM")
background_HCvWD <- NULL
for(i in 1:length(chroms_HCvWD)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_healthy_v_WD_DMRs/Reclass_healthy_v_WD_",chroms_HCvWD[i],
                                 "_background_DMRs.bed", sep=""), header = FALSE, sep = "\t")
        temp <- temp[2:length(temp[,1]),]
        background_HCvWD <- rbind(background_HCvWD, temp)
}

# Combining Background DMR Coverage into one table (Healthy Ctrl vs WD, Reclassification, Chr X, Y excluded)
chroms_HCvWD = c(paste("chr",1:22,sep=""), "chrM")
backgroundCov_HCvWD <- NULL
for(i in 1:length(chroms_HCvWD)){
        temp <- NULL
        temp <- read.delim(paste("Reclass_HCvsWD_background_DF5/Reclass_HCvsWD_DMRs_",chroms_HCvWD[i],
                                 "_background_DMR_coverage.txt", sep=""), header = TRUE, sep = "\t", stringsAsFactors=FALSE)
        backgroundCov_HCvWD <- rbind(backgroundCov_HCvWD, temp)
}
backgroundCov_HCvWD$end <- as.integer(backgroundCov_HCvWD$end)

# Healthy Control vs WD Analysis ###################

# Merge Methylation and Coverage with DMRs_HCvWD file
meth_HCvWD_info <- merge(meth_HCvWD, DMRs_HCvWD, by = c("chr", "start", "end"))
meth_HCvWD_info <- meth_HCvWD_info[,c(1:3,20:31,4:19)]

cov_HCvWD_info <- merge(cov_HCvWD, DMRs_HCvWD, by = c("chr", "start", "end"))
cov_HCvWD_info <- cov_HCvWD_info[,c(1:3,20:31,4:19)]

# Subset DMRs and get methdiff
meth <- meth_HCvWD_info[,16:31]
methavg <- rowMeans(meth, na.rm = TRUE)
methdiff <- meth - methavg
reads <- as.matrix(cov_HCvWD_info[,16:31])
table(is.na(reads)) # No NA values
#reads[is.na(reads)] <- 0 # replaced NA values with 0 (no reads)
minreads <- apply(reads, 1, min)
meanreads <- apply(reads, 1, mean)

samplesCtrl <- 6
samplesWD <- 10
precisionCtrl <- (1/samplesCtrl)*2
precisionWD <- (1/samplesWD)*2
precisionSum <- precisionCtrl + precisionWD
minDiff <- 0.1
coverage <- ceiling(precisionSum/minDiff)

# T-tests
DMRs <- 1:length(meth[,1])
ttest_HCvWD <- matrix(nrow = length(DMRs), ncol = 5)
for(i in DMRs){
        association <- NULL
        association <- t.test(meth[i,7:16], meth[i,1:6])
        ttest_HCvWD[i,1] <- association$estimate[1] - association$estimate[2]
        ttest_HCvWD[i,2] <- association$conf.int[1]
        ttest_HCvWD[i,3] <- association$conf.int[2]
        ttest_HCvWD[i,4] <- association$statistic
        ttest_HCvWD[i,5] <- association$p.value
}
colnames(ttest_HCvWD) <- c("meanDiff", "confIntL", "confIntR", "tstat", "pValue")
ttest_HCvWD <- as.data.frame(ttest_HCvWD)
ttest_HCvWD$pValueFDR <- p.adjust(ttest_HCvWD$pValue, "fdr")
ttest_HCvWD$pValueBonf <- p.adjust(ttest_HCvWD$pValue, "bonf")
ttest_HCvWD$minReads <- minreads
ttest_HCvWD$meanReads <- meanreads

ttest_HCvWD <- cbind(meth_HCvWD_info[,1:15], ttest_HCvWD[,2:length(ttest_HCvWD[1,])])
pValues <- ttest_HCvWD$pValue
meanDiffs <- ttest_HCvWD$meanDiff
ttest_HCvWD <- subset(ttest_HCvWD, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05)
#write.table(ttest_HCvWD, "Healthy Ctrl vs WD Liver DMR Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Subset HC vs WD Background
minreads <- apply(backgroundCov_HCvWD[,4:dim(backgroundCov_HCvWD)[2]], 1, min)

samplesCtrl <- 6
samplesWD <- 10
precisionCtrl <- (1/samplesCtrl)*2
precisionWD <- (1/samplesWD)*2
precisionSum <- precisionCtrl + precisionWD
minDiff <- 0.1
coverage <- ceiling(precisionSum/minDiff) # Need 6 reads in every sample for every DMR

plot(backgroundCov_HCvWD$start, background_HCvWD$start) # Same order

backgroundSub_HCvWD <- subset(background_HCvWD, minreads >= coverage) #573680 -> 279600

backgroundPrint_HCvWD <- backgroundSub_HCvWD[,c("chr", "start", "end")]
backgroundPrint_HCvWD <- backgroundPrint_HCvWD[order(backgroundPrint_HCvWD$chr, backgroundPrint_HCvWD$start),]
write.table(backgroundPrint_HCvWD, "Healthy_Ctrl_vs_WD_Liver_Subsetted_DMR_Background.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Healthy Control vs Disease Control Analysis ####################

# Merge Methylation and Coverage with DMRs_HCvWD file
meth_HCvDC_info <- merge(meth_HCvDC, DMRs_HCvDC, by = c("chr", "start", "end"))
meth_HCvDC_info <- meth_HCvDC_info[,c(1:3,15:26,4:14)]

cov_HCvDC_info <- merge(cov_HCvDC, DMRs_HCvDC, by = c("chr", "start", "end"))
cov_HCvDC_info <- cov_HCvDC_info[,c(1:3,15:26,4:14)]

# Subset DMRs and get methdiff
meth <- meth_HCvDC_info[,16:26]
methavg <- rowMeans(meth, na.rm = TRUE)
methdiff <- meth - methavg
reads <- as.matrix(cov_HCvDC_info[,16:26])
table(is.na(reads)) # No NA values
#reads[is.na(reads)] <- 0 # replaced NA values with 0 (no reads)
minreads <- apply(reads, 1, min)
meanreads <- apply(reads, 1, mean)

samplesCtrl <- 6
samplesDC <- 5
precisionCtrl <- (1/samplesCtrl)*2
precisionDC <- (1/samplesDC)*2
precisionSum <- precisionCtrl + precisionDC
minDiff <- 0.1
coverage <- ceiling(precisionSum/minDiff)

# T-tests
DMRs <- 1:length(meth[,1])
ttest_HCvDC <- matrix(nrow = length(DMRs), ncol = 5)
for(i in DMRs){
        association <- NULL
        association <- t.test(meth[i,7:11], meth[i,1:6])
        ttest_HCvDC[i,1] <- association$estimate[1] - association$estimate[2]
        ttest_HCvDC[i,2] <- association$conf.int[1]
        ttest_HCvDC[i,3] <- association$conf.int[2]
        ttest_HCvDC[i,4] <- association$statistic
        ttest_HCvDC[i,5] <- association$p.value
}
colnames(ttest_HCvDC) <- c("meanDiff", "confIntL", "confIntR", "tstat", "pValue")
ttest_HCvDC <- as.data.frame(ttest_HCvDC)
ttest_HCvDC$pValueFDR <- p.adjust(ttest_HCvDC$pValue, "fdr")
ttest_HCvDC$pValueBonf <- p.adjust(ttest_HCvDC$pValue, "bonf")
ttest_HCvDC$minReads <- minreads
ttest_HCvDC$meanReads <- meanreads

ttest_HCvDC <- cbind(meth_HCvDC_info[,1:15], ttest_HCvDC[,2:length(ttest_HCvDC[1,])])
pValues <- ttest_HCvDC$pValue
meanDiffs <- ttest_HCvDC$meanDiff
ttest_HCvDC <- subset(ttest_HCvDC, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05)
#write.table(ttest_HCvDC, "Healthy Ctrl vs Disease Ctrl Liver DMR Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# DC vs WD Analysis #####

# Merge Methylation and Coverage with DMRs_DCvWD file
meth_DCvWD_info <- merge(meth_DCvWD, DMRs_DCvWD, by = c("chr", "start", "end"))
meth_DCvWD_info <- meth_DCvWD_info[,c(1:3,19:30,4:18)]

cov_DCvWD_info <- merge(cov_DCvWD, DMRs_DCvWD, by = c("chr", "start", "end"))
cov_DCvWD_info <- cov_DCvWD_info[,c(1:3,19:30,4:18)]

# Subset DMRs and get methdiff
meth <- meth_DCvWD_info[,16:30]
methavg <- rowMeans(meth, na.rm = TRUE)
methdiff <- meth - methavg
reads <- as.matrix(cov_DCvWD_info[,16:30])
table(is.na(reads)) # No NA values
#reads[is.na(reads)] <- 0 # replaced NA values with 0 (no reads)
minreads <- apply(reads, 1, min)
meanreads <- apply(reads, 1, mean)

samplesCtrl <- 5
samplesWD <- 10
precisionCtrl <- (1/samplesCtrl)*2
precisionWD <- (1/samplesWD)*2
precisionSum <- precisionCtrl + precisionWD
minDiff <- 0.1
coverage <- ceiling(precisionSum/minDiff) # 7 reads per sample (35 reads DC, 70 reads WD)

# T-tests
DMRs <- 1:length(meth[,1])
ttest_DCvWD <- matrix(nrow = length(DMRs), ncol = 5)
for(i in DMRs){
        association <- NULL
        association <- t.test(meth[i,6:15], meth[i,1:5])
        ttest_DCvWD[i,1] <- association$estimate[1] - association$estimate[2]
        ttest_DCvWD[i,2] <- association$conf.int[1]
        ttest_DCvWD[i,3] <- association$conf.int[2]
        ttest_DCvWD[i,4] <- association$statistic
        ttest_DCvWD[i,5] <- association$p.value
}
colnames(ttest_DCvWD) <- c("meanDiff", "confIntL", "confIntR", "tstat", "pValue")
ttest_DCvWD <- as.data.frame(ttest_DCvWD)
ttest_DCvWD$pValueFDR <- p.adjust(ttest_DCvWD$pValue, "fdr")
ttest_DCvWD$pValueBonf <- p.adjust(ttest_DCvWD$pValue, "bonf")
ttest_DCvWD$minReads <- minreads
ttest_DCvWD$meanReads <- meanreads

ttest_DCvWD <- cbind(meth_DCvWD_info[,1:15], ttest_DCvWD[,2:length(ttest_DCvWD[1,])])
pValues <- ttest_DCvWD$pValue
meanDiffs <- ttest_DCvWD$meanDiff
ttest_DCvWD <- subset(ttest_DCvWD, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #15024 -> 3774
#write.table(ttest_DCvWD, "Tables/Disease Ctrl vs WD Liver DMR Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Overlap HCvsWD, HCvsDC, and DCvWD DMRs ####
ttest_HCvWD <- read.delim("Healthy Ctrl vs WD Liver DMR Stats.txt", sep="\t")
ttest_HCvDC <- read.delim("Healthy Ctrl vs Disease Ctrl Liver DMR Stats.txt", sep="\t")

ttest_HCvWD_hyper <- subset(ttest_HCvWD, meanDiff > 0)
ttest_HCvWD_hypo <- subset(ttest_HCvWD, meanDiff < 0)
ttest_HCvDC_hyper <- subset(ttest_HCvDC, meanDiff > 0)
ttest_HCvDC_hypo <- subset(ttest_HCvDC, meanDiff < 0)
ttest_DCvWD_hyper <- subset(ttest_DCvWD, meanDiff > 0)
ttest_DCvWD_hypo <- subset(ttest_DCvWD, meanDiff < 0)

GR_HCvWD <- GRanges(seqnames = ttest_HCvWD$chr, ranges=IRanges(start=ttest_HCvWD$start, end=ttest_HCvWD$end))
GR_HCvDC <- GRanges(seqnames = ttest_HCvDC$chr, ranges=IRanges(start=ttest_HCvDC$start, end=ttest_HCvDC$end))
GR_DCvWD <- GRanges(seqnames = ttest_DCvWD$chr, ranges=IRanges(start=ttest_DCvWD$start, end=ttest_DCvWD$end))
GR_HCvWD_hyper <- GRanges(seqnames = ttest_HCvWD_hyper$chr, ranges=IRanges(start=ttest_HCvWD_hyper$start, end=ttest_HCvWD_hyper$end))
GR_HCvWD_hypo <- GRanges(seqnames = ttest_HCvWD_hypo$chr, ranges=IRanges(start=ttest_HCvWD_hypo$start, end=ttest_HCvWD_hypo$end))
GR_HCvDC_hyper <- GRanges(seqnames = ttest_HCvDC_hyper$chr, ranges=IRanges(start=ttest_HCvDC_hyper$start, end=ttest_HCvDC_hyper$end))
GR_HCvDC_hypo <- GRanges(seqnames = ttest_HCvDC_hypo$chr, ranges=IRanges(start=ttest_HCvDC_hypo$start, end=ttest_HCvDC_hypo$end))
GR_DCvWD_hyper <- GRanges(seqnames = ttest_DCvWD_hyper$chr, ranges=IRanges(start=ttest_DCvWD_hyper$start, end=ttest_DCvWD_hyper$end))
GR_DCvWD_hypo <- GRanges(seqnames = ttest_DCvWD_hypo$chr, ranges=IRanges(start=ttest_DCvWD_hypo$start, end=ttest_DCvWD_hypo$end))

pdf(file="Figures/Subsetted Hyper DMR Overlap Venn.pdf", width=10, height=9, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_HCvDC_hyper, GR_HCvWD_hyper, GR_DCvWD_hyper), NameOfPeaks = c("HC_vs_WD_Hyper", "HC_vs_DC_Hyper", "DC_vs_WD_Hyper"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 240, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightpink", "lightblue", "lightgreen"), cat.pos = c(330,180,30), 
                        cat.dist = c(0.1, 0.05, 0.1), cat.col=c("white", "white", "white"), fontfamily="sans")
dev.off()

pdf(file="Figures/Subsetted Hypo DMR Overlap Venn.pdf", width=10, height=9, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_HCvWD_hypo, GR_HCvDC_hypo, GR_DCvWD_hypo), NameOfPeaks = c("HC_vs_WD_Hypo", "HC_vs_DC_Hypo", "DC_vs_WD_Hypo"),
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 180, margin = 0.02,
                        cat.cex = 2, cex = 2.5, fill = c("lightblue", "lightpink", "lightgreen"), cat.pos = c(15,345,180), 
                        cat.dist = c(0.02, 0.08, 0.04), cat.col=c("white", "white", "white"), fontfamily="sans")
dev.off()

GR_WD_Spec_hyper <- GR_HCvWD_hyper[HCvWD_overlaps_DCvWD_hyper & !HCvWD_overlaps_HCvDC_hyper,] #969 DMRs
GR_WD_Spec_hypo <- GR_HCvWD_hypo[HCvWD_overlaps_DCvWD_hypo & !HCvWD_overlaps_HCvDC_hypo,] #871 DMRs
GR_WD_Spec <- c(GR_WD_Spec_hyper, GR_WD_Spec_hypo) #1840 DMRs
isDisjoint(GR_WD_Spec) #TRUE (not overlapping)

DMRs_WD_Spec <- as.data.frame(GR_WD_Spec)
DMRs_WD_Spec_bed <- DMRs_WD_Spec[,c("seqnames", "start", "end")]
DMRs_WD_Spec_bed <- DMRs_WD_Spec_bed[order(DMRs_WD_Spec_bed$seqnames, DMRs_WD_Spec_bed$start),]
write.table(DMRs_WD_Spec_bed, "WD_Specific_DMRs.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

DMRs_WD_Spec_hyper <- as.data.frame(GR_WD_Spec_hyper)
DMRs_WD_Spec_hyper_bed <- DMRs_WD_Spec_hyper[,c("seqnames", "start", "end")]
DMRs_WD_Spec_hyper_bed <- DMRs_WD_Spec_hyper_bed[order(DMRs_WD_Spec_hyper_bed$seqnames, DMRs_WD_Spec_hyper_bed$start),]
write.table(DMRs_WD_Spec_hyper_bed, "WD_Specific_hyper_DMRs.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

DMRs_WD_Spec_hypo <- as.data.frame(GR_WD_Spec_hypo)
DMRs_WD_Spec_hypo_bed <- DMRs_WD_Spec_hypo[,c("seqnames", "start", "end")]
DMRs_WD_Spec_hypo_bed <- DMRs_WD_Spec_hypo_bed[order(DMRs_WD_Spec_hypo_bed$seqnames, DMRs_WD_Spec_hypo_bed$start),]
write.table(DMRs_WD_Spec_hypo_bed, "WD Specific_hypo_DMRs.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

ttest_WD_Spec <- merge(ttest_HCvWD, DMRs_WD_Spec, by.x=c("chr", "start", "end"), by.y=c("seqnames", "start", "end"), all.x=FALSE)

# GREAT WD Specific DMRs ####
# Get background file
Background <- read.delim("Healthy_Ctrl_vs_WD_Liver_Subsetted_DMR_Background.bed", sep="\t", header=FALSE)
colnames(Background) <- c("chr", "start", "end")
Background$chr <- as.character(Background$chr)
Background$start <- as.integer(Background$start)
Background$end <- as.integer(Background$end)
Background <- Background[order(Background$chr, Background$start),]
GR_Background <- GRanges(seqnames = Background$chr, ranges = IRanges(start = Background$start, end = Background$end))

# Get chain file
url <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
download.file(url, basename(url))
gunzip(basename(url))
chain <- import.chain("hg38ToHg19.over.chain")

# liftOver from hg38 to hg19
seqlevelsStyle(GR_WD_Spec) <- "UCSC"
seqlevelsStyle(GR_WD_Spec_hyper) <- "UCSC"
seqlevelsStyle(GR_WD_Spec_hypo) <- "UCSC"
seqlevelsStyle(GR_Background) <- "UCSC"
GR_WD_Spec_hg19 <- unlist(liftOver(GR_WD_Spec, chain))
GR_WD_Spec_hyper_hg19 <- unlist(liftOver(GR_WD_Spec_hyper, chain))
GR_WD_Spec_hypo_hg19 <- unlist(liftOver(GR_WD_Spec_hypo, chain))
GR_Background_hg19 <- unlist(liftOver(GR_Background, chain))

# If regions are overlapping, make them non-overlapping (disjoint)
if(!isDisjoint(GR_WD_Spec_hg19)){GR_WD_Spec_hg19 <- disjoin(GR_WD_Spec_hg19)} 
if(!isDisjoint(GR_WD_Spec_hyper_hg19)){GR_WD_Spec_hyper_hg19 <- disjoin(GR_WD_Spec_hyper_hg19)} 
if(!isDisjoint(GR_WD_Spec_hypo_hg19)){GR_WD_Spec_hypo_hg19 <- disjoin(GR_WD_Spec_hypo_hg19)} 
if(!isDisjoint(GR_Background_hg19)){GR_Background_hg19 <- disjoin(GR_Background_hg19)} 

# Reduce background regions to < 1M (Not Needed)
# Background must have < 1M regions for GREAT, this merges regions closer than 400 bases, change if needed
#if(length(GR_Background_hg19) > 1000000){GR_Background_hg19 <- reduce(GR_Background_hg19, min.gapwidth=400)} # New length: 937954

# Redefine so DMRs are subset of background regions (slow)
GR_WD_Spec_hg19 <- redefineUserSets(GR_WD_Spec_hg19, GR_Background_hg19)
GR_WD_Spec_hyper_hg19 <- redefineUserSets(GR_WD_Spec_hyper_hg19, GR_Background_hg19)
GR_WD_Spec_hypo_hg19 <- redefineUserSets(GR_WD_Spec_hypo_hg19, GR_Background_hg19)

# 5. Make bed format and write files
# All DMRs
chroms <- c(paste("chr",1:22,sep=""), "chrM") # Excluded chrX,Y
allDMRs_hg19 <- unlist(GRangesList(GR_WD_Spec_hg19))
allDMRs_hg19 <- as.data.frame(allDMRs_hg19)[,c("seqnames", "start", "end")]
colnames(allDMRs_hg19) <- c("chr", "start", "end")
allDMRs_hg19$chr <- as.character(allDMRs_hg19$chr)
allDMRs_hg19 <- allDMRs_hg19[order(allDMRs_hg19$chr, allDMRs_hg19$start),]
allDMRs_hg19 <- unique(subset(allDMRs_hg19, chr %in% chroms))

# Hyper DMRs
hyperDMRs_hg19 <- unlist(GRangesList(GR_WD_Spec_hyper_hg19))
hyperDMRs_hg19 <- as.data.frame(hyperDMRs_hg19)[,c("seqnames", "start", "end")]
colnames(hyperDMRs_hg19) <- c("chr", "start", "end")
hyperDMRs_hg19$chr <- as.character(hyperDMRs_hg19$chr)
hyperDMRs_hg19 <- hyperDMRs_hg19[order(hyperDMRs_hg19$chr, hyperDMRs_hg19$start),]
hyperDMRs_hg19 <- unique(subset(hyperDMRs_hg19, chr %in% chroms))

# Hypo DMRs
hypoDMRs_hg19 <- unlist(GRangesList(GR_WD_Spec_hypo_hg19))
hypoDMRs_hg19 <- as.data.frame(hypoDMRs_hg19)[,c("seqnames", "start", "end")]
colnames(hypoDMRs_hg19) <- c("chr", "start", "end")
hypoDMRs_hg19$chr <- as.character(hypoDMRs_hg19$chr)
hypoDMRs_hg19 <- hypoDMRs_hg19[order(hypoDMRs_hg19$chr, hypoDMRs_hg19$start),]
hypoDMRs_hg19 <- unique(subset(hypoDMRs_hg19, chr %in% chroms))

# Background
Background_hg19 <- as.data.frame(GR_Background_hg19)[,c("seqnames", "start", "end")]
colnames(Background_hg19) <- c("chr", "start", "end")
Background_hg19$chr <- as.character(Background_hg19$chr)
Background_hg19 <- Background_hg19[order(Background_hg19$chr, Background_hg19$start),]
Background_hg19 <- unique(subset(Background_hg19, chr %in% chroms))

# Check that DMRs are a subset of background
table(allDMRs_hg19$chr %in% Background_hg19$chr & allDMRs_hg19$start %in% Background_hg19$start & allDMRs_hg19$end %in% Background_hg19$end) # All TRUE
table(hyperDMRs_hg19$chr %in% Background_hg19$chr & hyperDMRs_hg19$start %in% Background_hg19$start & hyperDMRs_hg19$end %in% Background_hg19$end) # All TRUE
table(hypoDMRs_hg19$chr %in% Background_hg19$chr & hypoDMRs_hg19$start %in% Background_hg19$start & hypoDMRs_hg19$end %in% Background_hg19$end) # All TRUE

# Write bed files
write.table(allDMRs_hg19, "WD_Specific_allDMRs_hg19_for_GREAT.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(hyperDMRs_hg19, "WD_Specific_hyperDMRs_hg19_for_GREAT.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(hypoDMRs_hg19, "WD_Specific_hypoDMRs_hg19_for_GREAT.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
write.table(Background_hg19, "WD_Specific_Background_hg19_for_GREAT.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

