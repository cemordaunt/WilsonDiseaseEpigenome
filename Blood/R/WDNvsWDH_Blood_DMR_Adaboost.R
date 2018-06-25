# WDN vs WDH Blood DMR Methylation for Adaboost Machine Learning ####
# Charles Mordaunt
# 6/13/18

# Working Directory ####
setwd()

# Packages ####
library(GenomicRanges)
library(ggdendro)
library(ggplot2)
library(reshape2)
library(grid)
library(scales)
library(plyr)
library(reshape)

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

# Heatmap with Pheno Data Functions
mydplot_pheno <- function(ddata, row=!col, col=!row, labels=col) {
        # plots a dendrogram
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
        # plots a legend
        # from http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

ggheatmap.show_pheno <- function(L, widths=c(0.02,0.79,0.17,0.02), heights=c(0.02,0.17,0.04,0.75,0.02)){
        # plots the heatmap
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
                                               legend.position = c(0.6, -0.2)))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=2))
        grid.draw(legend)
        upViewport(0)
        
        ## add pheno legend
        phenoLegend <- g_legend_pheno(L$phenoData +
                                              theme(legend.title = element_text(size=16), 
                                                    legend.text = element_text(size = 15),
                                                    legend.direction = "vertical",
                                                    legend.position = c(0.917, 0.92),
                                                    legend.background = element_blank()))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=3))
        grid.draw(phenoLegend)
        upViewport(0)
}

ggheatmap2_pheno <- function(x, phenoData, name="Diagnosis", breaks=c("HC", "WD"), values=c("HC"="#3366CC", "WD"="#FF3366"), hm.colours=c("#0000FF", "Black", "#FF0000"), 
                             my.values=c(0,0.5,1), low=min(x), high=max(x)) {
        # makes the heatmap
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
                theme(plot.margin = unit(c(0,-0.5,0,-0.5), "lines"),
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
                scale_color_manual(name=name, breaks = breaks, values = values) +
                scale_fill_manual(name=name, breaks = breaks, values = values)
        ret <- list(col=col.plot,row=row.plot,centre=centre.plot, phenoData=phenoData.plot)
        invisible(ret)
}

# WDNvWDH DMR Methylation Data ####
WDNvWDH_DMRs <- read.delim("UCSC Tracks/WD_Blood_2_WDNvWDH_DMRs.bed", sep="\t", header=FALSE)
colnames(WDNvWDH_DMRs) <- c("chr", "start", "end")
WDNvWDH_DMRs <- WDNvWDH_DMRs[order(WDNvWDH_DMRs$chr, WDNvWDH_DMRs$start),]
GR_WDNvWDH_DMRs <- GRanges(seqnames=WDNvWDH_DMRs$chr, ranges=IRanges(start=WDNvWDH_DMRs$start, end=WDNvWDH_DMRs$end))

# Blood 2 WDNvWDH DMR Methylation (Training Set)
NvH_blood2_meth <- combineFiles(chroms=paste("chr",1:22, sep=""), prefix="Tables/WDNvWDH_DMR_Blood2_Methyl/WDNvWDH_",
                            suffix="_silver_DMR_methylation.txt")
NvH_blood2_meth <- NvH_blood2_meth[order(NvH_blood2_meth$chr, NvH_blood2_meth$start),]
GR_NvH_meth <- GRanges(seqnames=NvH_blood2_meth$chr, ranges=IRanges(start=NvH_blood2_meth$start, end=NvH_blood2_meth$end))
NvH_blood2_meth <- subset(NvH_blood2_meth, GR_NvH_meth %over% GR_WDNvWDH_DMRs)
table(WDNvWDH_DMRs$start == NvH_blood2_meth$start) #All true

NvH_blood2_cov <- combineFiles(chroms=paste("chr",1:22, sep=""), prefix="Tables/WDNvWDH_DMR_Blood2_Methyl/WDNvWDH_",
                           suffix="_silver_DMR_coverage.txt")
NvH_blood2_cov <- NvH_blood2_cov[order(NvH_blood2_cov$chr, NvH_blood2_cov$start),]
NvH_blood2_cov <- subset(NvH_blood2_cov, GR_NvH_meth %over% GR_WDNvWDH_DMRs)
table(WDNvWDH_DMRs$start == NvH_blood2_cov$start) #All true

NvH_blood2_mincov <- apply(NvH_blood2_cov[4:ncol(NvH_blood2_cov)], 1, min)
summary(NvH_blood2_mincov)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.00    5.00    8.00   11.02   13.00  146.00 All DMRs have at least 4 reads per sample

# Blood 1 WDNvWDH DMR Methylation (Test Set)
NvH_blood1_meth <- combineFiles(chroms=paste("chr",1:22, sep=""), prefix="Tables/WDNvWDH_DMR_blood1_Methyl/WDNvWDH_DMR_Blood1_Methyl_",
                                suffix="_DMR_methylation.txt")
NvH_blood1_meth <- NvH_blood1_meth[order(NvH_blood1_meth$chr, NvH_blood1_meth$start),]

NvH_blood1_cov <- combineFiles(chroms=paste("chr",1:22, sep=""), prefix="Tables/WDNvWDH_DMR_blood1_Methyl/WDNvWDH_DMR_Blood1_Methyl_",
                               suffix="_DMR_coverage.txt")
NvH_blood1_cov <- NvH_blood1_cov[order(NvH_blood1_cov$chr, NvH_blood1_cov$start),]

NvH_blood1_mincov <- apply(NvH_blood1_cov[4:ncol(NvH_blood1_cov)], 1, min)
summary(NvH_blood1_mincov)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  1.000   2.000   4.000   6.619   8.000  72.000     363 

NvH_blood1_meth <- subset(NvH_blood1_meth, NvH_blood1_mincov >= 1) # All samples must have at least 1 read in each DMR
GR_NvH_blood1_meth <- GRanges(seqnames=NvH_blood1_meth$chr, ranges=IRanges(start=NvH_blood1_meth$start, end=NvH_blood1_meth$end))
write.table(NvH_blood1_meth, "Tables/WDN vs WDH Blood Test Sample Methylation.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Subset so NvH Blood 1 and 2 Methylation is on same DMRs
NvH_blood2_meth <- subset(NvH_blood2_meth, GR_WDNvWDH_DMRs %over% GR_NvH_blood1_meth)
table(NvH_blood1_meth$start == NvH_blood2_meth$start) #All TRUE, 2142 DMRs
write.table(NvH_blood2_meth, "Tables/WDN vs WDH Blood Training Sample Methylation.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# WDN vs WDH Samples ####
# Blood 2 (Training Set)
blood2_samples <- read.csv("Tables/WD Blood 2 Sample Info.csv", header=TRUE, stringsAsFactors = FALSE)
blood2_samples <- blood2_samples[match(colnames(NvH_blood2_meth[4:ncol(NvH_blood2_meth)]), blood2_samples$Sequencing.ID),]
blood2_samples$Diagnosis[blood2_samples$Diagnosis == "Wilson's Disease"] <- "Wilson Disease"
blood2_samples$WD_Phenotype <- NA
blood2_samples$WD_Phenotype[blood2_samples$Neurologic.form == 1] <- "Neurologic"
blood2_samples$WD_Phenotype[blood2_samples$Hepatic.form == 1] <- "Hepatic"
blood2_samples <- subset(blood2_samples, Diagnosis == "Wilson Disease", select=c(Sequencing.ID, Diagnosis, WD_Phenotype))
write.table(blood2_samples, "Tables/WDN vs WDH Training Samples.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Blood 1 (Test Set)
blood1_samples <- read.csv("Tables/WD Blood 1 Sample Info.csv", header=TRUE, stringsAsFactors = FALSE)
blood1_samples <- subset(blood1_samples, Diagnosis == "Wilson Disease")
blood1_samples <- blood1_samples[match(colnames(NvH_blood1_meth[4:ncol(NvH_blood1_meth)]), blood1_samples$Sequencing.ID),]
hep <- c("VMNS007A", "VMNS009A", "VMNS006A", "VMNS005A", "VMNS008A")
neu <- c("VMNS007C", "VMNS008C", "VMNS005C", "VMNS009C", "VMNS006C")
blood1_samples$WD_Phenotype <- NA
blood1_samples$WD_Phenotype[blood1_samples$Sequencing.ID %in% hep] <- "Hepatic"
blood1_samples$WD_Phenotype[blood1_samples$Sequencing.ID %in% neu] <- "Neurologic"
write.table(blood1_samples, "Tables/WDN vs WDH Test Samples.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Adaboost Classifier Sample Dotplot ####
sample_stats <- read.delim("Tables/Adaboost All DMR Sample Predictions.txt", sep="\t", stringsAsFactors = FALSE, header=TRUE)
sample_stats$WD_Phenotype <- factor(sample_stats$WD_Phenotype, levels=c("Neurologic", "Hepatic"), ordered=TRUE)
sample_stats$Classification <- factor(sample_stats$Classification, levels=c("Neurologic", "Hepatic"), ordered=TRUE)
means <- data.frame(WD_Phenotype=factor(c("Neurologic", "Hepatic"), levels=c("Neurologic", "Hepatic"), ordered=TRUE), 
                    Y_predict=c(mean(sample_stats$Y_predict[sample_stats$WD_Phenotype == "Neurologic"]),
                           mean(sample_stats$Y_predict[sample_stats$WD_Phenotype == "Hepatic"])))

g <- ggplot()
g + 
        geom_hline(aes(yintercept=0.5), size=1.25, lty=2) +
        geom_dotplot(data=sample_stats, aes(x=WD_Phenotype, y=Y_predict, fill=Classification, color=Classification), 
                     binaxis = "y", stackdir = "center", binwidth = 0.04, stackratio = 1.2, dotsize = 0.8) +
        geom_errorbar(data=means, aes(x=WD_Phenotype, ymin=Y_predict, ymax=Y_predict), width=0.75, color="black", size=2) +
        #annotate("text", x=c(1.5), y=c(65), label="*", size=14) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(0.17,0.9), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=24), 
              legend.text=element_text(size=20), legend.key.size=unit(1.5, "lines"), axis.text = element_text(color = "black", size=22), 
              legend.background = element_blank(), axis.title = element_text(size=24)) +
        scale_y_continuous(breaks=pretty_breaks(n=3)) +
        coord_cartesian(ylim=c(0,1)) +
        scale_fill_manual(breaks = c("Neurologic", "Hepatic"), 
                          values = c("Neurologic"="#3366CC", "Hepatic"="#FF3366")) +
        scale_color_manual(breaks = c("Neurologic", "Hepatic"), 
                           values = c("Neurologic"="#3366CC", "Hepatic"="#FF3366")) +
        xlab("WD Phenotype") +
        ylab("Classifier Value")
ggsave("Figures/WDN vs WDH Adaboost Classifier Sample Dotplot.png", dpi = 600, width = 8, height = 7, units = "in")

# Adaboost Important DMR Training Sample Methylation Heatmap ####
important_DMR_meth <- read.csv("Tables/WD_Ada YZ.csv", header=TRUE, stringsAsFactors = FALSE)
meth <- important_DMR_meth[,5:ncol(important_DMR_meth)]
blood2_samples <- read.delim("Tables/WDN vs WDH Training Samples.txt", sep="\t", stringsAsFactors = FALSE, header = TRUE)
phenoData <- blood2_samples[,c("Sequencing.ID", "WD_Phenotype")]
colnames(phenoData)[1] <- "Sample"
phenoData <- phenoData[match(colnames(meth), phenoData$Sample),]
table(colnames(meth) == phenoData$Sample) # All TRUE must be the same order
phenoData$Sample <- factor(phenoData$Sample, levels=unique(phenoData$Sample), ordered=TRUE)
phenoData$WD_Phenotype[phenoData$WD_Phenotype == "Neurologic"] <- "WDN"
phenoData$WD_Phenotype[phenoData$WD_Phenotype == "Hepatic"] <- "WDH"
phenoData$WD_Phenotype <- factor(phenoData$WD_Phenotype, levels=c("WDN", "WDH"), ordered=TRUE)

methavg <- rowMeans(meth, na.rm = TRUE)
methdiff <- meth - methavg # subtract mean methylation for each DMR across all samples
methdiff <- as.matrix(methdiff)*100 # Transform to 0-100% scale

methplot <- ggheatmap2_pheno(x=methdiff, phenoData=phenoData, name="", breaks=c("WDN", "WDH"), values=c("WDN"="#3366CC", "WDH"="#FF3366"))
pdf(file="Figures/WDN vs WDH Adaboost Important DMR Training Sample Methylation Heatmap.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show_pheno(methplot)
dev.off()
