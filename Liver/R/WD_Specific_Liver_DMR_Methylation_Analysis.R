# WD Specific DMR Methylation Analysis
# HC, DC, WD Liver Samples
# Charles Mordaunt
# 12/14/17

# Packages ####
library(ggdendro)
library(ggplot2)
library(reshape2)
library(reshape)
library(grid)
library(scales)
library(ggbiplot)
library(RColorBrewer)

# Functions ####
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
                theme(plot.margin = unit(c(0,1,0.5,1), "lines"))
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
                                               legend.position = c(0.49, -0.63)))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=2))
        grid.draw(legend)
        upViewport(0)
        
        ## add pheno legend
        phenoLegend <- g_legend_pheno(L$phenoData +
                                              theme(legend.title = element_blank(), 
                                                    legend.text = element_text(size = 15),
                                                    legend.direction = "vertical",
                                                    legend.position = c(0.923, 0.9),
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
                theme(plot.margin = unit(c(0,-0.7,0,-1), "lines"),
                      axis.text.x = element_blank())
        row.plot <- mydplot_pheno(row.dendro, row=TRUE, labels=FALSE) +
                theme(plot.margin = unit(c(0,1,0,0), "lines"))
        
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
                scale_color_manual(breaks = c("HC", "WD", "DC", "Male", "Female"), 
                                   values = c("HC"="#3366CC", "WD"="#FF3366", "DC"="#009933", "Male"="#FFFF33", "Female"="#FF6633")) +
                scale_fill_manual(breaks = c("HC", "WD", "DC", "Male", "Female"), 
                                   values = c("HC"="#3366CC", "WD"="#FF3366", "DC"="#009933", "Male"="#FFFF33", "Female"="#FF6633"))
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
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, color = "black"))
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
                theme(plot.margin = unit(c(0,1.5,0.5,1), "lines"))
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

# Data ####
# Sample Info
info <- read.csv("Tables/WGBS Liver Sample Info.csv", header=TRUE)
info <- info[,c("Sequencing.ID", "Original.ID", "Sample.ID", "Specimen.provider", "Group", "Sex", "Age", "BMI", "FibrosisStage", "Steatosis", "InflammationGrade")]
info$FigureLabel <- c(paste("HC", 1:6,sep=""), paste("WD", 1:10, sep=""), paste("DC", 1:5, sep=""))
write.table(info, "Tables/WGBS Liver Sample Info 2.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")

# Combining DMR Smoothed Methylation into one table (WD Specific DMRs, HC, DC, WD samples, Chr X, Y excluded)
chroms <- c(paste("chr",1:22,sep=""), "chrM")
meth <- NULL
for(i in 1:length(chroms)){
        temp <- NULL
        temp <- read.delim(paste("Tables/WD Specific Methylation/WD_specific_",chroms[i],
                                 "_DMR_methylation.txt", sep=""), header = TRUE, sep = "\t")
        meth <- rbind(meth, temp)
}
meth$chr <- as.character(meth$chr)

# Combining DMR Coverage into one table (WD Specific DMRs, HC, DC, WD samples, Chr X, Y excluded)
cov <- NULL
for(i in 1:length(chroms)){
        temp <- NULL
        temp <- read.delim(paste("Tables/WD Specific Methylation/WD_specific_",chroms[i],
                                 "_DMR_coverage.txt", sep=""), header = TRUE, sep = "\t")
        cov <- rbind(cov, temp)
}
cov$chr <- as.character(cov$chr)
rm(temp, i, chroms)

# Coverage Summary ####
covStats <- as.data.frame(cbind(cov[,1:3], t(apply(cov[,4:dim(cov)[2]], 1, summary))))
hist(covStats$Min., breaks=200)
hist(covStats$Mean, breaks=200)
table(covStats$Min. < 6) #18 DMRs with at least 1 sample below coverage threshold
covStatsLow <- subset(covStats, Min. < 6)
hist(covStatsLow$Min., breaks=10)
hist(covStatsLow$Mean, breaks=10)
covLow <- subset(cov, covStats$Min. < 6)
# DC Samples have low coverage, but comparable to WD and HC samples
# No DC Samples have coverage < 3
# Keep DMRs in
rm(covStats, covStatsLow, covLow, cov)

# Methylation T-tests ####
DMRs <- 1:length(meth[,1])
HCsamples <- 4:9
WDsamples <- 10:19
DCsamples <- 20:24
ttest <- matrix(nrow = length(DMRs), ncol = 15)
for(i in DMRs){
        # WD vs HC
        association <- t.test(meth[i,WDsamples], meth[i,HCsamples])
        ttest[i,1] <- association$estimate[1] - association$estimate[2]
        ttest[i,2] <- association$conf.int[1]
        ttest[i,3] <- association$conf.int[2]
        ttest[i,4] <- association$statistic
        ttest[i,5] <- association$p.value
        # DC vs HC
        association <- t.test(meth[i,DCsamples], meth[i,HCsamples])
        ttest[i,6] <- association$estimate[1] - association$estimate[2]
        ttest[i,7] <- association$conf.int[1]
        ttest[i,8] <- association$conf.int[2]
        ttest[i,9] <- association$statistic
        ttest[i,10] <- association$p.value
        # WD vs DC
        association <- t.test(meth[i,WDsamples], meth[i,DCsamples])
        ttest[i,11] <- association$estimate[1] - association$estimate[2]
        ttest[i,12] <- association$conf.int[1]
        ttest[i,13] <- association$conf.int[2]
        ttest[i,14] <- association$statistic
        ttest[i,15] <- association$p.value
}
colnames(ttest) <- c("WDvsHC_meanDiff", "WDvsHC_confIntL", "WDvsHC_confIntR", "WDvsHC_tstat", "WDvsHC_pValue",
                     "DCvsHC_meanDiff", "DCvsHC_confIntL", "DCvsHC_confIntR", "DCvsHC_tstat", "DCvsHC_pValue",
                     "WDvsDC_meanDiff", "WDvsDC_confIntL", "WDvsDC_confIntR", "WDvsDC_tstat", "WDvsDC_pValue")
ttest <- as.data.frame(ttest)
ttest$WDvsHC_pValuebonf <- p.adjust(ttest$WDvsHC_pValue, "bonf")
ttest$DCvsHC_pValuebonf <- p.adjust(ttest$DCvsHC_pValue, "bonf")
ttest$WDvsDC_pValuebonf <- p.adjust(ttest$WDvsDC_pValue, "bonf")
ttest <- cbind(meth[,1:3], ttest)
write.table(ttest, "Tables/WD Specific DMR Stats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Methylation Heatmap with phenodata ####
methavg <- rowMeans(meth[,4:dim(meth)[2]], na.rm = TRUE)
methdiff <- meth[,4:dim(meth)[2]] - methavg
methdiff <- as.matrix(methdiff)
methdiff <- methdiff*100

phenoData <- info[,c("Sequencing.ID", "Group", "Sex")]
colnames(phenoData)[1] <- "Sample"
table(phenoData$Sample == colnames(methdiff)) #All TRUE

phenoData$Group <- as.character(phenoData$Group)
phenoData$Group[phenoData$Group == "Healthy Control"] <- "HC"
phenoData$Group[phenoData$Group == "Disease Control"] <- "DC"
phenoData$Group <- factor(phenoData$Group, levels=c("HC", "WD", "DC"), ordered=TRUE)

phenoData$Sex <- as.character(phenoData$Sex)
phenoData$Sex <- factor(phenoData$Sex, levels=c("Male", "Female"), ordered=TRUE)

methplot <- ggheatmap2_pheno(x = methdiff, phenoData = phenoData, low=min(methdiff), high=max(methdiff),
                       hm.colours = c("#0000FF", "#0000FF", "Black", "#FF0000", "#FF0000"), my.values = c(0,0,0.5,1,1))
pdf(file="Figures/WD Specific DMRs Heatmap pheno.pdf", width=8, height=8, onefile = FALSE)
ggheatmap.show_pheno(methplot, widths=c(0.02,0.81,0.15,0.02), heights=c(0.02,0.15,0.06,0.75,0.02))
dev.off()

# PCA ####
data <- t(as.matrix(meth[,4:dim(meth)[2]]))
rownames(data) <- c(paste("HC", 1:6,sep=""), paste("WD", 1:10, sep=""), paste("DC", 1:5, sep=""))
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
plot(data.pca, type = "l")
summary(data.pca)
# PC1 61.44%
# PC2 4.11%

# PCA Points
data <- t(as.matrix(meth[,4:dim(meth)[2]]))
rownames(data) <- c(paste("HC", 1:6,sep=""), paste("WD", 1:10, sep=""), paste("DC", 1:5, sep=""))
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 

diagnosis <- c(rep("HC", 6), rep("WD", 10), rep("DC", 5))
g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = diagnosis, ellipse = TRUE, circle = FALSE, ellipse.prob = 0.95,
              var.axes = FALSE, labels = NULL, varname.abbrev = FALSE, choices = 1:2, labels.size=6)
g + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(0.9, 0.93), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank()) +
        coord_cartesian(xlim = c(-60, 90), ylim = c(-75,75)) +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        xlab("PC1 (61% of Variance)") +
        ylab("PC2 (4% of Variance)") +
        scale_color_manual(breaks = c("HC", "WD", "DC"), 
                           values = c("HC"="#3366CC", "WD"="#FF3366", "DC"="#009933")) +
        geom_point(aes(color = diagnosis), size=3)
ggsave("Figures/WD Specific DMR PCA Plot points.png", dpi = 600, width = 8, height = 8, units = "in")

# Plot PC1 vs Fibrosis Stage
fibrosisVsPC1 <- lm(as.numeric(info$FibrosisStage) ~ data.pca$x[,"PC1"])
summary(fibrosisVsPC1) #p=2.5E-4
cor(x=data.pca$x[,"PC1"], y=as.numeric(info$FibrosisStage)) #r=0.7176
g <- ggplot()
g + 
        geom_smooth(aes(x=data.pca$x[,"PC1"], y=info$FibrosisStage), method="lm") +
        geom_point(aes(x=data.pca$x[,"PC1"], y=info$FibrosisStage, color=diagnosis), size=3) +
        annotate("text", x=-35, y=3.5, label="r = 0.72\np = 2.5E-4", hjust=0, size=7) +
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(0.2, 0.92), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank()) +
        coord_cartesian(xlim = c(-35, 55), ylim = c(-0.5,4.5)) +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        xlab("DMR Methylation PC1") +
        ylab("Fibrosis Stage") +
        scale_color_manual(breaks = c("Healthy Control", "Wilson's Disease", "Disease Control"), 
                           values = c("Healthy Control"="#3366CC", "Wilson's Disease"="#FF3366", "Disease Control"="#009933"))
ggsave("Figures/WD Specific DMR PC1 vs Fibrosis with Diagnosis.png", dpi = 600, width = 8, height = 7, units = "in")

# Covariate Association Analysis #### 
methData <- meth[,4:dim(meth)[2]]
covAssoc <- matrix(nrow = dim(methData)[1], ncol = 41)
covAssoc[,1:3] <- as.matrix(meth[,c("chr", "start", "end")])
table(colnames(methData) == as.character(info$Sequencing.ID)) #All True

for(i in 1:dim(methData)[1]){
        # Sex
        ttest_sex <- NULL
        ttest_sex <- t.test(methData[i,info$Sequencing.ID[info$Sex == "Male"]], 
                            methData[i,info$Sequencing.ID[info$Sex == "Female"]])
        covAssoc[i,4] <- ttest_sex$estimate[1] - ttest_sex$estimate[2] #meandiff
        covAssoc[i,5] <- ttest_sex$conf.int[1] #conf.int L
        covAssoc[i,6] <- ttest_sex$conf.int[2] #conf.int R
        covAssoc[i,7] <- ttest_sex$statistic #tstat
        covAssoc[i,8] <- ttest_sex$p.value #pvalue
        rm(ttest_sex)
        
        # Specimen.provider
        aov_Specimen.provider <- NULL
        tukey_Specimen.provider <- NULL
        aov_Specimen.provider <- aov(as.numeric(methData[i,]) ~ info$Specimen.provider)
        tukey_Specimen.provider <- TukeyHSD(aov_Specimen.provider)
        covAssoc[i,9] <- summary(aov_Specimen.provider)[[1]][[4]][1] #Fstat
        covAssoc[i,10] <- summary(aov_Specimen.provider)[[1]][[5]][1] #Pvalue
        covAssoc[i,11] <- tukey_Specimen.provider$`info$Specimen.provider`['Germany-CPMC','diff'] # Germany vs CPMC diff
        covAssoc[i,12] <- tukey_Specimen.provider$`info$Specimen.provider`['Germany-CPMC','p adj'] # Germany vs CPMC pvalue
        covAssoc[i,13] <- tukey_Specimen.provider$`info$Specimen.provider`['UCD GI BioBank-CPMC','diff'] # UCD vs CPMC diff
        covAssoc[i,14] <- tukey_Specimen.provider$`info$Specimen.provider`['UCD GI BioBank-CPMC','p adj'] # UCD vs CPMC pvalue
        covAssoc[i,15] <- tukey_Specimen.provider$`info$Specimen.provider`['UCD GI BioBank-Germany','diff'] # UCD vs Germany diff
        covAssoc[i,16] <- tukey_Specimen.provider$`info$Specimen.provider`['UCD GI BioBank-Germany','p adj'] # UCD vs Germany pvalue
        
        # Age
        lm_age <- NULL
        sum_age <- NULL
        lm_age <- lm(as.numeric(methData[i,]) ~ info$Age)
        sum_age <- summary(lm_age)
        covAssoc[i,17] <- sum_age$coefficients['info$Age','Estimate'] #Effect
        covAssoc[i,18] <- sum_age$coefficients['info$Age','Std. Error'] #Std Error
        covAssoc[i,19] <- sum_age$r.squared #Rsquared
        covAssoc[i,20] <- sum_age$coefficients['info$Age','t value'] #tstat
        covAssoc[i,21] <- sum_age$coefficients['info$Age', 'Pr(>|t|)'] #pvalue
        
        # BMI
        lm_BMI <- NULL
        sum_BMI <- NULL
        lm_BMI <- lm(as.numeric(methData[i,]) ~ info$BMI)
        sum_BMI <- summary(lm_BMI)
        covAssoc[i,22] <- sum_BMI$coefficients['info$BMI','Estimate'] #Effect
        covAssoc[i,23] <- sum_BMI$coefficients['info$BMI','Std. Error'] #Std Error
        covAssoc[i,24] <- sum_BMI$r.squared #Rsquared
        covAssoc[i,25] <- sum_BMI$coefficients['info$BMI','t value'] #tstat
        covAssoc[i,26] <- sum_BMI$coefficients['info$BMI', 'Pr(>|t|)'] #pvalue
        
        # FibrosisStage
        lm_FibrosisStage <- NULL
        sum_FibrosisStage <- NULL
        lm_FibrosisStage <- lm(as.numeric(methData[i,]) ~ info$FibrosisStage)
        sum_FibrosisStage <- summary(lm_FibrosisStage)
        covAssoc[i,27] <- sum_FibrosisStage$coefficients['info$FibrosisStage','Estimate'] #Effect
        covAssoc[i,28] <- sum_FibrosisStage$coefficients['info$FibrosisStage','Std. Error'] #Std Error
        covAssoc[i,29] <- sum_FibrosisStage$r.squared #Rsquared
        covAssoc[i,30] <- sum_FibrosisStage$coefficients['info$FibrosisStage','t value'] #tstat
        covAssoc[i,31] <- sum_FibrosisStage$coefficients['info$FibrosisStage', 'Pr(>|t|)'] #pvalue
        
        # Steatosis
        lm_Steatosis <- NULL
        sum_Steatosis <- NULL
        lm_Steatosis <- lm(as.numeric(methData[i,]) ~ info$Steatosis)
        sum_Steatosis <- summary(lm_Steatosis)
        covAssoc[i,32] <- sum_Steatosis$coefficients['info$Steatosis','Estimate'] #Effect
        covAssoc[i,33] <- sum_Steatosis$coefficients['info$Steatosis','Std. Error'] #Std Error
        covAssoc[i,34] <- sum_Steatosis$r.squared #Rsquared
        covAssoc[i,35] <- sum_Steatosis$coefficients['info$Steatosis','t value'] #tstat
        covAssoc[i,36] <- sum_Steatosis$coefficients['info$Steatosis', 'Pr(>|t|)'] #pvalue
        
        # InflammationGrade
        lm_InflammationGrade <- NULL
        sum_InflammationGrade <- NULL
        lm_InflammationGrade <- lm(as.numeric(methData[i,]) ~ info$InflammationGrade)
        sum_InflammationGrade <- summary(lm_InflammationGrade)
        covAssoc[i,37] <- sum_InflammationGrade$coefficients['info$InflammationGrade','Estimate'] #Effect
        covAssoc[i,38] <- sum_InflammationGrade$coefficients['info$InflammationGrade','Std. Error'] #Std Error
        covAssoc[i,39] <- sum_InflammationGrade$r.squared #Rsquared
        covAssoc[i,40] <- sum_InflammationGrade$coefficients['info$InflammationGrade','t value'] #tstat
        covAssoc[i,41] <- sum_InflammationGrade$coefficients['info$InflammationGrade', 'Pr(>|t|)'] #pvalue
}

covAssoc <- as.data.frame(covAssoc, stringsAsFactors = FALSE)
colnames(covAssoc) <- c("chr", "start", "end", 
                        "sex_estimate", "sex_confIntL", "sex_confIntR", "sex_tstat", "sex_pvalue", 
                        "Specimen.provider_fstat", "Specimen.provider_pvalue", "Specimen.provider_GerCPMC_estimate", "Specimen.provider_GerCPMC_pvalue", 
                        "Specimen.provider_UCDCPMC_estimate", "Specimen.provider_UCDCPMC_pvalue", "Specimen.provider_UCDGer_estimate", "Specimen.provider_UCDGer_pvalue",
                        "age_estimate", "age_stdError", "age_Rsquared", "age_tstat", "age_pvalue",
                        "BMI_estimate", "BMI_stdError", "BMI_Rsquared", "BMI_tstat", "BMI_pvalue",
                        "FibrosisStage_estimate", "FibrosisStage_stdError", "FibrosisStage_Rsquared", "FibrosisStage_tstat", "FibrosisStage_pvalue",
                        "Steatosis_estimate", "Steatosis_stdError", "Steatosis_Rsquared", "Steatosis_tstat", "Steatosis_pvalue",
                        "InflammationGrade_estimate", "InflammationGrade_stdError", "InflammationGrade_Rsquared", "InflammationGrade_tstat", "InflammationGrade_pvalue")

covAssoc$start <- as.integer(covAssoc$start)
covAssoc$end <- as.integer(covAssoc$end)
covAssoc[,4:41] <- apply(covAssoc[,4:41], 2, function(x) as.numeric(x))

covAssoc$sex_pAdj <- p.adjust(covAssoc$sex_pvalue, "bonf")
covAssoc$Specimen.provider_pAdj <- p.adjust(covAssoc$Specimen.provider_pvalue, "bonf")
covAssoc$Specimen.provider_GerCPMC_pAdj <- p.adjust(covAssoc$Specimen.provider_GerCPMC_pvalue, "bonf")
covAssoc$Specimen.provider_UCDCPMC_pAdj <- p.adjust(covAssoc$Specimen.provider_UCDCPMC_pvalue, "bonf")
covAssoc$Specimen.provider_UCDGer_pAdj <- p.adjust(covAssoc$Specimen.provider_UCDGer_pvalue, "bonf")
covAssoc$age_pAdj <- p.adjust(covAssoc$age_pvalue, "bonf")
covAssoc$BMI_pAdj <- p.adjust(covAssoc$BMI_pvalue, "bonf")
covAssoc$FibrosisStage_pAdj <- p.adjust(covAssoc$FibrosisStage_pvalue, "bonf")
covAssoc$Steatosis_pAdj <- p.adjust(covAssoc$Steatosis_pvalue, "bonf")
covAssoc$InflammationGrade_pAdj <- p.adjust(covAssoc$InflammationGrade_pvalue, "bonf")

table(covAssoc$sex_pAdj < 0.05) #0/1839
table(covAssoc$Specimen.provider_pAdj < 0.05) #137/1839
table(covAssoc$Specimen.provider_GerCPMC_pAdj < 0.05) #7/1839
table(covAssoc$Specimen.provider_UCDCPMC_pAdj < 0.05) #4/1839
table(covAssoc$Specimen.provider_UCDGer_pAdj < 0.05) #151/1839
table(covAssoc$age_pAdj < 0.05) #0/1839
table(covAssoc$BMI_pAdj < 0.05) #0/1839
table(covAssoc$FibrosisStage_pAdj < 0.05) #23/1839
table(covAssoc$Steatosis_pAdj < 0.05) #0/1839
table(covAssoc$InflammationGrade_pAdj < 0.05) #0/1839

covAssoc <- cbind(covAssoc, ttest[,4:21])
covAssoc <- covAssoc[,c(1:3,52:69,22:26,48,27:31,49,32:36,50,37:41,51,4:8,42,17:21,47,9,10,43,11,12,44,13,14,45,15,16,46)]
write.table(covAssoc, "WD Specific DMRs Diagnosis and Covariate Associations.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Heatmap sorted by WDvsHC -log(p-value)
pvals <- covAssoc[,c("WDvsHC_pValue", "WDvsDC_pValue", "DCvsHC_pValue", "BMI_pvalue", "FibrosisStage_pvalue", "Steatosis_pvalue",
                     "InflammationGrade_pvalue", "sex_pvalue", "age_pvalue")]
logPvals <- -log10(pvals)
logPvals <- as.matrix(logPvals)
colnames(logPvals) <- c("WDvsHC", "WDvsDC", "DCvsHC", "BMI", "FibrosisStage", "Steatosis", "InflammationGrade", "Sex", "Age")

logPvals_sort <- logPvals[order(logPvals[,"WDvsHC"]),]
plot <- ggheatmap2_sort(x = logPvals_sort, custom.label = c("WDvsHC", "WDvsDC", "DCvsHC", "BMI", "Fibrosis\nStage", "Steatosis", 
                                                            "Inflammation\nGrade", "Sex", "Age"), hm.colours = c("Black", "#FF0000"), 
                        my.values = c(0,1))
pdf(file="Figures/WD Specific Liver DMRs Covariate Association Heatmap sorted.pdf", width=6, height=7, onefile = FALSE)
ggheatmap.show_sort(plot, col.width = 0.25, row.width = 0.1)
dev.off()
