# WD Specific Blood DMR Methylation Analysis ####
# Charles Mordaunt
# 2/13/18

# Packages ####
library(ggdendro)
library(ggplot2)
library(reshape2)
library(grid)
library(scales)
library(plyr)
library(reshape)
library(ggbiplot)

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
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, color = "black"))
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
                scale_color_manual(breaks = c("HC", "WDN", "WDH", "NAFLD", "PSC", "Male", "Female"), 
                                   values = c("HC"="#3366CC", "WDN"="#FF3366", "WDH"="#FF66CC", "NAFLD"="#009933",
                                              "PSC"="#9933CC", "Male"="#FFFF33", "Female"="#FF6633")) +
                scale_fill_manual(breaks = c("HC", "WDN", "WDH", "NAFLD", "PSC", "Male", "Female"), 
                                  values = c("HC"="#3366CC", "WDN"="#FF3366", "WDH"="#FF66CC", "NAFLD"="#009933",
                                             "PSC"="#9933CC", "Male"="#FFFF33", "Female"="#FF6633"))
        ret <- list(col=col.plot,row=row.plot,centre=centre.plot, phenoData=phenoData.plot)
        invisible(ret)
}

ggheatmap2_pheno_cont <- function(x, phenoData, hm.colours=my.colours, my.values, low, high) {
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
                scale_fill_gradientn(name="Scale", colours = hm.colours, values = my.values, na.value = "black") +
                scale_color_gradientn(name="Scale", colours = hm.colours, values = my.values, na.value = "black") 
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
                theme(plot.margin = unit(c(0,0.35,0.5,0.1), "lines"))
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

# Load Data ####
# Methylation
meth <- combineFiles(chroms=paste("chr",1:22, sep=""), prefix="Tables/WD Specific Blood DMR Blood Methylation/WD_Specific_Blood_DMR_Blood_Methyl_",
                     suffix="_DMR_methylation.txt")
meth <- meth[order(meth$chr, meth$start),]

# DMRs
DMRs <- read.delim(file="Tables/WD Specific DMR Stats with annotation.txt", sep="\t", stringsAsFactors = FALSE)
DMRs <- DMRs[order(DMRs$chr, DMRs$start, DMRs$distance),]
DMRs <- DMRs[!duplicated(DMRs[,1:3]),]

table(DMRs$start == meth$start) # All True

# Samples
samples <- read.csv("Sample Info/WD Blood 2 Sample Info.csv", header=TRUE, stringsAsFactors = FALSE)

# Heatmap ####
meth <- meth[,4:ncol(meth)]
phenoData <- samples[,c("Sequencing.ID", "Diagnosis", "Sex")]
colnames(phenoData)[1] <- "Sample"
phenoData <- phenoData[match(colnames(meth), phenoData$Sample),]
table(colnames(meth) == phenoData$Sample) # All TRUE
phenoData$Sample <- factor(phenoData$Sample, levels=unique(phenoData$Sample), ordered=TRUE)
phenoData$Diagnosis[phenoData$Diagnosis == "Healthy Control"] <- "HC"
phenoData$Diagnosis[phenoData$Diagnosis == "Wilson's Disease"] <- "WD"
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels=c("HC", "WD", "NAFLD", "PSC"), ordered=TRUE)
phenoData$Sex <- factor(phenoData$Sex, levels=c("Male", "Female"), ordered=TRUE)

methavg <- rowMeans(meth, na.rm = TRUE)
methdiff <- meth - methavg
methdiff <- as.matrix(methdiff)
methplot <- ggheatmap2_pheno(x = methdiff, phenoData=phenoData, hm.colours = c("#0000FF", "#0000FF", "Black", "#FF0000", "#FF0000"), 
                             my.values = c(0,0,0.5,1,1), low = -0.5, high = 0.5)
pdf(file="Figures/WD Blood Specific DMRs All Samples Heatmap phenoData.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show_pheno(methplot)
dev.off()

# PCA Plot ####
data <- t(as.matrix(meth))
diagnosis <- c(rep("HC", 12), rep("WD", 40), rep("NAFLD", 10), 
               rep("PSC", 10))
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
plot(data.pca, type = "l")
summary(data.pca)
# PC1 20.3%
# PC2 3.7%

g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = diagnosis, ellipse = TRUE, circle = FALSE, 
              var.axes = FALSE, varname.abbrev = FALSE, choices = 1:2,ellipse.prob = 0.95)
g + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'vertical', legend.position = c(0.88, 0.91), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), legend.text=element_text(size=20)) +
        coord_cartesian(xlim = c(-17, 16), ylim = c(-17,16)) +
        xlab("PC1 (20% of Variance)") +
        ylab("PC2 (4% of Variance)") +
        scale_color_manual(breaks = c("HC", "WD", "NAFLD", "PSC"), 
                           values = c("HC"="#3366CC", "WD"="#FF3366", "NAFLD"="#009933",
                           "PSC"="#9933CC")) +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        geom_point(aes(color = diagnosis), size=2.5)
ggsave("Figures/WD Specific DMRs All Samples PCA plot Specific DC.png", dpi = 600, width = 8, height = 8, units = "in")

# Compare covariates to DMR methylation ####
samples$Group <- factor(samples$Group, levels=c("Healthy Control", "Disease Control", "Wilson's Disease"), ordered=TRUE)
samples$Diagnosis <- factor(samples$Diagnosis, levels=c("Healthy Control", "NAFLD", "PSC", "Wilson's Disease"), ordered=TRUE)
samples$Sex <- factor(samples$Sex, levels=c("Male", "Female"), ordered=TRUE)
samples <- samples[,c("Sequencing.ID", "Group", "Diagnosis", "Sex", "Age", "Age.at.onset", "Height..cm.", "Weight..kg.", "BMI", 
                      "Neurologic.form", "Hepatic.form")]
colnames(samples) <- c("Sample", "Group", "Diagnosis", "Sex", "Age", "OnsetAge", "Height", "Weight", "BMI", "NeuroWD", "HepaticWD")
meth <- combineFiles(chroms=paste("chr",1:22, sep=""), prefix="Tables/WD Specific Blood DMR Blood Methylation/WD_Specific_Blood_DMR_Blood_Methyl_",
                     suffix="_DMR_methylation.txt")
meth <- meth[order(meth$chr, meth$start),]
methData <- meth[,4:dim(meth)[2]]
covAssoc <- matrix(nrow = dim(methData)[1], ncol = 88)
covAssoc[,1:3] <- as.matrix(meth[,c("chr", "start", "end")])
samples <- samples[match(colnames(methData), samples$Sample),]
table(colnames(methData) == as.character(samples$Sample)) #All True

for(i in 1:dim(methData)[1]){
        # WDvsHC
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$Group == "Wilson's Disease"]], 
                        methData[i,samples$Sample[samples$Group == "Healthy Control"]])
        covAssoc[i,4] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,5] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,6] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,7] <- ttest$statistic #tstat
        covAssoc[i,8] <- ttest$p.value #pvalue
        
        # WDvsDC
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$Group == "Wilson's Disease"]], 
                        methData[i,samples$Sample[samples$Group == "Disease Control"]])
        covAssoc[i,9] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,10] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,11] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,12] <- ttest$statistic #tstat
        covAssoc[i,13] <- ttest$p.value #pvalue
        
        # DCvsHC
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$Group == "Disease Control"]], 
                        methData[i,samples$Sample[samples$Group == "Healthy Control"]])
        covAssoc[i,14] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,15] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,16] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,17] <- ttest$statistic #tstat
        covAssoc[i,18] <- ttest$p.value #pvalue
        
        # WDvsNAFLD
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$Diagnosis == "Wilson's Disease"]], 
                        methData[i,samples$Sample[samples$Diagnosis == "NAFLD"]])
        covAssoc[i,19] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,20] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,21] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,22] <- ttest$statistic #tstat
        covAssoc[i,23] <- ttest$p.value #pvalue
        
        # WDvsPSC
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$Diagnosis == "Wilson's Disease"]], 
                        methData[i,samples$Sample[samples$Diagnosis == "PSC"]])
        covAssoc[i,24] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,25] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,26] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,27] <- ttest$statistic #tstat
        covAssoc[i,28] <- ttest$p.value #pvalue
        
        # NAFLDvsHC
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$Diagnosis == "NAFLD"]], 
                        methData[i,samples$Sample[samples$Diagnosis == "Healthy Control"]])
        covAssoc[i,29] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,30] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,31] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,32] <- ttest$statistic #tstat
        covAssoc[i,33] <- ttest$p.value #pvalue
        
        # PSCvsHC
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$Diagnosis == "PSC"]], 
                        methData[i,samples$Sample[samples$Diagnosis == "Healthy Control"]])
        covAssoc[i,34] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,35] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,36] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,37] <- ttest$statistic #tstat
        covAssoc[i,38] <- ttest$p.value #pvalue
        
        # PSCvsNAFLD
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$Diagnosis == "PSC"]], 
                        methData[i,samples$Sample[samples$Diagnosis == "NAFLD"]])
        covAssoc[i,39] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,40] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,41] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,42] <- ttest$statistic #tstat
        covAssoc[i,43] <- ttest$p.value #pvalue
        
        # FvsM
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$Sex == "Female"]], 
                        methData[i,samples$Sample[samples$Sex == "Male"]])
        covAssoc[i,44] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,45] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,46] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,47] <- ttest$statistic #tstat
        covAssoc[i,48] <- ttest$p.value #pvalue
        
        # Age
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(methData[i,]) ~ samples$Age)
        sum <- summary(lm)
        covAssoc[i,49] <- sum$coefficients['samples$Age','Estimate'] #Effect
        covAssoc[i,50] <- sum$coefficients['samples$Age','Std. Error'] #Std Error
        covAssoc[i,51] <- sum$r.squared #Rsquared
        covAssoc[i,52] <- sum$coefficients['samples$Age','t value'] #tstat
        covAssoc[i,53] <- sum$coefficients['samples$Age', 'Pr(>|t|)'] #pvalue
        
        # OnsetAge
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(methData[i,]) ~ samples$OnsetAge)
        sum <- summary(lm)
        covAssoc[i,54] <- sum$coefficients['samples$OnsetAge','Estimate'] #Effect
        covAssoc[i,55] <- sum$coefficients['samples$OnsetAge','Std. Error'] #Std Error
        covAssoc[i,56] <- sum$r.squared #Rsquared
        covAssoc[i,57] <- sum$coefficients['samples$OnsetAge','t value'] #tstat
        covAssoc[i,58] <- sum$coefficients['samples$OnsetAge', 'Pr(>|t|)'] #pvalue
        
        # Height
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(methData[i,]) ~ samples$Height)
        sum <- summary(lm)
        covAssoc[i,59] <- sum$coefficients['samples$Height','Estimate'] #Effect
        covAssoc[i,60] <- sum$coefficients['samples$Height','Std. Error'] #Std Error
        covAssoc[i,61] <- sum$r.squared #Rsquared
        covAssoc[i,62] <- sum$coefficients['samples$Height','t value'] #tstat
        covAssoc[i,63] <- sum$coefficients['samples$Height', 'Pr(>|t|)'] #pvalue
        
        # Weight
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(methData[i,]) ~ samples$Weight)
        sum <- summary(lm)
        covAssoc[i,64] <- sum$coefficients['samples$Weight','Estimate'] #Effect
        covAssoc[i,65] <- sum$coefficients['samples$Weight','Std. Error'] #Std Error
        covAssoc[i,66] <- sum$r.squared #Rsquared
        covAssoc[i,67] <- sum$coefficients['samples$Weight','t value'] #tstat
        covAssoc[i,68] <- sum$coefficients['samples$Weight', 'Pr(>|t|)'] #pvalue
        
        # BMI
        lm <- NULL
        sum <- NULL
        lm <- lm(as.numeric(methData[i,]) ~ samples$BMI)
        sum <- summary(lm)
        covAssoc[i,69] <- sum$coefficients['samples$BMI','Estimate'] #Effect
        covAssoc[i,70] <- sum$coefficients['samples$BMI','Std. Error'] #Std Error
        covAssoc[i,71] <- sum$r.squared #Rsquared
        covAssoc[i,72] <- sum$coefficients['samples$BMI','t value'] #tstat
        covAssoc[i,73] <- sum$coefficients['samples$BMI', 'Pr(>|t|)'] #pvalue
        
        # Neuro vs Hepatic
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$NeuroWD == 1 & !is.na(samples$NeuroWD)]], 
                        methData[i,samples$Sample[samples$HepaticWD == 1 & !is.na(samples$HepaticWD)]])
        covAssoc[i,74] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,75] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,76] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,77] <- ttest$statistic #tstat
        covAssoc[i,78] <- ttest$p.value #pvalue
        
        # Neuro vs HC
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$NeuroWD == 1 & !is.na(samples$NeuroWD)]], 
                        methData[i,samples$Sample[samples$Group == "Healthy Control"]])
        covAssoc[i,79] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,80] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,81] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,82] <- ttest$statistic #tstat
        covAssoc[i,83] <- ttest$p.value #pvalue
        
        # Hepatic vs HC
        ttest <- NULL
        ttest <- t.test(methData[i,samples$Sample[samples$HepaticWD == 1 & !is.na(samples$HepaticWD)]], 
                        methData[i,samples$Sample[samples$Group == "Healthy Control"]])
        covAssoc[i,84] <- ttest$estimate[1] - ttest$estimate[2] #meandiff
        covAssoc[i,85] <- ttest$conf.int[1] #conf.int L
        covAssoc[i,86] <- ttest$conf.int[2] #conf.int R
        covAssoc[i,87] <- ttest$statistic #tstat
        covAssoc[i,88] <- ttest$p.value #pvalue
}

covAssoc <- as.data.frame(covAssoc, stringsAsFactors = FALSE)
newCols <- c("chr", "start", "end", paste(rep(c("WDvsHC", "WDvsDC", "DCvsHC", "WDvsNAFLD", "WDvsPSC", "NAFLDvsHC", "PSCvsHC", "PSCvsNAFLD", "FvsM"), each=5),
                                          c("meanDiff", "confIntL", "confIntR", "tstat", "pvalue"), sep="_"),
             paste(rep(c("Age", "OnsetAge", "Height", "Weight", "BMI"), each=5),
                   c("Effect", "StdErr", "Rsquared", "tstat", "pvalue"), sep="_"),
             paste(rep(c("WDNvsWDH", "WDNvsHC", "WDHvsHC"), each=5),
                   c("meanDiff", "confIntL", "confIntR", "tstat", "pvalue"), sep="_"))
colnames(covAssoc) <- newCols
covAssoc$start <- as.integer(covAssoc$start)
covAssoc$end <- as.integer(covAssoc$end)
covAssoc[,4:dim(covAssoc)[2]] <- apply(covAssoc[,4:dim(covAssoc)[2]], 2, function(x) as.numeric(x))

# Adj p-values
covAssoc$WDvsHC_pAdj <- p.adjust(covAssoc$WDvsHC_pvalue, "bonf")
covAssoc$WDvsDC_pAdj <- p.adjust(covAssoc$WDvsDC_pvalue, "bonf")
covAssoc$DCvsHC_pAdj <- p.adjust(covAssoc$DCvsHC_pvalue, "bonf")
covAssoc$WDvsNAFLD_pAdj <- p.adjust(covAssoc$WDvsNAFLD_pvalue, "bonf")
covAssoc$WDvsPSC_pAdj <- p.adjust(covAssoc$WDvsPSC_pvalue, "bonf")
covAssoc$NAFLDvsHC_pAdj <- p.adjust(covAssoc$NAFLDvsHC_pvalue, "bonf")
covAssoc$PSCvsHC_pAdj <- p.adjust(covAssoc$PSCvsHC_pvalue, "bonf")
covAssoc$PSCvsNAFLD_pAdj <- p.adjust(covAssoc$PSCvsNAFLD_pvalue, "bonf")
covAssoc$FvsM_pAdj <- p.adjust(covAssoc$FvsM_pvalue, "bonf")
covAssoc$Age_pAdj <- p.adjust(covAssoc$Age_pvalue, "bonf")
covAssoc$OnsetAge_pAdj <- p.adjust(covAssoc$OnsetAge_pvalue, "bonf")
covAssoc$Height_pAdj <- p.adjust(covAssoc$Height_pvalue, "bonf")
covAssoc$Weight_pAdj <- p.adjust(covAssoc$Weight_pvalue, "bonf")
covAssoc$BMI_pAdj <- p.adjust(covAssoc$BMI_pvalue, "bonf")
covAssoc$WDNvsWDH_pAdj <- p.adjust(covAssoc$WDNvsWDH_pvalue, "bonf")
covAssoc$WDNvsHC_pAdj <- p.adjust(covAssoc$WDNvsHC_pvalue, "bonf")
covAssoc$WDHvsHC_pAdj <- p.adjust(covAssoc$WDHvsHC_pvalue, "bonf")

write.table(covAssoc, file="Tables/WD Specific Blood DMRs Covariate Associations.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Significant DMRs
table(covAssoc$WDvsHC_pAdj < 0.05)      # 3
table(covAssoc$WDvsDC_pAdj < 0.05)      # 9
table(covAssoc$DCvsHC_pAdj < 0.05)      # 0
table(covAssoc$WDvsNAFLD_pAdj < 0.05)   # 3
table(covAssoc$WDvsPSC_pAdj < 0.05)     # 3
table(covAssoc$NAFLDvsHC_pAdj < 0.05)   # 0
table(covAssoc$PSCvsHC_pAdj < 0.05)     # 0
table(covAssoc$PSCvsNAFLD_pAdj < 0.05)  # 0
table(covAssoc$FvsM_pAdj < 0.05)        # 0
table(covAssoc$Age_pAdj < 0.05)         # 1
table(covAssoc$OnsetAge_pAdj < 0.05)    # 0
table(covAssoc$Height_pAdj < 0.05)      # 0
table(covAssoc$Weight_pAdj < 0.05)      # 0
table(covAssoc$BMI_pAdj < 0.05)         # 0
table(covAssoc$WDNvsWDH_pAdj < 0.05)    # 0
table(covAssoc$WDNvsHC_pAdj < 0.05)     # 2
table(covAssoc$WDHvsHC_pAdj < 0.05)     # 0

# Heatmap Sorted by WDvsHC logp-value ####
pvals <- covAssoc[,c("WDvsHC_pvalue", "WDvsDC_pvalue", "DCvsHC_pvalue", "WDvsNAFLD_pvalue", "WDvsPSC_pvalue", "NAFLDvsHC_pvalue",
                     "PSCvsHC_pvalue", "PSCvsNAFLD_pvalue", "FvsM_pvalue", "Age_pvalue", "OnsetAge_pvalue", "Height_pvalue", 
                     "Weight_pvalue", "BMI_pvalue", "WDNvsWDH_pvalue", "WDNvsHC_pvalue", "WDHvsHC_pvalue")]
logPvals <- -log10(pvals)
logPvals <- as.matrix(logPvals)
newCols <- c("WDvsHC", "WDvsDC", "DCvsHC", "WDvsNAFLD", "WDvsPSC", "NAFLDvsHC",
                        "PSCvsHC", "PSCvsNAFLD", "FvsM", "Age", "OnsetAge", "Height", 
                        "Weight", "BMI", "WDNvsWDH", "WDNvsHC", "WDHvsHC")
colnames(logPvals) <- newCols
logPvals_sort <- logPvals[order(logPvals[,"WDvsHC"]),]
plot <- ggheatmap2_sort(x = logPvals_sort, custom.label = newCols, hm.colours = c("Black", "#FF0000"), my.values = c(0,1))
pdf(file="Figures/WD Blood DMRs Covariate Association Heatmap sorted.pdf", width=6, height=7, onefile = FALSE)
ggheatmap.show_sort(plot, col.width = 0.25, row.width = 0.1)
dev.off()
