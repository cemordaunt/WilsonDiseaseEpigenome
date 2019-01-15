# GREAT Regulatory Domains hg38 ####
# Charles Mordaunt
# 8/23/18

# Packages ####
library(dplyr)
library(GenomicRanges)

# Data ####
genes <- read.delim("Tables/ncbiRefSeqCurated.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)
chrom_size <- read.delim("Tables/chromsizes_hg38.txt", sep="\t", stringsAsFactors = FALSE, header=FALSE)
colnames(chrom_size) <- c("chr", "width")

# Get Gene Regulatory Domains ####
# Download Curated RefSeq genes from UCSC Table Browser
# Columns: chrom, strand, txStart, txEnd, name2, locusLinkId
# Filter for one isoform for each gene, take isoform with most upstream TSS and the longest, by strand
# Assign proximal domains: 5kb upstream and 1kb downstream of TSS, can be overlapping proximal domains of other genes, by strand
# Assign distal domains: 1Mb upstream and 1Mb downstream, not overlapping with proximal domains from other genes, by strand

# Format table ####
colnames(genes) <- c("chr", "strand", "start", "end", "name", "entrezID")
genes <- genes[,c("chr", "start", "end", "strand", "name", "entrezID")]
sapply(genes, function(x) table(is.na(x))) # All FALSE
chroms <- paste("chr", c(1:22, "X", "Y", "M"), sep="")
genes <- subset(genes, chr %in% chroms) %>% unique
genes <- genes[order(genes$chr, genes$start, genes$end),]

# Filter unique isoforms ####
gene_list <- genes$name %>% unlist %>% as.character %>% unique
genes_ftr <- NULL
for(i in 1:length(gene_list)){
        temp <- NULL
        temp <- subset(genes, name == gene_list[i])
        if("chrX" %in% temp$chr & "chrY" %in% temp$chr){ # Accounts for 26 homologs on X and Y, filter separately
                tempX <- NULL
                tempX <- subset(temp, chr == "chrX")
                tempY <- NULL
                tempY <- subset(temp, chr == "chrY")
                if(tempX$strand[1] == "+"){
                        tempX <- subset(tempX, start == min(tempX$start)) 
                        tempX <- subset(tempX, end == max(tempX$end))
                } else {
                        tempX <- subset(tempX, end == max(tempX$end)) 
                        tempX <- subset(tempX, start == min(tempX$start))
                }
                if(tempY$strand[1] == "+"){
                        tempY <- subset(tempY, start == min(tempY$start)) 
                        tempY <- subset(tempY, end == max(tempY$end))
                } else {
                        tempY <- subset(tempY, end == max(tempY$end)) 
                        tempY <- subset(tempY, start == min(tempY$start))
                }
                genes_ftr <- rbind(genes_ftr, tempX, tempY)
        } else {
                if(temp$strand[1] == "+"){
                        temp <- subset(temp, start == min(temp$start)) 
                        temp <- subset(temp, end == max(temp$end))
                } else {
                        temp <- subset(temp, end == max(temp$end)) 
                        temp <- subset(temp, start == min(temp$start))
                }
                genes_ftr <- rbind(genes_ftr, temp)
        }
}
rm(temp, tempX, tempY, i)

# Assign proximal domains, -5kb to +1kb, different for each strand ####
proximal <- NULL
for(i in 1:nrow(genes_ftr)){
        temp <- NULL
        temp <- genes_ftr[i,]
        chrom_width <- NULL
        chrom_width <- chrom_size$width[chrom_size$chr == temp$chr]
        if(temp$strand == "+"){
                temp$proximal_start <- max(temp$start - 5000, 0)
                temp$proximal_end <- min(temp$start + 1000, chrom_width)
        } else {
                temp$proximal_start <- max(temp$end - 1000, 0)
                temp$proximal_end <- min(temp$end + 5000, chrom_width)
        }
        proximal <- rbind(proximal, temp)
}
rm(temp, i, chrom_width)

# Assign distal domains, -1Mb to +1Mb can't overlap with adjacent proximal domains ####
proximal <- proximal[order(proximal$chr, proximal$proximal_start, proximal$proximal_end),]
distal <- NULL
for(i in 1:nrow(proximal)){
        if(i %% 500 == 0){cat(i, "\t")}
        temp <- NULL
        prev_region <- NULL            
        next_region <- NULL
        chrom_width <- NULL
        temp <- proximal[i,]
        chrom_width <- chrom_size$width[chrom_size$chr == temp$chr]
        if(i == 1){
                next_region <- proximal[i+1,]
                if(temp$strand == "+"){
                        temp$distal_start <- min(max(temp$start - 1000000, 0), temp$proximal_start)
                        temp$distal_end <- max(min(temp$start + 1000000, chrom_width, next_region$proximal_start), temp$proximal_end)
                }
                if(temp$strand == "-"){
                        temp$distal_start <- min(max(temp$end - 1000000, 0), temp$proximal_start)
                        temp$distal_end <- max(min(temp$end + 1000000, chrom_width, next_region$proximal_start), temp$proximal_end)
                }
        }
        if(i == nrow(proximal)){
                prev_region <- proximal[i-1,]
                if(temp$strand == "+"){
                        temp$distal_start <- min(max(temp$start - 1000000, 0, prev_region$proximal_end), temp$proximal_start)
                        temp$distal_end <- max(min(temp$start + 1000000, chrom_width), temp$proximal_end)
                }
                if(temp$strand == "-"){
                        temp$distal_start <- min(max(temp$end - 1000000, 0, prev_region$proximal_end), temp$proximal_start)
                        temp$distal_end <- max(min(temp$end + 1000000, chrom_width), temp$proximal_end)
                }
        }
        if(i > 1 & i < nrow(proximal)){
                prev_region <- proximal[i-1,]
                next_region <- proximal[i+1,]
                if(!temp$chr == prev_region$chr){
                        if(temp$strand == "+"){
                                temp$distal_start <- min(max(temp$start - 1000000, 0), temp$proximal_start)
                                temp$distal_end <- max(min(temp$start + 1000000, chrom_width, next_region$proximal_start), temp$proximal_end)
                        }
                        if(temp$strand == "-"){
                                temp$distal_start <- min(max(temp$end - 1000000, 0), temp$proximal_start)
                                temp$distal_end <- max(min(temp$end + 1000000, chrom_width, next_region$proximal_start), temp$proximal_end)
                        }
                }
                if(!temp$chr == next_region$chr){
                        if(temp$strand == "+"){
                                temp$distal_start <- min(max(temp$start - 1000000, 0, prev_region$proximal_end), temp$proximal_start)
                                temp$distal_end <- max(min(temp$start + 1000000, chrom_width), temp$proximal_end)
                        }
                        if(temp$strand == "-"){
                                temp$distal_start <- min(max(temp$end - 1000000, 0, prev_region$proximal_end), temp$proximal_start)
                                temp$distal_end <- max(min(temp$end + 1000000, chrom_width), temp$proximal_end)
                        }
                }
                if(temp$chr == prev_region$chr & temp$chr == next_region$chr){
                        if(temp$strand == "+"){
                                temp$distal_start <- min(max(temp$start - 1000000, 0, prev_region$proximal_end), temp$proximal_start)
                                temp$distal_end <- max(min(temp$start + 1000000, chrom_width, next_region$proximal_start), temp$proximal_end)
                        }
                        if(temp$strand == "-"){
                                temp$distal_start <- min(max(temp$end - 1000000, 0, prev_region$proximal_end), temp$proximal_start)
                                temp$distal_end <- max(min(temp$end + 1000000, chrom_width, next_region$proximal_start), temp$proximal_end)
                        }
                }
        }
        distal <- rbind(distal, temp)
}
regDomains <- distal
colnames(regDomains) <- c("gene_chr", "gene_start", "gene_end", "gene_strand", "gene_name", "gene_entrezID", "proximal_start",
                          "proximal_end", "distal_start", "distal_end")
write.table(regDomains, "Tables/Regulatory domains hg38.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

bed <- distal[,c("chr", "distal_start", "distal_end", "name")]
bed$score <- rep(0, nrow(bed))
bed <- cbind(bed, distal[,c("strand", "proximal_start", "proximal_end")])
write.table(bed, "Tables/Regulatory domains hg38.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
