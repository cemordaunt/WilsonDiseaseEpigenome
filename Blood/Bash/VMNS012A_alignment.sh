#!/bin/bash
#
#SBATCH --workdir /share/lasallelab/Charles/CM_WGBS_WD_Blood_2/
#SBATCH --mem=100000
#SBATCH --time=1-0:00
#SBATCH --partition=gc,gc64,gc128,gc256,gc512
#SBATCH --mail-type=END                     
#SBATCH --mail-user=cemordaunt@ucdavis.edu  

PATH=$PATH:/share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/
module load bowtie/1.1.1
module load samtools/0.1.19
module load sratoolkit/2.4.2-3
module load bedtools2/2.25.0
export PYTHONPATH=/share/lasallelab/pysam/lib/python2.7/site-packages/

#Complete Run in hg38 for VMNS012A
gunzip -c raw_sequences/VMNS12A*fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > raw_sequences/VMNS012A_filtered.fq
gzip raw_sequences/VMNS012A_filtered.fq
perl /share/lasallelab/programs/perl_script/adapter_split.pl raw_sequences/VMNS012A_filtered.fq.gz raw_sequences/VMNS012A_noadap.fq.gz raw_sequences/VMNS012A_withadap.fq.gz
perl /share/lasallelab/programs/perl_script/adapter_trimmer.pl raw_sequences/VMNS012A_withadap.fq.gz raw_sequences/VMNS012A_trimmed.fq.gz 45 10
mkdir VMNS012A
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 90 -m 3 -f bam -g /share/lasallelab/genomes/hg38/hg38.fa -d /share/lasallelab/genomes/hg38/BSseek2_refgen/ -i raw_sequences/VMNS012A_noadap.fq.gz -o VMNS012A/VMNS012A_noadap.bam
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 80 -m 2 -f bam -g /share/lasallelab/genomes/hg38/hg38.fa -d /share/lasallelab/genomes/hg38/BSseek2_refgen/ -i raw_sequences/VMNS012A_trimmed.fq.gz -o VMNS012A/VMNS012A_trimmed.bam
samtools sort VMNS012A/VMNS012A_noadap.bam VMNS012A/VMNS012A_noadap_sorted
samtools sort VMNS012A/VMNS012A_trimmed.bam VMNS012A/VMNS012A_trimmed_sorted
samtools merge VMNS012A/VMNS012A.bam VMNS012A/VMNS012A_noadap_sorted.bam VMNS012A/VMNS012A_trimmed_sorted.bam
samtools view VMNS012A/VMNS012A.bam > VMNS012A/VMNS012A.sam
mkdir VMNS012A/tmp
perl /share/lasallelab/programs/perl_script/SAMsorted_to_permeth.pl VMNS012A/VMNS012A.sam VMNS012A/tmp/PerMeth_VMNS012A VMNS012A hg38 CG combined 1
mkdir VMNS012A/PerMeth_VMNS012A
perl /share/lasallelab/programs/perl_script/gbcompliance.pl hg38 VMNS012A/tmp/PerMeth_VMNS012A_ VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_ VMNS012A VMNS012A
rm -r VMNS012A/tmp
mkdir VMNS012A/NoCGI_PerMeth_VMNS012A
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr1.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr1.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr2.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr2.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr3.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr3.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr4.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr4.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr5.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr5.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr6.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr6.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr7.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr7.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr8.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr8.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr9.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr9.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr10.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr10.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr11.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr11.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr12.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr12.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr13.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr13.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr14.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr14.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr15.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr15.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr16.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr16.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr17.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr17.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr18.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr18.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr19.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr19.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr20.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr20.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr21.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr21.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chr22.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chr22.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chrX.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chrX.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chrY.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chrY.bed
bedtools subtract -a VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_chrM.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMNS012A/NoCGI_PerMeth_VMNS012A/NoCGI_PerMeth_VMNS012A_chrM.bed
perl /share/lasallelab/programs/perl_script/PerMeth_to_DSSformat.pl hg38 VMNS012A/PerMeth_VMNS012A/PerMeth_VMNS012A_ DSS_files/VMNS012A_