#!/bin/bash
#
#SBATCH --workdir /share/lasallelab/Dorothy/DK_WGBS_HumanLiver/
#SBATCH --mem=60000					# memory per CPU
#SBATCH --mail-type=END                     # notifications for job done & fail
#SBATCH --mail-user=dakieffer@ucdavis.edu  # send-to address
#SBATCH --time=2-0

PATH=$PATH:/share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/
module load bowtie/1.1.1
module load samtools/0.1.19
module load sratoolkit/2.4.2-3
module load bedtools2/2.25.0
export PYTHONPATH=/share/lasallelab/pysam/lib/python2.7/site-packages/

#Complete Run in hg38 for VMDK1A
gunzip -c raw_sequences/VMDK1A_S2_L002_R1_001.fastq.gz | grep -A 3 '^@.* [^:]*:N:[^:]*:' |   grep -v "^--$" > raw_sequences/VMDK1A_filtered.fq
gzip raw_sequences/VMDK1A_filtered.fq
perl /share/lasallelab/programs/perl_script/adapter_split.pl raw_sequences/VMDK1A_filtered.fq.gz raw_sequences/VMDK1A_noadap.fq.gz raw_sequences/VMDK1A_withadap.fq.gz
perl /share/lasallelab/programs/perl_script/adapter_trimmer.pl raw_sequences/VMDK1A_withadap.fq.gz raw_sequences/VMDK1A_trimmed.fq.gz 45 10
mkdir VMDK1A
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 90 -m 3 -f bam -g /share/lasallelab/genomes/hg38/hg38.fa -d /share/lasallelab/genomes/hg38/BSseek2_refgen/ -i raw_sequences/VMDK1A_noadap.fq.gz -o VMDK1A/VMDK1A_noadap.bam
python /share/lasallelab/programs/tuxedo/BSseeker2-master_v2.0.8/bs_seeker2-align.py --bt-p 12 -e 80 -m 2 -f bam -g /share/lasallelab/genomes/hg38/hg38.fa -d /share/lasallelab/genomes/hg38/BSseek2_refgen/ -i raw_sequences/VMDK1A_trimmed.fq.gz -o VMDK1A/VMDK1A_trimmed.bam
samtools sort VMDK1A/VMDK1A_noadap.bam VMDK1A/VMDK1A_noadap_sorted
samtools sort VMDK1A/VMDK1A_trimmed.bam VMDK1A/VMDK1A_trimmed_sorted
samtools merge VMDK1A/VMDK1A.bam VMDK1A/VMDK1A_noadap_sorted.bam VMDK1A/VMDK1A_trimmed_sorted.bam
samtools view VMDK1A/VMDK1A.bam > VMDK1A/VMDK1A.sam
mkdir VMDK1A/tmp
perl /share/lasallelab/programs/perl_script/SAMsorted_to_permeth.pl VMDK1A/VMDK1A.sam VMDK1A/tmp/PerMeth_VMDK1A VMDK1A hg38 CG combined 1
mkdir VMDK1A/PerMeth_VMDK1A
perl /share/lasallelab/programs/perl_script/gbcompliance.pl hg38 VMDK1A/tmp/PerMeth_VMDK1A_ VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_ VMDK1A VMDK1A
rm -r VMDK1A/tmp
mkdir VMDK1A/NoCGI_Permeth_VMDK1A
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr1.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr1.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr2.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr2.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr3.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr3.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr4.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr4.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr5.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr5.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr6.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr6.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr7.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr7.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr8.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr8.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr9.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr9.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr10.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr10.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr11.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr11.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr12.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr12.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr13.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr13.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr14.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr14.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr15.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr15.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr16.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr16.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr17.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr17.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr18.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr18.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr19.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr19.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr20.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr20.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr21.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr21.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chr22.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chr22.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chrX.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chrX.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chrY.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chrY.bed
bedtools subtract -a VMDK1A/PerMeth_VMDK1A/PerMeth_VMDK1A_chrM.bed -b /share/lasallelab/genomes/hg38/GTF/hg38_genome_CGI.bed > VMDK1A/NoCGI_Permeth_VMDK1A/NoCGI_Permeth_VMDK1A_chrM.bed
