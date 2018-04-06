#!/bin/bash
#
#SBATCH --workdir /share/lasallelab/Charles/CM_WGBS_WD_Blood_2
#SBATCH --mem=100000						# total memory
#SBATCH --time=1-0							# time (day-hr)
#SBATCH -n 2								# cores
#SBATCH --mail-type=END                     # notifications for job done & fail
#SBATCH --mail-user=cemordaunt@ucdavis.edu  # send-to address
#SBATCH --partition=gc,gc64,gc128,gc256,gc512
