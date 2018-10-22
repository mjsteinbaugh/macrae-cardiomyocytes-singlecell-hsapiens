#!/bin/bash

#SBATCH --job-name=DE_lrt
#SBATCH --partition=medium             # Partition name
#SBATCH --time=1-23:59                 # Runtime in D-HH:MM format
#SBATCH --nodes=1                      # Number of nodes (keep at 1)
#SBATCH --ntasks=1                     # Number of tasks per node (keep at 1)
#SBATCH --cpus-per-task=8              # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem-per-cpu=16G               # Memory needed per CPU
#SBATCH --error=jobid_%j.err           # File to which STDERR will be written, including job ID
#SBATCH --output=jobid_%j.out          # File to which STDOUT will be written, including job ID
#SBATCH --mail-type=ALL                # Type of email notification (BEGIN, END, FAIL, ALL)

module load gcc/6.2.0
module load gsl/2.3
module load hdf5
module load R/3.5.1

Rscript /home/vb83/PIs/calum_macrae/macrae_scRNASeq_cardiomyocytes/cca_analysis/DE_res_0.6_by_cluster_lrt_variation_cell_cycle.R
