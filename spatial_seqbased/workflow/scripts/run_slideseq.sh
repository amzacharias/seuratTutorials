#!/bin/bash
#SBATCH --job-name=run_slideseq
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=16amz1@queensu.ca
#SBATCH --cpus-per-task=10
#SBATCH --mem=50MB  # Job memory request
#SBATCH --time=0-12:00:00  # Day-Hours-Minutes-Seconds
#SBATCH --output=run_slideseq.out
#SBATCH --error=run_slideseq.err
# Title: Run the slide_seq.R script
# Author: Amanda Zacharias
# Date: 2025-11-03
# Email: 16amz1@queensu.ca
#-------------------------------------------------
# Notes -------------------------------------------
#
# Code -------------------------------------------
echo Job started at "$(date +%T)"
# Options

# Dependencies
module load StdEnv/2023 r/4.5.0
# Variables

# Body
Rscript slide_seq.R

echo Job ended at "$(date +%T)"