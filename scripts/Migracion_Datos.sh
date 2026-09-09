#!/bin/bash
#SBATCH --job-name=Backup_Mig_Datos
#SBATCH --partition=defq
#SBATCH --output=logs/Backup_Mig_Datos_%A_%a.out
#SBATCH --error=logs/Backup_Mig_Datos_%A_%a.err
#SBATCH --array=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=biol_alessandro@yahoo.com

mkdir -p ../../backup/Alessandro/

tar --exclude='data/results_Discrete/outputs_LL/independent_loci/neutros' \
    --exclude='data/results_Discrete/outputs_slim/independent_loci/D_FULL_seleccion_m1_*.csv' \
    -czvf ../../backup/Alessandro/data_filtered.tar.gz ./data