#!/bin/bash
#SBATCH --job-name=Global_Sel
#SBATCH --partition=defq
#SBATCH --output=logs/Global_Plot_%A_%a.out
#SBATCH --error=logs/Global_Plot_%A_%a.err
#SBATCH --array=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00

module load r/4.1.3
Rscript --vanilla /mnt/data/dortega/hlopezh/Inf_MigSelection/scripts/Discrete_Space_Inference/Plot_Global_Sel_D_Log_LL.R