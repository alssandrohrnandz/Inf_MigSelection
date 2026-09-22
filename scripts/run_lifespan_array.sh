#!/bin/bash
#SBATCH --job-name=lifespan_arr
#SBATCH --output=logs/lifespan_%A_%a.out
#SBATCH --error=logs/lifespan_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --partition=normal
#SBATCH --array=1-20            # 20 chunks en paralelo
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tu@correo

module purge
module load r/4.1.3
mkdir -p logs

Rscript procesar_lifespan.R "$SLURM_ARRAY_TASK_ID" 20

echo "[$(date)] Chunk $SLURM_ARRAY_TASK_ID terminado"