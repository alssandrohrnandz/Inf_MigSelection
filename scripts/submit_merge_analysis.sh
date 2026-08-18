#!/bin/bash
#SBATCH --job-name=merge_chr_R
#SBATCH --output=logs/merge_chr_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G       # Incrementamos la memoria porque R carga el cromosoma completo
#SBATCH --time=1:00:00
#SBATCH --array=1-22    # Array para lanzar los 22 cromosomas simultáneamente

# Cargar el módulo de R. 
# (Verifica con 'module spider R' o 'module avail R' la versión exacta de tu clúster)
module load r/4.1.3

# Variable del número de cromosoma provista por el Slurm Array
CHR=$SLURM_ARRAY_TASK_ID

# Crear carpeta de resultados si no existe
mkdir -p results

echo "Iniciando análisis del Cromosoma $CHR usando R..."

# Ejecutar el script de R pasando los parámetros en orden estricto
Rscript scripts/merge_and_filter_by_chr.R \
    $CHR \
    AADNA_data/Population_info.txt \
    "chr_frequencies/1240k_freq_chr{CHR}.frq.strat" \
    AADNA_data/ventajosos.txt \
    "results/unificado" \
    "AADNA_data/CADD_Data/Stats_CADD_Ancestrales_0_5.csv"

echo "Proceso del Cromosoma $CHR finalizado."