#!/bin/bash
#SBATCH --job-name=Ventajosos
#SBATCH --output=logs/Ventajosos_%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1:00:00
#SBATCH --array=1  # Ejecutar para los cromosomas autómicos 1 al 22

module load admixtools/8.0.2 plink/1.9

mkdir -p chr_frequencies

echo "Procesando Cromosoma $CHR"

# 1. Extraer solo el cromosoma CHR y generar frecuencias estratificadas

plink --bfile AADNA_data/1240k_bin \
      --extract AADNA_data/CADD_Data/ventajosos.txt \
      --freq \
      --within AADNA_data/Populations_Cleaned.pop \
      --out chr_frequencies/1240k_freq_chr_ventajosos

echo "Cromosoma $CHR completado."
