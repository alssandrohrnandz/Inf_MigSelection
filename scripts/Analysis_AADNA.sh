#!/bin/bash
#SBATCH --job-name=Sim_Mig_Wave
#SBATCH --partition=defq
#SBATCH --output=logs/AADNA_%A_%a.out
#SBATCH --error=logs/AADNA_%A_%a.err
#SBATCH --array=1            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=biol_alessandro@yahoo.com

set -e # Detiene el script inmediatamente si ocurre algún error

#=== LOADING MODULES ===#
module load plink/1.9
module load admixtools/8.0.2

mkdir -p AADNA_data

# 1. Descarga de datos en paralelo
if [ ! -f "AADNA_data/v66.1240K.aadr.PUB.tgeno" ]; then
    echo "Archivos principales no encontrados. Descargando en paralelo..."
    
    # Descargas simultáneas en background
    curl -L -o AADNA_data/v66.1240K.aadr.PUB.snp "https://dataverse.harvard.edu/api/access/datafile/13664260" &
    curl -L -o AADNA_data/v66.1240K.aadr.PUB.ind  "https://dataverse.harvard.edu/api/access/datafile/13663698" &
    curl -L -o AADNA_data/v66.1240K.aadr.PUB.tgeno  "https://dataverse.harvard.edu/api/access/datafile/13664080" &
    curl -L -o AADNA_data/v66.1240K.aadr.PUB.anno "https://dataverse.harvard.edu/api/access/datafile/13663706" &
    
    wait # Espera a que los 4 subprocesos de curl terminen
    echo "Descarga del dataset AADR completada."
else
    echo "Los datos crudos ya están presentes en AADNA_data/."
fi

# 2. Conversión a formato PLINK
if [ ! -f "AADNA_data/1240k.ped" ]; then
    echo "Ejecutando convertf..."
    convertf -p par.PACKEDPED.PED
else
    echo "El archivo 1240k.ped ya existe en AADNA_data/."
fi

# 3. Filtrado de genotipos y cálculo de frecuencias
if [ ! -f "AADNA_data/1240k_bin.bed" ]; then
    echo "Filtrando individuos con PLINK y calculando frecuencias..."
    # Corregido: --make-bed en lugar de --maked-bed
    plink --file AADNA_data/1240k --make-bed --out AADNA_data/1240k_bin
    
    # Usamos los archivos binarios recién creados para sacar las frecuencias
    #plink --bfile AADNA_data/1240k_filtered --freq --out AADNA_data/1240k_freq
else
    echo "Los archivos filtrados (1240k_bin.bed) y las frecuencias ya existen en AADNA_data/."
fi
