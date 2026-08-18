#!/bin/bash
#SBATCH --job-name=PDE_spatial_model
#SBATCH --output=logs/PDE_model_chr_%a.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G          # Reservamos 32GB; la grilla espacial temporal y las PDEs consumen bastante RAM
#SBATCH --time=12:00:00    # Tiempo holgado, la función ode.2D iterada toma su tiempo
#SBATCH --array=1-22       # Array del 1 al 22 para paralelizar por cromosoma

# Cargar el módulo de R correspondiente en tu clúster
module load r/4.1.3

# ========================================== #
# CONFIGURACIÓN DE VARIABLES
# ========================================== #
CHR=$SLURM_ARRAY_TASK_ID

# Cambia la palabra 'selection' por 'neutral' si vas a correr el modelo base
TIPO_ANALISIS="neutral" 
MODEL_NAME="D_FULL_${TIPO_ANALISIS}"
TASK_ID="Run_Eurasia_CHR${CHR}"
OUTPUT_DIR="results_PDE"
if [ "$TIPO_ANALISIS" == "selection" ]; then
    FREQ_FILE="results/unificado_CHR${CHR}_${TIPO_ANALISIS}.csv"
else
    FREQ_FILE="results/unificado_CHR${CHR}_${TIPO_ANALISIS}_CADD.csv"
fi
# Rutas de los archivos
# Asumiendo que el script anterior generó 'unificado_CHR1_selection.csv', etc.

POP_FILE="AADNA_data/Population_info.txt"
SNP_FILE="AADNA_data/v66.1240K.aadr.PUB.snp"

# ========================================== #
# EJECUCIÓN
# ========================================== #
# Crear directorios si no existen
mkdir -p $OUTPUT_DIR
mkdir -p logs

echo "=================================================="
echo "Iniciando modelo PDE para Cromosoma $CHR"
echo "Archivo de entrada: $FREQ_FILE"
echo "=================================================="

# Pasar los 6 argumentos posicionales exactos al script de R
Rscript scripts/spatial_pde_analysis.R \
    "$FREQ_FILE" \
    "$POP_FILE" \
    "$SNP_FILE" \
    "$TASK_ID" \
    "$MODEL_NAME" \
    "$OUTPUT_DIR"

echo "Procesamiento del Cromosoma $CHR finalizado."
