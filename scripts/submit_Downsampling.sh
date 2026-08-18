#!/bin/bash
#SBATCH --job-name=Downsampling_AADNA
#SBATCH --partition=defq
#SBATCH --output=logs/Downsampling_%A_%a.out
#SBATCH --error=logs/Downsampling_%A_%a.err
#SBATCH --array=1-280            
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00

# =========================================================================
# CONFIGURACIÓN DEL ENTORNO
# =========================================================================

# Cargar el módulo de R (Ajusta la versión según lo que uses en el clúster)
module load r/4.1.3  

# =========================================================================
# DEFINICIÓN DE RUTAS Y VARIABLES GLOBALES
# =========================================================================

# Directorio base del proyecto
BASE_DIR="/mnt/data/dortega/hlopezh/Inf_MigSelection"

# Directorios de entrada y salida
SLIM_DIR="${BASE_DIR}/data/results_Discrete/outputs_slim"
OUT_DIR="${BASE_DIR}/results_validations"
R_SCRIPT="scripts/DownSampling.R"
EMPIRICAL_FILE="${BASE_DIR}/Eurasia_Populations_Periods_Grid.csv"

# Crear directorios si no existen (solo el nodo maestro de la tarea lo hará)
mkdir -p ${OUT_DIR}
mkdir -p logs

# =========================================================================
# PARÁMETROS ESPECÍFICOS DE LA TAREA (JOB ARRAY)
# =========================================================================

TASK_ID=$SLURM_ARRAY_TASK_ID

# Define el tipo de modelo que vas a correr (Cambiar a "D_FULL_selection" cuando toque)
MODEL_NAME="D_FULL_seleccion"

# Construcción dinámica del nombre del archivo SLiM basado en el TASK_ID
SLIM_FILE="${SLIM_DIR}/v2_${MODEL_NAME}_m2_${TASK_ID}.csv"

# =========================================================================
# EJECUCIÓN
# =========================================================================

echo "=========================================================="
echo " Iniciando Job Array ID : $SLURM_ARRAY_JOB_ID"
echo " Task ID Actual         : $TASK_ID"
echo " Modelo a evaluar       : $MODEL_NAME"
echo " Archivo SLiM objetivo  : $SLIM_FILE"
echo "=========================================================="

# Verificamos que el archivo de SLiM exista antes de lanzar R para evitar fallos fantasma
if [ -f "$SLIM_FILE" ]; then
    
    Rscript ${R_SCRIPT} \
        "${TASK_ID}" \
        "${MODEL_NAME}" \
        "${SLIM_FILE}" \
        "${EMPIRICAL_FILE}" \
        "${OUT_DIR}"
        
    echo "Task $TASK_ID completada con éxito."
    
else
    echo "ERROR Crítico: No se encontró el archivo $SLIM_FILE"
    exit 1
fi