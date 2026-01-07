#!/bin/bash
#SBATCH --job-name=CompositeLL
#SBATCH --partition=defq            
#SBATCH --output=logs/job_%A_%a.out 
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-10                
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1           
#SBATCH --mem=8G
#SBATCH --time=24:00:00

# 1. Start time
START_TIME=$(date +%s)
echo "Iniciando Job Array ID: $SLURM_ARRAY_TASK_ID en $(hostname)"
# TODO: ME QUEDE AQUÍ. HAY QUE REVISAR ESTA PARTE
# FIXME: REVISAR LAS RUTAS Y CONSULTAR EN GEMINI EL CODIGO RESTANTE PARA EDITARLO AQUI.

DIR_BASE=

# 2. Loading all files into PREFIJO
FILES_TO_PROCESS=(
    "C_FULL_seleccion_m2"
    "C_FULL_neutros_m1"
    "C_aDNA_scattered"
    "D_FULL_seleccion_m2"
    "D_FULL_neutros_m1"
    "D_aDNA_scattered"
)

echo "Procesando archivos del Job ID: ${TASK_ID} ..."

for PREFIJO in "${FILES_TO_PROCESS[@]}"; do
    INPUT_FILE=

# Cargar modulo (ajusta segun tu servidor, a veces es R/4.x.x)
module load r/4.1.3

# === 3. Crear directorios de salida ===
# Es buena practica crearlos antes, pero esto asegura que existan
mkdir -p "${DIR_BASE}/data/subsets"
mkdir -p "${DIR_BASE}/data/outputs"
mkdir -p "${DIR_BASE}/logs"

# === 4. Variables del Job actual ===
TASK_ID=$SLURM_ARRAY_TASK_ID
SUBSET_FILE="${DIR_BASE}/data/subsets/alelos_subset_${TASK_ID}.txt"
OUT_FILE="${DIR_BASE}/data/outputs/resultados_${TASK_ID}.txt"

# === 5. Avoiding extra work ===
#if [ -f "${OUT_FILE}" ]; then
#    echo "El archivo ${OUT_FILE} ya existe. Saltando este paso."
 #   exit 0
#fi

# === 6. Muestreo Aleatorio ===
# IMPORTANTE: Esto asume que mutaciones_ID.txt NO tiene encabezado.
# Si tiene encabezado, `shuf` lo mezclará como un dato más.
echo "Generando subset de 1000 alelos..."
shuf -n 1000 "${ALLELE_LIST}" > "${SUBSET_FILE}"

# === 7. Ejecutar R ===
# Pasamos argumentos: 1=FreqFile, 2=SubsetFile, 3=TaskID
echo "Ejecutando script de R..."
Rscript --vanilla "${R_SCRIPT}" "${FREQ_FILE}" "${SUBSET_FILE}" "${TASK_ID}"

# === 8. Finalizar y reportar tiempo ===
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Job ${TASK_ID} finalizado correctamente."
echo "Duracion: $DURATION segundos (~$(($DURATION / 60)) min)"
