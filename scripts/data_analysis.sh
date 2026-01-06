#!/bin/bash
#SBATCH --job-name=CompositeLL
#SBATCH --partition=defq            # Asegurate que esta particion existe (la vimos en tu sinfo)
#SBATCH --output=logs/job_%A_%a.out # %A=ID Master, %a=ID Array (ej. job_1234_1.out)
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-10                # 10 trabajos en paralelo
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1           # R rara vez usa mas de 1 CPU a menos que paralelos internamente
#SBATCH --mem=8G
#SBATCH --time=24:00:00

# === 1. Marcar tiempo de inicio ===
START_TIME=$(date +%s)
echo "Iniciando Job Array ID: $SLURM_ARRAY_TASK_ID en $(hostname)"

DIR_BASE="/mnt/data/dortega/hlopezh/Inf_MigSelection/"

ALLELE_LIST="${DIR_BASE}/data/mutaciones_ID.txt"  # Asumiendo que moviste esto a data/
FREQ_FILE="${DIR_BASE}/data/frecuencias_mutaciones_m1_10_600Gen.csv"
R_SCRIPT="${DIR_BASE}/scripts/inference/Imagenes_Likelihood_frecuencias_mutaciones_m1_10.R"

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

# === 5. Evitar retrabajo ===
if [ -f "${OUT_FILE}" ]; then
    echo "El archivo ${OUT_FILE} ya existe. Saltando este paso."
    exit 0
fi

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
