#!/bin/bash
#SBATCH --job-name=Plot_Mig_Sel
#SBATCH --partition=defq
#SBATCH --output=logs/plot_%A_%a.out
#SBATCH --error=logs/plot_%A_%a.err
#SBATCH --array=1,11,21,31,41,51,61,71,81,91,101,111,121,131,141,151,161,171,181,191,201,211,221,231,241
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

module load r/4.1.3

# Definimos variables
MIG_VALUES=(0.1 0.01 0.001 0.0001 0.00001)
SEL_VALUES=(0.5 0.25 0.1 0.05 0.01)
REPLICAS_PER_VAL=10

# === LÓGICA DE ÍNDICES CORREGIDA (2 DIMENSIONES) ===

# 1. Calculamos el índice global del grupo (0 a 24)
IDX_GLOBAL=$(( ($SLURM_ARRAY_TASK_ID - 1) / $REPLICAS_PER_VAL ))

# 2. Obtenemos el tamaño del array de selección para hacer la división
NUM_SEL=${#SEL_VALUES[@]}

# 3. Matemáticas para separar Migración y Selección
# Migración cambia cada 5 bloques (división entera)
IDX_MIG=$(( $IDX_GLOBAL / $NUM_SEL ))
# Selección cambia en cada bloque (módulo)
IDX_SEL=$(( $IDX_GLOBAL % $NUM_SEL ))

# 4. Asignamos los valores correctos
CURRENT_MIG=${MIG_VALUES[$IDX_MIG]}
CURRENT_SEL=${SEL_VALUES[$IDX_SEL]}

# === LÓGICA DE RANGO (REGEX) ===
START_ID=$SLURM_ARRAY_TASK_ID
END_ID=$(( START_ID + REPLICAS_PER_VAL - 1 ))

# Generar patrón: (1|2|...|10)
SEQ_IDS=$(seq -s "|" $START_ID $END_ID)
TASK_ID_PATTERN="(${SEQ_IDS})"

echo "DEBUG INFO:"
echo "Job Array ID: $SLURM_ARRAY_TASK_ID"
echo " -> Global Index: $IDX_GLOBAL"
echo " -> Migración [$IDX_MIG]: $CURRENT_MIG"
echo " -> Selección [$IDX_SEL]: $CURRENT_SEL"
echo " -> Rango IDs: $START_ID - $END_ID"

# Validación de seguridad
if [ -z "$CURRENT_MIG" ] || [ -z "$CURRENT_SEL" ]; then
    echo "Error: Indices fuera de rango."
    exit 1
fi

DIR_BASE="/mnt/data/dortega/hlopezh/Inf_MigSelection"
FILES_TO_PROCESS=()
MODO=${1:-ambos}

# --- Selección de Archivos ---

if [[ "$MODO" == "continuo" || "$MODO" == "ambos" ]]; then
    FILES_TO_PROCESS+=(
        "Analysis_C_FULL_neutros_m1"
        # "Analysis_C_FULL_seleccion_m2"
    )
fi

if [[ "$MODO" == "discreto" || "$MODO" == "ambos" ]]; then
    FILES_TO_PROCESS+=(
        "Analysis_D_FULL_seleccion_m2" # CORRECCIÓN 2: Cambiado de C_ a D_
    )
fi

if [ ${#FILES_TO_PROCESS[@]} -eq 0 ]; then
    echo "Error: Modo desconocido '$MODO'. Usa: continuo, discreto o ambos."
    exit 1
fi

echo "--> Iniciando Plot en R..."

for PREFIJO in "${FILES_TO_PROCESS[@]}"; do
    
    # 1. Definir rutas dinámicamente
    if [[ "$PREFIJO" == "Analysis_C_"* ]]; then
        BASE_PATH_TYPE="results_Continuous"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Discrete_Space_Inference/Plot_Likelihood_Sel.R"
    else
        BASE_PATH_TYPE="results_Discrete"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Discrete_Space_Inference/Plot_Likelihood_Sel.R"
    fi

    INPUT_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_LL"
    OUTPUT_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/figures"
    mkdir -p "${OUTPUT_DIR}"

    # 2. Verificar existencia con el primer ID del grupo
    count=$(ls "${INPUT_DIR}/${PREFIJO}_${START_ID}_"* 2>/dev/null | wc -l)
    
    if [ "$count" -gt 0 ]; then
        echo "--> Ejecutando R para ${PREFIJO} (Mig: $CURRENT_MIG)..."
        
        Rscript --vanilla "${SCRIPT_R_PATH}" \
            "${INPUT_DIR}" \
            "${OUTPUT_DIR}" \
            "${TASK_ID_PATTERN}" \
            "${CURRENT_MIG}" \
            "${PREFIJO}" \
            "${CURRENT_SEL}"
        
    else
        echo "ALERTA: No se encontraron archivos iniciales (ID ${START_ID}) para ${PREFIJO}."
    fi
    
done