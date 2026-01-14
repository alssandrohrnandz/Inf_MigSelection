#!/bin/bash
#SBATCH --job-name=Plot_Mig_Wave
#SBATCH --partition=defq
#SBATCH --output=logs/plot_%A_%a.out
#SBATCH --error=logs/plot_%A_%a.err
#SBATCH --array=1,51,101,151,201        # Inicios de cada grupo (1 Job por valor de Migración)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

module load r/4.1.3
MIG_VALUES=(0.1 0.01 0.001 0.0001 0.00001)
REPLICAS_PER_VAL=50

# Calcular índices para saber qué Valor de Migración es este
IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) / $REPLICAS_PER_VAL ))
CURRENT_MIG=${MIG_VALUES[$IDX]}

# === NUEVA LÓGICA DE RANGO ===
START_ID=$SLURM_ARRAY_TASK_ID
END_ID=$(( START_ID + REPLICAS_PER_VAL - 1 ))

# Generamos una cadena Regex con todos los IDs del grupo: "(1|2|...|50)"
# Esto le dirá a R: "Busca cualquier archivo que tenga el ID 1, O el 2, O el 3..."
SEQ_IDS=$(seq -s "|" $START_ID $END_ID)
TASK_ID_PATTERN="(${SEQ_IDS})"

echo "DEBUG INFO:"
echo "Job Array ID: $SLURM_ARRAY_TASK_ID"
echo " -> Valor Migración (m): $CURRENT_MIG"
echo " -> Procesando réplicas del $START_ID al $END_ID"
echo " -> Patrón generado para R: ${TASK_ID_PATTERN:0:20}..." # Imprimimos solo el inicio para no llenar el log

DIR_BASE="/mnt/data/dortega/hlopezh/Inf_MigSelection"
FILES_TO_PROCESS=()
MODO=${1:-ambos}

# ... (Bloques de selección de archivos MODO continuo/discreto IGUAL QUE ANTES) ...
if [[ "$MODO" == "continuo" || "$MODO" == "ambos" ]]; then
    FILES_TO_PROCESS+=(
        "Analysis_C_FULL_seleccion_m2"
        "Analysis_C_FULL_neutros_m1"
        "Analysis_C_aDNA_scattered_neutros_m1"
        "Analysis_C_aDNA_scattered_seleccion_m2"
    )
fi

if [[ "$MODO" == "discreto" || "$MODO" == "ambos" ]]; then
    FILES_TO_PROCESS+=(
        "Analysis_D_FULL_neutros_m1"
    )
fi

if [ ${#FILES_TO_PROCESS[@]} -eq 0 ]; then
    echo "Error: Modo desconocido '$MODO'. Usa: continuo, discreto o ambos."
    exit 1
fi

echo "--> Iniciando Plot en R (Modo Batch Agrupado)..."

for PREFIJO in "${FILES_TO_PROCESS[@]}"; do
    
    # 1. Definir rutas
    if [[ "$PREFIJO" == "Analysis_C_"* ]]; then
        BASE_PATH_TYPE="results_Continuous"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Continuous_Space_Inference/Plot_Likelihood.R"
    else
        BASE_PATH_TYPE="results_Discrete"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Discrete_Space_Inference/Plot_Likelihood.R"
    fi

    INPUT_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_LL"
    OUTPUT_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/figures"
    mkdir -p "${OUTPUT_DIR}"

    # 2. Verificar existencia (Safety Check)
    # Buscamos si existe al menos el archivo de la primera réplica (START_ID)
    # Esto evita usar el regex complejo en 'ls' que podría fallar en bash.
    count=$(ls "${INPUT_DIR}/${PREFIJO}_${START_ID}_"* 2>/dev/null | wc -l)
    
    if [ "$count" -gt 0 ]; then
        echo "--> Ejecutando R para ${PREFIJO} (Grupo Migración $CURRENT_MIG)..."
        
        # 3. LLAMADA A R
        # Argumentos actualizados para soportar la nueva lógica:
        # 1: Input Dir
        # 2: Output Dir
        # 3: Patrón de IDs (Regex string) <-- AQUÍ VA LA SECUENCIA
        # 4: Valor de Migración (Para títulos y lógica)
        # 5: Prefijo
        
        Rscript --vanilla "${SCRIPT_R_PATH}" \
            "${INPUT_DIR}" \
            "${OUTPUT_DIR}" \
            "${TASK_ID_PATTERN}" \
            "${CURRENT_MIG}" \
            "${PREFIJO}"
        
    else
        echo "ALERTA: No se encontraron archivos iniciales (ID ${START_ID}) para ${PREFIJO}. Saltando."
    fi
    
done