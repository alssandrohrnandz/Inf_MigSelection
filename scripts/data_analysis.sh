#!/bin/bash
#SBATCH --job-name=Sim_Mig_Wave
#SBATCH --partition=defq
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-250                # 5 valores migración * 50 réplicas = 250
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00

# === 1. Configuración Inicial ===
module load r/4.1.3
module load slim/5.1

# === 2. Parameter Sweep Math ===
MIG_VALUES=(0.1 0.01 0.001 0.0001 0.00001)
REPLICAS_PER_VAL=50

#rm data/results_Discrete/outputs_LL/*.txt
#rm data/results_Continuous/outputs_LL/*.txt

# Borra todo lo que haya en outputs_LL para empezar fresco


# Calcular índices
IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) / $REPLICAS_PER_VAL ))
CURRENT_MIG=${MIG_VALUES[$IDX]}
REAL_REP=$(( ($SLURM_ARRAY_TASK_ID - 1) % $REPLICAS_PER_VAL + 1 ))

echo "DEBUG INFO:"
echo "Job ID Global: $SLURM_ARRAY_TASK_ID"
echo "  -> Valor Migración: $CURRENT_MIG"
echo "  -> Réplica #$REAL_REP del grupo"

START_TIME=$(date +%s)
TASK_ID=$SLURM_ARRAY_TASK_ID
echo "Iniciando Job ID: $TASK_ID en $(hostname)"

echo "--> Modo: $MODO | Acción: $ACCION"

DIR_BASE="/mnt/data/dortega/hlopezh/Inf_MigSelection"
# Argumento de entrada (continuo, discreto, ambos)
MODO=${1:-ambos}
# Argumento 2: Acción (completo, solo_analisis) - Por defecto: completo
ACCION=${2:-completo}

echo "--> Mode: $MODO | Action: $ACCION"

# Definiendo argumentos de SLiM
SLIM_ARGS="-d id_replica=$TASK_ID -d MIG=$CURRENT_MIG"

FILES_TO_PROCESS=()

# === 3. Ejecución de SLiM ===

# --- MODO CONTINUO ---
if [[ "$MODO" == "continuo" || "$MODO" == "ambos" ]]; then

    # Crear carpetas necesarias
    mkdir -p "${DIR_BASE}/data/results_Continuous/subsets"
    mkdir -p "${DIR_BASE}/data/results_Continuous/outputs_slim"
    mkdir -p "${DIR_BASE}/data/results_Continuous/outputs_LL"

    FILES_TO_PROCESS+=(
        "C_FULL_seleccion_m2"
        "C_FULL_neutros_m1"
        "C_aDNA_scattered_neutros_m1"
        "C_aDNA_scattered_seleccion_m2"
    )
    if [[ "$ACCION" != "solo_analisis" ]]; then
        echo "    Ejecutando SLiM: Continuous Space..."
        slim $SLIM_ARGS "${DIR_BASE}/scripts/Continuous_Space_Inference/Continuous_Space.slim"
    else
        echo "    SALTANDO SLiM (Continuous Space) - Se usarán archivos existentes."
    fi
fi

# --- MODO DISCRETO ---
if [[ "$MODO" == "discreto" || "$MODO" == "ambos" ]]; then

    # Crear carpetas necesarias
    mkdir -p "${DIR_BASE}/data/results_Discrete/subsets"
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_slim"
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_LL"
    
    FILES_TO_PROCESS+=(
        #"D_FULL_seleccion_m2"
        "D_FULL_neutros_m1" #TODO: QUITAR EL TEST
        #"D_aDNA_scattered_neutros_m1"
        #"D_aDNA_scattered_seleccion_m2"
    )
    #TODO: QUITAR EL TEST
    FILE_CHECK="${DIR_BASE}/data/results_Discrete/outputs_slim/D_FULL_neutros_m1_${TASK_ID}.csv"

    if [[ "$ACCION" == "solo_analisis" ]] || [[ -f "$FILE_CHECK" && -s "$FILE_CHECK" ]]; then
        
        echo "--> [SKIP] Saltando SLiM (Solicitado 'solo_analisis' o archivo ya existente)."
        
    else 
        
        echo "--> [RUN] Ejecutando SLiM: Discrete Space..."
        slim $SLIM_ARGS "${DIR_BASE}/scripts/Discrete_Space_Inference/Discrete_Space.slim"
        
    fi
fi

# Verificación de seguridad
if [ ${#FILES_TO_PROCESS[@]} -eq 0 ]; then
    echo "Error: Modo desconocido '$MODO'. Usa: continuo, discreto o ambos."
    exit 1
fi

# === 4. Análisis en R (Dinámico) ===
echo "--> Iniciando extracción y análisis en R..."

for PREFIJO in "${FILES_TO_PROCESS[@]}"; do
    # CORRECCIÓN 2: Determinar rutas dinámicamente según el prefijo del archivo
    # Si empieza con "C_", es Continuo. Si es "D_", es Discreto.
    if [[ "$PREFIJO" == "C_"* ]]; then
        BASE_PATH_TYPE="results_Continuous"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Continuous_Space_Inference/infLikelihood_mutations.R"
    else
        BASE_PATH_TYPE="results_Discrete"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Discrete_Space_Inference/infLikelihood_mutations.R"
    fi
    
    CURRENT_SLIM_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_slim"
    CURRENT_SUBSET_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/subsets"
    LL_OUTPUT="${DIR_BASE}"/data/"${BASE_PATH_TYPE}/outputs_LL"

    SLIM_OUTPUT="${CURRENT_SLIM_DIR}/${PREFIJO}_${TASK_ID}.csv"
    SUBSET_OUTPUT="${CURRENT_SUBSET_DIR}/subset_${PREFIJO}_${TASK_ID}.txt"
    
    # Verificación y Extracción (AWK)
    if [ -f "${SLIM_OUTPUT}" ]; then
        
        if [[ "$PREFIJO" == *"m1"* ]]; then
            echo "    [Subsampling] Seleccionando 1000 SNPs neutros al azar..."
            awk -F "," 'NR>1 {print $2}' "${SLIM_OUTPUT}" | sort | uniq | shuf | head -n 1000 > "${SUBSET_OUTPUT}"
        else
            echo "    [Full] Conservando todas las mutaciones bajo selección..."
            awk -F "," 'NR>1 {print $2}' "${SLIM_OUTPUT}" | sort | uniq > "${SUBSET_OUTPUT}"
        fi
        
        # Ejecutar R
        if [ -s "${SUBSET_OUTPUT}" ]; then
            echo "    [${PREFIJO}] Analizando en R..."
            
            # CORRECCIÓN 4: Cerrada la comilla al final y variables correctas
            Rscript --vanilla "${SCRIPT_R_PATH}" "${SLIM_OUTPUT}" "${SUBSET_OUTPUT}" "${TASK_ID}" "${PREFIJO}" "${LL_OUTPUT}"
            
        else
            echo "    ALERTA: El subset para ${PREFIJO} quedó vacío."
        fi
        
    else
        echo "    ALERTA: No se encontró la salida de SLiM: ${SLIM_OUTPUT}"
    fi
done

# === 5. Finalizar ===
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Job ${TASK_ID} finalizado."
echo "Duracion: $DURATION segundos (~$(($DURATION / 60)) min)"