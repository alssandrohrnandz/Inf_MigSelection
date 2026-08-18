#!/bin/bash
#SBATCH --job-name=Neu_Mig_Optim
#SBATCH --partition=defq
#SBATCH --output=logs/Neu_Optim_%A_%a.out
#SBATCH --error=logs/Neu_Optim_%A_%a.err
#SBATCH --array=1-700          
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20     # <--- CRÍTICO PARA PARALELIZAR EN R
#SBATCH --mem=32G              # <--- Aumentado para soportar múltiples procesos paralelos
#SBATCH --time=24:00:00

# === 1. Configuración Inicial ===
module load r/4.1.3
module load slim/5.1

# === 2. Parameter Sweep Math ===
MIG_VALUES=(0.0 0.000005 0.00001 0.00005  0.0001 0.0005 0.001 0.005 0.01 0.025 0.05 0.075 0.1 0.125)
SEL_VALUES=(0.0)
REPLICAS_PER_VAL=50

# Calcular índices
IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) / $REPLICAS_PER_VAL ))
NUM_SEL=${#SEL_VALUES[@]} 
IDX_MIG=$(( $IDX / $NUM_SEL )) 
IDX_SEL=$(( $IDX % $NUM_SEL ))
CURRENT_MIG=${MIG_VALUES[$IDX_MIG]}
CURRENT_SEL=${SEL_VALUES[$IDX_SEL]}

REAL_REP=$(( ($SLURM_ARRAY_TASK_ID - 1) % $REPLICAS_PER_VAL + 1 ))

# === DEBUG INFO ===
echo "Job ID: $SLURM_ARRAY_TASK_ID"
echo "  -> Migración [$IDX_MIG]: $CURRENT_MIG"
echo "  -> Selección [$IDX_SEL]: $CURRENT_SEL"
echo "  -> Réplica: $REAL_REP"

START_TIME=$(date +%s)
TASK_ID=$SLURM_ARRAY_TASK_ID
echo "Iniciando Job ID: $TASK_ID en $(hostname)"

if [ -z "$CURRENT_MIG" ] || [ -z "$CURRENT_SEL" ]; then
    echo "Error: Indices fuera de rango. Revisa --array vs Arrays de valores."
    exit 1
fi

MODO=${1:-ambos}
ACCION=${2:-completo}
echo "--> Mode: $MODO | Action: $ACCION"

DIR_BASE="/mnt/data/dortega/hlopezh/Inf_MigSelection"
SLIM_ARGS="-d id_replica=$TASK_ID -d MIG=$CURRENT_MIG -d Sel_V=$SEL_VALUES" 
FILES_TO_PROCESS=()

# === 3. Ejecución de SLiM ===
if [[ "$MODO" == "continuo" || "$MODO" == "ambos" ]]; then
    mkdir -p "${DIR_BASE}/data/results_Continuous/subsets"
    mkdir -p "${DIR_BASE}/data/results_Continuous/outputs_slim/independent_loci"
    mkdir -p "${DIR_BASE}/data/results_Continuous/outputs_LL"

    FILES_TO_PROCESS+=("C_FULL_seleccion_m2" "C_FULL_neutros_m1" "C_aDNA_scattered_neutros_m1" "C_aDNA_scattered_seleccion_m2")
    if [[ "$ACCION" != "solo_analisis" ]]; then
        echo "    Ejecutando SLiM: Continuous Space..."
    else
        echo "    SALTANDO SLiM (Continuous Space) - Se usarán archivos existentes."
    fi
fi

if [[ "$MODO" == "discreto" || "$MODO" == "ambos" ]]; then
    mkdir -p "${DIR_BASE}/data/results_Discrete/subsets"
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_slim/independent_loci"
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_LL/independent_loci"

    FILES_TO_PROCESS+=("D_FULL_neutros_m1")
    FILE_CHECK="${DIR_BASE}/data/results_Discrete/outputs_slim/independent_loci/D_FULL_neutros_m1_${TASK_ID}.csv" 

    if [[ "$ACCION" == "solo_analisis" ]] || [[ -f "$FILE_CHECK" && -s "$FILE_CHECK" ]]; then
        echo "--> [SKIP] Saltando SLiM (Solicitado 'solo_analisis' o archivo ya existente)."
    else 
        echo "--> [RUN] Ejecutando SLiM: Discrete Space..."
        slim $SLIM_ARGS "${DIR_BASE}/scripts/Discrete_Space_Inference/Discrete_Space_Sel.slim"
    fi
fi

if [ ${#FILES_TO_PROCESS[@]} -eq 0 ]; then
    echo "Error: Modo desconocido '$MODO'. Usa: continuo, discreto o ambos."
    exit 1
fi

# === 4. Análisis en R (Dinámico y Paralelo) ===
echo "--> Iniciando extracción y análisis en R..."

for PREFIJO in "${FILES_TO_PROCESS[@]}"; do
    if [[ "$PREFIJO" == "C_"* ]]; then
        BASE_PATH_TYPE="results_Continuous"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Continuous_Space_Inference/infLikelihood_mutations_optim.R"
    else
        BASE_PATH_TYPE="results_Discrete"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Discrete_Space_Inference/infLikelihood_mutations_optim.R"
    fi
    
    CURRENT_SLIM_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_slim/independent_loci"
    LL_OUTPUT="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_LL/independent_loci"
    SLIM_OUTPUT="${CURRENT_SLIM_DIR}/${PREFIJO}_${TASK_ID}.csv"
    
    if [ -f "${SLIM_OUTPUT}" ]; then
        Rscript --vanilla "${SCRIPT_R_PATH}" "${SLIM_OUTPUT}" "${TASK_ID}" "${PREFIJO}" "${LL_OUTPUT}"
    else
        echo "    ALERTA: No se encontró la salida de SLiM: ${SLIM_OUTPUT}"
    fi
done

# === 5. Finalizar ===
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Job ${TASK_ID} finalizado."
echo "Duracion: $DURATION segundos (~$(($DURATION / 60)) min)"