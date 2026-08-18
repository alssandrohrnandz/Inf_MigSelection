#!/bin/bash
#SBATCH --job-name=Sim_aDNA_Sweep
#SBATCH --partition=defq
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-280           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00

# === 1. Configuración Inicial ===
module load r/4.1.3
module load slim/5.1

# === 2. Parameter Sweep Math ===
MIG_VALUES=(0.0 0.01 0.05 0.1) 
SEL_VALUES=(0.1 0.05 0.01 0.005 0.001 0.0005 0.0001)
REPLICAS_PER_VAL=10

IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) / $REPLICAS_PER_VAL ))
NUM_SEL=${#SEL_VALUES[@]} 
IDX_MIG=$(( $IDX / $NUM_SEL )) 
IDX_SEL=$(( $IDX % $NUM_SEL ))
CURRENT_MIG=${MIG_VALUES[$IDX_MIG]}
CURRENT_SEL=${SEL_VALUES[$IDX_SEL]}

REAL_REP=$(( ($SLURM_ARRAY_TASK_ID - 1) % $REPLICAS_PER_VAL + 1 ))

echo "Job ID: $SLURM_ARRAY_TASK_ID"
echo "  -> Migración [$IDX_MIG]: $CURRENT_MIG"
echo "  -> Selección [$IDX_SEL]: $CURRENT_SEL"
echo "  -> Réplica: $REAL_REP"

START_TIME=$(date +%s)
TASK_ID=$SLURM_ARRAY_TASK_ID

DIR_BASE="/mnt/data/dortega/hlopezh/Inf_MigSelection"
MODO=${1:-discreto}  # Por defecto solo discreto para esta prueba
ACCION=${2:-solo_analisis} # Forzamos solo_analisis para no re-correr SLiM

echo "--> Modo: $MODO | Acción: aDNA Downsampling"

# Usamos los archivos FULL como base
FILES_TO_PROCESS=(
    "D_FULL_neutros_m1"
    "D_FULL_seleccion_m2" 
)

echo "--> Iniciando extracción y simulación de aDNA en R..."

for PREFIJO in "${FILES_TO_PROCESS[@]}"; do
    
    BASE_PATH_TYPE="results_Discrete"
    SCRIPT_R_PATH="${DIR_BASE}/scripts/Discrete_Space_Inference/infLikelihood_aDNA.R"
    
    CURRENT_SLIM_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_slim"
    CURRENT_SUBSET_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/subsets"
    LL_OUTPUT="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_LL"

    SLIM_OUTPUT="${CURRENT_SLIM_DIR}/${PREFIJO}_${TASK_ID}.csv"
    SUBSET_OUTPUT="${CURRENT_SUBSET_DIR}/subset_${PREFIJO}_${TASK_ID}.txt"
    
    if [ -f "${SLIM_OUTPUT}" ]; then

        if [[ "$PREFIJO" == *"m1"* ]]; then
            echo "    [Subsampling] Seleccionando N SNPs neutros al azar..."
            awk -F "," 'NR>1 {print $2}' "${SLIM_OUTPUT}" | sort | uniq  > "${SUBSET_OUTPUT}" #T
        else
            echo "    [Full] Conservando todas las mutaciones bajo selección..."
            awk -F "," 'NR>1 {print $2}' "${SLIM_OUTPUT}" | sort | uniq > "${SUBSET_OUTPUT}"
        fi


        if [ -s "${SUBSET_OUTPUT}" ]; then
            echo "    [${PREFIJO}] Aplicando filtro aDNA y analizando..."
            
            # NOTA EL SEXTO ARGUMENTO "TRUE": Esto activa el filtro aDNA en R
            Rscript --vanilla "${SCRIPT_R_PATH}" "${SLIM_OUTPUT}" "${SUBSET_OUTPUT}" "${TASK_ID}" "${PREFIJO}_aDNA" "${LL_OUTPUT}" "TRUE"
            
        else
            echo "    ALERTA: El subset para ${PREFIJO} quedó vacío."
        fi
    else
        echo "    ALERTA: No se encontró la salida de SLiM base: ${SLIM_OUTPUT}"
    fi
done

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Job ${TASK_ID} finalizado. Duracion: $DURATION segundos."