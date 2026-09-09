#!/bin/bash
#SBATCH --job-name=TEST_Loci
#SBATCH --partition=defq
#SBATCH --output=logs/Loci_%A_%a.out
#SBATCH --error=logs/Loci_%A_%a.err
#SBATCH --array=11-550        
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00

# === 1. Configuración Inicial ===
module load r/4.1.3
module load slim/5.1

DIR_BASE="/mnt/data/dortega/hlopezh/Inf_MigSelection"
# Argumento de entrada (continuo, discreto, ambos)
MODO=${1:-ambos}
# Argumento 2: Acción (completo, solo_analisis) - Por defecto: completo
ACCION=${2:-completo}
# Argumento 3: Tipo de alelo (neutros, seleccion)
MODEL=${3:-mixto} 

# CORRECCIÓN: Espacios obligatorios y uso de "neutros" (con 's')
if [ "${MODEL}" == "mixto" ]; then
    SEL_VALUES=(0.01)
else
    SEL_VALUES=(0.1 0.075 0.05 0.025 0.01 0.0075 0.0050 0.0025 0.001 0.00075 0.0005 0.00025 0.0001 0.0)
fi
# === 2. Parameter Sweep Math ===
#MIG_VALUES=(0.0 0.000005 0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.025 0.05 0.075 0.1 0.125)
MIG_VALUES=(0.0 0.0001 0.0005 0.001 0.005 0.01 0.025 0.05 0.075 0.1 0.125) #Solo para bajo seleccion, para reducir el tiempo de ejecucion. Se pueden agregar mas valores para mayor robustez
# Equivale a calcular el % de migrantes por generacion en cada deme
# 0.00001 = 0.01 migrantes de un deme de acuerdo a N*m
# 0.001 = 1 migrante 
# 0.1= 100 migrantes de un deme 


 ## TODO: Modificado para agregar variacion en la seleccion
REPLICAS_PER_VAL=50 # Modificado a 10 para reducir el tiempo de ejecucion, pero se pueden aumentar para mayor robustez

# Calcular índices
IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) / $REPLICAS_PER_VAL ))
NUM_SEL=${#SEL_VALUES[@]} #Valores de seleccion que tenemos
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

# Validación de seguridad 
if [ -z "$CURRENT_MIG" ] || [ -z "$CURRENT_SEL" ]; then
    echo "Error: Indices fuera de rango. Revisa --array vs Arrays de valores."
    exit 1
fi

echo "--> Mode: $MODO | Action: $ACCION"

# Definiendo argumentos de SLiM
SLIM_ARGS="-d id_replica=$TASK_ID -d MIG=$CURRENT_MIG -d Sel_V=$SEL_VALUES -d MODEL='$MODEL' "

FILES_TO_PROCESS=()

# === 3. Ejecución de SLiM ===

# --- MODO CONTINUO --- este modo esta en construccion, solo usaremos el "discreto" por ahora
if [[ "$MODO" == "continuo" || "$MODO" == "ambos" ]]; then

    # Crear carpetas necesarias
    mkdir -p "${DIR_BASE}/data/results_Continuous/subsets"
    mkdir -p "${DIR_BASE}/data/results_Continuous/outputs_slim/independent_loci"
    mkdir -p "${DIR_BASE}/data/results_Continuous/outputs_LL"

    FILES_TO_PROCESS+=(
        "C_FULL_${MODEL}_m2"
        "C_FULL_${MODEL}_m1"
        #"C_aDNA_scattered_neutros_m1"
        #"C_aDNA_scattered_seleccion_m2"
    )
    if [[ "$ACCION" != "solo_analisis" ]]; then
        echo "    Ejecutando SLiM: Continuous Space..."
        #slim $SLIM_ARGS "${DIR_BASE}/scripts/Continuous_Space_Inference/Continuous_Space.slim"
    else
        echo "    SALTANDO SLiM (Continuous Space) - Se usarán archivos existentes."
    fi
fi

# --- MODO DISCRETO ---
if [[ "$MODO" == "discreto" || "$MODO" == "ambos" ]]; then

    # Crear carpetas necesarias
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_slim/loci/"
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_LL/loci/"
    FILES_TO_PROCESS+=(
        "D_FULL_${MODEL}_m1"
    )
    FILE_CHECK="${DIR_BASE}/data/results_Discrete/outputs_slim/loci/D_FULL_${MODEL}_m1_${TASK_ID}.csv" 
    if [[ "$ACCION" == "solo_analisis" ]] || [[ -f "$FILE_CHECK" && -s "$FILE_CHECK" ]]; then
        
        echo "--> [SKIP] Saltando SLiM (Solicitado 'solo_analisis' o archivo ya existente)."
        
    else 

        # EJECUCIÓN DE SLiM

        echo "--> [RUN] Ejecutando SLiM: Discrete Space..."
        slim $SLIM_ARGS "${DIR_BASE}/scripts/Discrete_Space_Inference/Discrete_Space_Loci.slim"
        
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

    if [[ "$PREFIJO" == "C_"* ]]; then
        BASE_PATH_TYPE="results_Continuous"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Continuous_Space_Inference/infLikelihood_mutations.R"
    else
        BASE_PATH_TYPE="results_Discrete"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Discrete_Space_Inference/infLikelihood_mutations_binomial.R"
    fi
    
    CURRENT_SLIM_DIR="${DIR_BASE}"/data/"${BASE_PATH_TYPE}/outputs_slim/loci"
    LL_OUTPUT="${DIR_BASE}"/data/"${BASE_PATH_TYPE}/outputs_LL/loci"

    SLIM_OUTPUT="${CURRENT_SLIM_DIR}/${PREFIJO}_${TASK_ID}.csv"
    
    # Verificación y Extracción (AWK)
    if [ -f "${SLIM_OUTPUT}" ]; then
        
        #CORRECCIÓN 4: Cerrada la comilla al final y variables correctas
        Rscript --vanilla "${SCRIPT_R_PATH}" "${SLIM_OUTPUT}" "${TASK_ID}" "${PREFIJO}" "${LL_OUTPUT}" "${CURRENT_MIG}" "${CURRENT_SEL}"
        
    else
        echo "    ALERTA: No se encontró la salida de SLiM: ${SLIM_OUTPUT}"
    fi
done

# === 5. Finalizar ===
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Job ${TASK_ID} finalizado."
echo "Duracion: $DURATION segundos (~$(($DURATION / 60)) min)"