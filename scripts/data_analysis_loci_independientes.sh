#!/bin/bash
#SBATCH --job-name=Sel_Mig
#SBATCH --partition=defq
#SBATCH --output=logs/IndLoci_%A_%a.out
#SBATCH --error=logs/IndLoci_%A_%a.err
#SBATCH --array=1-3200
#SBATCH --cpus-per-task=4
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
MODEL=${3:-neutros} 

# CORRECCIÓN: Espacios obligatorios y uso de "neutros" (con 's')
if [ "${MODEL}" == "neutros" ]; then
    SEL_VALUES=(0.0)
else
    SEL_VALUES=(0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.0) 
fi
# === 2. Parameter Sweep Math ===
MIG_VALUES=(0.0 0.0001 0.0005 0.001 0.005 0.01 0.05 0.1)
#si son neutros ir de 351 a 400

# Valores de migracion para alelos neutros
# MIG_VALUES=(0.0001 0.0005 0.001 0.005 0.01 0.025 0.05 0.075 0.1 0.125)
# Task ID va de 351 a 400, con 50 réplicas por combinación de parámetros. Total 50 archivos
# 0.0001: Nm=0.1 control casi aislado
# 0.0005: Nm=0.5 control aislado
# 0.001: Nm=1 control aislado
# 0.005: Nm=5 control moderado
# 0.01: Nm=10 control moderado
# 0.1: Nm=100 control alto


REPLICAS_PER_VAL=50 

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
SLIM_ARGS="-d id_replica=$TASK_ID -d MIG=$CURRENT_MIG -d Sel_V=$CURRENT_SEL -d MODEL='$MODEL' "

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
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_slim/independent_loci/${MODEL}"
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_LL/independent_loci/${MODEL}"

    FILES_TO_PROCESS+=(
        "D_FULL_${MODEL}_m1"
    )
    FILE_CHECK="${DIR_BASE}/data/results_Discrete/outputs_slim/independent_loci/${MODEL}/D_FULL_${MODEL}_m1_${TASK_ID}.csv" 
    if [[ "$ACCION" == "solo_analisis" ]] || [[ -f "$FILE_CHECK" && -s "$FILE_CHECK" ]]; then
        
        echo "--> [SKIP] Saltando SLiM (Solicitado 'solo_analisis' o archivo ya existente)."
        
    else 
        
        T_START_SLIM=$(date +%s)

        echo "--> [RUN] Ejecutando SLiM: Discrete Space..."
        slim $SLIM_ARGS "${DIR_BASE}/scripts/Discrete_Space_Inference/Discrete_Space_Sel.slim"

        T_END_SLIM=$(date +%s)
        DUR_SLIM=$((T_END_SLIM - T_START_SLIM))
        echo "Duracion ejecución: $DUR_SLIM segundos (~$(($DUR_SLIM / 60)) min)"

    fi
fi

# Verificación de seguridad
if [ ${#FILES_TO_PROCESS[@]} -eq 0 ]; then
    echo "Error: Modo desconocido '$MODO'. Usa: continuo, discreto o ambos."
    exit 1
fi

# === 4. Análisis en R (paralelo por modo de filtro) ===
echo "--> Iniciando extracción y análisis en R (paralelo)..."

# Evita que cada proceso de R use varios hilos
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

MODOS_FILTRO=("all" "nuevos_400") # "vigentes_400") desaparecido

# Función que lanza el análisis de un prefijo y un modo
lanzar_R() {
    local SCRIPT="$1" SLIM_OUT="$2" PREF="$3" LL_OUT="$4" MODO_F="$5" LOG="$6"
    mkdir -p "$LL_OUT" "$(dirname "$LOG")"
    Rscript --vanilla "$SCRIPT" "$SLIM_OUT" "$TASK_ID" "$PREF" "$LL_OUT" \
        "$CURRENT_MIG" "$CURRENT_SEL" "$MODO_F" > "$LOG" 2>&1
}

PIDS=()
LOGS=()
T_START_R=$(date +%s)

for PREFIJO in "${FILES_TO_PROCESS[@]}"; do
    if [[ "$PREFIJO" == "C_"* ]]; then
        BASE_PATH_TYPE="results_Continuous"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Continuous_Space_Inference/infLikelihood_mutations.R"
    else
        BASE_PATH_TYPE="results_Discrete"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Discrete_Space_Inference/infLikelihood_mutations.R"
    fi

    CURRENT_SLIM_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_slim/independent_loci/${MODEL}"
    SLIM_OUTPUT="${CURRENT_SLIM_DIR}/${PREFIJO}_${TASK_ID}.csv"

    echo "Archivo a analizar: ${SLIM_OUTPUT}"

    if [ ! -f "${SLIM_OUTPUT}" ]; then
        echo "    ALERTA: No se encontró la salida de SLiM: ${SLIM_OUTPUT}"
        continue
    fi

    if [[ "$CURRENT_MIG" == "0.0" && "$CURRENT_SEL" == "0.0" ]]; then
        LL_OUTPUT_NEUTROS="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_LL/independent_loci/neutros"
        LL_OUTPUT="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_LL/independent_loci/${MODEL}"
        # Caso neutro: cada modo se lanza para el archivo de neutros y para el normal
        for MODO_F in "${MODOS_FILTRO[@]}"; do
            LOG="${DIR_BASE}/logs/R_${PREFIJO}_neutros_${TASK_ID}_${MODO_F}.log"
            lanzar_R "$SCRIPT_R_PATH" "$SLIM_OUTPUT" "D_FULL_neutros_m1" "$LL_OUTPUT_NEUTROS" "$MODO_F" "$LOG" &
            PIDS+=($!); LOGS+=("$LOG")

            LOG="${DIR_BASE}/logs/R_${PREFIJO}_${TASK_ID}_${MODO_F}.log"
            lanzar_R "$SCRIPT_R_PATH" "$SLIM_OUTPUT" "$PREFIJO" "$LL_OUTPUT" "$MODO_F" "$LOG" &
            PIDS+=($!); LOGS+=("$LOG")
        done
    else
        LL_OUTPUT="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_LL/independent_loci/${MODEL}"
        for MODO_F in "${MODOS_FILTRO[@]}"; do
            LOG="${DIR_BASE}/logs/R_${PREFIJO}_${TASK_ID}_${MODO_F}.log"
            lanzar_R "$SCRIPT_R_PATH" "$SLIM_OUTPUT" "$PREFIJO" "$LL_OUTPUT" "$MODO_F" "$LOG" &
            PIDS+=($!); LOGS+=("$LOG")
        done
    fi
done

# Esperar a todos los procesos y contar fallos
FALLOS=0
for i in "${!PIDS[@]}"; do
    if ! wait "${PIDS[$i]}"; then
        FALLOS=$((FALLOS + 1))
        echo "    FALLÓ: ${LOGS[$i]}"
    fi
done

DUR_R_SCRIPT=$(( $(date +%s) - T_START_R ))
echo "Análisis R terminado. Procesos: ${#PIDS[@]} | Fallos: $FALLOS | Duración: $DUR_R_SCRIPT s (~$((DUR_R_SCRIPT / 60)) min)"