#!/bin/bash
#SBATCH --job-name=Sim_Mig_Wave
#SBATCH --partition=defq
#SBATCH --output=logs/Sim_aDNA_%A_%a.out
#SBATCH --error=logs/Sim_aDNA_%A_%a.err
#SBATCH --array=1-280           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=100:00:00

# === 1. Configuración Inicial ===
module load r/4.1.3
module load slim/5.1

#TODO: Guardar todos los archivos .csv que se generen
# === 2. Parameter Sweep Math ===
MIG_VALUES=(0.0 0.01 0.05 0.1) #Valores de migracion que tenemos
SEL_VALUES=(0.1 0.05 0.01 0.005 0.001 0.0005 0.0001)
 ## TODO: Modificado para agregar variacion en la seleccion
REPLICAS_PER_VAL=10


IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) / $REPLICAS_PER_VAL ))
NUM_SEL=${#SEL_VALUES[@]} #Valores de seleccion que tenemos
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
echo "Iniciando Job ID: $TASK_ID en $(hostname)"


if [ -z "$CURRENT_MIG" ] || [ -z "$CURRENT_SEL" ]; then
    echo "Error: Indices fuera de rango. Revisa --array vs Arrays de valores."
    exit 1
fi

echo "--> Modo: $MODO | Acción: $ACCION"

DIR_BASE="/mnt/data/dortega/hlopezh/Inf_MigSelection"
# Argumento de entrada (continuo, discreto, ambos)
MODO=${1:-ambos}
# Argumento 2: Acción (completo, solo_analisis) - Por defecto: completo
ACCION=${2:-completo}

echo "--> Mode: $MODO | Action: $ACCION"


SLIM_ARGS="-d id_replica=$TASK_ID -d MIG=$CURRENT_MIG -d Sel_V=$CURRENT_SEL" 

FILES_TO_PROCESS=()



if [[ "$MODO" == "continuo" || "$MODO" == "ambos" ]]; then

    mkdir -p "${DIR_BASE}/data/results_Continuous/subsets"
    mkdir -p "${DIR_BASE}/data/results_Continuous/outputs_slim"
    mkdir -p "${DIR_BASE}/data/results_Continuous/outputs_LL"

    FILES_TO_PROCESS+=(
        "C_FULL_seleccion_m1" #modificado
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


if [[ "$MODO" == "discreto" || "$MODO" == "ambos" ]]; then


    mkdir -p "${DIR_BASE}/data/results_Discrete/subsets"
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_slim"
    mkdir -p "${DIR_BASE}/data/results_Discrete/outputs_LL"
    ##TODO: Mejorar para que se permita el análisis de archivo bajo seleccion y neutros
    FILES_TO_PROCESS+=(
        "v2_D_FULL_neutros_m1"
        "v2_D_FULL_seleccion_m2" 
    )
    FILE_CHECK="${DIR_BASE}/data/results_Discrete/outputs_slim/v2_D_FULL_neutros_m1_${TASK_ID}.csv" ##TODO: Aqui cambie el D_FULL_neutros_m1

    if [[ "$ACCION" == "solo_analisis" ]] || [[ -f "$FILE_CHECK" && -s "$FILE_CHECK" ]]; then
        
        echo "--> [SKIP] Saltando SLiM (Solicitado 'solo_analisis' o archivo ya existente)."
        
    else 
        
        echo "--> [RUN] Ejecutando SLiM: Discrete Space..."
        slim $SLIM_ARGS "${DIR_BASE}/scripts/Discrete_Space_Inference/slim_62x18.slim"
        
    fi
fi


if [ ${#FILES_TO_PROCESS[@]} -eq 0 ]; then
    echo "Error: Modo desconocido '$MODO'. Usa: continuo, discreto o ambos."
    exit 1
fi


echo "--> Iniciando extracción y análisis en R..."

for PREFIJO in "${FILES_TO_PROCESS[@]}"; do
    # CORRECCIÓN 2: Determinar rutas dinámicamente según el prefijo del archivo
    # Si empieza con "C_", es Continuo. Si es "D_", es Discreto.
    if [[ "$PREFIJO" == "C_"* ]]; then
        BASE_PATH_TYPE="results_Continuous"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/Continuous_Space_Inference/infLikelihood_mutations.R"
    else
        BASE_PATH_TYPE="results_Discrete"
        SCRIPT_R_PATH="${DIR_BASE}/scripts/spatial_pde_analysis_T_onset.R"
    fi

    
    CURRENT_SLIM_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/outputs_slim"
    CURRENT_SUBSET_DIR="${DIR_BASE}/data/${BASE_PATH_TYPE}/subsets"
    LL_OUTPUT="${DIR_BASE}"/data/"${BASE_PATH_TYPE}/outputs_LL"

    SLIM_OUTPUT="${CURRENT_SLIM_DIR}/${PREFIJO}_${TASK_ID}.csv"
    SUBSET_OUTPUT="${CURRENT_SUBSET_DIR}/subset_${PREFIJO}_${TASK_ID}.txt"
    
   if [[ "$PREFIJO" == "v2_"* ]]; then
        MODELO="SLiM_Data"
        # Usamos una variable vacía para mantener el mismo número de argumentos
        EXTRA_PATH="NONE"
    else
        MODELO="aDNA_Data"
        EXTRA_PATH="/mnt/data/dortega/hlopezh/Inf_MigSelection/AADNA_data/Population_info.txt"
    fi

# Ahora defines el array una sola vez, usando las variables que cambiaron arriba
    args=(
        --vanilla
        "${SCRIPT_R_PATH}"
        "${MODELO}"       
        "${EXTRA_PATH}"   
        "${SLIM_OUTPUT}"  
        "${SUBSET_OUTPUT}" 
        "${TASK_ID}"      
        "${PREFIJO}"      
        "${LL_OUTPUT}"    
    )


  
    if [ -f "${SLIM_OUTPUT}" ]; then
        
        if [[ "$PREFIJO" == *"m1"* ]]; then
            echo "    [Subsampling] Seleccionando N SNPs neutros al azar..."
            awk -F "," 'NR>1 {print $2}' "${SLIM_OUTPUT}" | sort | uniq  > "${SUBSET_OUTPUT}" #T
        else
            echo "    [Full] Conservando todas las mutaciones bajo selección..."
            awk -F "," 'NR>1 {print $2}' "${SLIM_OUTPUT}" | sort | uniq > "${SUBSET_OUTPUT}"
        fi
        

        if [ -s "${SUBSET_OUTPUT}" ]; then
            echo "    [${PREFIJO}] Analizando en R..."
            
            
            Rscript "${args[@]}"
            
        else
            echo "    ALERTA: El subset para ${PREFIJO} quedó vacío."
        fi
        
    else
        echo "    ALERTA: No se encontró la salida de SLiM: ${SLIM_OUTPUT}"
    fi
done


END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Job ${TASK_ID} finalizado."
echo "Duracion: $DURATION segundos (~$(($DURATION / 60)) min)"