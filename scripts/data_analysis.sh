#!/bin/bash
#SBATCH --job-name=Sim_And_Analyze
#SBATCH --partition=defq
#SBATCH --output=logs/job_%A_%a.out
#SBATCH --error=logs/job_%A_%a.err
#SBATCH --array=1-5                # Correrá 5 simulaciones independientes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
# === 1. Configuración Inicial ===
module load r/4.1.3
module load slim/5.1

START_TIME=$(date +%s)
TASK_ID=$SLURM_ARRAY_TASK_ID   # <--- Definimos esto PRIMERO para poder usarlo abajo
echo "Iniciando Job ID: $TASK_ID en $(hostname)"

DIR_BASE="/mnt/data/dortega/hlopezh/Inf_MigSelection"

# Definimos rutas de scripts y archivos base (Asegúrate que existen)
R_SCRIPT="${DIR_BASE}/scripts/inference/infLikelihood_mutations.R"

# Crear carpetas si no existen
mkdir -p "${DIR_BASE}/data/subsets"
mkdir -p "${DIR_BASE}/data/results_simulations" # Aseguramos donde SLiM guarda
mkdir -p "${DIR_BASE}/data/results_LL"
mkdir -p "${DIR_BASE}/logs"

# === 2. Ejecutar Simulaciones SLiM ===
echo "--> Ejecutando SLiM (Continuous y Discrete)..."

# Nota: Asegúrate que tus scripts de SLiM guarden los CSVs en: ${DIR_BASE}/data/results_simulations/
slim -d id_replica=$TASK_ID "${DIR_BASE}/scripts/slim/Continuous_Space.slim"
slim -d id_replica=$TASK_ID "${DIR_BASE}/scripts/slim/Discrete_Space.slim"

# === 3. Procesamiento y Análisis (Bucle Maestro) ===

# Lista de modelos a procesar
#Los prefijos son los modelos. Nos ayudará a entender qué estamos analizando
FILES_TO_PROCESS=(
    "C_FULL_seleccion_m2"
    "C_FULL_neutros_m1"
    "C_aDNA_scattered_neutros_m1"     
    "C_aDNA_scattered_seleccion_m2"    
    "D_FULL_seleccion_m2"
    "D_FULL_neutros_m1"
    "D_aDNA_scattered_neutros_m1"      
    "D_aDNA_scattered_seleccion_m2"    
)

echo "--> Iniciando extracción y análisis en R..."

for PREFIJO in "${FILES_TO_PROCESS[@]}"; do
    
    # A. Rutas Dinámicas
    # Archivo que acaba de salir de SLiM
    SLIM_OUTPUT="${DIR_BASE}/data/results_simulations/${PREFIJO}_${TASK_ID}.csv"
    # Archivo intermedio con solo los IDs limpios
    SUBSET_OUTPUT="${DIR_BASE}/data/subsets/subset_${PREFIJO}_${TASK_ID}.txt"
    
    # B. Verificación y Extracción (AWK)
    if [ -f "${SLIM_OUTPUT}" ]; then
        
        # Extraer IDs (columna 2, saltando encabezado)
        awk -F "," 'NR>1 {print $2}' "${SLIM_OUTPUT}" | sort | uniq > "${SUBSET_OUTPUT}"
        
        # C. Ejecutar R para ESTE modelo específico
        if [ -s "${SUBSET_OUTPUT}" ]; then
            echo "    [${PREFIJO}] Analizando en R..."
         
            Rscript --vanilla "${R_SCRIPT}" "${SLIM_OUTPUT}" "${SUBSET_OUTPUT}" "${TASK_ID}" "${PREFIJO}
            
        else
            echo "    ALERTA: El subset para ${PREFIJO} quedó vacío."
        fi
        
    else
        echo "    ALERTA: No se encontró la salida de SLiM: ${SLIM_OUTPUT}"
    fi
done

# === 4. Finalizar ===
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Job ${TASK_ID} finalizado."
echo "Duracion: $DURATION segundos (~$(($DURATION / 60)) min)"