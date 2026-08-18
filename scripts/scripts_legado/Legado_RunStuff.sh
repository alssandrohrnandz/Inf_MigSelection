#!/bin/bash
#SBATCH --job-name=ANA_Run
#SBATCH --partition=defq
#SBATCH --output=logs/ANA_%A_%a.out
#SBATCH --error=logs/ANA_%A_%a.err
#SBATCH --array=1           
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=biol_alessandro@yahoo.com

module load r/4.1.3

echo "Inicio del trabajo: $(date)"
start_time=$(date +%s)

# 1. Configuración de bloque
BLOCK_SIZE=1000
TASK_ID=$SGE_TASK_ID
START_LINE=$(( (TASK_ID - 1) * BLOCK_SIZE + 1 ))
END_LINE=$(( TASK_ID * BLOCK_SIZE ))

BIM_FILE="/mnt/Timina/dortega/hlopezh/data/data.bim"
FRQ_FILE="/mnt/Timina/dortega/hlopezh/data/frq_v4.txt"
DATA_FILE="OrderedCADDDatasets1.txt"
OUT_FILE="SNP/SNPData${TASK_ID}.txt"

echo "SGE_TASK_ID = $TASK_ID | Buscando SNPs del Chr1 (Bloque: $START_LINE a $END_LINE)"

rm -f $OUT_FILE

# ==========================================
# 2. Extracción Inteligente y Paginada con AWK
# ==========================================
# Extraemos el bloque y además guardamos el CHR y POS del primer SNP para el nombre del archivo
awk -v start="$START_LINE" -v end="$END_LINE" -v task="$TASK_ID" '
    BEGIN { first=1; count=0 }
    NR > 1 && $2 == 1 {
        count++
        if (count >= start && count <= end) {
            if (first == 1) {
                print $2, $3 > ("tmp_metadata_" task)
                first = 0
            }
            print $2 "\t" $3
        }
        if (count == end) { exit }
    }
' "$DATA_FILE" > tmp_patterns_${TASK_ID}.txt
# Validación de seguridad
if [ ! -s tmp_patterns_${TASK_ID}.txt ]; then
    echo "No hay más SNPs en el Chr1 para procesar. Finalizando."
    rm -f tmp_patterns_${TASK_ID}.txt
    exit 0
fi

# ==========================================
# 3. PROCESAMIENTO INDIVIDUAL (Renglón por Renglón)
# ==========================================
# Leemos el archivo línea por línea obteniendo explícitamente el CHR y la POS de cada SNP
while read -r CHR POS; do

    echo "-------------------------------------------------------"
    echo "Procesando SNP individual: Chr ${CHR} | Posición ${POS}"
    echo "-------------------------------------------------------"
    # 3.1 Crear el patrón específico para ESTE SNP
    # Genera una expresión regular estricta para evitar falsos positivos en el archivo .bim
    echo -e "${CHR} \t ${POS}" > tmp_single_pattern_${TASK_ID}.txt

    # 3.2 Buscar el rsID correspondiente en data.bim
    grep -f tmp_single_pattern_${TASK_ID}.txt $BIM_FILE | awk '{print $2}' > tmp_rsnumbers_${TASK_ID}.txt

    # Si encontramos el rsID, procedemos a extraer frecuencias y correr R
    if [ -s tmp_rsnumbers_${TASK_ID}.txt ]; then

        # Definir el nombre de salida único para este SNP específico
        OUT_FILE="SNP/SNPData_Chr${CHR}_Pos${POS}.txt"

        # 3.3 Extraer las frecuencias desde el archivo .frq
        grep -F -w -f tmp_rsnumbers_${TASK_ID}.txt $FRQ_FILE > $OUT_FILE

        # 3.4 Llamar al script de R pasándole los metadatos de este SNP específico
        Rscript --vanilla LikelihoodCalculations_August25_2022.R $TASK_ID $CHR $POS

    else
        echo "Advertencia: No se encontró un rsID para Chr ${CHR} Pos ${POS} en el archivo .bim"
    fi

done < tmp_patterns_${TASK_ID}.txt
# ==========================================
# 4. LIMPIEZA DE ARCHIVOS TEMPORALES del Job
# ==========================================
rm -f tmp_patterns_${TASK_ID}.txt tmp_single_pattern_${TASK_ID}.txt tmp_rsnumbers_${TASK_ID}.txt
echo "Job ${TASK_ID} finalizado con éxito para todo el bloque."

end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Fin del trabajo: $(date)"
echo "Duración del trabajo: $execution_time segundos"