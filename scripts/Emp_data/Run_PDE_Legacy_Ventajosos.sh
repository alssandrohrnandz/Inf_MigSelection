#!/bin/bash
#SBATCH --job-name=Ventajosos
#SBATCH --partition=defq
#SBATCH --output=logs/A_VENTAJOSOS_%A_%a.out
#SBATCH --error=logs/A_VENTAJOSOS_%A_%a.err
#SBATCH --array=1     
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=biol_alessandro@yahoo.com

module load r/4.1.3

TASK_ID=$SLURM_ARRAY_TASK_ID

# === DEBUG INFO ===
echo "Job ID: $SLURM_ARRAY_TASK_ID"

START_TIME=$(date +%s)
TASK_ID=$SLURM_ARRAY_TASK_ID
echo "Iniciando Job ID: $TASK_ID en $(hostname)"

mkdir -p results/Empirical_data

CHUNK_FILE="chr_frequencies/1240k_freq_chr_ventajosos.frq.strat"
SCRIPT_R_PATH="scripts/Emp_data/infLikelihood_mutations_binomial_ventajosos.R"
LL_OUTPUT="results/Empirical_data/"
POP_INFO="AADNA_data/Eurasia_Populations_Periods_Grid.csv"

Rscript --vanilla "${SCRIPT_R_PATH}" "${CHUNK_FILE}" "${TASK_ID}" "${LL_OUTPUT}" "${POP_INFO}"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "Job ${TASK_ID} finalizado."
echo "Duracion: $DURATION segundos (~$(($DURATION / 60)) min)"