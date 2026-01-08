# Proyecto de Selección y Migración

Este repositorio contiene un pipeline bioinformático basado en **SLiM** y **R** para simular la evolución de 100 subpoblaciones bajo distintos escenarios de selección y migración.

**Autor:** Héctor Alessandro López Hernández  
**Fecha:** Enero 2026

## 🚀 Orden de Ejecución (Pipeline)

El flujo de trabajo está automatizado en dos etapas principales:

### 1. Simulación e Inferencia
Ejecutar el script maestro en Bash:
`scripts/data_analysis.sh`

Este script se encarga de orquestar todo el proceso automáticamente:
1.  **Simulación:** Genera los archivos de salida de los modelos espaciales (`D_space.slim` y `C_space.slim`).
2.  **Inferencia:** Inmediatamente después de la simulación, procesa las frecuencias y datos espaciales para calcular la verosimilitud utilizando el script de R: `scripts/inference/infLikelihood_mutations.R`.

### 2. Visualización de Resultados
Para observar e interpretar los valores de verosimilitud calculados, ejecutar:
`scripts/inference/plotting/Plot_Likelihood.R`

* 📂 **Salida:** Los gráficos resultantes se guardarán automáticamente en la carpeta `figures/`.

---

# Selection and Migration Project

This repository hosts a bioinformatics pipeline using **SLiM** and **R** to simulate the evolution of 100 subpopulations under various selection and migration scenarios.

**Author:** Héctor Alessandro López Hernández  
**Date:** January 2026

## 🚀 Execution Workflow

The analysis pipeline is divided into two main stages:

### 1. Simulation & Inference
Run the master Bash script:
`scripts/data_analysis.sh`

This script automates the entire upstream process:
1.  **Simulation:** Generates output files from the spatial models (`D_space.slim` and `C_space.slim`).
2.  **Inference:** Immediately calculates likelihood values based on the spatial data and allele frequencies using the R script: `scripts/inference/infLikelihood_mutations.R`.

### 2. Visualization
To visualize the likelihood results, run the plotting script:
`scripts/inference/plotting/Plot_Likelihood.R`

* 📂 **Output:** The resulting plots will be saved in the `figures/` directory.

---