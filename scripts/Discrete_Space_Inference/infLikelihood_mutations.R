# ==============================================================================
# INFERENCIA ESPACIAL DE SELECCIÓN Y MIGRACIÓN (MODELO PDE 2D)
# ==============================================================================
# Autor:        Héctor Alessandro López Hernandez
# Fecha:        Mayo 2026
# Proyecto:     Inferencia de Selección y Migración a partir de Frecuencias Alélicas
# Script:       infLikelihood_mutations.R
#
# DESCRIPCIÓN:
#   Script para estimar parámetros de Difusión (D) y Selección (s) mediante
#   Grid Search sobre series temporales y espaciales de frecuencias alélicas.
#   Modela la dinámica de alelos mediante una PDE en 2D (Ecuación de Difusión-
#   Reacción) y evalúa la verosimilitud bajo tres modelos probabilísticos:
#     1. Binomial Clásica (Legacy)
#     2. Beta-Binomial (Ares)
#     3. Beta-Binomial con Spikes en bordes 0 y 1 (Spikes)
#
# USO / EJECUCIÓN:
#   Rscript infLikelihood_mutations.R <freq_file> <task_id> <model_name> <output_dir> <mig> <sel>
#   - Admite fallback a parámetros de prueba locales si no se pasan argumentos.
#
# MODELOS SOPORTADOS (<model_name>):
#   - "D_FULL_seleccion_m1" : Evalúa cuadrícula de Selección (s) y Migración/Difusión (D).
#   - "D_FULL_neutros_m1"   : Evalúa únicamente modelo Neutro (s = 0.0).
#
# ENTRADAS (INPUTS):
#   - CSV de frecuencias temporales/espaciales de SNPs (mínimo 7 u 8 columnas:
#     Generation, MutationID, [TypeMut], X, Y, Frequency, Count, Chr_Tot).
#
# SALIDAS (OUTPUTS):
#   - TRON_LEGACY_Grid_TaskID_<id>_<modelo>.txt  (Log-Likelihood Binomial)
#   - TRON_ARES_Grid_TaskID_<id>_<modelo>.txt    (Log-Likelihood Beta-Binomial)
#   - TRON_SPIKES_Grid_TaskID_<id>_<modelo>.txt  (Log-Likelihood Beta-Binomial c/ Spikes)
#
# DEPENDENCIAS / LIBRERÍAS:
#   - deSolve   : Resolución numérica de Ecuaciones Diferenciales (ode.2D, rk4).
#   - data.table: Manipulación de estructuras de datos masivas.
#   - VGAM      : Distribución Beta-Binomial (dbetabinom).
# ==============================================================================

library(deSolve)
library(data.table)
library(VGAM)

args <- commandArgs(trailingOnly = TRUE)
freq_file   <- args[1] #Archivo de las frecuencias
task_id     <- args[2]
model_name  <- args[3]
output_dir  <- args[4]
current_mig <- args[5]
current_sel <- args[6]

if(length(args) < 3) {
  print("Usando parámetros de prueba locales...")
  freq_file   <- "data/results_Discrete/outputs_slim/independent_loci/seleccion/D_FULL_seleccion_m1_1.csv"
  # freq_file   <- "/mnt/data/dortega/hlopezh/Inf_MigSelection/data/results_Discrete/outputs_slim/independent_loci/D_FULL_neutros_m1_3.csv"
  task_id     <- "1"
  # task_id <- "3"
  model_name  <- "D_FULL_selecion_m1"
  output_dir  <- getwd()
  current_mig <- "0.0"
  current_sel <- "0.01"
}

print(paste("Procesando archivo:", freq_file))
print(paste("Modelo:", model_name))
print(paste("Task ID:", task_id))
print(paste("Directorio de salida:", output_dir))
if (model_name == "D_FULL_seleccion_m1") {
  print(paste0("Evaluando solo modelo con mig (",current_mig,"). y selección (",current_sel,")."))
} else {
  print(paste0("Evaluando modelo neutro y mig (",current_mig,")."))
}

# stop("Fin del test")

# Función del modelo espacial
diffusion2D <- function(t, conc, par) {
  Conc <- matrix(nrow = n, ncol = n, data = conc)
  dConc <- Conc*(1-Conc)*(Conc*par[2]+par[3]*(1-2*Conc))
  Flux <- -par[1] * rbind(rep(0, n), (Conc[2:n,]-Conc[1:(n-1),]), rep(0, n))/dx
  dConc <- dConc - (Flux[2:(n+1),] - Flux[1:n,])/dx
  Flux <- -par[1] * cbind(rep(0, n), (Conc[,2:n]-Conc[,1:(n-1)]), rep(0, n))/dy
  dConc <- dConc - (Flux[,2:(n+1)]-Flux[,1:n])/dy
  return(list(as.vector(dConc)))
}

dbetabinom_spikes <- function(x, size, prob, rho, log = FALSE) {
  pi_0 <- (1 - prob) * rho
  pi_1 <- prob * rho
  
  # CORRECCIÓN 1: Usar pmax() en lugar de max() para evaluación vectorizada
  pi_mid <- pmax(1 - pi_0 - pi_1, 1e-10) 
  
  ll_bb <- dbetabinom(x = x, size = size, prob = prob, rho = rho, log = FALSE)
  dens <- pi_mid * ll_bb
  
  # CORRECCIÓN 2: Filtrar pi_0 y pi_1 con el mismo índice lógico que dens
  dens[x == 0] <- dens[x == 0] + pi_0[x == 0]
  dens[x == size] <- dens[x == size] + pi_1[x == size]
  
  dens <- pmax(dens, 1e-300)
  
  if (log) return(sum(log(dens))) else return(dens)
}

GRID_SIZE <- 10
dy <- dx <- 1
n <- GRID_SIZE
N_eff <- 1000

if (model_name == "D_FULL_neutros_m1") {
  MAX_PADDED_GENERATIONS <- 1 
} else {
  MAX_PADDED_GENERATIONS <- 1 
}

# Parámetros para el Grid Search
exponentes_D <- -6:0
valores_base <- c(1, 2.5, 5, 7.5)
DifussionValuesToCheck <- sort(unique(c(0, as.vector(outer(valores_base, 10^exponentes_D)), 0.125)))
DifussionValuesToCheck <- DifussionValuesToCheck[DifussionValuesToCheck <= 1]

if (model_name == "D_FULL_neutros_m1") {
  SelectionValuesToCheck <- c(0.0)
} else {
  bases <- c(1, 2.5, 5, 7.5)
  exponentes <- -5:0
  valores_pos <- as.vector(outer(bases, 10^exponentes))
  valores_pos <- valores_pos[valores_pos <= 1]
  #SelectionValuesToCheck <- sort(unique(c(-valores_pos, 0, valores_pos)))
  SelectionValuesToCheck <- sort(unique(c(0, valores_pos)))
  DifussionValuesToCheck <- current_mig <- as.numeric(current_mig)
}

# Esta parte esta en construccion
#if (model_name == "D_FULL_seleccion_m1") {
#  SelectionValuesToCheck <- current_sel <- as.numeric(current_sel)
#}
#if (length(SelectionValuesToCheck) == 1) {
#  FIXED_SELECTION <- TRUE
#} else {
#  FIXED_SELECTION <- FALSE
#}

# Lectura y preparación de datos
freq_data_raw <- read.csv(freq_file, header=TRUE, stringsAsFactors=FALSE)
num_cols <- ncol(freq_data_raw)

if (num_cols == 8) {
  colnames(freq_data_raw) <- c("Generation","MutationID","TypeMut","X","Y","Frequency","Count","Chr_Tot")
} else if (num_cols == 7) {
  colnames(freq_data_raw) <- c("Generation","MutationID","X","Y","Frequency","Count","Chr_Tot")
  freq_data_raw$TypeMut <- NA 
}
# TODO Tomar solo los alelos que hayan surgido dentro de las últimas generaciones
# Repetir Seleccion y Dispersion
freq_data <- freq_data_raw #[freq_data_raw$Generation >= 600, ]
freq_data$ChrOBS <- freq_data$Chr_Tot
tiempos_globales_muestreo <- sort(unique(freq_data$Generation))
snps_to_analyze <- unique(sort(freq_data$MutationID))
print(paste("Número de SNPs a analizar en 1000 generaciones:", length(snps_to_analyze)))
#all_results_legacy <- list()
#all_results_ares <- list()
#all_results_spikes <- list()
#param_grid <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)
#print(param_grid)

#library(dplyr)

# 1. Identificar los alelos (MutationID) que surgieron en o después de la gen 600
#snps_recientes <- freq_data_raw %>%
#  group_by(MutationID) %>%
#  summarise(Gen_Origen = min(Generation), .groups = "drop") %>%
#  filter(Gen_Origen >= 600) %>%
#  pull(MutationID)

# 2. Filtrar el dataframe para quedarse solo con los registros de esos alelos
#freq_data <- freq_data_raw %>%
#  filter(MutationID %in% snps_recientes)

# 3. Continuar con el resto de tu código
#freq_data$ChrOBS <- freq_data$Chr_Tot
#tiempos_globales_muestreo <- sort(unique(freq_data$Generation))
#snps_to_analyze <- unique(sort(freq_data$MutationID))

#print(paste("Número de SNPs surgidos a partir de la gen 600 a analizar:", length(snps_to_analyze)))

all_results_legacy <- list()
all_results_ares <- list()
all_results_spikes <- list()
param_grid <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)
print(param_grid)

# === BUCLE PRINCIPAL ===
for (snp_actual in snps_to_analyze) {
  
  df_snp <- freq_data[freq_data$MutationID == snp_actual, ]
  if(nrow(df_snp) == 0 || max(df_snp$Frequency) == 0) next
  
  df_snp <- df_snp[order(df_snp$Generation, df_snp$Y, df_snp$X), ]

  First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
  if(is.na(First_OcurrenceData)) next
  
  AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]

  max_gen_snp    <- max(df_snp$Generation)
  max_gen_global <- max(tiempos_globales_muestreo)
 
 # --- ZERO-PADDING ESPACIO-TEMPORAL COMPLETO ---
  
  # --- Determinar el tope temporal del análisis ---

  # Última generación con alelo presente (sobre datos empíricos, sin padding aún)
  gens_con_alelo <- df_snp$Generation[df_snp$Frequency > 0]
  last_positive_gen <- max(gens_con_alelo)

  gens_presentes <- sort(unique(df_snp$Generation))
  n_total <- n * n # Total de coordenadas posibles en tu grid (X * Y)

  # ¿El alelo está fijado en la última observación?
  freq_final <- max(df_snp$Frequency[df_snp$Generation == last_positive_gen])
  esta_fijado <- freq_final >= (1 - 1e-6)   # freq == 1 en la práctica

  gen_fijacion_global <- NA
  for (g in gens_presentes) {
    df_g <- df_snp[df_snp$Generation == g, ]
    # Contar en cuántos demes la frecuencia es prácticamente 1
    demes_fijos <- sum(df_g$Frequency >= (1 - 1e-6))
  
    if (demes_fijos == n_total) {
      gen_fijacion_global <- g
      break # Encontramos la primera generación donde se fijó globalmente
    }
  }

  if (!is.na(gen_fijacion_global)) {
    # Caso A: Alelo fijado globalmente.
    # El tope será exactamente esa primera generación de fijación total.
    gen_limite <- gen_fijacion_global
  
  } else {
    # Caso B: El alelo NO está fijado globalmente (se perdió o sigue segregando).
    # Evaluamos hasta UNA generación después de la última observación positiva.
    idx_last <- which(tiempos_globales_muestreo == last_positive_gen)
    if (length(idx_last) == 0) {
      # Última observación no está en el calendario global (caso raro)
      gen_limite <- last_positive_gen
    } else if (idx_last < length(tiempos_globales_muestreo)) {
      # Tomar un punto de muestreo extra para capturar el cero confirmatorio
      gen_limite <- tiempos_globales_muestreo[idx_last + 1]
    } else {
      # Se pierde en el último punto muestreado absoluto
      gen_limite <- last_positive_gen
    }
  }

  # Generaciones a evaluar: desde el origen hasta gen_limite, cubriendo huecos intermedios
  generaciones_totales <- tiempos_globales_muestreo[
    tiempos_globales_muestreo >= AlleleOriginAge &
    tiempos_globales_muestreo <= gen_limite]
  
  # 2. Crear el universo completo de coordenadas y tiempos posibles
  # Asume que 'n' es la dimensión de tu grid (ej. 10 para una metapoblación 10x10)
  grid_completo <- expand.grid(Generation = generaciones_totales, 
                               X = 1:n, 
                               Y = 1:n)
  grid_completo$MutationID <- snp_actual
  
  # 3. Unir (Merge) el grid completo con tus datos empíricos de SLiM
  df_snp_completo <- merge(grid_completo, df_snp, 
                           by = c("Generation", "X", "Y", "MutationID"), 
                           all.x = TRUE)
  
  # 4. Rellenar los verdaderos ceros (los que SLiM omitió por optimización)
  df_snp_completo$Frequency[is.na(df_snp_completo$Frequency)] <- 0
  df_snp_completo$Count[is.na(df_snp_completo$Count)] <- 0
  
  # 5. Rellenar la cantidad de cromosomas (Chr_Tot) para las zonas sin datos
  chr_tot_promedio <- max(round(mean(df_snp$Chr_Tot, na.rm = TRUE)), 10)
  df_snp_completo$Chr_Tot[is.na(df_snp_completo$Chr_Tot)] <- chr_tot_promedio
  
  # Opcional: limpiar NAs residuales de otras columnas si tu código las usa
  df_snp_completo$TypeMut[is.na(df_snp_completo$TypeMut)] <- NA
  df_snp_completo$ChrOBS[is.na(df_snp_completo$ChrOBS)] <- chr_tot_promedio
  
  # Sustituir y ordenar para que la PDE lo lea correctamente
  df_snp <- df_snp_completo[order(df_snp_completo$Generation,
                                 df_snp_completo$Y,
                                 df_snp_completo$X), ]
  
  # --- CONFIGURACIÓN DE LA PDE ---
  # Se resuelve generación por generación para mantener estabilidad numérica en rk4
  times_run <- seq(min(generaciones_totales), max(generaciones_totales), by = 1)
  
  # SEGURIDAD: ode.2D requiere al menos 2 puntos de tiempo distintos
  if(length(times_run) < 2) {
    cat("SNP", snp_actual, "omitido: Datos insuficientes para serie temporal (n < 2).\n")
    next
  }
  
  Conc0 <- matrix(0, nrow=n, ncol=n)
  origin_data <- df_snp[df_snp$Generation == AlleleOriginAge & df_snp$Frequency > 0, ]
  for(k in 1:nrow(origin_data)) {
    ox <- origin_data$X[k]; oy <- origin_data$Y[k]
    if(ox >= 1 && ox <= n && oy >= 1 && oy <= n) Conc0[ox, oy] <- origin_data$Frequency[k]
  }
  if(sum(Conc0) == 0) next
  
  # --- PRE-CALCULO DE INDICES (UNA sola vez por SNP, fuera del grid) ---

  time_idx_vec    <- match(df_snp$Generation, times_run)
  spatial_idx_vec <- (df_snp$Y - 1) * n + df_snp$X

  filas_validas <- !is.na(time_idx_vec)
  if (!all(filas_validas)) {
    cat("SNP", snp_actual, ": se omiten", sum(!filas_validas), "filas cuya Generation no está en times_run.\n")
  }

  counts_vec <- df_snp$Count[filas_validas]
  totals_vec <- df_snp$Chr_Tot[filas_validas]
  idx_matrix <- cbind(time_idx_vec[filas_validas], spatial_idx_vec[filas_validas])

  # rho no depende de D/s: se calcula una sola vez por SNP
  t_elapsed_vec <- pmax(df_snp$Generation[filas_validas] - AlleleOriginAge, 0.5)
  rho_ares_vec  <- pmax(1 - exp(-t_elapsed_vec / (2 * N_eff)), 1e-6)

  piso_minimo <- 1e-6

  ll_legacy_vals <- numeric(nrow(param_grid))
  ll_ares_vals <- numeric(nrow(param_grid))
  ll_spikes_vals <- numeric(nrow(param_grid))
  
  cat("Evaluando SNP:", snp_actual, "- Ventana temporal:", length(times_run), "generaciones.\n")
  
  for (i in 1:nrow(param_grid)) {
    D_curr <- param_grid$D[i]; s_curr <- param_grid$s[i]; d_curr <- 2 * s_curr
  
    ST3 <- ode.2D(y = Conc0, times = times_run, func = diffusion2D,
                parms = c(D_curr, d_curr, s_curr), dimens = c(n, n),
                method = "rk4", atol = 1e-7, rtol = 1e-7)
  
    ST3_mat <- as.matrix(ST3[,-1])

    if (any(!is.finite(ST3_mat))) {
      ll_legacy_vals[i] <- NA; ll_ares_vals[i] <- NA; ll_spikes_vals[i] <- NA
      next
    }
  
    pred_freq_raw <- ST3_mat[idx_matrix]
    pred_freq_vec <- pmin(pmax(pred_freq_raw, piso_minimo), 1 - piso_minimo)
  
    ll_legacy_vals[i] <- sum(dbinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, log = TRUE))
    ll_ares_vals[i]   <- sum(dbetabinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, rho = rho_ares_vec, log = TRUE))
    ll_spikes_vals[i] <- sum(dbetabinom_spikes(x = counts_vec, size = totals_vec, prob = pred_freq_vec, rho = rho_ares_vec, log = TRUE))
  }
  
  col_names <- paste0("D_", param_grid$D, "_s_", param_grid$s)
  all_results_legacy[[length(all_results_legacy) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_legacy_vals)), c("SNP", col_names))
  all_results_ares[[length(all_results_ares) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_ares_vals)), c("SNP", col_names))
  all_results_spikes[[length(all_results_spikes)+1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_spikes_vals)), c("SNP", col_names))
}

# --- GUARDADO ---
if(length(all_results_legacy) > 0) {
  df_legacy <- do.call(rbind, all_results_legacy)
  out_legacy <- file.path(output_dir, paste0("TRON_LEGACY_Grid_TaskID_", task_id, "_", model_name, ".txt"))
  write.table(df_legacy, file=out_legacy, row.names=FALSE, quote = FALSE, sep="\t")
  print(paste("Guardado Grid Legacy:", out_legacy))
} else {
  print("No se procesaron SNPs para Legacy.")
}

if(length(all_results_ares) > 0) {
  df_ares <- do.call(rbind, all_results_ares)
  out_ares <- file.path(output_dir, paste0("TRON_ARES_Grid_TaskID_", task_id, "_", model_name, ".txt"))
  write.table(df_ares, file=out_ares, row.names=FALSE, quote = FALSE, sep="\t")
  print(paste("Guardado Grid Ares:", out_ares))
} else {
  print("No se procesaron SNPs para Ares.")
}

if(length(all_results_spikes) > 0){
  df_spikes <- do.call(rbind, all_results_spikes)
  out_spikes <- file.path(output_dir, paste0("TRON_SPIKES_Grid_TaskID_", task_id,"_", model_name,".txt"))
  write.table(df_spikes, file=out_spikes, row.names=FALSE, quote= FALSE, sep="\t")
  print(paste("Guardado Grid Spikes:", out_spikes))
} else {
  print ("No se procesaro SNPs para Spikes")
}


print("Ejecución de Grid Search completada.")