library(deSolve)
library(dplyr)
library(tidyverse)
library(VGAM)

args <- commandArgs(trailingOnly = TRUE)
freq_file   <- args[1] #Archivo de las frecuencias
subset_file <- args[2] #Tiene el nombre de la mutaciones
task_id     <- args[3]
model_name  <- args[4]
output_dir  <- args[5]

if(length(args) < 4) {
  print("Usando parámetros de prueba locales...")
  freq_file   <- "data/results_Discrete/outputs_slim/D_FULL_seleccion_m2_200.csv"
  subset_file <- "data/results_Discrete/subsets/subset_D_FULL_seleccion_m2_200.txt"
  task_id     <- "TEST_GRID"
  model_name  <- "D_FULL"
  output_dir  <- getwd()
}

print(paste("Procesando archivo:", freq_file))
print(paste("Modelo:", model_name))
print(paste("Task ID:", task_id))

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

GRID_SIZE <- 10
dy <- dx <- 1
n <- GRID_SIZE
N_eff <- 2000

if (model_name == "D_FULL_neutros_m1") {
  MAX_PADDED_GENERATIONS <- 1 
} else {
  MAX_PADDED_GENERATIONS <- 2 
}

# Parámetros para el Grid Search
exponentes_D <- -3:0
valores_base <- c(1, 5)
DifussionValuesToCheck <- sort(unique(c(0, as.vector(outer(valores_base, 10^exponentes_D)))))
DifussionValuesToCheck <- DifussionValuesToCheck[DifussionValuesToCheck <= 1]

if (model_name == "D_FULL_neutros_m1") {
  SelectionValuesToCheck <- c(0.0)
} else {
  bases <- c(1, 2, 2.5, 5)
  exponentes <- -5:0
  valores_pos <- as.vector(outer(bases, 10^exponentes))
  valores_pos <- valores_pos[valores_pos <= 1]
  SelectionValuesToCheck <- sort(unique(c(-valores_pos, 0, valores_pos)))
}

# Lectura y preparación de datos
freq_data_raw <- read.csv(freq_file, header=TRUE, stringsAsFactors=FALSE)
num_cols <- ncol(freq_data_raw)

if (num_cols == 8) {
  colnames(freq_data_raw) <- c("Generation","MutationID","TypeMut","X","Y","Frequency","Count","Chr_Tot")
} else if (num_cols == 7) {
  colnames(freq_data_raw) <- c("Generation","MutationID","X","Y","Frequency","Count","Chr_Tot")
  freq_data_raw$TypeMut <- NA 
}

freq_data <- freq_data_raw
freq_data$ChrOBS <- freq_data$Chr_Tot
tiempos_globales_muestreo <- sort(unique(freq_data$Generation))

# Filtramos subset
if (file.exists(subset_file) && file.info(subset_file)$size > 0) {
  snps_subset <- readLines(subset_file)
  snps_to_analyze <- intersect(unique(freq_data$MutationID), as.integer(snps_subset))
  print(paste("Analizando", length(snps_to_analyze), "SNPs del subset."))
} else {
  snps_to_analyze <- sort(unique(freq_data$MutationID))
}

all_results_legacy <- list()
all_results_ares <- list()
param_grid <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)

# === BUCLE PRINCIPAL ===
for (snp_actual in snps_to_analyze) {
  
  df_snp <- freq_data[freq_data$MutationID == snp_actual, ]
  if(nrow(df_snp) == 0 || max(df_snp$Frequency) == 0) next
  
  First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
  if(is.na(First_OcurrenceData)) next
  
  AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]
  max_gen_snp <- max(df_snp$Generation)
  max_gen_global <- max(tiempos_globales_muestreo)
  
  # --- ZERO-PADDING CONDICIONAL ---
  # Solo aplicamos padding si el alelo desaparece ANTES del final de la simulación (500)
  if (max_gen_snp < max_gen_global) {
    gens_faltantes <- tiempos_globales_muestreo[tiempos_globales_muestreo > max_gen_snp]
    if(length(gens_faltantes) > MAX_PADDED_GENERATIONS) {
      gens_faltantes <- gens_faltantes[1:MAX_PADDED_GENERATIONS]
    }
    
    coordenadas_historicas <- unique(df_snp[, c("X", "Y")])
    chr_tot_promedio <- max(round(mean(df_snp$Chr_Tot, na.rm = TRUE)), 10)
    
    filas_ceros <- expand.grid(Generation = gens_faltantes, X = coordenadas_historicas$X, Y = coordenadas_historicas$Y)
    filas_ceros$MutationID <- snp_actual
    filas_ceros$Frequency <- 0
    filas_ceros$Count <- 0
    filas_ceros$Chr_Tot <- chr_tot_promedio
    filas_ceros$TypeMut <- NA                   
    filas_ceros$ChrOBS <- chr_tot_promedio      
    
    filas_ceros <- filas_ceros[, colnames(df_snp)]
    df_snp <- rbind(df_snp, filas_ceros)
    df_snp <- df_snp[order(df_snp$Generation, df_snp$Y, df_snp$X), ]
  }
  
  # --- CONFIGURACIÓN DE LA PDE ---
  times_run <- min(df_snp$Generation):max(df_snp$Generation)
  
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
  
  ll_legacy_vals <- numeric(nrow(param_grid))
  ll_ares_vals <- numeric(nrow(param_grid))
  
  cat("Evaluando SNP:", snp_actual, "- Ventana temporal:", length(times_run), "generaciones.\n")
  
  for (i in 1:nrow(param_grid)) {
    D_curr <- param_grid$D[i]; s_curr <- param_grid$s[i]; d_curr <- 2 * s_curr
    
    ST3 <- ode.2D(y = Conc0, times = times_run, func = diffusion2D,
                  parms = c(D_curr, d_curr, s_curr), dimens = c(n, n),
                  method = "rk4", atol = 1e-7, rtol = 1e-7)
    
    ST3_mat <- as.matrix(ST3[,-1])
    piso_minimo <- 1e-6 
    
    # Recolección vectorizada
    counts_vec <- df_snp$Count
    totals_vec <- df_snp$Chr_Tot
    pred_freq_vec <- numeric(nrow(df_snp))
    rho_ares_vec <- numeric(nrow(df_snp))
    
    for(j in 1:nrow(df_snp)) {
      t_abs <- df_snp$Generation[j]
      time_idx <- match(t_abs, times_run)
      spatial_idx <- (df_snp$Y[j] - 1) * n + df_snp$X[j] 
      
      pred_freq_raw <- ST3_mat[time_idx, spatial_idx]
      pred_freq <- max(min(pred_freq_raw, 1 - piso_minimo), piso_minimo)
      
      t_elapsed <- max(t_abs - AlleleOriginAge, 0.5)
      rho_val <- max(1 - exp(-t_elapsed / (2 * N_eff)), 1e-6)
      
      pred_freq_vec[j] <- pred_freq
      rho_ares_vec[j] <- rho_val
    }
    
    ll_legacy_vals[i] <- sum(dbinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, log = TRUE))
    ll_ares_vals[i] <- sum(dbetabinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, rho = rho_ares_vec, log = TRUE))
  }
  
  col_names <- paste0("D_", param_grid$D, "_s_", param_grid$s)
  all_results_legacy[[length(all_results_legacy) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_legacy_vals)), c("SNP", col_names))
  all_results_ares[[length(all_results_ares) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_ares_vals)), c("SNP", col_names))
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


print("Ejecución de Grid Search completada.")