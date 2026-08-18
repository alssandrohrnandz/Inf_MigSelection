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
  freq_file   <- "data/results_Discrete/outputs_slim/independent_loci/D_FULL_seleccion_m1_1.csv"
  # freq_file   <- "/mnt/data/dortega/hlopezh/Inf_MigSelection/data/results_Discrete/outputs_slim/independent_loci/D_FULL_neutros_m1_3.csv"
  task_id     <- "1"
  # task_id <- "3"
  model_name  <- "D_FULL_selecion_m1"
  output_dir  <- getwd()
}

print(paste("Procesando archivo:", freq_file))
print(paste("Modelo:", model_name))
print(paste("Task ID:", task_id))
if (model_name == "D_FULL_seleccion_m1") {
  print(paste0("Evaluando solo modelo con mig (",current_mig,"). y selección (",current_sel,")."))
} else {
  print(paste0("Evaluando modelo neutro y mig (",current_mig,")."))
}

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


if (model_name == "D_FULL_neutros_m1") {
  MAX_PADDED_GENERATIONS <- 1 
} else {
  MAX_PADDED_GENERATIONS <- 1 
}

# Parámetros para el Grid Search
exponentes_D <- -4:0
valores_base <- c(1, 2.5, 5)
DifussionValuesToCheck <- sort(unique(c(0, as.vector(outer(valores_base, 10^exponentes_D)), 0.125, 0.075)))
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

freq_data <- freq_data_raw[freq_data_raw$Generation >= 600, ]
freq_data$ChrOBS <- freq_data$Chr_Tot
tiempos_globales_muestreo <- sort(unique(freq_data$Generation))
snps_to_analyze <- unique(sort(freq_data$MutationID))
print(paste("Número de SNPs a analizar en 400 generaciones:", length(snps_to_analyze)))
all_results_legacy <- list()
all_results_ares <- list()
all_results_spikes <- list()
param_grid <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)

# ============================================================
# INFERENCIA DE D y s — SOLO MODELO BINOMIAL (simplificado)
# ============================================================

for (snp_actual in snps_to_analyze) {

  df_snp <- freq_data[freq_data$MutationID == snp_actual, ]
  if (nrow(df_snp) == 0 || max(df_snp$Frequency) == 0) next

  df_snp <- df_snp[order(df_snp$Generation, df_snp$Y, df_snp$X), ]

  First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
  if (is.na(First_OcurrenceData)) next

  AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]
  max_gen_global  <- max(tiempos_globales_muestreo)

  df_snp <- df_snp[df_snp$Generation >= AlleleOriginAge, ]
  max_gen_snp <- max(df_snp$Generation)

  # --- ZERO-PADDING CONDICIONAL (igual que antes) ---
  if (max_gen_snp < max_gen_global) {
    print(paste("SNP", snp_actual, "desaparece antes del final de la simulación. Aplicando zero-padding hasta generación", max_gen_global))
    gens_faltantes <- tiempos_globales_muestreo[tiempos_globales_muestreo > max_gen_snp]
    if (length(gens_faltantes) > MAX_PADDED_GENERATIONS) {
      gens_faltantes <- gens_faltantes[1:MAX_PADDED_GENERATIONS]
    }

    coordenadas_historicas <- unique(df_snp[, c("X", "Y")])
    chr_tot_promedio <- max(round(mean(df_snp$Chr_Tot, na.rm = TRUE)), 10)

    filas_ceros <- expand.grid(Generation = gens_faltantes, X = coordenadas_historicas$X, Y = coordenadas_historicas$Y)
    filas_ceros$MutationID <- snp_actual
    filas_ceros$Frequency  <- 0
    filas_ceros$Count      <- 0
    filas_ceros$Chr_Tot    <- chr_tot_promedio
    filas_ceros$TypeMut    <- NA
    filas_ceros$ChrOBS     <- chr_tot_promedio

    filas_ceros <- filas_ceros[, colnames(df_snp)]
    df_snp <- rbind(df_snp, filas_ceros)
    df_snp <- df_snp[order(df_snp$Generation, df_snp$Y, df_snp$X), ]
  }

  # --- CONFIGURACIÓN DE LA PDE ---
  times_run <- AlleleOriginAge:max(df_snp$Generation)

  if (length(times_run) < 2) {
    cat("SNP", snp_actual, "omitido: Datos insuficientes para serie temporal (n < 2).\n")
    next
  }

  Conc0 <- matrix(0, nrow = n, ncol = n)
  origin_data <- df_snp[df_snp$Generation == AlleleOriginAge & df_snp$Frequency > 0, ]
  for (k in 1:nrow(origin_data)) {
    ox <- origin_data$X[k]; oy <- origin_data$Y[k]
    if (ox >= 1 && ox <= n && oy >= 1 && oy <= n) Conc0[ox, oy] <- origin_data$Frequency[k]
  }
  if (sum(Conc0) == 0) next

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

  piso_minimo <- 1e-6

  ll_legacy_vals <- numeric(nrow(param_grid))

  cat("Evaluando SNP:", snp_actual, "- Ventana temporal:", length(times_run), "generaciones.\n")

  for (i in 1:nrow(param_grid)) {
    D_curr <- param_grid$D[i]; s_curr <- param_grid$s[i]; d_curr <- 2 * s_curr
    
    ST3 <- ode.2D(y = Conc0, times = times_run, func = diffusion2D,
                  parms = c(D_curr, d_curr, s_curr), dimens = c(n, n),
                  method = "rk4", atol = 1e-7, rtol = 1e-7)
    
    ST3_mat <- as.matrix(ST3[,-1])
    piso_minimo <- 1e-6 

    if (any(!is.finite(ST3_mat))) {
      ll_legacy_vals[i] <- NA
      next
    }

    pred_freq_raw <- ST3_mat[idx_matrix]
    pred_freq_vec <- pmin(pmax(pred_freq_raw, piso_minimo), 1 - piso_minimo)

    ll_legacy_vals[i] <- sum(dbinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, log = TRUE))
  }

  col_names <- paste0("D_", param_grid$D, "_s_", param_grid$s)
  all_results_legacy[[length(all_results_legacy) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_legacy_vals)), c("SNP", col_names))
}

# --- GUARDADO ---
if (length(all_results_legacy) > 0) {
  df_legacy <- do.call(rbind, all_results_legacy)
  out_legacy <- file.path(output_dir, paste0("P_TRON_LEGACY_Grid_TaskID_", task_id, "_", model_name, ".txt"))
  write.table(df_legacy, file = out_legacy, row.names = FALSE, quote = FALSE, sep = "\t")
  print(paste("Guardado Grid Legacy:", out_legacy))
} else {
  print("No se procesaron SNPs para Legacy.")
}

matriz_LL<-df_legacy[,-1]
composite_ll<-colSums(matriz_LL, na.rm = TRUE)
best_idx <- which.max(composite_ll)
best_params <- param_grid[best_idx, ]
print(paste("Mejores parámetros encontrados: D =", best_params$D, ", s =", best_params$s, "con log-likelihood =", composite_ll[best_idx]))