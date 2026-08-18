library(deSolve)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
freq_file   <- args[1] #Archivo de las frecuencias
task_id <- as.numeric(args[2])
output_dir  <- args[3]
pop_info  <- args[4]


#if(length(args) < 3) {
#  print("Usando parámetros de prueba locales...")
#  freq_file   <- "AADNA_data/frq_chunks/Master_Chunk_1.rds"
#  task_id     <- "1"
#  pop_info <- "AADNA_data/Eurasia_Populations_Periods_Grid.csv"
#  output_dir  <- "results/Empirical_data"
#}

print(paste("Procesando archivo:", freq_file))
print(paste("Task ID:", task_id))


# Función del modelo espacial
diffusion2D <- function(t, conc, par) {
  Conc <- matrix(nrow = n, ncol = n, data = conc)
  dConc <- matrix(0, nrow = n, ncol = n)
  Flux <- -par[1] * rbind(rep(0, n), (Conc[2:n, ] - Conc[1:(n-1), ]), rep(0, n)) / dx
  dConc <- dConc - (Flux[2:(n+1), ] - Flux[1:n, ]) / dx
  Flux <- -par[1] * cbind(rep(0, n), (Conc[, 2:n] - Conc[, 1:(n-1)]), rep(0, n)) / dy
  dConc <- dConc - (Flux[, 2:(n+1)] - Flux[, 1:n]) / dy
  
  return(list(as.vector(dConc)))
}

GRID_SIZE <- 10
dy <- dx <- 1
n <- GRID_SIZE

# Parámetros para el Grid Search
exponentes_D <- -4:0
valores_base <- c(1, 2.5, 5)
DifussionValuesToCheck <- sort(unique(c(0, as.vector(outer(valores_base, 10^exponentes_D)), 0.125, 0.075)))
DifussionValuesToCheck <- DifussionValuesToCheck[DifussionValuesToCheck <= 1]

#SelectionValuesToCheck <- sort(unique(c(-valores_pos, 0, valores_pos)))
SelectionValuesToCheck <- 0 #sort(unique(c(0, valores_pos)))


# Lectura y preparación de datos

freq_data_raw <- setDT(readRDS(freq_file))
colnames(freq_data_raw) <- c("Chr","Pos","MutationID","CLST","A1","A2","Frequency","Count","Chr_Tot")

freq_data_raw$Frequency <- as.numeric(freq_data_raw$Frequency)
freq_data_raw$Count <- as.numeric(freq_data_raw$Count)
freq_data_raw$Chr_Tot <- as.numeric(freq_data_raw$Chr_Tot)
freq_data_raw$A1<-as.numeric(freq_data_raw$A1)
freq_data_raw$A2<-as.numeric(freq_data_raw$A2)

medians <- fread(pop_info)
colnames(medians) <- c("CLST", "Mean_Date", "Lat", "Long", "N_Samples", "Generation", "X_un", "Y_un", "Period")
min_lat <- min(medians$Lat); max_lat <- max(medians$Lat)
min_long <- min(medians$Long); max_long <- max(medians$Long)
medians[, X := round(((Lat - min_lat)/(max_lat - min_lat)) * 9 + 1)]
medians[, Y := round(((Long - min_long)/(max_long - min_long)) * 9 + 1)]

pop_info_merged <- merge(freq_data_raw, medians[, .(CLST, X, Y, Generation)], by.x = "CLST", by.y = "CLST", all.x = TRUE)

pop_info_merged$ChrOBS <- pop_info_merged$Chr_Tot
freq_data <- pop_info_merged
param_grid <- expand.grid(D=DifussionValuesToCheck)#, s=SelectionValuesToCheck)

# Lista de SNPS para analizar
snps_to_analyze <- unique(sort(freq_data$MutationID))
tiempos_globales_muestreo <- sort(unique(freq_data$Generation))

all_results_legacy <- list()

# ============================================================
# INFERENCIA DE D y s — SOLO MODELO BINOMIAL (simplificado)
# ============================================================

for (snp_actual in snps_to_analyze) {

  df_snp <- freq_data[freq_data$MutationID == snp_actual, ]
  Pos_actual <- df_snp$Pos[1]
  if (nrow(df_snp) == 0 || max(df_snp$Frequency) == 0) next

  df_snp <- df_snp[order(df_snp$Generation, df_snp$Y, df_snp$X), ]

  First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
  if (is.na(First_OcurrenceData)) next

  AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]
  max_gen_snp <- max(df_snp$Generation)
  max_gen_global  <- max(tiempos_globales_muestreo)

  df_snp <- df_snp[df_snp$Generation >= AlleleOriginAge, ]
  

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

  if(sum(Conc0) == 0) next
  

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
    D_curr <- param_grid$D[i]#; s_curr <- param_grid$s[i]; d_curr <- 2 * s_curr
    
    ST3 <- ode.2D(y = Conc0, times = times_run, func = diffusion2D,
              parms = D_curr, dimens = c(n, n), # <--- Solo pasamos D_curr
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

  col_names <- paste0("D_", param_grid$D)#, "_s_", param_grid$s)
  all_results_legacy[[length(all_results_legacy) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual),t(ll_legacy_vals)), c("SNP", col_names))
  print(paste("-> SNP", snp_actual, "procesado. Mejores D:", param_grid$D[which.max(ll_legacy_vals)], "con log-likelihood =", max(ll_legacy_vals, na.rm = TRUE)))
}
#print("Probando hasta aquí…")
#stop("Fin del test")
# --- GUARDADO ---
if (length(all_results_legacy) > 0) {
  df_legacy <- do.call(rbind, all_results_legacy)
  out_legacy <- file.path(output_dir, paste0("Analysis_LL_", task_id, "_neutrales", ".txt"))
  write.table(df_legacy, file = out_legacy, row.names = FALSE, quote = FALSE, sep = "\t")
  print(paste("Guardado Grid:", out_legacy))
} else {
  print("No se procesaron SNPs para Legacy.")
}

#matriz_LL<-df_legacy[,-1]
#composite_ll<-colSums(matriz_LL, na.rm = TRUE)
#best_idx <- which.max(composite_ll)
#best_params <- param_grid[best_idx, ]
#print(paste("Mejores parámetros encontrados: D =", best_params$D, ", s =", best_params$s, "con log-likelihood =", composite_ll[best_idx]))