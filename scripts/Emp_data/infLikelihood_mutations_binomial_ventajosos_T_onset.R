library(deSolve)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
freq_file   <- args[1] #Archivo de las frecuencias
task_id <- as.numeric(args[2])
output_dir  <- args[3]
pop_info  <- args[4]



print("Usando parámetros de prueba locales...")
freq_file   <- "chr_frequencies/1240k_freq_chr_ventajosos.frq.strat"
task_id     <- "1"
pop_info <- "AADNA_data/Eurasia_Populations_Periods_Grid.csv"
output_dir  <- "results/Empirical_data"


print(paste("Procesando archivo:", freq_file))
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

# ============================================================
# PARAMETROS DEL GRID SEARCH ORIGINAL (fase de selección, usado
# como referencia / también disponible si se quiere reutilizar)
# ============================================================
exponentes_D <- -4:0
valores_base <- c(1, 2.5, 5)
DifussionValuesToCheck <- sort(unique(c(0, as.vector(outer(valores_base, 10^exponentes_D)), 0.125, 0.075)))
DifussionValuesToCheck <- DifussionValuesToCheck[DifussionValuesToCheck <= 1]

bases <- c(1, 2.5, 5, 7.5)
exponentes <- -5:0
valores_pos <- as.vector(outer(bases, 10^exponentes))
valores_pos <- valores_pos[valores_pos <= 1]
SelectionValuesToCheck <- sort(unique(c(0, valores_pos)))

# ============================================================
# GRID REDUCIDO PARA EL PERFIL DE T_onset
# (más grueso a propósito: se evalúa una vez por cada candidato
# T_onset, así que el costo total es candidatos x este grid)
# Ajustar aquí si luego se quiere afinar la resolución.
# ============================================================
DifussionValuesToCheck_Tonset <- sort(unique(c(0, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.25, 0.5, 1)))
SelectionValuesToCheck_Tonset <- sort(unique(c(0, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.25, 0.5, 1)))
param_grid_Tonset <- expand.grid(D = DifussionValuesToCheck_Tonset, s = SelectionValuesToCheck_Tonset)

# Para la fase NEUTRAL (s = 0 forzado) sólo se varía D
param_grid_neutral_Tonset <- data.frame(D = DifussionValuesToCheck_Tonset, s = 0)

piso_minimo <- 1e-6

# ============================================================
# Helper: corre la PDE para un set de parámetros (D, s) sobre una
# ventana temporal "times_run", partiendo de Conc0, y devuelve el
# array de concentraciones (tiempos x n x n) ya aplanado a matriz
# (tiempos x n*n), tal como en el script original.
# ============================================================
run_pde <- function(Conc0, times_run, D_curr, s_curr) {
  d_curr <- 2 * s_curr
  ST3 <- ode.2D(y = Conc0, times = times_run, func = diffusion2D,
                parms = c(D_curr, s_curr, d_curr), dimens = c(n, n),
                method = "rk4", atol = 1e-7, rtol = 1e-7)
  ST3_mat <- as.matrix(ST3[, -1])
  ST3_mat
}

# ============================================================
# Helper: dada una matriz de salida de la PDE (tiempos x n*n) y los
# índices pre-calculados de tiempo/espacio, calcula la log-verosimilitud
# binomial. Devuelve NA si la PDE produjo valores no finitos.
# ============================================================
loglik_from_pde <- function(ST3_mat, idx_matrix, counts_vec, totals_vec) {
  if (any(!is.finite(ST3_mat))) return(NA_real_)
  pred_freq_raw <- ST3_mat[idx_matrix]
  pred_freq_vec <- pmin(pmax(pred_freq_raw, piso_minimo), 1 - piso_minimo)
  sum(dbinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, log = TRUE))
}

# Lectura y preparación de datos

freq_data_raw <- fread(freq_file, header = TRUE)
colnames(freq_data_raw) <- c('Chrom_frq', 'MutationID', 'CLST', 'A1', 'A2', 'Frequency', 'Count', 'Chr_Tot')

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

# Lista de SNPS para analizar
snps_to_analyze <- unique(sort(freq_data$MutationID))
tiempos_globales_muestreo <- sort(unique(freq_data$Generation))

all_results_legacy <- list()      # resumen: mejor (D,s) bajo modelo neutral puro, como en el script original
all_results_Tonset_best <- list() # resumen: mejor T_onset por SNP (D,s en cada fase, LL total)
all_results_Tonset_profile <- list() # perfil COMPLETO: LL total para cada T_onset candidato (para graficar)

# ============================================================
# INFERENCIA DE D, s, y T_onset (punto de inicio de la selección)
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

  # --- CONFIGURACIÓN GENERAL DE LA VENTANA TEMPORAL ---
  times_run_full <- AlleleOriginAge:max(df_snp$Generation)

  if (length(times_run_full) < 2) {
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

  # --- PRE-CALCULO DE INDICES PARA TODA LA VENTANA (una sola vez) ---
  time_idx_vec_full    <- match(df_snp$Generation, times_run_full)
  spatial_idx_vec_full <- (df_snp$Y - 1) * n + df_snp$X

  filas_validas_full <- !is.na(time_idx_vec_full)
  if (!all(filas_validas_full)) {
    cat("SNP", snp_actual, ": se omiten", sum(!filas_validas_full), "filas cuya Generation no está en times_run_full.\n")
  }

  counts_vec_full <- df_snp$Count[filas_validas_full]
  totals_vec_full <- df_snp$Chr_Tot[filas_validas_full]
  idx_matrix_full  <- cbind(time_idx_vec_full[filas_validas_full], spatial_idx_vec_full[filas_validas_full])

  cat("Evaluando SNP:", snp_actual, "- Ventana temporal:", length(times_run_full), "generaciones.\n")

  # ============================================================
  # PASO 1: GRID SEARCH "LEGACY" — un único régimen (D,s) en toda
  # la ventana, igual que el script original. Se conserva por
  # compatibilidad y como referencia comparativa frente a T_onset.
  # ============================================================
  param_grid <- expand.grid(D = DifussionValuesToCheck, s = SelectionValuesToCheck)
  ll_legacy_vals <- numeric(nrow(param_grid))

  for (i in 1:nrow(param_grid)) {
    ST3_mat <- run_pde(Conc0, times_run_full, param_grid$D[i], param_grid$s[i])
    ll_legacy_vals[i] <- loglik_from_pde(ST3_mat, idx_matrix_full, counts_vec_full, totals_vec_full)
  }

  col_names <- paste0("D_", param_grid$D, "_s_", param_grid$s)
  all_results_legacy[[length(all_results_legacy) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_legacy_vals)), c("SNP", col_names))

  # ============================================================
  # PASO 2: PERFIL DE VEROSIMILITUD SOBRE T_onset
  #
  # Candidatos T_onset = generaciones observadas para ESTE SNP
  # (densidad de datos propia del SNP), incluyendo el caso límite
  # T_onset = AlleleOriginAge (selección desde el inicio, sin fase
  # neutral previa). Se excluye el último tiempo observado, porque
  # ahí no quedaría fase de selección con datos para ajustar.
  # ============================================================
  generaciones_snp <- sort(unique(df_snp$Generation))
  T_onset_candidatos <- generaciones_snp[generaciones_snp < max(generaciones_snp)]

  if (length(T_onset_candidatos) == 0) {
    cat("SNP", snp_actual, ": sin candidatos válidos de T_onset, se omite el perfil.\n")
    next
  }

  perfil_Tonset <- vector("list", length(T_onset_candidatos))

  for (ti in seq_along(T_onset_candidatos)) {

    T_onset_cand <- T_onset_candidatos[ti]

    # --- FASE NEUTRAL: desde AlleleOriginAge hasta T_onset_cand, s=0 ---
    if (T_onset_cand > AlleleOriginAge) {

      times_run_neutral <- AlleleOriginAge:T_onset_cand
      time_idx_vec_n    <- match(df_snp$Generation[df_snp$Generation <= T_onset_cand], times_run_neutral)
      rows_neutral       <- df_snp$Generation <= T_onset_cand
      spatial_idx_vec_n <- (df_snp$Y[rows_neutral] - 1) * n + df_snp$X[rows_neutral]
      filas_validas_n   <- !is.na(time_idx_vec_n)

      counts_vec_n <- df_snp$Count[rows_neutral][filas_validas_n]
      totals_vec_n <- df_snp$Chr_Tot[rows_neutral][filas_validas_n]
      idx_matrix_n <- cbind(time_idx_vec_n[filas_validas_n], spatial_idx_vec_n[filas_validas_n])

      ll_neutral_vals <- numeric(nrow(param_grid_neutral_Tonset))
      ST3_mats_neutral <- vector("list", nrow(param_grid_neutral_Tonset))

      for (i in 1:nrow(param_grid_neutral_Tonset)) {
        ST3_mat <- run_pde(Conc0, times_run_neutral, param_grid_neutral_Tonset$D[i], 0)
        ST3_mats_neutral[[i]] <- ST3_mat
        ll_neutral_vals[i] <- loglik_from_pde(ST3_mat, idx_matrix_n, counts_vec_n, totals_vec_n)
      }

      best_neutral_i <- which.max(ll_neutral_vals)
      if (length(best_neutral_i) == 0 || is.na(ll_neutral_vals[best_neutral_i])) {
        # Toda la fase neutral fue inválida (PDE no finita) para este candidato
        perfil_Tonset[[ti]] <- data.frame(
          SNP = snp_actual, T_onset = T_onset_cand,
          D_neutral = NA, LL_neutral = NA,
          D_seleccion = NA, s_seleccion = NA, LL_seleccion = NA,
          LL_total = NA
        )
        next
      }

      D_neutral_best  <- param_grid_neutral_Tonset$D[best_neutral_i]
      LL_neutral_best <- ll_neutral_vals[best_neutral_i]

      # Estado final de la fase neutral = condición inicial de la fase de selección
      Conc_at_Tonset <- matrix(ST3_mats_neutral[[best_neutral_i]][nrow(ST3_mats_neutral[[best_neutral_i]]), ], nrow = n, ncol = n)

    } else {
      # Caso límite: T_onset == AlleleOriginAge -> sin fase neutral previa
      D_neutral_best  <- NA
      LL_neutral_best <- 0
      Conc_at_Tonset  <- Conc0
    }

    # --- FASE DE SELECCIÓN: desde T_onset_cand hasta el final ---
    times_run_seleccion <- T_onset_cand:max(df_snp$Generation)

    if (length(times_run_seleccion) < 2) {
      perfil_Tonset[[ti]] <- data.frame(
        SNP = snp_actual, T_onset = T_onset_cand,
        D_neutral = D_neutral_best, LL_neutral = LL_neutral_best,
        D_seleccion = NA, s_seleccion = NA, LL_seleccion = NA,
        LL_total = NA
      )
      next
    }

    rows_seleccion <- df_snp$Generation >= T_onset_cand
    time_idx_vec_s    <- match(df_snp$Generation[rows_seleccion], times_run_seleccion)
    spatial_idx_vec_s <- (df_snp$Y[rows_seleccion] - 1) * n + df_snp$X[rows_seleccion]
    filas_validas_s   <- !is.na(time_idx_vec_s)

    counts_vec_s <- df_snp$Count[rows_seleccion][filas_validas_s]
    totals_vec_s <- df_snp$Chr_Tot[rows_seleccion][filas_validas_s]
    idx_matrix_s <- cbind(time_idx_vec_s[filas_validas_s], spatial_idx_vec_s[filas_validas_s])

    # Grid search reducido (D,s) para la fase de selección
    ll_seleccion_vals <- numeric(nrow(param_grid_Tonset))
    for (i in 1:nrow(param_grid_Tonset)) {
      ST3_mat <- run_pde(Conc_at_Tonset, times_run_seleccion, param_grid_Tonset$D[i], param_grid_Tonset$s[i])
      ll_seleccion_vals[i] <- loglik_from_pde(ST3_mat, idx_matrix_s, counts_vec_s, totals_vec_s)
    }

    best_sel_i <- which.max(ll_seleccion_vals)

    if (length(best_sel_i) == 0 || is.na(ll_seleccion_vals[best_sel_i])) {
      perfil_Tonset[[ti]] <- data.frame(
        SNP = snp_actual, T_onset = T_onset_cand,
        D_neutral = D_neutral_best, LL_neutral = LL_neutral_best,
        D_seleccion = NA, s_seleccion = NA, LL_seleccion = NA,
        LL_total = NA
      )
      next
    }

    # --- Refinamiento continuo con L-BFGS-B alrededor del óptimo del grid reducido ---
    D0 <- param_grid_Tonset$D[best_sel_i]
    s0 <- param_grid_Tonset$s[best_sel_i]

    neg_ll_seleccion <- function(par) {
      D_try <- par[1]; s_try <- par[2]
      ST3_mat <- run_pde(Conc_at_Tonset, times_run_seleccion, D_try, s_try)
      ll <- loglik_from_pde(ST3_mat, idx_matrix_s, counts_vec_s, totals_vec_s)
      if (is.na(ll) || !is.finite(ll)) return(1e10)
      -ll
    }

    opt_res <- tryCatch({
      optim(par = c(max(D0, 1e-6), s0), fn = neg_ll_seleccion, method = "L-BFGS-B",
            lower = c(0, 0), upper = c(1, 1))
    }, error = function(e) NULL)

    if (!is.null(opt_res) && is.finite(opt_res$value)) {
      D_sel_best  <- opt_res$par[1]
      s_sel_best  <- opt_res$par[2]
      LL_sel_best <- -opt_res$value
    } else {
      D_sel_best  <- D0
      s_sel_best  <- s0
      LL_sel_best <- ll_seleccion_vals[best_sel_i]
    }

    LL_total <- LL_neutral_best + LL_sel_best

    perfil_Tonset[[ti]] <- data.frame(
      SNP = snp_actual, T_onset = T_onset_cand,
      D_neutral = D_neutral_best, LL_neutral = LL_neutral_best,
      D_seleccion = D_sel_best, s_seleccion = s_sel_best, LL_seleccion = LL_sel_best,
      LL_total = LL_total
    )
  }

  perfil_df <- do.call(rbind, perfil_Tonset)
  all_results_Tonset_profile[[length(all_results_Tonset_profile) + 1]] <- perfil_df

  # Mejor T_onset = el que maximiza la verosimilitud TOTAL del perfil
  if (all(is.na(perfil_df$LL_total))) {
    cat("SNP", snp_actual, ": no se pudo estimar un T_onset válido (todas las fases fallaron).\n")
    next
  }

  best_row <- perfil_df[which.max(perfil_df$LL_total), ]
  all_results_Tonset_best[[length(all_results_Tonset_best) + 1]] <- best_row

  cat("  -> SNP", snp_actual, ": T_onset estimado =", best_row$T_onset,
      "| D_neutral =", best_row$D_neutral,
      "| D_sel =", round(best_row$D_seleccion, 5),
      "| s_sel =", round(best_row$s_seleccion, 5),
      "| LL_total =", round(best_row$LL_total, 3), "\n")
}

# --- GUARDADO ---

# 1) Grid LEGACY (igual que el script original, para comparación)
if (length(all_results_legacy) > 0) {
  df_legacy <- do.call(rbind, all_results_legacy)
  out_legacy <- file.path(output_dir, paste0("Analysis_LL_Ventajosos_", task_id, "_neutrales", ".txt"))
  write.table(df_legacy, file = out_legacy, row.names = FALSE, quote = FALSE, sep = "\t")
  print(paste("Guardado Grid:", out_legacy))
} else {
  print("No se procesaron SNPs para Legacy.")
}

# 2) Mejor T_onset por SNP (resumen, una fila por SNP)
if (length(all_results_Tonset_best) > 0) {
  df_Tonset_best <- do.call(rbind, all_results_Tonset_best)
  out_Tonset_best <- file.path(output_dir, paste0("Analysis_Tonset_BestEstimate_", task_id, ".txt"))
  write.table(df_Tonset_best, file = out_Tonset_best, row.names = FALSE, quote = FALSE, sep = "\t")
  print(paste("Guardado mejor T_onset por SNP:", out_Tonset_best))
} else {
  print("No se procesaron SNPs para T_onset.")
}

# 3) Perfil COMPLETO de verosimilitud vs T_onset (formato largo, ideal para graficar)
#    Columnas: SNP, T_onset, D_neutral, LL_neutral, D_seleccion, s_seleccion, LL_seleccion, LL_total
#    -> para cada SNP se puede graficar LL_total (eje Y) vs T_onset (eje X)
if (length(all_results_Tonset_profile) > 0) {
  df_Tonset_profile <- do.call(rbind, all_results_Tonset_profile)
  out_Tonset_profile <- file.path(output_dir, paste0("Analysis_Tonset_Profile_", task_id, ".txt"))
  write.table(df_Tonset_profile, file = out_Tonset_profile, row.names = FALSE, quote = FALSE, sep = "\t")
  print(paste("Guardado perfil de verosimilitud T_onset:", out_Tonset_profile))
} else {
  print("No se generó perfil de T_onset.")
}