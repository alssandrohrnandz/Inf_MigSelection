library(deSolve)
library(data.table)
library(VGAM)

args <- commandArgs(trailingOnly = TRUE)

freq_file   <- args[3]
pop_file    <- args[2]
snp_file    <- "AADNA_data/v66.1240K.aadr.PUB.snp"
task_id     <- args[5]
model_name  <- args[6]
output_dir  <- args[7]





if(length(args) < 5) {
  freq_file   <- "results/unificado_CHR1_selection.csv"
  pop_file    <- "Population_info.txt"
  snp_file    <- "AADNA_data/v66.1240K.aadr.PUB.snp"
  task_id     <- "TEST_OPTIM"
  model_name  <- "D_FULL_selection"
  output_dir  <- getwd()
}

cat("=========================================\n")
cat("Iniciando Pipeline con Optimización de 's'\n")
cat("Modelo:", model_name, "| Task ID:", task_id, "\n")
cat("=========================================\n")

# ==================== FUNCIÓN DE DETECCIÓN DE TIPO DE DATOS ====================
detect_data_type <- function(freq_file, model_name) {
  # Detección 1: Por nombre de archivo (prioridad alta)
  if(grepl("^v2_D_FULL|FULL_[a-z]+_m[12]_", basename(freq_file), ignore.case = TRUE)) {
    cat("✓ Detectado: SIMULACIÓN SLiM\n")
    return("SLIM")
  }
  
  # Detección 2: Por estructura de columnas (fallback)
  tryCatch({
    sample_data <- fread(freq_file, nrows = 1, colClasses = "character")
    col_names <- tolower(names(sample_data))
    
    # Columnas esperadas en SLiM: Generation, MutationID, X, Y, Frequency, AlleleCount, Chr_Tot
    slim_cols <- c("generation", "mutationid", "x", "y", "frequency", "allelecount", "chr_tot")
    slim_detected <- sum(slim_cols %in% col_names) >= 5
    
    # Columnas esperadas en aDNA: CLST, Media_Lat, Media_Long, Media_YBP, MAC, NCHROBS
    adna_cols <- c("clst", "media_lat", "media_long", "media_ybp", "mac", "nchrobs")
    adna_detected <- sum(adna_cols %in% col_names) >= 4
    
    if(slim_detected) {
      cat("✓ Detectado: SIMULACIÓN SLiM (por estructura de columnas)\n")
      return("SLIM")
    } else if(adna_detected) {
      cat("✓ Detectado: DATOS REALES aDNA\n")
      return("aDNA")
    }
  }, error = function(e) {
    cat("⚠ No se pudo detectar por estructura; asumiendo aDNA\n")
  })
  
  cat("✓ Detectado: DATOS REALES aDNA (default)\n")
  return("aDNA")
}

# ==================== PROCESAMIENTO DE DATOS SEGÚN TIPO ====================
process_frequency_data <- function(freq_file, pop_file, data_type, nx, ny) {
  
  if(data_type == "SLIM") {
    cat("\n--- PROCESANDO DATOS SLiM ---\n")
    
    # Carga directa de SLiM
    freq_data <- fread(freq_file)
    
    # Validar columnas esperadas
    required_cols <- c("Generation", "MutationID", "X", "Y", "Frequency", "AlleleCount", "Chr_Tot")
    missing_cols <- setdiff(required_cols, names(freq_data))
    if(length(missing_cols) > 0) {
      stop(paste("Columnas faltantes en datos SLiM:", paste(missing_cols, collapse=", ")))
    }
    
    # Mapeo directo: SLiM ya tiene estructura de grid
    grid_data <- freq_data[, .(
      MutationID = as.character(MutationID),
      Generation = as.integer(Generation),
      X = as.integer(X),
      Y = as.integer(Y),
      Frequency = as.numeric(Frequency),
      Count = as.numeric(AlleleCount),
      Chr_Tot = as.numeric(Chr_Tot),
      ChrOBS = as.numeric(Chr_Tot)
    )]
    
    cat("✓ Datos SLiM cargados:", nrow(grid_data), "registros\n")
    cat("  Rango Generaciones:", min(grid_data$Generation), "-", max(grid_data$Generation), "\n")
    cat("  SNPs únicos:", length(unique(grid_data$MutationID)), "\n")
    cat("  Grid: X =", min(grid_data$X), "-", max(grid_data$X), 
        "| Y =", min(grid_data$Y), "-", max(grid_data$Y), "\n")
    
    return(grid_data)
    
  } else if(data_type == "aDNA") {
    cat("\n--- PROCESANDO DATOS aDNA REALES ---\n")
    
    # Carga de datos reales con mapeo geográfico-genético
    freq_data <- fread(freq_file)
    pop_info <- fread(pop_file)
    setnames(pop_info, "Group_Master_ID", "CLST", skip_absent = TRUE)
    
    merged_df <- merge(freq_data, pop_info[, .(CLST, Media_Lat, Media_Long)], 
                       by = "CLST", all.x = TRUE)
    merged_df[, c("Media_Lat.y", "Media_Long.y") := NULL]
    
    setnames(merged_df, 
             old = c("Media_Lat.x", "Media_Long.x"), 
             new = c("Media_Lat", "Media_Long"))
    
    # Filtrar por región Eurasia
    eurasia_df <- merged_df[Media_Lat >= 30 & Media_Lat <= 75 & 
                            Media_Long >= -15 & Media_Long <= 55]
    
    eurasia_df[, Media_Long := as.numeric(Media_Long)]
    eurasia_df[, Media_Lat := as.numeric(Media_Lat)]
    eurasia_df[, X := as.integer(cut(Media_Long, 
                                      breaks = seq(-15, 55, length.out = nx + 1), 
                                      include.lowest = TRUE))]
    eurasia_df[, Y := as.integer(cut(Media_Lat, 
                                      breaks = seq(30, 75, length.out = ny + 1), 
                                      include.lowest = TRUE))]
    
    # Convertir YBP a Generation (27 años por generación)
    eurasia_df[, Generation := floor((11700 - Media_YBP) / 27) + 1]
    eurasia_df <- eurasia_df[Generation >= 1]
    
    # Agregar por SNP-Generation-GridCell
    grid_data <- eurasia_df[, .(
      Frequency = sum(MAC, na.rm=TRUE) / sum(NCHROBS, na.rm=TRUE), 
      Count = sum(MAC, na.rm=TRUE), 
      Chr_Tot = sum(NCHROBS, na.rm=TRUE)
    ), by = .(SNP, Generation, X, Y)]
    
    setnames(grid_data, "SNP", "MutationID")
    grid_data$ChrOBS <- grid_data$Chr_Tot
    
    cat("✓ Datos aDNA cargados y mapeados\n")
    cat("  Registros procesados:", nrow(eurasia_df), "\n")
    cat("  Registros en grid final:", nrow(grid_data), "\n")
    cat("  Rango Generaciones:", min(grid_data$Generation), "-", max(grid_data$Generation), "\n")
    cat("  SNPs únicos:", length(unique(grid_data$MutationID)), "\n")
    
    return(grid_data)
  }
}

# ==================== FUNCIÓN PDE ====================
diffusion2D <- function(t, conc, par) {
  Conc <- matrix(nrow = nx, ncol = ny, data = conc)
  D_eff <- par[1]
  
  if (t >= par[4]) {
    d_eff <- par[2]
    s_eff <- par[3]
  } else {
    d_eff <- 0
    s_eff <- 0
  }
  
  dConc <- Conc*(1-Conc)*(Conc*d_eff + s_eff*(1-2*Conc))
  FluxX <- -D_eff * rbind(rep(0, ny), (Conc[2:nx,] - Conc[1:(nx-1),]), rep(0, ny)) / dx
  dConc <- dConc - (FluxX[2:(nx+1),] - FluxX[1:nx,]) / dx
  FluxY <- -D_eff * cbind(rep(0, nx), (Conc[,2:ny] - Conc[,1:(ny-1)]), rep(0, nx)) / dy
  dConc <- dConc - (FluxY[,2:(ny+1)] - FluxY[,1:ny]) / dy
  
  return(list(as.vector(dConc)))
}

# ==================== FUNCIONES DE VEROSIMILITUD ====================
dbetabinom_spikes <- function(x, size, prob, rho, log = FALSE) {
  pi_0 <- (1 - prob) * rho; pi_1 <- prob * rho
  pi_mid <- max(1 - pi_0 - pi_1, 1e-10)
  ll_bb <- dbetabinom(x = x, size = size, prob = prob, rho = rho, log = FALSE)
  dens <- pi_mid * ll_bb
  dens[x == 0] <- dens[x == 0] + pi_0
  dens[x == size] <- dens[x == size] + pi_1
  dens <- pmax(dens, 1e-300)
  if (log) return(log(dens)) else return(dens)
}

# ==================== CONFIGURACIÓN PRINCIPAL ====================
nx <- 31; ny <- 9; dx <- 1; dy <- 1; N_eff <- 1000

# Detectar tipo de datos
data_type <- detect_data_type(freq_file, model_name)
is_neutral <- grepl("neutros", model_name, ignore.case = TRUE)
MAX_PADDED_GENERATIONS <- ifelse(is_neutral, 1, 2)

# Grid estático para D
exponentes_D <- -3:0; valores_base <- c(1, 5)
DifussionValuesToCheck <- sort(unique(c(0, as.vector(outer(valores_base, 10^exponentes_D)))))
DifussionValuesToCheck <- DifussionValuesToCheck[DifussionValuesToCheck <= 1]

# ==================== PROCESAMIENTO DE DATOS ====================
grid_data <- process_frequency_data(freq_file, pop_file, data_type, nx, ny)

tiempos_globales_muestreo <- sort(unique(grid_data$Generation))
snps_to_analyze <- sort(unique(grid_data$MutationID))

# ==================== CARGA DE INFORMACIÓN SNP ====================
if(data_type == "SLIM") {
  # Para SLiM, creamos un snp_info dummy (opcional, ya tenemos MutationID)
  cat("\n--- Procesando metadatos SNP ---\n")
  snp_info <- data.table(SNP = unique(grid_data$MutationID))
  snp_info[, CHR := 1]  # SLiM no tiene info de cromosoma, asignamos default
  snp_info[, POS := seq_len(.N)]
  cat("✓ SNP info creada para", nrow(snp_info), "mutaciones\n")
} else {
  # Para aDNA, cargar archivo SNP real
  cat("\n--- Cargando metadatos SNP reales ---\n")
  snp_info <- fread(snp_file, select = c(1, 2, 4), col.names = c("SNP", "CHR", "POS"))
  cat("✓ SNP info cargada:", nrow(snp_info), "SNPs\n")
}

# ==================== OPTIMIZACIÓN DE PARÁMETROS ====================
all_results_legacy <- list(); all_results_ares <- list(); all_results_spikes <- list()

cat("\n--- Iniciando optimización de 's' para", length(snps_to_analyze), "SNPs ---\n\n")

for (snp_idx in seq_along(snps_to_analyze)) {
  snp_actual <- snps_to_analyze[snp_idx]
  
  df_snp <- grid_data[MutationID == snp_actual]
  df_snp <- df_snp[!is.na(Frequency) & !is.nan(Frequency)]
  if(nrow(df_snp) == 0 || max(df_snp$Frequency, na.rm=TRUE) == 0) next
  
  df_snp <- df_snp[order(Generation)]
  AlleleOriginAge <- df_snp$Generation[which(df_snp$Frequency > 0)[1]]
  max_gen_snp <- max(df_snp$Generation)
  
  # Padding de generaciones posteriores (para datos incompletos)
  if (max_gen_snp < max(tiempos_globales_muestreo)) {
    gens_faltantes <- tiempos_globales_muestreo[tiempos_globales_muestreo > max_gen_snp]
    if(length(gens_faltantes) > MAX_PADDED_GENERATIONS) 
      gens_faltantes <- gens_faltantes[1:MAX_PADDED_GENERATIONS]
    filas_ceros <- as.data.table(expand.grid(
      Generation = gens_faltantes, 
      X = unique(df_snp$X), 
      Y = unique(df_snp$Y)))
    filas_ceros[, `:=`(
      MutationID = snp_actual, 
      Frequency = 0, Count = 0, 
      Chr_Tot = max(round(mean(df_snp$Chr_Tot, na.rm = T)), 10), 
      ChrOBS = max(round(mean(df_snp$Chr_Tot, na.rm = T)), 10)
    )]
    df_snp <- rbind(df_snp, filas_ceros, use.names=TRUE)[order(Generation, Y, X)]
  }
  
  times_run <- min(df_snp$Generation):max(df_snp$Generation)
  if(length(times_run) < 2) next
  
  # Condición inicial
  Conc0 <- matrix(0, nrow=nx, ncol=ny)
  origin_data <- df_snp[Generation == AlleleOriginAge & Frequency > 0]
  for(k in 1:nrow(origin_data)) {
    if(origin_data$X[k] >= 1 && origin_data$X[k] <= nx && 
       origin_data$Y[k] >= 1 && origin_data$Y[k] <= ny) {
      Conc0[origin_data$X[k], origin_data$Y[k]] <- origin_data$Frequency[k]
    }
  }
  
  tiempos_onset_probar <- if(is_neutral) AlleleOriginAge else 
    round(seq(AlleleOriginAge, max_gen_snp, length.out = 4))
  
  local_grid <- expand.grid(D = DifussionValuesToCheck, T_onset = tiempos_onset_probar)
  
  counts_vec <- df_snp$Count; totals_vec <- df_snp$Chr_Tot
  rho_ares_vec <- numeric(nrow(df_snp))
  for(j in 1:nrow(df_snp)) 
    rho_ares_vec[j] <- max(1 - exp(-max(df_snp$Generation[j] - AlleleOriginAge, 0.5) / 
                                    (2 * N_eff)), 1e-6)
  
  cat("[", snp_idx, "/", length(snps_to_analyze), "] Optimizando SNP:", snp_actual, "\n")
  
  # Estructuras para resultados
  res_leg <- list(); res_ares <- list(); res_spikes <- list()
  
  for (i in 1:nrow(local_grid)) {
    D_curr <- local_grid$D[i]; T_onset_curr <- local_grid$T_onset[i]
    
    # Función objetivo para optimizar 's'
    eval_s_optim <- function(s_val, model_type="LEGACY") {
      d_val <- 2 * s_val
      
      # Solucionar PDE
      ST3 <- tryCatch({
        ode.2D(y = Conc0, times = times_run, func = diffusion2D, 
               parms = c(D_curr, d_val, s_val, T_onset_curr), 
               dimens = c(nx, ny), method = "rk4", atol = 1e-6, rtol = 1e-6)
      }, error = function(e) return(NULL))
      
      if(is.null(ST3)) return(1e9)
      
      ST3_mat <- as.matrix(ST3[,-1])
      pred_freq_vec <- numeric(nrow(df_snp))
      for(j in 1:nrow(df_snp)) {
        time_idx <- match(df_snp$Generation[j], times_run)
        spatial_idx <- (df_snp$Y[j] - 1) * nx + df_snp$X[j] 
        pred_freq_vec[j] <- max(min(ST3_mat[time_idx, spatial_idx], 1 - 1e-6), 1e-6)
      }
      
      ll <- switch(model_type,
                   "LEGACY" = sum(dbinom(x = counts_vec, size = totals_vec, 
                                        prob = pred_freq_vec, log = TRUE)),
                   "ARES"   = sum(dbetabinom(x = counts_vec, size = totals_vec, 
                                            prob = pred_freq_vec, rho = rho_ares_vec, log = TRUE)),
                   "SPIKES" = sum(dbetabinom_spikes(x = counts_vec, size = totals_vec, 
                                                   prob = pred_freq_vec, rho = rho_ares_vec, log = TRUE)))
      
      return(ifelse(is.finite(ll), -ll, 1e9))
    }
    
    if (is_neutral) {
      ll_leg <- -eval_s_optim(0, "LEGACY")
      ll_ar <- -eval_s_optim(0, "ARES")
      ll_spk <- -eval_s_optim(0, "SPIKES")
      res_leg[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=0, LL=ll_leg, SE_s=0)
      res_ares[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=0, LL=ll_ar, SE_s=0)
      res_spikes[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=0, LL=ll_spk, SE_s=0)
    } else {
      opt_leg <- optim(par=0.01, fn=eval_s_optim, model_type="LEGACY", 
                       method="L-BFGS-B", lower=-1, upper=1, hessian=TRUE)
      se_leg <- ifelse(opt_leg$hessian[1,1] > 0, sqrt(1/opt_leg$hessian[1,1]), NA)
      res_leg[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, 
                                 s_opt=opt_leg$par, LL=-opt_leg$value, SE_s=se_leg)
      
      opt_ar <- optim(par=0.01, fn=eval_s_optim, model_type="ARES", 
                      method="L-BFGS-B", lower=-1, upper=1, hessian=TRUE)
      se_ar <- ifelse(opt_ar$hessian[1,1] > 0, sqrt(1/opt_ar$hessian[1,1]), NA)
      res_ares[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, 
                                  s_opt=opt_ar$par, LL=-opt_ar$value, SE_s=se_ar)
      
      opt_spk <- optim(par=0.01, fn=eval_s_optim, model_type="SPIKES", 
                       method="L-BFGS-B", lower=-1, upper=1, hessian=TRUE)
      se_spk <- ifelse(opt_spk$hessian[1,1] > 0, sqrt(1/opt_spk$hessian[1,1]), NA)
      res_spikes[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, 
                                    s_opt=opt_spk$par, LL=-opt_spk$value, SE_s=se_spk)
    }
  }
  
  # Consolidar resultados
  dt_leg <- rbindlist(res_leg); dt_leg[, SNP := snp_actual]
  dt_ares <- rbindlist(res_ares); dt_ares[, SNP := snp_actual]
  dt_spikes <- rbindlist(res_spikes); dt_spikes[, SNP := snp_actual]
  
  all_results_legacy[[length(all_results_legacy)+1]] <- dt_leg
  all_results_ares[[length(all_results_ares)+1]] <- dt_ares
  all_results_spikes[[length(all_results_spikes)+1]] <- dt_spikes
}

# ==================== GUARDADO DE RESULTADOS ====================
cat("\n--- Guardando resultados finales ---\n")

save_results <- function(results_list, model_type, snp_info_dt, data_type) {
  if(length(results_list) == 0) {
    cat("⚠ No hay resultados para", model_type, "\n")
    return(NULL)
  }
  
  dt <- rbindlist(results_list)
  
  # Merge con info SNP
  if(data_type == "SLIM") {
    setnames(dt, "SNP", "SNP_temp")
    dt_final <- merge(snp_info_dt, dt, by.x="SNP", by.y="SNP_temp", all.y=TRUE)
  } else {
    dt_final <- merge(snp_info_dt, dt, by = "SNP")
  }
  
  # Calcular IC 95%
  dt_final[, CI_95_Lower := s_opt - (1.96 * SE_s)]
  dt_final[, CI_95_Upper := s_opt + (1.96 * SE_s)]
  
  out_name <- file.path(output_dir, 
                        paste0("OPTIM_", model_type, "_ALL_TaskID_", task_id, 
                               "_", model_name, ".txt"))
  fwrite(dt_final, file = out_name, sep = "\t", quote = FALSE)
  cat("✓ Guardado:", out_name, "\n")
}

save_results(all_results_legacy, "LEGACY", snp_info, data_type)
save_results(all_results_ares, "ARES", snp_info, data_type)
save_results(all_results_spikes, "SPIKES", snp_info, data_type)

cat("\n=========================================\n")
cat("✓ EJECUCIÓN COMPLETADA\n")
cat("Tipo de datos:", data_type, "\n")
cat("SNPs analizados:", length(snps_to_analyze), "\n")
cat("=========================================\n")
