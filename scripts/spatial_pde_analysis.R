library(deSolve)
library(data.table)
library(VGAM)

args <- commandArgs(trailingOnly = TRUE)
freq_file   <- args[1]
pop_file    <- args[2]
snp_file    <- args[3]
task_id     <- args[4]
model_name  <- args[5]
output_dir  <- args[6]

if (length(args) < 6) {
  freq_file  <- "results/unificado_CHR1_selection.csv"
  pop_file   <- "Population_info.txt"
  snp_file   <- "AADNA_data/v66.1240K.aadr.PUB.snp"
  task_id    <- "TEST_OPTIM"
  model_name <- "D_FULL_selection"
  output_dir <- getwd()
}

cat("=========================================\n")
cat("Iniciando Pipeline con Optimización de 's'\n")
cat("Modelo:", model_name, "| Task ID:", task_id, "\n")
cat("=========================================\n")

# ── Definir out_files PRIMERO, luego limpiar ──────────────────────────────────
out_files <- list(
  LEGACY = file.path(output_dir, paste0("OPTIM_LEGACY_ALL_TaskID_", task_id, "_", model_name, ".txt")),
  ARES   = file.path(output_dir, paste0("OPTIM_ARES_ALL_TaskID_",   task_id, "_", model_name, ".txt")),
  SPIKES = file.path(output_dir, paste0("OPTIM_SPIKES_ALL_TaskID_", task_id, "_", model_name, ".txt"))
)

for (f in out_files) {
  if (file.exists(f)) file.remove(f)
}

# ── Cargar metadatos de SNPs para el merge final ──────────────────────────────
snp_info <- fread(snp_file, select = c(1, 2, 4), col.names = c("SNP", "CHR", "POS"))

# ── Función de escritura incremental (una vez por SNP) ────────────────────────
finalize_snp <- function(res_list, model_type, snp_id) {
  dt <- rbindlist(res_list)
  dt[, SNP := snp_id]
  dt[, CI_95_Lower := s_opt - 1.96 * SE_s]
  dt[, CI_95_Upper := s_opt + 1.96 * SE_s]
  
  # Agregar CHR y POS
  dt <- merge(snp_info[SNP == snp_id], dt, by = "SNP")
  
  target <- out_files[[model_type]]
  fwrite(dt, file = target, sep = "\t", quote = FALSE,
         append = TRUE,
         col.names = !file.exists(target))
}

# ── Parámetros globales ───────────────────────────────────────────────────────
nx <- 31; ny <- 9; dx <- 1; dy <- 1; N_eff <- 2000

diffusion2D <- function(t, conc, par) {
  Conc  <- matrix(nrow = nx, ncol = ny, data = conc)
  D_eff <- par[1]; d_eff <- par[2]; s_eff <- par[3]
  dConc <- Conc * (1 - Conc) * (Conc * d_eff + s_eff * (1 - 2 * Conc))
  FluxX <- -D_eff * rbind(rep(0, ny), (Conc[2:nx,] - Conc[1:(nx-1),]), rep(0, ny)) / dx
  dConc <- dConc - (FluxX[2:(nx+1),] - FluxX[1:nx,]) / dx
  FluxY <- -D_eff * cbind(rep(0, nx), (Conc[,2:ny] - Conc[,1:(ny-1)]), rep(0, nx)) / dy
  dConc <- dConc - (FluxY[,2:(ny+1)] - FluxY[,1:ny]) / dy
  return(list(as.vector(dConc)))
}

dbetabinom_spikes <- function(x, size, prob, rho, log = FALSE) {
  pi_0   <- (1 - prob) * rho
  pi_1   <- prob * rho
  pi_mid <- max(1 - pi_0 - pi_1, 1e-10)
  ll_bb  <- dbetabinom(x = x, size = size, prob = prob, rho = rho, log = FALSE)
  dens   <- pi_mid * ll_bb
  dens[x == 0]    <- dens[x == 0]    + pi_0
  dens[x == size] <- dens[x == size] + pi_1
  dens   <- pmax(dens, 1e-300)
  if (log) return(log(dens)) else return(dens)
}

is_neutral <- grepl("neutral", model_name, ignore.case = TRUE)
MAX_PADDED_GENERATIONS <- ifelse(is_neutral, 1, 2)

exponentes_D <- -3:0; valores_base <- c(1, 5)
DifussionValuesToCheck <- sort(unique(c(0, as.vector(outer(valores_base, 10^exponentes_D)))))
DifussionValuesToCheck <- DifussionValuesToCheck[DifussionValuesToCheck <= 1]

# ── Carga y mapeo de datos ────────────────────────────────────────────────────
freq_data <- fread(freq_file)
pop_info  <- fread(pop_file)
setnames(pop_info, "Group_Master_ID", "CLST", skip_absent = TRUE)

merged_df <- merge(freq_data, pop_info[, .(CLST, Media_Lat, Media_Long)], by = "CLST", all.x = TRUE)

if ("Media_Lat.y" %in% names(merged_df)) {
  merged_df[, c("Media_Lat.y", "Media_Long.y") := NULL]
  setnames(merged_df, old = c("Media_Lat.x", "Media_Long.x"), new = c("Media_Lat", "Media_Long"))
}

eurasia_df <- merged_df[Media_Lat >= 30 & Media_Lat <= 75 & Media_Long >= -15 & Media_Long <= 140]
eurasia_df[, Media_Long := as.numeric(Media_Long)]
eurasia_df[, Media_Lat  := as.numeric(Media_Lat)]
eurasia_df[, X := as.integer(cut(Media_Long, breaks = seq(-15, 140, length.out = nx + 1), include.lowest = TRUE))]
eurasia_df[, Y := as.integer(cut(Media_Lat,  breaks = seq(30,   75, length.out = ny + 1), include.lowest = TRUE))]
eurasia_df[, Generation := floor((11700 - Media_YBP) / 27) + 1]
eurasia_df <- eurasia_df[Generation >= 1]

grid_data <- eurasia_df[, .(
  Frequency = sum(MAC, na.rm = TRUE) / sum(NCHROBS, na.rm = TRUE),
  Count     = sum(MAC, na.rm = TRUE),
  Chr_Tot   = sum(NCHROBS, na.rm = TRUE)
), by = .(SNP, Generation, X, Y)]

setnames(grid_data, "SNP", "MutationID")
grid_data$ChrOBS <- grid_data$Chr_Tot
tiempos_globales_muestreo <- sort(unique(grid_data$Generation))
snps_to_analyze <- sort(unique(grid_data$MutationID))

# ── Loop principal por SNP ────────────────────────────────────────────────────
for (snp_actual in snps_to_analyze) {
  
  df_snp <- grid_data[MutationID == snp_actual]
  df_snp <- df_snp[!is.na(Frequency) & !is.nan(Frequency)]
  if (nrow(df_snp) == 0 || max(df_snp$Frequency) == 0) next
  
  df_snp <- df_snp[order(Generation)]
  AlleleOriginAge <- df_snp$Generation[which(df_snp$Frequency > 0)[1]]
  max_gen_snp     <- max(df_snp$Generation)
  
  if (max_gen_snp < max(tiempos_globales_muestreo)) {
    gens_faltantes <- tiempos_globales_muestreo[tiempos_globales_muestreo > max_gen_snp]
    if (length(gens_faltantes) > MAX_PADDED_GENERATIONS)
      gens_faltantes <- gens_faltantes[1:MAX_PADDED_GENERATIONS]
    filas_ceros <- as.data.table(expand.grid(
      Generation = gens_faltantes, X = unique(df_snp$X), Y = unique(df_snp$Y)))
    filas_ceros[, `:=`(MutationID = snp_actual, Frequency = 0, Count = 0,
                       Chr_Tot = max(round(mean(df_snp$Chr_Tot, na.rm = TRUE)), 10),
                       ChrOBS  = max(round(mean(df_snp$Chr_Tot, na.rm = TRUE)), 10))]
    df_snp <- rbind(df_snp, filas_ceros, use.names = TRUE)[order(Generation, Y, X)]
  }
  
  times_run <- min(df_snp$Generation):max(df_snp$Generation)
  if (length(times_run) < 2) next
  
  Conc0 <- matrix(0, nrow = nx, ncol = ny)
  origin_data <- df_snp[Generation == AlleleOriginAge & Frequency > 0]
  for (k in 1:nrow(origin_data))
    if (origin_data$X[k] >= 1 && origin_data$Y[k] >= 1)
      Conc0[origin_data$X[k], origin_data$Y[k]] <- origin_data$Frequency[k]
  
  counts_vec  <- df_snp$Count
  totals_vec  <- df_snp$Chr_Tot
  rho_ares_vec <- vapply(1:nrow(df_snp), function(j) {
    max(1 - exp(-max(df_snp$Generation[j] - AlleleOriginAge, 0.5) / (2 * N_eff)), 1e-6)
  }, numeric(1))
  
  cat("Procesando SNP:", snp_actual, "\n")
  
  res_leg <- list(); res_ares <- list(); res_spikes <- list()
  
  for (i in seq_len(nrow(expand.grid(D = DifussionValuesToCheck)))) {
    D_curr       <- DifussionValuesToCheck[i]
    T_onset_curr <- NA
    
    # ── Función objetivo interna ──────────────────────────────────────
    eval_s_optim <- function(s_val, model_type = "LEGACY") {
      h     <- 0.5   # modelo aditivo estándar
      d_val <- h * s_val
      ST3 <- tryCatch(
        ode.2D(y = Conc0, times = times_run, func = diffusion2D,
               parms = c(D_curr, d_val, s_val),
               dimens = c(nx, ny), method = "rk4", atol = 1e-6, rtol = 1e-6),
        error = function(e) NULL)
      if (is.null(ST3)) return(1e9)
      ST3_mat <- as.matrix(ST3[, -1])
      pred_freq_vec <- vapply(1:nrow(df_snp), function(j) {
        time_idx    <- match(df_snp$Generation[j], times_run)
        spatial_idx <- (df_snp$Y[j] - 1) * nx + df_snp$X[j]
        max(min(ST3_mat[time_idx, spatial_idx], 1 - 1e-6), 1e-6)
      }, numeric(1))
      ll <- switch(model_type,
        "LEGACY" = sum(dbinom(counts_vec, totals_vec, pred_freq_vec, log = TRUE)),
        "ARES"   = sum(dbetabinom(counts_vec, totals_vec, pred_freq_vec, rho_ares_vec, log = TRUE)),
        "SPIKES" = sum(dbetabinom_spikes(counts_vec, totals_vec, pred_freq_vec, rho_ares_vec, log = TRUE)))
      return(ifelse(is.finite(ll), -ll, 1e9))
    }
    
    # ── Rama neutral: PDE una sola vez ───────────────────────────────
    if (is_neutral) {
      ST3 <- tryCatch(
        ode.2D(y = Conc0, times = times_run, func = diffusion2D,
               parms = c(D_curr, 0, 0),
               dimens = c(nx, ny), method = "rk4", atol = 1e-6, rtol = 1e-6),
        error = function(e) NULL)
      if (is.null(ST3)) {
        res_leg[[i]]    <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=0, LL=-1e9, SE_s=0)
        res_ares[[i]]   <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=0, LL=-1e9, SE_s=0)
        res_spikes[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=0, LL=-1e9, SE_s=0)
        next
      }
      ST3_mat <- as.matrix(ST3[, -1])
      pred_freq_vec <- vapply(1:nrow(df_snp), function(j) {
        time_idx    <- match(df_snp$Generation[j], times_run)
        spatial_idx <- (df_snp$Y[j] - 1) * nx + df_snp$X[j]
        max(min(ST3_mat[time_idx, spatial_idx], 1 - 1e-6), 1e-6)
      }, numeric(1))
      ll_leg <- sum(dbinom(counts_vec, totals_vec, pred_freq_vec, log = TRUE))
      ll_ar  <- sum(dbetabinom(counts_vec, totals_vec, pred_freq_vec, rho_ares_vec, log = TRUE))
      ll_spk <- sum(dbetabinom_spikes(counts_vec, totals_vec, pred_freq_vec, rho_ares_vec, log = TRUE))
      res_leg[[i]]    <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=0, LL=ll_leg, SE_s=0)
      res_ares[[i]]   <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=0, LL=ll_ar,  SE_s=0)
      res_spikes[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=0, LL=ll_spk, SE_s=0)
      
    # ── Rama selección: optimizar s ──────────────────────────────────
    } else {
      opt_leg <- optim(0.01, eval_s_optim, model_type="LEGACY",
                       method="L-BFGS-B", lower=-1, upper=1, hessian=TRUE)
      res_leg[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr,
                                 s_opt=opt_leg$par, LL=-opt_leg$value,
                                 SE_s=ifelse(opt_leg$hessian[1,1]>0, sqrt(1/opt_leg$hessian[1,1]), NA))
      
      opt_ar <- optim(0.01, eval_s_optim, model_type="ARES",
                      method="L-BFGS-B", lower=-1, upper=1, hessian=TRUE)
      res_ares[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr,
                                  s_opt=opt_ar$par, LL=-opt_ar$value,
                                  SE_s=ifelse(opt_ar$hessian[1,1]>0, sqrt(1/opt_ar$hessian[1,1]), NA))
      
      opt_spk <- optim(0.01, eval_s_optim, model_type="SPIKES",
                       method="L-BFGS-B", lower=-1, upper=1, hessian=TRUE)
      res_spikes[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr,
                                    s_opt=opt_spk$par, LL=-opt_spk$value,
                                    SE_s=ifelse(opt_spk$hessian[1,1]>0, sqrt(1/opt_spk$hessian[1,1]), NA))
    }
  }
  
  # ── Escribir a disco y liberar memoria ───────────────────────────────────
  finalize_snp(res_leg,    "LEGACY", snp_actual)
  finalize_snp(res_ares,   "ARES",   snp_actual)
  finalize_snp(res_spikes, "SPIKES", snp_actual)
  
  rm(res_leg, res_ares, res_spikes, df_snp,
     counts_vec, totals_vec, rho_ares_vec, Conc0)
  gc(verbose = FALSE)
}

cat("Ejecución completada.\n")