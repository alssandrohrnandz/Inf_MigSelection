library(deSolve)
library(data.table)
library(VGAM)

# =========================================================================
# 1. PARSEO DE ARGUMENTOS SLURM
# =========================================================================
args <- commandArgs(trailingOnly = TRUE)

task_id       <- if(length(args) >= 1) args[1] else "1"
model_name    <- if(length(args) >= 2) args[2] else "D_FULL_neutros"
slim_file     <- if(length(args) >= 3) args[3] else paste0("/mnt/data/dortega/hlopezh/Inf_MigSelection/data/results_Discrete/outputs_slim/v2_D_FULL_neutros_m1_", task_id, ".csv")
empirical_csv <- if(length(args) >= 4) args[4] else "Eurasia_Populations_Periods_Grid.csv"
output_dir    <- if(length(args) >= 5) args[5] else getwd()
prop_data     <- if(length(args) >= 6) as.numeric(args[6]) else 0.60 

is_neutral <- grepl("neutros", model_name, ignore.case = TRUE)
nx <- 31; ny <- 9; dx <- 1; dy <- 1; N_eff <- 1000

# =========================================================================
# 2. DEFINICIÓN DE LA PDE (A prueba de NA en el caso neutro)
# =========================================================================
diffusion2D <- function(t, conc, par) {
  Conc <- matrix(nrow = nx, ncol = ny, data = conc)
  
  D_eff <- par[1]
  
  # Si T_onset es NA (modelo neutro) o aún no llegamos a T_onset, s=0
  if (!is.na(par[4]) && t >= par[4]) {
    d_eff <- par[2]
    s_eff <- par[3]
  } else {
    d_eff <- 0
    s_eff <- 0
  }
  
  dConc <- Conc * (1 - Conc) * (Conc * d_eff + s_eff * (1 - 2 * Conc))
  
  FluxX <- -D_eff * rbind(rep(0, ny), (Conc[2:nx,] - Conc[1:(nx-1),]), rep(0, ny)) / dx
  dConc <- dConc - (FluxX[2:(nx+1),] - FluxX[1:nx,]) / dx
  
  FluxY <- -D_eff * cbind(rep(0, nx), (Conc[,2:ny] - Conc[,1:(ny-1)]), rep(0, nx)) / dy
  dConc <- dConc - (FluxY[,2:(ny+1)] - FluxY[,1:ny]) / dy
  
  return(list(as.vector(dConc)))
}

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

# =========================================================================
# 3. CARGA DE DATOS Y PERFILADO EMPÍRICO
# =========================================================================
emp_grid <- fread(empirical_csv)
emp_summary <- emp_grid[, .(ChrOBS = sum(N_Samples) * 2), by = .(Generation, X = X, Y = Y)]

gen_emp_totales <- length(unique(emp_summary$Generation))
loc_promedio_emp <- round(nrow(emp_summary) / gen_emp_totales)
chr_promedio_emp <- round(mean(emp_summary$ChrOBS, na.rm = TRUE))

slim_data <- fread(slim_file)
slim_data[, X := as.integer(X)]
slim_data[, Y := as.integer(Y)]
slim_data[, Generation := as.integer(Generation)]

# Filtro de Viabilidad SLiM
slim_summary <- slim_data[, .(Max_Freq = max(Frequency), Lifespan = max(Generation) - min(Generation)), by = MutationID]
snps_viables <- slim_summary[Max_Freq > 0.01 & Lifespan >= 50]$MutationID

master_leg <- list(); master_ares <- list(); master_spikes <- list()

# =========================================================================
# 4. FUNCIÓN CENTRAL DE VEROSIMILITUD (Limpia y reutilizable)
# =========================================================================
# Esta función simplemente corre la PDE y devuelve la verosimilitud (sin optimizar nada)
compute_LL <- function(D_val, s_val, Tonset_val, model_type, Conc0, times_run, counts_vec, totals_vec, rho_ares_vec) {
  
  if (any(is.na(counts_vec)) || any(is.na(totals_vec))) {
    cat("   [ERROR] Los vectores de conteo tienen NAs.\n")
    return(-1e9)
  }

  h <- 2  
  d_val <- h * s_val
  
  ST3 <- tryCatch({
    ode.2D(y = Conc0, times = times_run, func = diffusion2D, 
           parms = c(D_val, d_val, s_val, Tonset_val), 
           dimens = c(nx, ny), method = "lsodes", 
           atol = 1e-6, rtol = 1e-6, 
           lrw = 5e5) # Asignación expandida de memoria interna
  }, error = function(e) {
    cat("   [ERROR ODE] Falla del solver:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(ST3)) return(-1e9)
  if (nrow(ST3) != length(times_run)) return(-1e9)
  
  ST3_mat <- as.matrix(ST3[, -1])
  
  if (any(is.nan(ST3_mat))) {
    return(-1e9)
  }
  
  pred_freq_vec <- vapply(1:length(counts_vec), function(j) {
    time_idx <- match(df_snp$Generation[j], times_run)
    spatial_idx <- (df_snp$Y[j] - 1) * nx + df_snp$X[j]
    max(min(ST3_mat[time_idx, spatial_idx], 1 - 1e-6), 1e-6)
  }, numeric(1))
  
  ll <- switch(model_type,
               "LEGACY" = sum(dbinom(counts_vec, totals_vec, pred_freq_vec, log = TRUE)),
               "ARES"   = sum(dbetabinom(counts_vec, totals_vec, pred_freq_vec, rho_ares_vec, log = TRUE)),
               "SPIKES" = sum(dbetabinom_spikes(counts_vec, totals_vec, pred_freq_vec, rho_ares_vec, log = TRUE)))
  
  if (!is.finite(ll)) return(-1e9)
  
  return(ll)
}

# =========================================================================
# 5. BUCLE DE SUBMUESTREO Y OPTIMIZACIÓN DUAL
# =========================================================================
for (snp_actual in snps_viables) {
  
  dt_sub <- slim_data[MutationID == snp_actual]

    if (sum(dt_sub$Frequency) < 0.001) {
    cat("   [!] Alelo demasiado raro para inferencia, saltando...\n")
    next
}

  Real_T_onset <- if(!is_neutral) min(dt_sub[Frequency > 0]$Generation) else NA
  todas_las_gen_slim <- unique(dt_sub[Frequency > 0]$Generation)
  
  set.seed(as.integer(task_id) + as.integer(snp_actual))
  
  # A & B. Submuestreo Estocástico
  n_gen_a_tomar <- max(round(length(todas_las_gen_slim) * prop_data), 5)
  gen_seleccionadas <- sort(sample(todas_las_gen_slim, size = n_gen_a_tomar))
  dt_sub <- dt_sub[Generation %in% gen_seleccionadas]
  
  # B. Submuestreo Espacial Estocástico (Priorizando presencias)
  n_loc_a_tomar <- max(round(loc_promedio_emp * prop_data), 1)
  
  dt_espacial_list <- list()
  for (g in gen_seleccionadas) {
    locs_disponibles <- dt_sub[Generation == g]
    
    if (nrow(locs_disponibles) > 0) {
      # ¿Cuántos lugares necesitamos en total para esta generación?
      target_n <- min(n_loc_a_tomar, nrow(locs_disponibles))
      
      # Separar localidades en "positivas" (tienen el alelo) y "negativas" (aún no llega)
      locs_con_alelo <- locs_disponibles[Frequency > 0]
      locs_sin_alelo <- locs_disponibles[Frequency == 0]
      
      n_con_alelo <- nrow(locs_con_alelo)
      
      if (n_con_alelo >= target_n) {
        # CASO 1: Hay suficientes lugares con el alelo. Tomamos el target solo de estos.
        locs_finales <- locs_con_alelo[sample(.N, target_n)]
      } else {
        # CASO 2: Hay pocos lugares con el alelo. Los tomamos TODOS obligatoriamente.
        # Luego, rellenamos el resto del target con lugares azarosos sin el alelo (para imitar que no ha llegado).
        faltantes <- target_n - n_con_alelo
        
        if (faltantes > 0 && nrow(locs_sin_alelo) > 0) {
          faltantes <- min(faltantes, nrow(locs_sin_alelo)) # Seguro contra pedir más de los que hay
          locs_extra <- locs_sin_alelo[sample(.N, faltantes)]
          locs_finales <- rbindlist(list(locs_con_alelo, locs_extra))
        } else {
          locs_finales <- locs_con_alelo
        }
      }
      
      dt_espacial_list[[as.character(g)]] <- locs_finales
    }
  }
  df_snp <- rbindlist(dt_espacial_list)
  
  if(sum(df_snp$Frequency) == 0) {
  cat("   [!] Frecuencia cero en toda la muestra. Saltando alelo...\n")
  next
  }
  # C. Downsampling Binomial
  df_snp[, ChrOBS := rpois(.N, chr_promedio_emp)]
  df_snp[ChrOBS < 2, ChrOBS := 2]
  df_snp[, Count := rbinom(.N, size = ChrOBS, prob = Frequency)]
  df_snp[, Freq_Obs := Count / ChrOBS]
  
  # Preparar tensores
  df_snp <- df_snp[order(Generation)]
  generaciones_positivas <- df_snp[Freq_Obs > 0]$Generation
  if (length(unique(generaciones_positivas)) < 5) next
  
  AlleleOriginAge <- min(generaciones_positivas)
  df_snp <- df_snp[Generation >= AlleleOriginAge]
  times_run <- min(df_snp$Generation):max(df_snp$Generation)
  
  Conc0 <- matrix(0, nrow = nx, ncol = ny)
  origin_data <- df_snp[Generation == AlleleOriginAge & Freq_Obs > 0]
  for(k in 1:nrow(origin_data)) {
    if(origin_data$X[k] >= 1 && origin_data$X[k] <= nx && origin_data$Y[k] >= 1 && origin_data$Y[k] <= ny) {
      Conc0[origin_data$X[k], origin_data$Y[k]] <- origin_data$Freq_Obs[k]
    }
  }
  
  # Valores D a probar
  exponentes_D <- -3:0; valores_base <- c(1, 5)
  DifussionValuesToCheck <- sort(unique(c(0.0001, as.vector(outer(valores_base, 10^exponentes_D)))))
  DifussionValuesToCheck <- DifussionValuesToCheck[DifussionValuesToCheck <= 1]
  
  counts_vec <- df_snp$Count
  totals_vec <- df_snp$ChrOBS
  rho_ares_vec <- vapply(1:nrow(df_snp), function(j) max(1 - exp(-max(df_snp$Generation[j] - AlleleOriginAge, 0.5) / (2 * N_eff)), 1e-6), numeric(1))
  
  snp_leg <- list(); snp_ares <- list(); snp_spikes <- list()
  
  # -------------------------------------------------------------------------
  # RAMA 1: FLUJO NEUTRO (Solo iterar D, sin buscar s ni T_onset)
  # -------------------------------------------------------------------------
  if (is_neutral) {
    for (i in 1:length(DifussionValuesToCheck)) {
      D_curr <- DifussionValuesToCheck[i]
      
      ll_leg <- compute_LL(D_curr, 0, NA, "LEGACY", Conc0, times_run, counts_vec, totals_vec, rho_ares_vec)
      ll_ar  <- compute_LL(D_curr, 0, NA, "ARES", Conc0, times_run, counts_vec, totals_vec, rho_ares_vec)
      ll_spk <- compute_LL(D_curr, 0, NA, "SPIKES", Conc0, times_run, counts_vec, totals_vec, rho_ares_vec)
      
      snp_leg[[i]]    <- data.table(D=D_curr, Tonset=NA, s_opt=0, LL=ll_leg, SE_s=NA)
      snp_ares[[i]]   <- data.table(D=D_curr, Tonset=NA, s_opt=0, LL=ll_ar,  SE_s=NA)
      snp_spikes[[i]] <- data.table(D=D_curr, Tonset=NA, s_opt=0, LL=ll_spk, SE_s=NA)
    }
  } 
  # -------------------------------------------------------------------------
  # RAMA 2: FLUJO DE SELECCIÓN (Iterar D y T_onset, optimizando s)
  # -------------------------------------------------------------------------
  else {
    T_onset_Values <- unique(round(seq(Real_T_onset, AlleleOriginAge, length.out = 4)))
    # o asegurar al menos 2 puntos:
    if(length(T_onset_Values) < 2) T_onset_Values <- c(Real_T_onset, AlleleOriginAge)
    local_grid <- expand.grid(D = DifussionValuesToCheck, T_onset = T_onset_Values)
    
    for (i in 1:nrow(local_grid)) {
      D_curr <- local_grid$D[i]
      T_onset_curr <- local_grid$T_onset[i]
      
      opt_leg <- optim(0.01, function(s) -compute_LL(D_curr, s, T_onset_curr, "LEGACY", Conc0, times_run, counts_vec, totals_vec, rho_ares_vec), method="L-BFGS-B", lower=-1, upper=1, hessian=TRUE)
      snp_leg[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=opt_leg$par, LL=-opt_leg$value, SE_s=ifelse(opt_leg$hessian[1,1]>0, sqrt(1/opt_leg$hessian[1,1]), NA))
      
      opt_ar <- optim(0.01, function(s) -compute_LL(D_curr, s, T_onset_curr, "ARES", Conc0, times_run, counts_vec, totals_vec, rho_ares_vec), method="L-BFGS-B", lower=-1, upper=1, hessian=TRUE)
      snp_ares[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=opt_ar$par, LL=-opt_ar$value, SE_s=ifelse(opt_ar$hessian[1,1]>0, sqrt(1/opt_ar$hessian[1,1]), NA))
      
      opt_spk <- optim(0.01, function(s) -compute_LL(D_curr, s, T_onset_curr, "SPIKES", Conc0, times_run, counts_vec, totals_vec, rho_ares_vec), method="L-BFGS-B", lower=-1, upper=1, hessian=TRUE)
      snp_spikes[[i]] <- data.table(D=D_curr, Tonset=T_onset_curr, s_opt=opt_spk$par, LL=-opt_spk$value, SE_s=ifelse(opt_spk$hessian[1,1]>0, sqrt(1/opt_spk$hessian[1,1]), NA))
    }
  }
  
  bind_and_append <- function(snp_list, master) {
    dt <- rbindlist(snp_list)
    dt[, `:=`(SNP = snp_actual, Prop_Data = prop_data, True_Tonset = Real_T_onset)]
    master[[length(master) + 1]] <- dt
    return(master)
  }
  
  master_leg <- bind_and_append(snp_leg, master_leg)
  master_ares <- bind_and_append(snp_ares, master_ares)
  master_spikes <- bind_and_append(snp_spikes, master_spikes)
}

# =========================================================================
# 6. ESCRITURA FINAL
# =========================================================================
finalize_sim <- function(master_list, model_type) {
  if (length(master_list) == 0) return()
  dt <- rbindlist(master_list)
  dt[, CI_95_Lower := s_opt - 1.96 * SE_s]
  dt[, CI_95_Upper := s_opt + 1.96 * SE_s]
  
  out_name <- file.path(output_dir, paste0("SIM_OPTIM_", model_type, "_Task_", task_id, "_Prop_", prop_data, "_", model_name, ".txt"))
  fwrite(dt, file = out_name, sep = "\t", quote = FALSE)
}

finalize_sim(master_leg, "LEGACY")
finalize_sim(master_ares, "ARES")
finalize_sim(master_spikes, "SPIKES")