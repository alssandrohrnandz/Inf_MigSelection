library(deSolve)
library(dplyr)
library(tidyverse)
library(VGAM)

args <- commandArgs(trailingOnly = TRUE)
freq_file   <- args[1] 
subset_file <- args[2] 
task_id     <- args[3]
model_name  <- args[4]
output_dir  <- args[5]
simulate_aDNA <- as.logical(args[6]) 

if(is.na(simulate_aDNA)) simulate_aDNA <- FALSE

print(paste("Procesando archivo:", freq_file))
print(paste("Modelo:", model_name, "| aDNA Mode:", simulate_aDNA))

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
  pi_mid <- 1 - pi_0 - pi_1
  pi_mid <- max(pi_mid, 1e-10)
  
  ll_bb <- dbetabinom(x = x, size = size, prob = prob, rho = rho, log = FALSE)
  dens <- pi_mid * ll_bb
  dens[x == 0]    <- dens[x == 0]    + pi_0
  dens[x == size] <- dens[x == size] + pi_1
  dens <- pmax(dens, 1e-300)
  if (log) return(log(dens)) else return(dens)
}

GRID_SIZE <- 10
dy <- dx <- 1
n <- GRID_SIZE
N_eff <- 1000

exponentes_D <- -3:0
valores_base <- c(1, 5)
DifussionValuesToCheck <- sort(unique(c(0, as.vector(outer(valores_base, 10^exponentes_D)))))
DifussionValuesToCheck <- DifussionValuesToCheck[DifussionValuesToCheck <= 1]

if (grepl("neutros", model_name)) {
  SelectionValuesToCheck <- c(0.0)
} else {
  bases <- c(1, 2, 2.5, 5)
  exponentes <- -5:0
  valores_pos <- as.vector(outer(bases, 10^exponentes))
  valores_pos <- valores_pos[valores_pos <= 1]
  SelectionValuesToCheck <- sort(unique(c(-valores_pos, 0, valores_pos)))
}

# PRE-DEFINIMOS LA REJILLA Y NOMBRES DE COLUMNA PARA GARANTIZAR LA ESTRUCTURA DEL OUTPUT
param_grid <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)
col_names <- paste0("D_", param_grid$D, "_s_", param_grid$s)

# ==========================================
# LECTURA DE DATOS
# ==========================================
freq_data_raw <- read.csv(freq_file, header=TRUE, stringsAsFactors=FALSE)
if (ncol(freq_data_raw) == 8) {
  colnames(freq_data_raw) <- c("Generation","MutationID","TypeMut","X","Y","Frequency","Count","Chr_Tot")
} else {
  colnames(freq_data_raw) <- c("Generation","MutationID","X","Y","Frequency","Count","Chr_Tot")
  freq_data_raw$TypeMut <- NA 
}

tiempos_globales_muestreo <- sort(unique(freq_data_raw$Generation))
gen_max_global <- max(tiempos_globales_muestreo) # Típicamente 500

# ==========================================
# NUEVO: FILTRO ESTRICTO DE SUPERVIVENCIA Y DURACIÓN
# ==========================================
alelos_fuertes <- freq_data_raw %>%
  group_by(MutationID) %>%
  summarise(
    Gen_Max = max(Generation),
    Gen_Min = min(Generation),
    Duracion = Gen_Max - Gen_Min,
    .groups = "drop"
  ) %>%
  # Condición: Llega al final de la simulación Y dura al menos 10 generaciones
  filter(Gen_Max == gen_max_global & Duracion >= 10) %>%
  pull(MutationID)

snps_subset <- readLines(subset_file)
# Intersectamos el subset deseado con los alelos que pasaron el filtro estricto
snps_to_analyze <- intersect(intersect(unique(freq_data_raw$MutationID), as.integer(snps_subset)), alelos_fuertes)

cat(sprintf("\nAlelos que cumplen criterios de supervivencia (>=10 gens y llegan a %d): %d\n", gen_max_global, length(snps_to_analyze)))

# ==========================================
# CREACIÓN DE LA MÁSCARA FÓSIL (YACIMIENTOS SESGADOS)
# ==========================================
if (simulate_aDNA) {
  mask_dir <- file.path(dirname(output_dir), "aDNA_Masks")
  if (!dir.exists(mask_dir)) dir.create(mask_dir, recursive = TRUE)
  mask_file <- file.path(mask_dir, paste0("aDNA_Global_Mask_", task_id, ".csv"))
  
  if (!file.exists(mask_file)) {
    set.seed(42 + as.integer(task_id)) 
    gen_min <- min(tiempos_globales_muestreo)
    rango_tiempo <- gen_max_global - gen_min
    
    # 1.1 Definimos el "esfuerzo de muestreo"
    total_yacimientos <- 50 
    inflexion <- gen_min + (rango_tiempo * 0.6) 
    
    # a) Forzar un muestreo base en el periodo antiguo (ej. 15% de los sitios asegurados)
    n_tempranos <- max(3, floor(total_yacimientos * 0.15)) # Garantiza al menos 3 sitios, en este caso serán ~7
    gens_tempranas <- tiempos_globales_muestreo[tiempos_globales_muestreo < inflexion]
    yacimientos_tempranos <- sample(gens_tempranas, size = n_tempranos, replace = TRUE)
    
    # b) El resto sigue la curva del "boom" Neolítico (Suavizamos la pendiente a 0.03)
    n_tardios <- total_yacimientos - n_tempranos
    probabilidades_tiempo <- 1 / (1 + exp(-0.03 * (tiempos_globales_muestreo - inflexion)))
    probabilidades_tiempo <- probabilidades_tiempo / sum(probabilidades_tiempo)
    
    yacimientos_tardios <- sample(tiempos_globales_muestreo, size = n_tardios, replace = TRUE, prob = probabilidades_tiempo)
    
    # c) Unimos ambos tiempos
    tiempos_yacimientos <- c(yacimientos_tempranos, yacimientos_tardios)
    
    # 1.4 Construimos el Data Frame de la Máscara
    patron_muestreo <- data.frame(
      Generation = tiempos_yacimientos,
      X = sample(1:GRID_SIZE, total_yacimientos, replace = TRUE),
      Y = sample(1:GRID_SIZE, total_yacimientos, replace = TRUE),
      # Nuevo tamaño de muestra
      Chr_Tot_Fijo = sample(2:20, total_yacimientos, replace = TRUE) 
    ) %>%
      group_by(Generation, X, Y) %>%
      summarise(Chr_Tot_Fijo = sum(Chr_Tot_Fijo), .groups = 'drop')
    
    write.csv(patron_muestreo, mask_file, row.names = FALSE)
  } else {
    patron_muestreo <- read.csv(mask_file)
  }
}

all_results_classic <- list()
all_results_bb <- list()
all_results_bws <- list()

# === BUCLE PRINCIPAL ===
for (snp_actual in snps_to_analyze) {
  
  df_snp_real <- freq_data_raw[freq_data_raw$MutationID == snp_actual, ]
  
  if (simulate_aDNA) {
    df_snp <- patron_muestreo %>%
      left_join(df_snp_real, by = c("Generation", "X", "Y")) %>%
      mutate(
        MutationID = snp_actual,
        Real_Frequency = replace_na(Frequency, 0), 
        Chr_Tot = Chr_Tot_Fijo
      )
    df_snp$Count <- rbinom(nrow(df_snp), size = df_snp$Chr_Tot, prob = df_snp$Real_Frequency)
    df_snp$Frequency <- df_snp$Count / df_snp$Chr_Tot
    df_snp <- df_snp %>% select(Generation, MutationID, X, Y, Frequency, Real_Frequency, Count, Chr_Tot)
    
    # Filtro post-máscara: ¿Sobrevivió suficiente información fósil?
    if(sum(df_snp$Count) == 0 || length(unique(df_snp$Generation[df_snp$Frequency > 0])) < 2) next
    
  } else {
    df_snp <- df_snp_real
    df_snp$Real_Frequency <- df_snp$Frequency 
  }

  First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
  if(is.na(First_OcurrenceData)) next
  AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]
  
  times_run <- min(df_snp$Generation):max(df_snp$Generation)
  if(length(times_run) < 2) next
  
  Conc0 <- matrix(0, nrow=n, ncol=n)
  origin_data <- df_snp[df_snp$Generation == AlleleOriginAge & df_snp$Frequency > 0, ]
  for(k in 1:nrow(origin_data)) {
    ox <- origin_data$X[k]; oy <- origin_data$Y[k]
    if(ox >= 1 && ox <= n && oy >= 1 && oy <= n) Conc0[ox, oy] <- origin_data$Frequency[k]
  }
  if(sum(Conc0) == 0) next
  
  ll_classic_vals <- numeric(nrow(param_grid))
  ll_bb_vals <- numeric(nrow(param_grid))
  ll_bws_vals <- numeric(nrow(param_grid))
  
  for (i in 1:nrow(param_grid)) {
    D_curr <- param_grid$D[i]; s_curr <- param_grid$s[i]; d_curr <- 2 * s_curr
    
    ST3 <- ode.2D(y = Conc0, times = times_run, func = diffusion2D,
                  parms = c(D_curr, d_curr, s_curr), dimens = c(n, n),
                  method = "rk4", atol = 1e-7, rtol = 1e-7)
    ST3_mat <- as.matrix(ST3[,-1])
    
    counts_vec <- df_snp$Count
    totals_vec <- df_snp$Chr_Tot
    pred_freq_vec <- numeric(nrow(df_snp))
    rho_ares_vec <- numeric(nrow(df_snp))
    
    for(j in 1:nrow(df_snp)) {
      t_abs <- df_snp$Generation[j]
      time_idx <- match(t_abs, times_run)
      spatial_idx <- (df_snp$Y[j] - 1) * n + df_snp$X[j] 
      
      pred_freq <- max(min(ST3_mat[time_idx, spatial_idx], 1 - 1e-6), 1e-6)
      t_elapsed <- max(t_abs - AlleleOriginAge, 0.5)
      rho_val <- max(1 - exp(-t_elapsed / (2 * N_eff)), 1e-6)
      
      pred_freq_vec[j] <- pred_freq
      rho_ares_vec[j] <- rho_val
    }
    
    ll_classic_vals[i] <- sum(dbinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, log = TRUE))
    ll_bb_vals[i]      <- sum(dbetabinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, rho = rho_ares_vec, log = TRUE))
    ll_bws_vals[i]     <- sum(dbetabinom_spikes(x= counts_vec, size = totals_vec, prob = pred_freq_vec, rho= rho_ares_vec, log= TRUE))
  }
  
  all_results_classic[[length(all_results_classic) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_classic_vals)), c("SNP", col_names))
  all_results_bb[[length(all_results_bb) + 1]]           <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_bb_vals)), c("SNP", col_names))
  all_results_bws[[length(all_results_bws) + 1]]         <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_bws_vals)), c("SNP", col_names))
}

# ==========================================
# GUARDADO INTELIGENTE (MANEJO DE FALTA DE DATOS)
# ==========================================
guardar_resultados <- function(lista_res, prefix) {
  archivo_salida <- file.path(output_dir, paste0(prefix, "_Grid_TaskID_", task_id, "_", model_name, ".txt"))
  
  if(length(lista_res) > 0) {
    df_out <- do.call(rbind, lista_res)
  } else {
    # Si no hubo datos (extinción prematura o máscara aDNA destruyó todo), creamos una fila dummy
    # con el SNP "INSUFFICIENT_DATA" y puros NAs, manteniendo los headers correctos
    df_out <- data.frame(SNP = "INSUFFICIENT_DATA")
    for (cn in col_names) df_out[[cn]] <- NA
    cat("    -> No sobrevivieron SNPs para", prefix, ". Generando archivo dummy con NAs.\n")
  }
  
  write.table(df_out, file = archivo_salida, row.names = FALSE, quote = FALSE, sep = "\t")
}

guardar_resultados(all_results_classic, "Classic")
guardar_resultados(all_results_bb, "B-B")
guardar_resultados(all_results_bws, "BwS")

print("Ejecución de Likelihood completada exitosamente.")