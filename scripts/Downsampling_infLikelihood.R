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
  DifussionValuesToCheck <- current_mig <- as.numeric(current_mig)
}

# Esta parte esta en construccion
if (model_name == "D_FULL_seleccion_m1") {
  SelectionValuesToCheck <- current_sel <- as.numeric(current_sel)
}
if (length(SelectionValuesToCheck) == 1) {
  FIXED_SELECTION <- TRUE
} else {
  FIXED_SELECTION <- FALSE
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

# 2. GENERACIÓN DE ESCENARIOS DE DOWNSAMPLING ---------------------------------
generar_escenarios <- function(df_sim, df_adna) {
  escenarios <- list()
  escenarios[["Completo_100"]] <- df_sim
  
  # A. Downsampling Temporal (Elimina el X% de las generaciones observadas)
  generaciones_unicas <- unique(df_sim$Generation)
  for (frac in c(0.75, 0.50, 0.25)) {
    gens_muestra <- sample(generaciones_unicas, size = floor(length(generaciones_unicas) * frac))
    escenarios[[paste0("Temporal_", frac * 100)]] <- df_sim %>% filter(Generation %in% gens_muestra)
  }
  
  # B. Downsampling Global (Elimina el X% de las filas aleatoriamente, lagunas espaciotemporales)
  for (frac in c(0.75, 0.50, 0.25)) {
    escenarios[[paste0("Global_", frac * 100)]] <- df_sim %>% sample_frac(frac)
  }
  
  # C. Máscara Empírica Realista (Solo mantiene coordenadas espaciotemporales que existen en aDNA)
  # Usamos un semi_join para quedarnos con los datos de SLiM que coinciden en Generation, X, y Y con la matriz de aDNA
  escenarios[["Empirico_aDNA"]] <- df_sim %>%
    semi_join(df_adna, by = c("Generation", "X", "Y"))
  
  return(escenarios)
}

lista_datos_evaluacion <- generar_escenarios(sim_data_raw, adna_metadata)

# 3. FUNCIÓN DE INFERENCIA DE VEROSIMILITUD (Basado en tu script) -------------
# Envolvemos tu código en una función para iterar sobre los escenarios
correr_inferencia <- function(freq_data, param_grid, snps_to_analyze, n, tiempos_globales_muestreo, N_eff, MAX_PADDED_GENERATIONS) {
  
  all_results_legacy <- list()
  all_results_ares <- list()
  all_results_spikes <- list()
  
  for (snp_actual in snps_to_analyze) {
    
    df_snp <- freq_data[freq_data$MutationID == snp_actual, ]
    if(nrow(df_snp) == 0 || max(df_snp$Frequency) == 0) next
    
    First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
    if(is.na(First_OcurrenceData)) next
    
    AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]
    max_gen_snp <- max(df_snp$Generation)
    max_gen_global <- max(tiempos_globales_muestreo)
    
    # --- ZERO-PADDING CONDICIONAL ---
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
      
      # Asegurar compatibilidad de columnas antes de hacer rbind
      columnas_comunes <- intersect(colnames(df_snp), colnames(filas_ceros))
      df_snp <- rbind(df_snp[, columnas_comunes], filas_ceros[, columnas_comunes])
      df_snp <- df_snp[order(df_snp$Generation, df_snp$Y, df_snp$X), ]
    }
    
    # --- CONFIGURACIÓN DE LA PDE ---
    times_run <- min(df_snp$Generation):max(df_snp$Generation)
    
    if(length(times_run) < 2) {
      next # Saltamos sin imprimir para no saturar la consola en iteraciones masivas
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
    ll_spikes_vals <- numeric(nrow(param_grid))
    
    for (i in 1:nrow(param_grid)) {
      D_curr <- param_grid$D[i]; s_curr <- param_grid$s[i]; d_curr <- 2 * s_curr
      
      ST3 <- ode.2D(y = Conc0, times = times_run, func = diffusion2D,
                    parms = c(D_curr, d_curr, s_curr), dimens = c(n, n),
                    method = "rk4", atol = 1e-7, rtol = 1e-7)
      
      ST3_mat <- as.matrix(ST3[,-1])
      piso_minimo <- 1e-6 
      
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
      
      # Verosimilitudes
      ll_legacy_vals[i] <- sum(dbinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, log = TRUE))
      ll_ares_vals[i] <- sum(dbetabinom(x = counts_vec, size = totals_vec, prob = pred_freq_vec, rho = rho_ares_vec, log = TRUE))
      # ll_spikes_vals[i] <- sum(dbetabinom_spikes(x= counts_vec, size = totals_vec, prob = pred_freq_vec, rho= rho_ares_vec, log= TRUE)) # Comentado si dbetabinom_spikes es una función custom no estándar de VGAM
    }
    
    col_names <- paste0("D_", param_grid$D, "_s_", param_grid$s)
    all_results_legacy[[length(all_results_legacy) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_legacy_vals)), c("SNP", col_names))
    all_results_ares[[length(all_results_ares) + 1]] <- setNames(cbind(data.frame(SNP = snp_actual), t(ll_ares_vals)), c("SNP", col_names))
  }
  
  # Combinar y retornar
  return(list(
    Legacy = bind_rows(all_results_legacy),
    Ares = bind_rows(all_results_ares)
  ))
}

# 4. EJECUCIÓN DEL ANÁLISIS SOBRE TODOS LOS ESCENARIOS ------------------------
# NOTA: Asegúrate de tener definidos param_grid, n, tiempos_globales_muestreo, N_eff, MAX_PADDED_GENERATIONS y diffusion2D en tu entorno global.

resultados_completos <- list()

for (nombre_escenario in names(lista_datos_evaluacion)) {
  cat("\n======================================================\n")
  cat("Evaluando Escenario de Muestreo:", nombre_escenario, "\n")
  cat("Observaciones disponibles:", nrow(lista_datos_evaluacion[[nombre_escenario]]), "\n")
  cat("======================================================\n")
  
  # Extraemos los datos del escenario actual
  datos_escenario <- lista_datos_evaluacion[[nombre_escenario]]
  
  # Corremos la PDE
  resultados_escenario <- correr_inferencia(
    freq_data = datos_escenario,
    param_grid = param_grid, 
    snps_to_analyze = snps_to_analyze, 
    n = n, 
    tiempos_globales_muestreo = tiempos_globales_muestreo, 
    N_eff = N_eff, 
    MAX_PADDED_GENERATIONS = MAX_PADDED_GENERATIONS
  )
  
  # Guardamos los perfiles de verosimilitud (Log-Likelihood surfaces) en RAM
  resultados_completos[[nombre_escenario]] <- resultados_escenario
}

# 5. POST-PROCESAMIENTO: Comparación de Máxima Verosimilitud ------------------
# Aquí puedes extraer el MLE (Maximum Likelihood Estimate) para cada escenario
# para graficar cómo se desplazan las estimaciones de D y s a medida que pierdes datos.