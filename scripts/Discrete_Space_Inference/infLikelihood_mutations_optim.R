library(deSolve)
library(data.table)
library(VGAM)
library(parallel)

# =========================================================================
# 1. CONFIGURACIÓN DEL CLÚSTER HPCC
# =========================================================================
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
if (num_cores < 1) num_cores <- 1
cat("Iniciando L-BFGS-B en paralelo con", num_cores, "núcleos.\n")

args <- commandArgs(trailingOnly = TRUE)
freq_file   <- if(length(args) >= 1) args[1] else "data/results_Discrete/outputs_slim/D_FULL_neutros_m1_1.csv"
task_id     <- if(length(args) >= 2) args[2] else "1"
model_name  <- if(length(args) >= 3) args[3] else "D_FULL_neutros_m1"
output_dir  <- if(length(args) >= 4) args[4] else getwd()

# Variables espaciales globales
GRID_SIZE <- 10; n <- GRID_SIZE; dx <- 1; dy <- 1; N_eff <- 1000
MAX_PADDED_GENERATIONS <- ifelse(model_name == "D_FULL_neutros_m1", 1, 2)

# Funciones Matemáticas
diffusion2D <- function(t, conc, par) {
  Conc <- matrix(nrow = n, ncol = n, data = conc)
  dConc <- Conc*(1-Conc)*(Conc*par[2]+par[3]*(1-2*Conc))
  Flux <- -par[1] * rbind(rep(0, n), (Conc[2:n,]-Conc[1:(n-1),]), rep(0, n))/dx
  dConc <- dConc - (Flux[2:(n+1),] - Flux[1:n,])/dx
  Flux <- -par[1] * cbind(rep(0, n), (Conc[,2:n]-Conc[,1:(n-1)]), rep(0, n))/dy
  dConc <- dConc - (Flux[,2:(n+1)]-Flux[,1:n])/dy
  return(list(as.vector(dConc)))
}


# =========================================================================
# 2. CARGA Y EMPAQUETADO DE DATOS EMPÍRICOS
# =========================================================================
freq_data_raw <- read.csv(freq_file, header=TRUE, stringsAsFactors=FALSE)
if (ncol(freq_data_raw) == 8) {
  colnames(freq_data_raw) <- c("Generation","MutationID","TypeMut","X","Y","Frequency","Count","Chr_Tot")
} else {
  colnames(freq_data_raw) <- c("Generation","MutationID","X","Y","Frequency","Count","Chr_Tot")
  freq_data_raw$TypeMut <- NA 
}

freq_data <- freq_data_raw[freq_data_raw$Generation >= 600, ]
freq_data$ChrOBS <- freq_data$Chr_Tot
tiempos_globales_muestreo <- sort(unique(freq_data$Generation))
snps_to_analyze <- unique(sort(freq_data$MutationID))

# Pre-procesar SNPs para evitar hacer esto en cada iteración del optimizador
lista_snps_empaquetados <- list()

for (snp in snps_to_analyze) {
  df_snp <- freq_data[freq_data$MutationID == snp, ]
  if(nrow(df_snp) == 0 || max(df_snp$Frequency) == 0) next
  First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
  if(is.na(First_OcurrenceData)) next
  
  AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]
  max_gen_snp <- max(df_snp$Generation)
  max_gen_global <- max(tiempos_globales_muestreo)
  
  # Zero-padding
  if (max_gen_snp < max_gen_global) {
    gens_faltantes <- tiempos_globales_muestreo[tiempos_globales_muestreo > max_gen_snp]
    if(length(gens_faltantes) > MAX_PADDED_GENERATIONS) gens_faltantes <- gens_faltantes[1:MAX_PADDED_GENERATIONS]
    coordenadas_historicas <- unique(df_snp[, c("X", "Y")])
    chr_tot_promedio <- max(round(mean(df_snp$Chr_Tot, na.rm = TRUE)), 10)
    filas_ceros <- expand.grid(Generation = gens_faltantes, X = coordenadas_historicas$X, Y = coordenadas_historicas$Y)
    filas_ceros$MutationID <- snp; filas_ceros$Frequency <- 0; filas_ceros$Count <- 0
    filas_ceros$Chr_Tot <- chr_tot_promedio; filas_ceros$TypeMut <- NA; filas_ceros$ChrOBS <- chr_tot_promedio      
    filas_ceros <- filas_ceros[, colnames(df_snp)]
    df_snp <- rbind(df_snp, filas_ceros)
    df_snp <- df_snp[order(df_snp$Generation, df_snp$Y, df_snp$X), ]
  }
  
  times_run <- min(df_snp$Generation):max(df_snp$Generation)
  if(length(times_run) < 2) next
  
  Conc0 <- matrix(0, nrow=n, ncol=n)
  origin_data <- df_snp[df_snp$Generation == AlleleOriginAge & df_snp$Frequency > 0, ]
  for(k in 1:nrow(origin_data)) {
    ox <- origin_data$X[k]; oy <- origin_data$Y[k]
    # Corrección de coordenadas: Ajustamos límite superior para evitar out of bounds
    if(ox >= 1 && ox <= n && oy >= 1 && oy <= n) Conc0[ox, oy] <- origin_data$Frequency[k]
  }
  if(sum(Conc0) == 0) next
  
  # Calcular rho y guardar datos pre-procesados
  rho_ares_vec <- numeric(nrow(df_snp))
  for(j in 1:nrow(df_snp)) {
    t_abs <- df_snp$Generation[j]
    t_elapsed <- max(t_abs - AlleleOriginAge, 0.5)
    rho_ares_vec[j] <- max(1 - exp(-t_elapsed / (2 * N_eff)), 1e-6)
  }
  
  lista_snps_empaquetados[[length(lista_snps_empaquetados) + 1]] <- list(
    id = snp, data = df_snp, Conc0 = Conc0, times_run = times_run, 
    AlleleOriginAge = AlleleOriginAge, rho_ares_vec = rho_ares_vec
  )
}

cat("Total de SNPs empaquetados viables:", length(lista_snps_empaquetados), "\n")

# =========================================================================
# 3. FUNCIÓN OBJETIVO PARA L-BFGS-B (Evaluación Conjunta)
# =========================================================================
funcion_objetivo_conjunta <- function(parametros) {
  D_propuesto <- parametros[1]
  s_propuesto <- parametros[2]
  d_propuesto <- 2*s_propuesto # Relación fija entre d y s para reducir la dimensionalidad del espacio de búsqueda 
  
  # Evaluación paralela de todos los SNPs para el D propuesto
  ll_lista <- mclapply(lista_snps_empaquetados, function(paquete) {
    
    ST3 <- tryCatch({
      ode.2D(y = paquete$Conc0, times = paquete$times_run, func = diffusion2D,
            parms = c(D_propuesto, d_propuesto, s_propuesto), dimens = c(n, n),
            method = "rk4", atol = 1e-6, rtol = 1e-6)
    }, error = function(e) NULL)
    
    if (is.null(ST3) || nrow(ST3) != length(paquete$times_run)) return(NA)
    
    ST3_mat <- as.matrix(ST3[,-1])
    if(any(is.nan(ST3_mat))) return(NA)
    
    piso_minimo <- 1e-6 
    pred_freq_vec <- numeric(nrow(paquete$data))
    
    for(j in 1:nrow(paquete$data)) {
      t_abs <- paquete$data$Generation[j]
      time_idx <- match(t_abs, paquete$times_run)
      
      # Filtro anti-caídas espaciales (Spatial bounds check)
      x_idx <- paquete$data$X[j]
      y_idx <- paquete$data$Y[j]
      if(x_idx < 1 || x_idx > n || y_idx < 1 || y_idx > n) next 
      
      spatial_idx <- (y_idx - 1) * n + x_idx 
      pred_freq_raw <- ST3_mat[time_idx, spatial_idx]
      pred_freq_vec[j] <- max(min(pred_freq_raw, 1 - piso_minimo), piso_minimo)
    }
    
    # Calcular Verosimilitud (Usamos Binomial normal)
    ll_binom <- sum(dbinom(x = paquete$data$Count, size = paquete$data$Chr_Tot, prob = pred_freq_vec, log = TRUE))
    return(ll_binom)
    
  }, mc.cores = num_cores)
  
  # Limpiar fallos numéricos antes de sumar (Resuelve el problema del -5.29e11)
  ll_vec <- unlist(ll_lista)
  ll_vec <- ll_vec[is.finite(ll_vec)] 
  
  composite_ll <- sum(ll_vec)
  
  # Mostrar el progreso en los logs de SLURM
  cat("L-BFGS-B evaluando -> D:", D_propuesto, "| s:", s_propuesto, "==> LL Conjunta:", composite_ll, "\n")
  
  # Optim minimiza, así que devolvemos el valor negativo
  # Si todo falló, devolvemos una penalización gigante
  if(length(ll_vec) == 0 || is.na(composite_ll)) return(1e9)
  
  return(-composite_ll)
}

# =========================================================================
# =========================================================================
# 4. EJECUCIÓN DEL OPTIMIZADOR
# =========================================================================
cat("\nIniciando optimización L-BFGS-B...\n")

opt_result <- optim(
  par = c(0.01, 0.0), 
  fn = funcion_objetivo_conjunta, 
  method = "L-BFGS-B", 
  lower = c(0.001, -1.0), 
  upper = c(0.05, 1.0), 
  control = list(
    maxit = 20,
    trace = 1,     
    REPORT = 1,
    factr = 1e7    
  ),
  hessian = TRUE
)

# --- CÁLCULO DE ERRORES ESTÁNDAR A PARTIR DEL HESSIANO ---

# 1. Extraer la Matriz Hessiana
H <- opt_result$hessian

# [⚠️ ADVERTENCIA DE ESCALADO]: Si dejaste el 'factor_escala = 1e6' en tu función 
# objetivo para ayudar al optimizador, debes deshacer el escalado en el Hessiano
# antes de invertirlo. Si quitaste el factor de escala, borra la siguiente línea:
# H <- H * 1e6 

# 2. Calcular la Matriz de Covarianza invirtiendo el Hessiano
# Usamos tryCatch porque si la superficie es plana o está en un límite (ej. D=0.001 exacto), 
# la matriz no será invertible (singular) y R arrojaría un error deteniendo el script.
cov_matrix <- tryCatch({
  solve(H)
}, error = function(e) {
  cat("  [Aviso] Matriz Hessiana singular. No se pueden calcular Errores Estándar.\n")
  matrix(NA, nrow = 2, ncol = 2)
})

# 3. Extraer el Error Estándar (Raíz cuadrada de la varianza en la diagonal)
# Se verifica que el valor exista y sea positivo antes de aplicar la raíz.
se_D <- ifelse(!is.na(cov_matrix[1,1]) & cov_matrix[1,1] > 0, sqrt(cov_matrix[1,1]), NA)
se_s <- ifelse(!is.na(cov_matrix[2,2]) & cov_matrix[2,2] > 0, sqrt(cov_matrix[2,2]), NA)

# --- ENSAMBLAJE DE RESULTADOS ---
df_resultado <- data.frame(
  TaskID = task_id,
  Model = model_name,
  D_Optimo = opt_result$par[1],
  s_Optimo = opt_result$par[2],
  SE_D = se_D,       # <--- Error Estándar de D
  SE_s = se_s,       # <--- Error Estándar de s
  Max_Composite_LL = -opt_result$value,
  Convergence_Code = opt_result$convergence
)

# =========================================================================
# 5. GUARDADO DE RESULTADOS
# =========================================================================
out_file <- file.path(output_dir, paste0("MLE_LBFGSB_TaskID_", task_id, "_", model_name, ".txt"))
write.table(df_resultado, file = out_file, row.names = FALSE, quote = FALSE, sep = "\t")

cat("\nOptimización completada. Resultados guardados en:", out_file, "\n")

