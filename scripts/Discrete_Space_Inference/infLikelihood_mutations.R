library(deSolve)
library(rootSolve)
library(dplyr)
library(tidyverse)
library(VGAM)

# === 1. Lectura de Argumentos ===
args <- commandArgs(trailingOnly = TRUE)
freq_file   <- args[1]
subset_file <- args[2]
task_id     <- args[3]
model_name  <- args[4]

if(length(args) < 4) {
  stop("Faltan argumentos. Se requieren: freq_file, subset_file, task_id, model_name")
}

print(paste("Procesando archivo:", freq_file))
print(paste("Modelo:", model_name))
print(paste("Task ID:", task_id))

# === 2. Configuración del Modelo ===
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
d <- 0.5
MAX_STEPS <- 10000
MIN_GENERATIONS <- 3
N_eff <- 1000

# Parámetros de búsqueda
DifussionValuesToCheck <- seq(0.0000001,0.000005,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)
SelectionValuesToCheck <- 0.0 #seq(-0.05, 0.05, by=0.01)

# === 3. Lectura de Datos Inteligente ===

# Leemos con header=TRUE porque SLiM lo escribe
freq_data_raw <- read.csv(freq_file, header=TRUE, stringsAsFactors=FALSE)

# Detectamos estructura de columnas dinámicamente
num_cols <- ncol(freq_data_raw)

if (num_cols == 8) {
    # Caso aDNA (Tiene MutType)
    # SLiM header: Generation,MutationID,MutType,X,Y,Frequency,AlleleCount,Chr_Tot
    colnames(freq_data_raw) <- c("Generation","MutationID","TypeMut","X","Y","Frequency","Count","Chr_Tot")
} else if (num_cols == 7) {
    # Caso FULL (No tiene MutType)
    # SLiM header: Generation,MutationID,X,Y,Frequency,AlleleCount,Chr_Tot
    colnames(freq_data_raw) <- c("Generation","MutationID","X","Y","Frequency","Count","Chr_Tot")
    # Creamos una columna dummy para que el código no falle si llamas a TypeMut luego
    freq_data_raw$TypeMut <- NA 
} else {
    stop(paste("El archivo tiene un número inesperado de columnas:", num_cols))
}

freq_data <- freq_data_raw

# === 4. Filtrado por Subset (Uso del archivo de Bash) ===
if (file.exists(subset_file) && file.info(subset_file)$size > 0) {
    snps_subset <- readLines(subset_file)
    # Filtramos para analizar solo lo que Bash nos dijo
    # Convertimos a integer por seguridad
    snps_to_analyze <- intersect(unique(freq_data$MutationID), as.integer(snps_subset))
    print(paste("Analizando", length(snps_to_analyze), "SNPs indicados en el subset."))
} else {
    print("Subset vacío o inexistente, analizando todos los SNPs del archivo.")
    snps_to_analyze <- sort(unique(freq_data$MutationID))
}

# Cálculos auxiliares si no vinieran de SLiM (SLiM ya da Count y Chr_Tot, pero esto asegura enteros)
freq_data$AlleleCount <- round(freq_data$Frequency * freq_data$Chr_Tot) 
freq_data$ChrOBS <- freq_data$Chr_Tot

# === 5. Bucle Principal ===

for (snp_actual in snps_to_analyze) {
  
  df_snp <- freq_data[freq_data$MutationID == snp_actual, ]
  
  # Validaciones iniciales
  if(nrow(df_snp) == 0) next
  if(max(df_snp$Frequency) == 0) next
  
  First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
  
  if(is.na(First_OcurrenceData)) next
  
  AlleleOriginLat <- df_snp$Y[First_OcurrenceData]
  AlleleOriginLong <- df_snp$X[First_OcurrenceData]
  AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]

  unique_generations <- sort(unique(df_snp$Generation))
  if(length(unique_generations) < MIN_GENERATIONS) {
    # warning(paste("SNP", snp_actual, "pocas generaciones. Skip."))
    next 
  }
  
  # Matriz Inicial
  Conc0 <- matrix(0, nrow=n, ncol=n)
  
  # Llenamos condiciones iniciales usando todos los puntos donde aparece por primera vez
  origin_data <- df_snp[df_snp$Generation == AlleleOriginAge & df_snp$Frequency > 0, ]
  
  for(k in 1:nrow(origin_data)) {
    ox <- origin_data$X[k]
    oy <- origin_data$Y[k]
    # Validación de límites de array por si acaso
    if(ox >= 1 && ox <= n && oy >= 1 && oy <= n) {
        Conc0[ox, oy] <- origin_data$Frequency[k]
    }
  }

  if(sum(Conc0) == 0) next

  # Grid Search
  results <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)
  results$LL <- NA
  
  TimesToTest <- sort(unique(df_snp$Generation))
  times_run <- min(TimesToTest):max(TimesToTest)

  for (i in 1:nrow(results)) {
    
    D_val <- results$D[i]
    s_val <- results$s[i]
    pars <- c(D_val, d, s_val) 
    
    # ODE
    ST3 <- ode.2D(
      y      = Conc0,
      times  = times_run,
      func   = diffusion2D,
      parms  = pars,
      dimens = c(n, n),
      method = rkMethod("rk45ck"),
      atol   = 1e-6,
      rtol   = 1e-6,
      maxsteps = 1e5
    )
    
    ST3_mat <- as.matrix(ST3[,-1])
    
    # Likelihood Calculation
    ll <- 0
    
    for(j in 1:nrow(df_snp)) {
      t_abs <- df_snp$Generation[j]
      time_idx <- match(t_abs, times_run)
      
      if(is.na(time_idx)) next 
      
      t_elapsed <- t_abs - AlleleOriginAge
      if (t_elapsed < 1) t_elapsed <- 0.5 
      
      xg <- df_snp$X[j]
      yg <- df_snp$Y[j]
      spatial_idx <- (xg - 1) * n + yg
      
      # Protección contra índices fuera de rango
      if(spatial_idx < 1 || spatial_idx > ncol(ST3_mat)) next

      pred_freq <- ST3_mat[time_idx, spatial_idx]
      
      piso_minimo <- 1 / (2 * N_eff)
      pred_freq <- max(min(pred_freq, 1 - piso_minimo), piso_minimo)
      
      rho_val <- 1 - exp(-t_elapsed / (2 * N_eff))
      rho_val <- max(rho_val, 1e-6)
      
      ll <- ll + dbetabinom(
        x    = df_snp$AlleleCount[j],
        size = df_snp$ChrOBS[j],
        prob = pred_freq,
        rho  = rho_val,
        log  = TRUE
      )
    }
    results$LL[i] <- ll
  }

  # Guardar Resultados
  # Construimos la ruta de salida basada en la entrada para no usar "../" ciego
  # Asumimos estructura: data/results_LL/
  # El script de bash definio resultados en data/results_LL
  
  # Extraemos el directorio base del archivo de entrada
  input_dir <- dirname(freq_file) 
  # Subimos un nivel y entramos a results_LL (asumiendo estructura estandar)
  # Si input es .../data/results_simulations/archivo.csv -> .../data/results_LL/
  output_dir <- file.path(dirname(input_dir), "results_LL")
  
  # Si no existe, usamos el actual
  if(!dir.exists(output_dir)) output_dir <- "."
  
  output_filename <- paste0("Analysis_", model_name, "_SNP_", snp_actual, ".txt")
  output_path <- file.path(output_dir, output_filename)
  
  best_params <- results[which.max(results$LL),]
  # print(paste("Guardando:", output_path, "Max LL:", best_params$LL))
  
  write.table(results, file=output_path, row.names=FALSE, quote = FALSE)
}
