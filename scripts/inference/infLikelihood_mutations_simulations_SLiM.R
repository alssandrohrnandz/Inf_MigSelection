library(deSolve)
library(rootSolve)
library(dplyr)
library(tidyverse)
library(VGAM)

# ==========================================
# CONFIGURACIÓN DE RUTAS Y ARCHIVOS
# ==========================================
# Directorio donde están los CSV generados por SLiM
input_dir <- "/Users/alessandrohernandez/Documents/Posgrado/Inf_MigSelection/Inference_SLiM/results_simulations/" 
# Directorio donde quieres guardar los resultados de Verosimilitud
output_dir <- "/Users/alessandrohernandez/Documents/Posgrado/Inf_MigSelection/Inference_SLiM/results_LL/"

# Aseguramos que el directorio de salida exista
if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }

setwd(input_dir) # Establecemos directorio de trabajo en los inputs

# Lista de los 6 archivos exactos que mostraste en la imagen
files_to_analyze <- c(
  "C_aDNA_scattered_rep1.csv",
  "D_aDNA_scattered_rep1.csv",
  "FULL_C_neutros_m1_rep1.csv",
  "FULL_C_seleccion_m2_rep1.csv",
  "FULL_D_neutros_m1_rep1.csv",
  "FULL_D_seleccion_m2_rep1.csv"
)

# ==========================================
# DEFINICIÓN DEL MODELO (Igual que antes)
# ==========================================

# Parámetros Globales del Grid
GRID_SIZE <- 10
dy <- dx <- 1
n <- GRID_SIZE
d <- 0.5    # Dominancia
N_eff <- 1000 # Tamaño efectivo local (para la deriva)

# Modelo de difusión (EDP)
diffusion2D <- function(t, conc, par) {
  # par[1] = D, par[2] = d, par[3] = s
  Conc <- matrix(nrow = n, ncol = n, data = conc)
  
  # Término de Reacción (Selección)
  # asumiendo gamma = p(1-p)(pd + s(1-2p)) ajustado a tus parms
  dConc <- Conc*(1-Conc)*(Conc*par[2]+par[3]*(1-2*Conc))
  
  # Término de Difusión en X
  Flux_X <- -par[1] * rbind(rep(0, n), (Conc[2:n,]-Conc[1:(n-1),]), rep(0, n))/dx
  dConc <- dConc - (Flux_X[2:(n+1),] - Flux_X[1:n,])/dx
  
  # Término de Difusión en Y
  Flux_Y <- -par[1] * cbind(rep(0, n), (Conc[,2:n]-Conc[,1:(n-1)]), rep(0, n))/dy
  dConc <- dConc - (Flux_Y[,2:(n+1)]-Flux_Y[,1:n])/dy
  
  return(list(as.vector(dConc)))
}

# Espacio de Búsqueda (Grid Search)
# Puedes ajustar esto si quieres "Francotirador" o búsqueda amplia
DifussionValuesToCheck <- c(0.01, 0.02, 0.05) # Ejemplo: probando varios D
SelectionValuesToCheck <- seq(0.0, 0.1, by=0.01)  


# ==========================================
# BUCLE MAESTRO: ITERAR SOBRE LOS 6 ARCHIVOS
# ==========================================

for (current_file in files_to_analyze) {
  
  print(paste("----------------------------------------------------------------"))
  print(paste("ANALIZANDO ARCHIVO:", current_file))
  print(paste("----------------------------------------------------------------"))
  
  full_path <- file.path(input_dir, current_file)
  
  # 1. LECTURA DE DATOS
  # IMPORTANTE: header=TRUE porque tu SLiM nuevo pone encabezados
  if(!file.exists(full_path)) {
    warning(paste("El archivo", current_file, "no existe. Saltando..."))
    next
  }
  
  freq_data <- read.csv(full_path, header=TRUE)
  
  # 2. ESTANDARIZACIÓN DE COLUMNAS
  # Los archivos SLiM tienen: Generation, MutationID, [MutType opcional], X, Y, Frequency, AlleleCount, Chr_Tot
  # Necesitamos mapear Chr_Tot a ChrOBS para el script de inferencia
  
  # Renombrar Chr_Tot a ChrOBS para mantener consistencia
  if("Chr_Tot" %in% colnames(freq_data)) {
    freq_data <- freq_data %>% rename(ChrOBS = Chr_Tot)
  }
  
  # Asegurar que AlleleCount existe (SLiM lo genera, pero por seguridad)
  if(!"AlleleCount" %in% colnames(freq_data)) {
    # Si por alguna razón no está, lo calculamos (backup)
    freq_data$AlleleCount <- round(freq_data$Frequency * freq_data$ChrOBS)
  }
  
  # 3. FILTRADO DE SNPS
  # Obtenemos la lista única de Mutaciones en este archivo
  snps_to_analyze <- freq_data %>% pull(MutationID) %>% sort() %>% unique()
  
  print(paste("SNPs encontrados en este archivo:", length(snps_to_analyze)))
  
  # ==========================================
  # BUCLE INTERNO: ITERAR SOBRE SNPS
  # ==========================================
  
  for (snp_actual in snps_to_analyze) {
    
    # Filtrar datos para el SNP actual
    df_snp <- freq_data %>% filter(MutationID == snp_actual)
    
    # Identificar Origen (Primera vez que Frequency > 0)
    First_Occur <- df_snp %>% filter(Frequency > 0) %>% slice(1)
    
    if(nrow(First_Occur) == 0) {
      warning(paste("SNP", snp_actual, "nunca tiene frecuencia > 0. Se salta."))
      next
    }
    
    AlleleOriginLat <- First_Occur$Y
    AlleleOriginLong <- First_Occur$X
    AlleleOriginAge <- First_Occur$Generation
    
    # Validar número de generaciones
    unique_generations <- sort(unique(df_snp$Generation))
    MIN_GENERATIONS <- 3
    if(length(unique_generations) < MIN_GENERATIONS) {
      # Opcional: imprimir aviso silencioso para no llenar la consola si son muchos
      # warning(paste("SNP", snp_actual, "pocas generaciones. Se omite."))
      next 
    }
    
    # Configurar Tiempos de Simulación
    TimesToTest <- sort(unique(df_snp$Generation))
    
    # Cortamos la simulación para que empiece en el origen del alelo
    # (Aunque times_run debe ser continuo para la ODE)
    times_run <- AlleleOriginAge:max(TimesToTest)
    
    # Matriz Inicial (Conc0)
    Conc0 <- matrix(0, nrow=n, ncol=n)
    
    # Llenar condiciones iniciales con los datos observados en el tiempo de origen
    origin_data <- df_snp %>% filter(Generation == AlleleOriginAge & Frequency > 0)
    
    valid_start <- FALSE
    for(k in 1:nrow(origin_data)) {
      ox <- origin_data$X[k]
      oy <- origin_data$Y[k]
      # Validación de límites de grid (por seguridad)
      if(ox >= 1 && ox <= n && oy >= 1 && oy <= n) {
        Conc0[ox, oy] <- origin_data$Frequency[k]
        valid_start <- TRUE
      }
    }
    
    if(!valid_start || sum(Conc0) == 0) {
      next
    }
    
    # ---- GRID SEARCH (D y s) ----
    results <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)
    results$LL <- NA
    
    for (i in 1:nrow(results)) {
      
      D_val <- results$D[i]
      s_val <- results$s[i]
      pars <- c(D_val, d, s_val) 
      
      # 1. Ejecutar ODE
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
      
      ST3_mat <- as.matrix(ST3[,-1])  # matriz resultados (filas=tiempo, cols=celdas)
      
      # 2. Likelihood Beta-Binomial
      ll <- 0
      
      for(j in 1:nrow(df_snp)) {
        
        t_abs <- df_snp$Generation[j]
        
        # Solo evaluamos si el tiempo está dentro de la simulación
        if(t_abs < AlleleOriginAge) next
        
        time_idx <- match(t_abs, times_run)
        if(is.na(time_idx)) next 
        
        # Tiempo relativo para deriva
        t_elapsed <- t_abs - AlleleOriginAge
        if (t_elapsed < 1) t_elapsed <- 0.5 
        
        xg <- df_snp$X[j]
        yg <- df_snp$Y[j]
        
        # Indexación espacial (vectorización por columna es estándar en deSolve/R)
        # SLiM output: X (1..10), Y (1..10). 
        # Fórmula: (Fila - 1) * nCols + Columna  -> Si ode.2D vectoriza por filas?
        # En R 'as.vector(matrix)' va por columnas.
        # Si Conc[x,y], entonces idx = (y-1)*n + x ó (x-1)*n + y?
        # Probemos estándar: (Columna - 1) * n + Fila. Asumiendo X=Fila, Y=Columna en la matriz
        spatial_idx <- (yg - 1) * n + xg 
        
        pred_freq <- ST3_mat[time_idx, spatial_idx]
        
        # Pisos mínimos
        piso_minimo <- 1 / (2 * N_eff)
        pred_freq <- max(min(pred_freq, 1 - piso_minimo), piso_minimo)
        
        # Rho (Deriva)
        rho_val <- 1 - exp(-t_elapsed / (2 * N_eff))
        rho_val <- max(rho_val, 1e-6)
        
        # Sumar Log-Likelihood
        ll <- ll + VGAM::dbetabinom(
          x    = df_snp$AlleleCount[j],
          size = df_snp$ChrOBS[j],
          prob = pred_freq,
          rho  = rho_val,
          log  = TRUE
        )
      }
      
      results$LL[i] <- ll
    }
    
    # Guardar el mejor resultado para consola
    best_params <- results[which.max(results$LL),]
    print(paste("File:", current_file, "| SNP:", snp_actual, "| Best D:", best_params$D, "| Best s:", best_params$s, "| Max LL:", round(best_params$LL, 2)))
    
    # ---- GUARDAR ARCHIVO DE SALIDA ----
    # Nombre dinámico: ANALISIS + NombreArchivoOriginal + SNP + .txt
    
    # Limpiamos la extensión .csv del nombre original para que quede limpio
    clean_filename <- tools::file_path_sans_ext(current_file)
    
    output_filename <- paste0("LL_Analysis_", clean_filename, "_SNP_", snp_actual, ".txt")
    full_output_path <- file.path(output_dir, output_filename)
    
    write.table(results, file=full_output_path, row.names=FALSE, quote = FALSE)
    
  } # Fin Bucle SNP
} # Fin Bucle Archivos

print("¡Análisis de todos los archivos completado!")