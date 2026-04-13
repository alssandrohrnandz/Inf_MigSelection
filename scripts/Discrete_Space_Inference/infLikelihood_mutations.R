library(deSolve)
library(rootSolve)
library(dplyr)
library(tidyverse)
library(VGAM)

args <- commandArgs(trailingOnly = TRUE)
freq_file   <- args[1]
subset_file <- args[2]
task_id     <- args[3]
model_name  <- args[4]
output_dir  <- args[5]

#Prueba
#freq_file   <- "data/results_Discrete/outputs_slim/D_FULL_seleccion_m1_1.csv"
#subset_file <- "data/results_Discrete/subsets/subset_D_FULL_seleccion_m1_1.txt"
#task_id     <- 1
#model_name  <- "D_FULL_seleccion_m1"
#output_dir  <- "data/results_Discrete/outputs_LL"

if(length(args) < 4) {
  stop("Faltan argumentos. Se requieren: freq_file, subset_file, task_id, model_name")
}

print(paste("Procesando archivo:", freq_file))
print(paste("Modelo:", model_name))
print(paste("Task ID:", task_id))

# the modern model
diffusion2D <- function(t, conc, par) {
  Conc <- matrix(nrow = n, ncol = n, data = conc)
  dConc <- Conc*(1-Conc)*(Conc*par[2]+par[3]*(1-2*Conc))
  Flux <- -par[1] * rbind(rep(0, n), (Conc[2:n,]-Conc[1:(n-1),]), rep(0, n))/dx
  dConc <- dConc - (Flux[2:(n+1),] - Flux[1:n,])/dx
  Flux <- -par[1] * cbind(rep(0, n), (Conc[,2:n]-Conc[,1:(n-1)]), rep(0, n))/dy
  dConc <- dConc - (Flux[,2:(n+1)]-Flux[,1:n])/dy
  return(list(as.vector(dConc)))
}

GRID_SIZE <- 5
dy <- dx <- 1
n <- GRID_SIZE

s<-0.0
d <- 2*s # Parametrización

MAX_STEPS <- 10000
MIN_GENERATIONS <- 1
N_eff <- 10000

# parámetros a buscar TODO:EDITAR ESTO PORQUE PUEDE ESTAR MAL
DifussionValuesToCheck <- c(0.0000000001) #sort(unique(c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)))
SelectionValuesToCheck <- sort(unique(c(
  -0.0001,-0.0005,-0.00001,-0.00005,-0.000001,-0.000005,
  0.0001,0.0005,0.00001,0.00005,0.000001,0.000005,
  -0.5,-0.1,-0.05,-0.01,-0.005,-0.001,0,0.001,0.005,0.05,0.01,0.5,0.1,0.25,0.025,-0.025,
  -0.25,-1,1,0.002,0.02,0.025,0.0025
  )
  )
  )
print(DifussionValuesToCheck)
print(SelectionValuesToCheck)

# lectura del archivo
freq_data_raw <- read.csv(freq_file, header=TRUE, stringsAsFactors=FALSE)

# numero de columnas por si acaso para datos empiricos
num_cols <- ncol(freq_data_raw)
###TODO: ESTA PARTE HAY QUE EDITARLA PARA QUE ENTREN LOS DATOS EMPIRICOS
###NO SE TE VAYA A OLVIDAR POR FAVOR 
if (num_cols == 8) {
    # SLiM header: Generation,MutationID,MutType,X,Y,Frequency,AlleleCount,Chr_Tot
    #AADNA header: Time, SNP, Lat, Long, Freq, AlleleCount, Chr_Tot
    #TODO: ESTA PARTE HAY QUE MODIFICARLA PARA QUE ENTREN LOS DATOS EMPIRICOS SIN VASELINA
    colnames(freq_data_raw) <- c("Generation","MutationID","TypeMut","X","Y","Frequency","Count","Chr_Tot")
} else if (num_cols == 7) {
    colnames(freq_data_raw) <- c("Generation","MutationID","X","Y","Frequency","Count","Chr_Tot") #<- editar aqui paa datos empiricos

    freq_data_raw$TypeMut <- NA 
} else {
    stop(paste("El archivo tiene un número inesperado de columnas:", num_cols))
}

freq_data <- freq_data_raw

# filtramos subset
if (file.exists(subset_file) && file.info(subset_file)$size > 0) {
    snps_subset <- readLines(subset_file)
    snps_to_analyze <- intersect(unique(freq_data$MutationID), as.integer(snps_subset))
    print(paste("Analizando", length(snps_to_analyze), "SNPs indicados en el subset."))
} else {
    print("Subset vacío o inexistente, analizando todos los SNPs del archivo.")
    snps_to_analyze <- sort(unique(freq_data$MutationID))
}

# calculos auxiliares si no vinieran de SLiM (SLiM ya da Count y Chr_Tot, pero esto asegura enteros)
#freq_data$AlleleCount <- round(freq_data$Frequency * freq_data$Chr_Tot) 
freq_data$ChrOBS <- freq_data$Chr_Tot

# === 5. Bucle Principal ===

# Inicializar lista ANTES del loop para acumular los resultados de todos los SNPs
all_results_list <- list() 

for (snp_actual in snps_to_analyze) {
  
  df_snp <- freq_data[freq_data$MutationID == snp_actual, ]
  
  if(nrow(df_snp) == 0) next
  if(max(df_snp$Frequency) == 0) next
  
  First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
  
  if(is.na(First_OcurrenceData)) next
  
  AlleleOriginLat <- df_snp$Y[First_OcurrenceData]
  AlleleOriginLong <- df_snp$X[First_OcurrenceData]
  AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]

  unique_generations <- sort(unique(df_snp$Generation))
  
  # --- MODIFICACIÓN 1: Detectar si solo hay 1 generación ---
  FlagExtraRun <- 0
  if (length(unique_generations) == 1) {
    FlagExtraRun <- 1
  } else if (length(unique_generations) < MIN_GENERATIONS) {
    next # Salta SNPs que tienen más de 1 pero menos que MIN_GENERATIONS
  }
  
  # Matriz Inicial
  Conc0 <- matrix(0, nrow=n, ncol=n)
  
  origin_data <- df_snp[df_snp$Generation == AlleleOriginAge & df_snp$Frequency > 0, ]
  
  for(k in 1:nrow(origin_data)) {
    ox <- origin_data$X[k]
    oy <- origin_data$Y[k]
    if(ox >= 1 && ox <= n && oy >= 1 && oy <= n) {
        Conc0[ox, oy] <- origin_data$Frequency[k]
    }
  }

  if(sum(Conc0) == 0) next

  results <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)
  results$LL <- NA
  
  TimesToTest <- sort(unique(df_snp$Generation))
  
  # Si solo hay 1 generación, añadimos una generación extra para evaluar ceros
  if (FlagExtraRun == 1) {
     # Agregamos 1 generación extra al final
     TimesToTest <- c(TimesToTest, max(TimesToTest) + 1)
  }
  
  times_run <- min(TimesToTest):max(TimesToTest)
  
  for (i in 1:nrow(results)) {
    
    D_val <- results$D[i]
    s_val <- results$s[i]
    d <- 2*s_val
    pars <- c(D_val, d, s_val) 

    ST3 <- ode.2D(
      y      = Conc0,
      times  = times_run + 1,
      func   = diffusion2D,
      parms  = pars,
      dimens = c(n, n),
      method = rkMethod("rk45ck"),
      atol   = 1e-10, 
      rtol   = 1e-10,
      maxsteps = 1e5
    )
    
    ST3_mat <- as.matrix(ST3[,-1])
    ll <- 0
    
    # 1. Cálculo de Likelihood para datos reales
    for(j in 1:nrow(df_snp)) {
      obs_freq <- df_snp$Frequency[j]
      t_abs <- df_snp$Generation[j]
      time_idx <- match(t_abs, times_run)
      
      t_elapsed <- t_abs - AlleleOriginAge
      if (t_elapsed < 1) t_elapsed <- 0.5 
      rho_val <- 1 - exp(-t_elapsed / (2 * N_eff))
      rho_val <- max(rho_val, 1e-6)

      xg <- df_snp$X[j]
      yg <- df_snp$Y[j]
      spatial_idx <- (yg - 1) * n + xg 
      
      if(spatial_idx < 1 || spatial_idx > ncol(ST3_mat)) next

      pred_freq_raw <- ST3_mat[time_idx, spatial_idx]
      piso_minimo <- 1 / 1000000 
      pred_freq <- max(min(pred_freq_raw, 1 - piso_minimo), piso_minimo)
      
      ll <- ll + dbetabinom(
        x    = df_snp$Count[j],
        size = df_snp$Chr_Tot[j],
        prob = pred_freq,
        rho  = rho_val,
        log  = TRUE
      )
    }
    
    # 2. Cálculo de Likelihood para la generación extra (Matriz = 0)
    if (FlagExtraRun == 1){
      # Tomamos el tiempo extra que agregamos
      t_abs_extra <- max(TimesToTest)
      time_idx_extra <- match(t_abs_extra, times_run)
      
      t_elapsed <- t_abs_extra - AlleleOriginAge
      if (t_elapsed < 1) t_elapsed <- 0.5 
      rho_val <- 1 - exp(-t_elapsed / (2 * N_eff))
      rho_val <- max(rho_val, 1e-6)

      # Calculamos un tamaño de muestra representativo para esta población en "0"
      rep_chr_tot <- round(mean(df_snp$Chr_Tot, na.rm=TRUE))
      if(is.nan(rep_chr_tot) || rep_chr_tot < 1) rep_chr_tot <- 10 
      
      for (x_temp in 1:n){
         for (y_temp in 1:n){
            spatial_idx <- (y_temp - 1) * n + x_temp 
            
            if(spatial_idx < 1 || spatial_idx > ncol(ST3_mat)) next

            pred_freq_raw <- ST3_mat[time_idx_extra, spatial_idx]
            piso_minimo <- 1 / 1000000
            pred_freq <- max(min(pred_freq_raw, 1 - piso_minimo), piso_minimo)
            
            # Penaliza duramente si el modelo predice > 0
            ll <- ll + dbetabinom(
              x    = 0,
              size = rep_chr_tot, 
              prob = pred_freq,
              rho  = rho_val,
              log  = TRUE
            )
         }
      }
    }
    
    # Guardamos el LL final de esta combinación de parámetros
    results$LL[i] <- ll
  }

  # --- MODIFICACIÓN 2: Reestructurar los resultados (Filas = SNP, Columnas = Parámetros) ---
  
  # Crear nombres de columnas descriptivos, ej: "D_0.01_s_0.05"
  col_names <- paste0("D_", results$D, "_s_", results$s)
  
  # Crear una fila nueva para este SNP
  snp_row <- data.frame(SNP = snp_actual, stringsAsFactors = FALSE)
  
  # Asignar los valores de LL como columnas
  snp_row[col_names] <- results$LL
  
  # Guardarlo en la lista maestra
  all_results_list[[length(all_results_list) + 1]] <- snp_row
  print(paste("Procesado SNP:", snp_actual))
}

# --- MODIFICACIÓN 3: Guardar TODO en un solo archivo por task_id ---
if(length(all_results_list) > 0) {
  # Unimos todas las filas
  final_results <- do.call(rbind, all_results_list)
  
  output_filename <- paste0("Analysis_", model_name, "_TaskID_", task_id, "_All_SNPs.txt")
  output_path <- file.path(output_dir, output_filename)
  
  print(paste("Guardando archivo acumulado:", output_path))
  # Se guarda separado por tabulaciones (.tsv pero con extension .txt)
  write.table(final_results, file=output_path, row.names=FALSE, quote = FALSE, sep="\t")
} else {
  print("No se procesaron SNPs que cumplieran los criterios.")
}