install.packages("bbmle")
library(deSolve)
library(rootSolve)
library(dplyr)
library(tidyverse)
library(VGAM)
library(bbmle)

args <- commandArgs(trailingOnly = TRUE)
freq_file   <- "C:/Users/miel_/Documents/Doctorado/input/D_FULL_neutros_m1_1.csv" #args[1]
subset_file <- "C:/Users/miel_/Documents/Doctorado/input/subset_D_FULL_neutros_m1_1.txt" #args[2]
task_id     <- 1 #args[3]
model_name <- "Prueba_Holoceno"
output_dir  <- "C:/Users/miel_/Documents/Doctorado/output" #args[5]


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

GRID_SIZE <- 10
dy <- dx <- 1
n <- GRID_SIZE

s<-0.0
d <- 2*s # Parametrización

MAX_STEPS <- 10000
MIN_GENERATIONS <- 2
N_eff <- 1000

# parámetros a buscar TODO:EDITAR ESTO PORQUE PUEDE ESTAR MAL
DifussionValuesToCheck <- c(0.0000000001) #sort(unique(c(0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1)))
SelectionValuesToCheck <- sort(unique(c(-0.5,-0.1,-0.05,-0.01,-0.005,-0.001,0,0.001,0.005,0.05,0.01,0.5,0.1,0.25,0.025,-0.025,-0.25,-1,1,0.002,0.02,0.025,0.0025)))


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
snp_actual<-932
prueba_manual <- objective_function(D_val = 0.01, s_val = 0.01)
print(prueba_manual)
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
  
  TimesToTest <- sort(unique(df_snp$Generation))
  times_run <- min(TimesToTest):max(TimesToTest)
  
  # 1. Definir la función objetivo (AFUERA de la ejecución del optimizador)
  # Notar que D_val y s_val ahora son argumentos directos para que mle2 los entienda
  objective_function <- function(D_val, s_val, df_snp_data = df_snp, init_Conc = Conc0, trun = times_run, grid_n = n, Ne = N_eff) {
    
    # Límites físicos para evitar que el optimizador colapse con valores imposibles
    if (D_val <= 0 || s_val < -1 || s_val > 1) return(1e10)
    
    d_val <- 2 * s_val
    
    # Ejecutar simulación
    ST3 <- ode.2D(
      y = init_Conc, times = trun, func = diffusion2D,
      parms = c(D_val, d_val, s_val), dimens = c(grid_n, grid_n),
      method = "rk4", atol = 1e-7, rtol = 1e-7
    )
    
    ST3_mat <- as.matrix(ST3[,-1])
    ll <- 0
    piso_minimo <- 1e-6
    
    for(j in 1:nrow(df_snp_data)) {
      t_abs <- df_snp_data$Generation[j]
      time_idx <- match(t_abs, trun)
      
      # Cálculo de deriva (rho)
      t_elapsed <- max(t_abs - min(df_snp_data$Generation), 0.5)
      rho_val <- 1 - exp(-t_elapsed / (2 * Ne))
      rho_val <- max(rho_val, 1e-6)
      
      xg <- df_snp_data$X[j]; yg <- df_snp_data$Y[j]
      spatial_idx <- (yg - 1) * grid_n + xg
      
      pred_freq_raw <- ST3_mat[time_idx, spatial_idx]
      pred_freq <- max(min(pred_freq_raw, 1 - piso_minimo), piso_minimo)
      
      ll <- ll + dbetabinom(
        x = df_snp_data$Count[j], size = df_snp_data$Chr_Tot[j],
        prob = pred_freq, rho = rho_val, log = TRUE
      )
    }
    return(-ll) 
  }
  
  # 2. Ejecutar el optimizador (Asegúrate de tener library(bbmle) al inicio de tu script)
  start_params <- list(D_val = 0.01, s_val = 0.01)
  
  # tryCatch evita que el script completo muera si un SNP no converge
  # ... (resto del código igual) ...
  fit_result <- tryCatch({
    mle2(
      minuslogl = objective_function, 
      start = start_params,
      method = "L-BFGS-B", 
      lower = c(D_val = 0.0000000001, s_val = -0.5),
      upper = c(D_val = 0.0000000009, s_val = 0.5),
      control = list(maxit = 300)
    )
  }, error = function(e) {
    # ¡AQUÍ ESTÁ EL CAMBIO! Agregamos e$message para ver el error real
    cat("Error optimizando SNP:", snp_actual, "- Razón:", e$message, "\n")
    return(NULL)
  })
  
  # 3. Guardar Resultados de manera limpia
  if (!is.null(fit_result)) {
    best_pars <- coef(fit_result)
    # mle2 minimiza el negativo, así que lo volvemos positivo para lectura humana
    max_ll <- -logLik(fit_result) 
    
    # Armamos el dataframe de salida directamente con el mejor valor
    results_out <- data.frame(
      SNP = snp_actual,
      D = best_pars["D_val"],
      s = best_pars["s_val"],
      LL = max_ll
    )
    
    output_filename <- paste0("Analysis_", model_name,"_" ,task_id,"_SNP_", snp_actual, ".txt")
    output_path <- file.path(output_dir, output_filename)
    
    print(paste("Guardando:", output_path, Best,"Max LL:", round(max_ll, 2)))
    write.table(results_out, file=output_path, row.names=FALSE, quote = FALSE)
  }
} # Fin del bucle for principal
