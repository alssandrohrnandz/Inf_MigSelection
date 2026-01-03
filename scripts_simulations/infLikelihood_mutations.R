library(deSolve)
library(rootSolve)
library(dplyr)
library(tidyverse)
library(VGAM)

args <- commandArgs(trailingOnly = TRUE)
freq_file <- "/Users/alessandrohernandez/Documents/Posgrado/SimulacionesSLIM/Discrete/dataset_completo_10.csv" #args[1]
subset_file <-"../archivos_descargados/alelos_subset_1.txt" #args[2]
task_id <- "SELECTION"#as.integer(args[3])
setwd("/Users/alessandrohernandez/Documents/Posgrado/Doctorado/Discrete/scripts")
directorio<-getwd()


print(freq_file)
print(subset_file)
print(task_id)
print(directorio)
# Modelo de difusión
diffusion2D <- function(t, conc, par) {
  Conc <- matrix(nrow = n, ncol = n, data = conc)
  dConc <- Conc*(1-Conc)*(Conc*par[2]+par[3]*(1-2*Conc))
  Flux <- -par[1] * rbind(rep(0, n), (Conc[2:n,]-Conc[1:(n-1),]), rep(0, n))/dx
  dConc <- dConc - (Flux[2:(n+1),] - Flux[1:n,])/dx
  Flux <- -par[1] * cbind(rep(0, n), (Conc[,2:n]-Conc[,1:(n-1)]), rep(0, n))/dy
  dConc <- dConc - (Flux[,2:(n+1)]-Flux[,1:n])/dy
  return(list(as.vector(dConc)))
}
  
  # Parámetros
GRID_SIZE <- 10
dy <- dx <- 1
n <- GRID_SIZE
d <- 0.5
s <- 0.0
MAX_STEPS <- 10000
MIN_GENERATIONS <- 3
N_eff <- 1000

# Francotirador V2 (Para m = 0.005)
# Si tu máximo fue 0.03, el real probablemente ande por 0.05 - 0.08
# Estrategia Francotirador (m = 0.005)
# Francotirador V2 (Para m = 0.005)
# Si tu máximo fue 0.03, el real probablemente ande por 0.05 - 0.08
#DifussionValuesToCheck <- c(
#  0,
#  seq(0.0001, 0.001, by = 0.0001),
#  seq(0.002,0.01,by=0.001),
#  seq(0.02,0.1,by=0.01)
#)
DifussionValuesToCheck <- 0.01
SelectionValuesToCheck <- seq(-0.05,0.05, by=0.01)  


freq_data <- read.csv(freq_file, header=FALSE, sep=',')
colnames(freq_data)<-c("Generation","MutationID","TypeMut","X","Y","Frequency","Count")
head(freq_data)
snps_to_analyze <- freq_data[,2] %>% sort() %>% unique()
## snps_to_analyze <- freq_data[freq_data$TypeMut==1,] %>% pull(2) %>% sort() %>% unique()
#freq_data <- read.csv("frecuencias_mutaciones_m1_10_600Gen.csv", header=TRUE)
#snps_to_analyze <- readLines("subsets/alelos_subset_1.txt.")
#freq_data_subset <- freq_data[freq_data$SNP_ID %in% snps_to_analyze, ]
freq_data$AlleleCount<-round(freq_data$Frequency*freq_data$Count*1)
freq_data$ChrOBS<-round(freq_data$Count*2)


for (snp_actual in snps_to_analyze) {
  df_snp <- freq_data[freq_data$MutationID == snp_actual, ]
  First_OcurrenceData <- which(df_snp$Frequency > 0)[1]
  AlleleOriginLat <- df_snp$Y[First_OcurrenceData]
  AlleleOriginLong <- df_snp$X[First_OcurrenceData]
  AlleleOriginAge <- df_snp$Generation[First_OcurrenceData]

  #print(paste("Allele origin age:", AlleleOriginAge))
  #print(paste("Allele origin lat:", AlleleOriginLat))
  #print(paste("Allele origin long:", AlleleOriginLong))
  origin_locations <- df_snp %>% filter(Generation == AlleleOriginAge & Frequency > 0) %>% distinct(X, Y)
 
 #if (nrow(origin_locations) > 1) {
  #warning(paste("SNP", snp_actual, "se originó en", nrow(origin_locations), "ubicaciones (X,Y) en Gen.", AlleleOriginAge, ". Se omite."))
  #next
  #}
  origin_data <- df_snp[df_snp$Generation == AlleleOriginAge & df_snp$Frequency > 0, ]

  unique_generations <- sort(unique(df_snp$Generation))
  if(length(unique_generations) < MIN_GENERATIONS) {
    warning(paste("SNP", snp_actual, 
                  "tiene solo", length(unique_generations),
                  "generaciones. Se omite."))
    next  # Salta este SNP y continúa con el siguiente
  }
  
  TimesToTest<-unique(sort(df_snp$Generation))
  start_time_index <- which(TimesToTest == AlleleOriginAge)
  times_simulacion <- TimesToTest[start_time_index:length(TimesToTest)]
  Conc0 <- matrix(0, nrow=n, ncol=n)

  origin_idx <- which(df_snp$Generation == AlleleOriginAge &
                      df_snp$X == AlleleOriginLong &
                      df_snp$Y == AlleleOriginLat)
  
  for(k in 1:nrow(origin_data)) {
    ox <- origin_data$X[k]
    oy <- origin_data$Y[k]
    # Usamos la frecuencia observada en cada punto
    Conc0[ox, oy] <- origin_data$Frequency[k]
  }
  #if(length(origin_idx) > 0) {
  #  Origin_X_Grid <- df_snp$X[origin_idx[1]]
  #  Origin_Y_Grid <- df_snp$Y[origin_idx[1]]
  #  Conc0 <- matrix(0, nrow=n, ncol=n)
  #  Conc0[Origin_X_Grid, Origin_Y_Grid] <- df_snp$Frequency[origin_idx[1]]
  #} else {
  #  stop("Error: No se encontró el punto de origen del alelo para Conc0.")
  #
  #}
  if(sum(Conc0) == 0) {
    warning(paste("SNP", snp_actual, "tiene Conc0 vacía. Se salta."))
    next
  }
  results <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)
  #print(str(times_simulacion))
  results$LL <- NA



  for (i in 1:nrow(results)) {
    
    D <- results$D[i]
    s_val <- results$s[i]
    # Asumo que 'd' está definido fuera de este bucle en tu entorno global
    pars <- c(D, d, s_val) 
    
    # ---- 1. Construir tiempos EXACTAMENTE como en el original ----
    TimesToTest <- sort(unique(df_snp$Generation))
    times_run <- min(TimesToTest):max(TimesToTest)
    
    # ---- 2. Ejecutar ODE con RK45 ----
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
    
    # Convertimos ST3 en matriz para acceso rápido
    ST3_mat <- as.matrix(ST3[,-1])  # quitamos columna "time"
    
    # ---- 3. Likelihood (Beta-Binomial) ----
    ll <- 0
    
    for(j in 1:nrow(df_snp)) {
      
      # Tiempo absoluto (para buscar en la matriz de difusión)
      t_abs <- df_snp$Generation[j]
      time_idx <- match(t_abs, times_run)
      
      if(is.na(time_idx)) next 
      
      # --- CORRECCIÓN CRÍTICA: TIEMPO RELATIVO PARA DERIVA ---
      # Calculamos cuántas generaciones han pasado desde el origen del alelo
      t_elapsed <- t_abs - AlleleOriginAge
      
      # Si estamos en el mismo momento del origen, t_elapsed es 0. 
      # Ajustamos a un mínimo pequeño para evitar rho=0 si N es grande
      if (t_elapsed < 1) t_elapsed <- 0.5 
      
      # Indexación espacial correcta
      xg <- df_snp$X[j]
      yg <- df_snp$Y[j]
      spatial_idx <- (xg - 1) * n + yg
      
      pred_freq <- ST3_mat[time_idx, spatial_idx]
      # Evitar 0 y 1 absolutos
      # El piso no debería ser menor a la inversa del tamaño poblacional muestreado.
      piso_minimo <- 1 / (2 * N_eff) # Ej: 1/2000 = 0.0005
      pred_freq <- max(min(pred_freq, 1 - piso_minimo), piso_minimo)
      
      # --- CÁLCULO DE RHO (DRIFT) ---
      # Usamos t_elapsed en lugar de t_abs
      rho_val <- 1 - exp(-t_elapsed / (2 * N_eff))
      
      # Ajuste de seguridad
      rho_val <- max(rho_val, 1e-6)
      
      # --- LIKELIHOOD BETA-BINOMIAL ---
      ll <- ll + dbetabinom(
        x    = df_snp$AlleleCount[j],
        size = df_snp$ChrOBS[j],
        prob = pred_freq,
        rho  = rho_val,
        log  = TRUE
      )
    }
    
    results$LL[i] <- ll
    #print(paste("D:", D, "s:", s_val, "LL:", ll))
  }

  best_params <- results[which.max(results$LL),]
  print(paste("Best D:", best_params$D, "Best s:", best_params$s, "Max LL:", best_params$LL, "Alelo:",snp_actual ))
  output_file <- paste0("../LL/SelectionDiffusion_Analysis_",snp_actual,"_m_",task_id,".txt", sep="")
  write.table(results, file=output_file, row.names=FALSE, quote = FALSE)
}


