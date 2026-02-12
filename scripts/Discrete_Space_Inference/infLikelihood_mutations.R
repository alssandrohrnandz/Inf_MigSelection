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
d <- 0.5
s<-0.0
MAX_STEPS <- 10000
MIN_GENERATIONS <- 2
N_eff <- 1000

# parámetros a buscar TODO:EDITAR ESTO PORQUE PUEDE ESTAR MAL
DifussionValuesToCheck <- sort(unique(c(0.0000001, 0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1)))
SelectionValuesToCheck <- sort(unique(c(-0.5,-0.1,-0.05,-0.01,0,0.05,0.01,0.5,0.1,0.25,0.025,-0.025,-0.25,-1,1)))


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

  results <- expand.grid(D=DifussionValuesToCheck, s=SelectionValuesToCheck)
  results$LL <- NA
  
  TimesToTest <- sort(unique(df_snp$Generation))
  times_run <- min(TimesToTest):max(TimesToTest)

  for (i in 1:nrow(results)) {
    
    D_val <- results$D[i]
    s_val <- results$s[i]
    pars <- c(D_val, d, s_val) 

    ST3 <- ode.2D(
      y      = Conc0,
      times  = times_run,
      func   = diffusion2D,
      parms  = pars,
      dimens = c(n, n),
      method = rkMethod("rk45ck"),
      atol   = 1e-10,
      rtol   = 1e-10,
      maxsteps = 1e5
    )
    
    #elimnar la columnas "time"
    ST3_mat <- as.matrix(ST3[,-1])
    
    ll <- 0
    

    # if(i == 1) print(paste("Max freq en matriz:", max(ST3_mat)))

    for(j in 1:nrow(df_snp)) {
      t_abs <- df_snp$Generation[j]
      time_idx <- match(t_abs, times_run)
      
      #if(is.na(time_idx)) next 
      
      # Cálculo de deriva
      t_elapsed <- t_abs - AlleleOriginAge
      if (t_elapsed < 1) t_elapsed <- 0.5 
      rho_val <- 1 - exp(-t_elapsed / (2 * N_eff))
      rho_val <- max(rho_val, 1e-6)

      xg <- df_snp$X[j]
      yg <- df_snp$Y[j]
      
      spatial_idx <- (yg - 1) * n + xg #indexacion por ver 
      

      if(spatial_idx < 1 || spatial_idx > ncol(ST3_mat)) next

      pred_freq_raw <- ST3_mat[time_idx, spatial_idx]
      
      # Piso mínimo (Evita log(0))
      piso_minimo <- 1e-6
      pred_freq <- max(min(pred_freq_raw, 1 - piso_minimo), piso_minimo)
      
     
      # Si pred_freq siempre es 0.0005, el modelo "no ve" el alelo en esa coordenada
      # if(i==1 && j < 5) print(paste("Gen:", t_abs, "X:", xg, "Y:", yg, "Raw:", pred_freq_raw, "Final:", pred_freq))

      ll <- ll + dbetabinom(
        x    = df_snp$Count[j],
        size = df_snp$Chr_Tot[j],
        prob = pred_freq,
        rho  = rho_val,
        log  = TRUE
      )
    }
    results$LL[i] <- ll
  }

  # Guardar Resultados <- segun chatgtp pero hay que modificar el nombre
  
  output_filename <- paste0("Analysis_", model_name,"_" ,task_id,"_SNP_", snp_actual, ".txt")
  output_path <- file.path(output_dir, output_filename)
  
  best_params <- results[which.max(results$LL),]
  print(paste("Guardando:", output_path, "Max LL:", best_params$LL))
  
  write.table(results, file=output_path, row.names=FALSE, quote = FALSE)
}
