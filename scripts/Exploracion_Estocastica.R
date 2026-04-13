#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(tidyr)

# === 1. CAPTURA DE ARGUMENTOS ===
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 6) stop("Faltan argumentos. Se requieren 7.")

freq_file   <- args[1]
subset_file <- args[2]
task_id     <- args[3]
prefijo     <- args[4]
output_dir  <- args[5]
true_s      <- as.numeric(args[6])

#archivos_esperados <- c(freq_file, subset_file, task_id, prefijo, output_dir, true_s)
#num_esperados <- length(argumentos_esperados)

# === 2. FUNCIONES TEÓRICAS ===
# Fórmula clásica para simular la trayectoria esperada determinista
simulate_freq <- function(p0, s, h, generations) {
  p <- numeric(generations)
  p[1] <- p0
  for (t in 2:generations) {
    p_prev <- p[t-1] 
    numerador <- (p_prev^2 * (1 + s)) + (p_prev * (1 - p_prev) * (1 + s * h))
    denominador <- (p_prev^2 * (1 + s)) + (2 * p_prev * (1 - p_prev) * (1 + s * h)) + ((1 - p_prev)^2)
    p[t] <- numerador / denominador
  }
  return(p)
}

# === 3. LECTURA Y LIMPIEZA DE DATOS ===
cat(sprintf("\n======================================================\n"))
cat(sprintf(" INICIANDO ANÁLISIS DE DERIVA Y SELECCIÓN (Tarea %s)\n", task_id))
cat(sprintf(" Selección Simulada Esperada (True s): %f\n", true_s))
cat(sprintf("======================================================\n"))

df_raw <- read.csv(freq_file, header=TRUE, stringsAsFactors=FALSE)

stats_alelos <- df_raw %>%
  group_by(MutationID) %>%
  summarise(
    Gen_Origen = min(Generation),
    Gen_Fin = max(Generation),
    Lifespan = max(Generation) - min(Generation) + 1,
    Max_Freq = max(Frequency),
    Fixed = ifelse(max(Frequency) >= 0.99, TRUE, FALSE),
    .groups = "drop"
  )

total_muts <- nrow(stats_alelos)

# === 4. ESTADÍSTICAS DE SUPERVIVENCIA Y EXTINCIÓN ===
extintos <- stats_alelos %>% filter(Fixed == FALSE)

singletons <- sum(extintos$Lifespan == 1)
perdidos_5 <- sum(extintos$Lifespan < 5)
perdidos_50 <- sum(extintos$Lifespan < 50)
perdidos_100 <- sum(extintos$Lifespan < 100)
perdidos_150 <- sum(extintos$Lifespan < 150)
fijados <- sum(stats_alelos$Fixed == TRUE)

cat(sprintf("\n--- DEMOGRAFÍA DE MUTACIONES ---\n"))
cat(sprintf("Total de mutaciones nacidas: %d\n", total_muts))
cat(sprintf("Singletons (Duraron 1 gen): %d (%.2f%%)\n", singletons, (singletons/total_muts)*100))
cat(sprintf("Extintos en < 5 gens: %d (%.2f%%)\n", perdidos_5, (perdidos_5/total_muts)*100))
cat(sprintf("Extintos en < 50 gens: %d (%.2f%%)\n", perdidos_50, (perdidos_50/total_muts)*100))
cat(sprintf("Extintos en < 100 gens: %d (%.2f%%)\n", perdidos_100, (perdidos_100/total_muts)*100))
cat(sprintf("Extintos en < 150 gens: %d (%.2f%%)\n", perdidos_150, (perdidos_150/total_muts)*100))
cat(sprintf("Alelos FIJADOS (>= 99%%): %d (%.2f%%)\n", fijados, (fijados/total_muts)*100))

# === 5. INFERENCIA DEL COEFICIENTE DE SELECCIÓN (s) ===
candidatos <- stats_alelos %>% filter(Max_Freq > 0.01)
alelos_inferidos <- data.frame(MutationID = integer(), s_inferido = numeric(), error_abs = numeric())

cat(sprintf("\n--- INFERENCIA DE SELECCIÓN ---\n"))
cat(sprintf("Alelos viables para inferencia (alcanzaron >1%% frec): %d\n", nrow(candidatos)))

for(mut_id in candidatos$MutationID) {
  df_mut <- df_raw %>% 
    filter(MutationID == mut_id, Frequency > 0.01, Frequency < 0.9)
  
  if(nrow(df_mut) >= 4) { 
    df_mut$Logit <- log(df_mut$Frequency / (1 - df_mut$Frequency))
    modelo <- lm(Logit ~ Generation, data = df_mut)
    pendiente <- coef(modelo)["Generation"]
    s_inf <- as.numeric(2 * pendiente) 
    
    alelos_inferidos <- rbind(alelos_inferidos, data.frame(
      MutationID = mut_id,
      s_inferido = s_inf,
      error_abs = abs(s_inf - true_s)
    ))
  }
}

if(nrow(alelos_inferidos) > 0) {
  s_promedio <- mean(alelos_inferidos$s_inferido, na.rm = TRUE)
  cat(sprintf("Promedio de 's' inferido en supervivientes: %f\n", s_promedio))
  cat(sprintf("Diferencia contra 's' real: %f\n", s_promedio - true_s))
  
  # === 6. EXTRACCIÓN DEL TOP 10 Y GRAFICACIÓN ===
  top10 <- alelos_inferidos %>% 
    arrange(error_abs) %>% 
    head(10)
  
  cat(sprintf("\nGraficando el Top %d de alelos más cercanos a la teoría...\n", nrow(top10)))
  
  df_plot <- df_raw %>% 
    filter(MutationID %in% top10$MutationID) %>%
    group_by(MutationID) %>%
    mutate(Gen_Relativa = Generation - min(Generation) + 1) %>% 
    ungroup()
  
  # --- CÁLCULO DE LA CURVA TEÓRICA ---
  max_gen_plot <- max(df_plot$Gen_Relativa)
  p0_val <- 1 / 20000 # 1 mutante en una población diploide de N=1000
  
  df_teorico <- data.frame(
    Gen_Relativa = 1:max_gen_plot,
    Frequency = simulate_freq(p0 = p0_val, s = true_s, h = 0.5, generations = max_gen_plot)
  )
  
  df_plot$MutationID <- as.factor(df_plot$MutationID)
  
  # --- GRÁFICA COMBINADA ---
  p <- ggplot() +
    # Líneas empíricas de SLiM (Top 10)
    geom_line(data = df_plot, aes(x = Gen_Relativa, y = Frequency, color = MutationID), linewidth = 1.2, alpha = 0.7) +
    # Línea Teórica Maestra (Negra, gruesa, punteada)
    geom_line(data = df_teorico, aes(x = Gen_Relativa, y = Frequency), color = "black", linewidth = 2, linetype = "dashed") +
    labs(
      title = paste("Trayectorias Top 10 vs Modelo Determinista (s =", true_s, ")"),
      subtitle = paste("Línea punteada: Expectativa teórica | s Inferido Promedio:", round(mean(top10$s_inferido), 4)),
      x = "Generaciones Relativas (Desde el origen)",
      y = "Frecuencia Alélica"
    ) +
    theme_bw(base_size = 18) +
    theme(legend.position = "right")
  
  plot_filename <- file.path(output_dir, paste0("Top10_Trayectorias_Task_", task_id, "_Sel_", true_s, ".png"))
  ggsave(plot_filename, plot = p, width = 12, height = 8, dpi = 300)
  cat(sprintf("Gráfica guardada en: %s\n", plot_filename))
  
} else {
  cat(sprintf("ADVERTENCIA: Ningún alelo sobrevivió lo suficiente para inferir 's'.\n"))
}
cat(sprintf("======================================================\n\n"))