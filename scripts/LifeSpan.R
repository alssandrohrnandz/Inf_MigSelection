#!/usr/bin/env Rscript
# ============================================================
# master_lifespan_analysis.R
# Pipeline completo para análisis de Lifespan de SLiM
# 1. Lee resúmenes (summary) y genera gráfica de medias + Tukey HSD.
# 2. Lee mutaciones individuales (per_mutation) y genera gráfica de proporciones.
# Textos y etiquetas configurados en inglés para publicación científica.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(multcompView)
  library(scales)
})

# ============================================================
# 0. CONFIGURACIÓN GLOBAL
# ============================================================
dir_datos    <- "."  # Directorio con los CSVs procesados
dir_salida   <- "publication_figures"
if (!dir.exists(dir_salida)) dir.create(dir_salida, recursive = TRUE)

SEL_VALUES <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.0)
MIG_VALUES <- c(0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)
REPLICAS_PER_VAL <- 50

# Tema base para publicación (blanco y negro, texto formal)
pub_theme <- theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle    = element_text(size = 10, color = "grey30", hjust = 0),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.8),
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.x      = element_text(angle = 45, hjust = 1, color = "black", size = 9),
    axis.text.y      = element_text(color = "black"),
    axis.title       = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold")
  )

cat(sprintf("[%s] INICIANDO PROCESAMIENTO...\n", Sys.time()))

# ============================================================
# PARTE 1: ANÁLISIS DE MEDIAS Y TUKEY HSD (SUMMARY DATA)
# ============================================================
cat("\n--- PARTE 1: Average Lifespan & Statistical Analysis ---\n")

archivos_res <- list.files(dir_datos, pattern = "slim_lifespan_summary_chunk.*\\.csv$", full.names = TRUE)
if (length(archivos_res) == 0 && file.exists("slim_lifespan_summary.csv")) archivos_res <- "slim_lifespan_summary.csv"

df_sum <- rbindlist(lapply(archivos_res, fread), use.names = TRUE, fill = TRUE)

# Reconstruir s y m
df_sum[, idx     := (TaskID - 1) %/% REPLICAS_PER_VAL]
df_sum[, idx_sel := idx %% length(SEL_VALUES)]
df_sum[, idx_mig := idx %/% length(SEL_VALUES)]
df_sum[, s_val := SEL_VALUES[idx_sel + 1]]
df_sum[, m_val := MIG_VALUES[idx_mig + 1]]

# Formatear factor s (s=0 como "0" para evitar NA)
s_unicos <- sort(unique(df_sum$s_val))
s_labels_clean <- ifelse(s_unicos == 0, "0", gsub("e-0", "e-", format(s_unicos, scientific = TRUE)))
df_sum[, s_factor := factor(s_val, levels = s_unicos, labels = s_labels_clean)]
df_sum[, Mig_Label := factor(paste0("True m: ", m_val), levels = paste0("True m: ", sort(MIG_VALUES)))]

# --- Cálculo de Tukey HSD y CLD ---
calcular_cld <- function(sub_dt) {
  # Nombres temporales para evitar fallo por guiones en multcompLetters
  s_niveles <- levels(sub_dt$s_factor)
  s_map_temp <- setNames(paste0("s", seq_along(s_niveles)), s_niveles)
  
  sub_dt_temp <- copy(sub_dt)
  sub_dt_temp[, s_temp := factor(s_map_temp[as.character(s_factor)], levels = unname(s_map_temp))]
  
  modelo <- aov(Mean_Lifespan ~ s_temp, data = sub_dt_temp)
  tukey  <- TukeyHSD(modelo, "s_temp")
  p_vals <- tukey$s_temp[, "p adj"]
  
  cld_temp <- multcompLetters(p_vals)$Letters
  names(cld_temp) <- names(s_map_temp)[match(names(cld_temp), s_map_temp)]
  
  res_df <- sub_dt[, .(Max_Y = max(Mean_Lifespan, na.rm = TRUE) + (max(Mean_Lifespan, na.rm=TRUE)*0.08)), by = s_factor]
  res_df[, Letter := cld_temp[as.character(s_factor)]]
  return(res_df)
}

df_cld <- df_sum[, calcular_cld(.SD), by = Mig_Label]
df_means <- df_sum[, .(Grand_Mean = mean(Mean_Lifespan, na.rm = TRUE)), by = .(s_factor, Mig_Label)]

# --- Gráfico de Publicación 1 ---
p_stats <- ggplot(df_sum, aes(x = s_factor, y = Mean_Lifespan)) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 1.0, color = "black") +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black", linewidth = 0.5, width = 0.4) +
  geom_line(data = df_means, aes(x = s_factor, y = Grand_Mean, group = 1), color = "black", linewidth = 0.8) +
  geom_point(data = df_means, aes(x = s_factor, y = Grand_Mean), color = "black", size = 2.0) +
  geom_text(data = df_cld, aes(x = s_factor, y = Max_Y, label = Letter), 
            size = 3.2, fontface = "bold", color = "#B22222") + # Rojo oscuro para contraste
  facet_wrap(~ Mig_Label, ncol = 4) +
  labs(
    title    = "Average Allele Lifespan by Selection Coefficient",
    subtitle = "Faded points represent individual simulation replicates. Solid lines connect global means.\nLetters above boxes indicate Tukey HSD groups (p < 0.05).",
    x        = "Simulated Selection Coefficient (True s)",
    y        = "Mean Lifespan (Generations)"
  ) +
  pub_theme

out_fig1 <- file.path(dir_salida, "Fig1_AverageLifespan_Tukey.png")
ggsave(out_fig1, plot = p_stats, width = 12, height = 6.5, dpi = 300)
cat(" -> Gráfico 1 guardado en:", out_fig1, "\n")


# ============================================================
# PARTE 2: PROPORCIÓN DE ESTATUS (PER_MUTATION DATA)
# ============================================================
cat("\n--- PARTE 2: Mutation Status Proportions ---\n")

archivos_pm <- list.files(dir_datos, pattern = "slim_lifespan_per_mutation_chunk.*\\.csv$", full.names = TRUE)
if (length(archivos_pm) == 0 && file.exists("slim_lifespan_per_mutation.csv")) archivos_pm <- "slim_lifespan_per_mutation.csv"

if (length(archivos_pm) > 0) {
  df_pm <- rbindlist(lapply(archivos_pm, fread), use.names = TRUE, fill = TRUE)
  
  # Reconstruir s_val usando el mismo mecanismo numérico
  df_pm[, idx     := (TaskID - 1) %/% REPLICAS_PER_VAL]
  df_pm[, idx_sel := idx %% length(SEL_VALUES)]
  df_pm[, s_val   := SEL_VALUES[idx_sel + 1]]
  
  # Usar los mismos labels limpios del sumario
  df_pm[, s_factor := factor(s_val, levels = s_unicos, labels = s_labels_clean)]
  df_pm[, Status := factor(Status, levels = c("Fixed", "Lost", "Segregating"))]
  
  # Calcular proporciones
  df_prop <- df_pm[, .(N = .N), by = .(s_factor, Status)]
  df_prop[, Total := sum(N), by = s_factor]
  df_prop[, Proportion := N / Total]
  
  status_colors <- c("Fixed" = "#1b9e77", "Lost" = "#d95f02", "Segregating" = "#7570b3") # ColorBrewer Dark2
  
  # --- Gráfico de Publicación 2 ---
  p_prop <- ggplot(df_prop, aes(x = s_factor, y = Proportion, fill = Status)) +
    geom_bar(stat = "identity", position = "fill", width = 0.75, color = "black", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(Proportion > 0.04, sprintf("%.1f%%", Proportion * 100), "")),
      position = position_fill(vjust = 0.5),
      color = "white", fontface = "bold", size = 3
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = status_colors, name = "Mutation Status") +
    labs(
      title    = "Proportion of Mutation Outcomes by Selection Coefficient",
      x        = "Simulated Selection Coefficient (True s)",
      y        = "Relative Proportion"
    ) +
    pub_theme
  
  out_fig2 <- file.path(dir_salida, "Fig2_MutationStatus_Proportions.png")
  ggsave(out_fig2, plot = p_prop, width = 8, height = 5, dpi = 300)
  cat(" -> Gráfico 2 guardado en:", out_fig2, "\n")
  
} else {
  cat(" -> ADVERTENCIA: No se encontraron archivos per_mutation. Saltando Gráfico 2.\n")
}

cat(sprintf("\n[%s] PROCESO FINALIZADO CON ÉXITO.\n", Sys.time()))