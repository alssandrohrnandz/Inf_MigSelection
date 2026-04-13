library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

# === 1. LECTURA Y ETIQUETADO DE ARCHIVOS ===
directorio <- "data/results_Discrete/outputs_LL"
rutas_archivos <- list.files(path = directorio, pattern = "\\.txt$", full.names = TRUE)

# Función para procesar un solo archivo
procesar_archivo <- function(ruta) {
  nombre_archivo <- basename(ruta)
  
  # Leer el archivo (usando read_tsv por las tabulaciones)
  df <- read_tsv(ruta, show_col_types = FALSE)
  
  # Extraer metadatos del nombre del archivo
  modelo <- ifelse(str_detect(nombre_archivo, "TRON_ARES"), "TRON_ARES", "TRON_LEGACY")
  tipo <- ifelse(str_detect(nombre_archivo, "neutros"), "Neutro", "Seleccion")
  task_id <- str_extract(nombre_archivo, "(?<=TaskID_)\\d+")
  
  # Agregar metadatos como columnas nuevas
  df %>%
    mutate(
      Model = modelo,
      Type = tipo,
      TaskID = as.integer(task_id),
      File = nombre_archivo
    )
}

# Aplicar la función a todos los archivos y unirlos en un solo Master DF
df_raw <- map_dfr(rutas_archivos, procesar_archivo)

# === 2. TRANSFORMACIÓN, LIMPIEZA DE NAs Y COMPOSITE LIKELIHOOD ===

# 1. Transformamos a formato largo y limpiamos
df_long <- df_raw %>%
  pivot_longer(
    cols = starts_with("D_"), 
    names_to = "Grid_Param", 
    values_to = "Likelihood"
  ) %>%
  drop_na(Likelihood, SNP) %>% 
  mutate(
    Inferred_D = as.numeric(str_extract(Grid_Param, "(?<=D_)[0-9\\.]+")),
    Inferred_s = as.numeric(str_extract(Grid_Param, "(?<=s_)[-0-9\\.]+"))
  )

# 2. CALCULAR EL COMPOSITE LIKELIHOOD (Sumar todos los SNPs por cada punto del Grid)
df_composite <- df_long %>%
  # Agrupamos por Archivo (TaskID/Model/Type) y por el Punto del Grid (Inferred_D, Inferred_s)
  group_by(Model, Type, TaskID, Inferred_D, Inferred_s) %>%
  summarise(
    Composite_LL = sum(Likelihood),  # La suma de los log-likelihoods de todos los SNPs
    N_SNPs_Usados = n(),             # Informativo: cuántos SNPs válidos se sumaron
    .groups = "drop"
  )

# 3. EXTRAER EL MLE GLOBAL POR ARCHIVO
df_mle <- df_composite %>%
  # Ahora solo agrupamos por Archivo
  group_by(Model, Type, TaskID) %>%
  # Y extraemos LA fila donde la suma (Composite_LL) alcanzó su máximo absoluto
  slice_max(order_by = Composite_LL, n = 1, with_ties = FALSE) %>%
  ungroup()

mig_values <- c(0.0, 0.01, 0.05, 0.1)
sel_values <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
replicas_per_val <- 10

df_mle <- df_mle %>%
  mutate(
    # Replicamos la lógica: IDX = (TaskID - 1) / REPLICAS
    idx = floor((TaskID - 1) / replicas_per_val),
    
    # R indexa desde 1, así que sumamos 1 a los índices de Bash
    idx_mig = floor(idx / length(sel_values)) + 1,
    idx_sel = (idx %% length(sel_values)) + 1,
    
    # Asignamos los valores reales a nuevas columnas
    True_Mig = mig_values[idx_mig],
    True_Sel = sel_values[idx_sel]
  ) %>%
  # Limpiamos las columnas temporales que ya no necesitamos
  select(-idx, -idx_mig, -idx_sel)

df_mle <- df_mle %>%
  mutate(
    # 1. Escalar el valor inferido de la EDP al modelo Clásico (Wright-Fisher)
    Inferred_s_scaled = Inferred_s * 2,
    
    # 2. Calcular los errores usando el valor escalado
    Error_D = Inferred_D - True_Mig,
    Abs_Error_D = abs(Error_D),
    
    Error_s = Inferred_s_scaled - True_Sel,
    Abs_Error_s = abs(Error_s)
  )

# === 3. SEPARACIÓN EN LOS DF SOLICITADOS ===
# Separamos por Neutro y Selección
df_neutros <- df_mle %>% filter(Type == "Neutro")
df_seleccion <- df_mle %>% filter(Type == "Seleccion")

# Si además quieres separarlos por modelo en variables distintas:
df_neutros_ARES <- df_neutros %>% filter(Model == "TRON_ARES")
df_neutros_LEGACY <- df_neutros %>% filter(Model == "TRON_LEGACY")
df_seleccion_ARES <- df_seleccion %>% filter(Model == "TRON_ARES")
df_seleccion_LEGACY <- df_seleccion %>% filter(Model == "TRON_LEGACY")

# 2. Generar una tabla resumen (MAE) para el manuscrito
summary_table <- df_mle %>%
  group_by(Type, Model, True_Mig, True_Sel) %>%
  summarise(
    Mean_Inferred_D = mean(Inferred_D, na.rm = TRUE),
    MAE_D = mean(Abs_Error_D, na.rm = TRUE),
    Mean_Inferred_s = mean(Inferred_s, na.rm = TRUE),
    MAE_s = mean(Abs_Error_s, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_table, "Estimator_Performance_Summary.csv", row.names = FALSE)

print("¡Procesamiento finalizado! El Master DF de MLE está listo.")

plot_neutros_recovery <- ggplot(df_neutros, aes(x = True_Mig, y = Inferred_D, color = Model)) +
  
  geom_jitter(width = 0.003, height = 0, size = 2.5, alpha = 0.7) +
  
  # LÍNEA ROJA DE RECUPERACIÓN PERFECTA (y = x)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1.2) +
  
  facet_wrap(~ Model) +
  scale_x_continuous(breaks = mig_values) +
  scale_color_manual(values = c("TRON_ARES" = "#F8766D", "TRON_LEGACY" = "#00BFC4")) +
  labs(
    title = "Parameter Recovery: Diffusion (D) for Neutral Alleles",
    subtitle = "Red dashed line indicates perfect inference (y = x)",
    x = "Simulated Diffusion Rate (True D)",
    y = "Estimated Diffusion Rate (Inferred D)"
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "none", 
    panel.grid.minor = element_blank(), 
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold")
  )

ggsave("Publication_Recovery_Neutros_Facet.png", plot = plot_neutros_recovery, width = 12, height = 6, dpi = 300)

df_seleccion <- df_mle %>% filter(Type == "Seleccion")


# === 2. PARAMETER RECOVERY PLOT (CORREGIDO) ===
linea_perfecta <- data.frame(
  True_Sel = c(1e-4, 1e-1),
  Inferred_s_scaled = c(1e-4, 1e-1) # Usamos el nombre de la nueva columna
)

plot_seleccion_recovery <- ggplot(df_seleccion, aes(x = True_Sel, y = Inferred_s_scaled, color = Model)) +
  
  geom_point(size = 2.5, alpha = 0.6) +
  
  geom_line(data = linea_perfecta, aes(x = True_Sel, y = Inferred_s_scaled), 
            color = "red", linetype = "dashed", linewidth = 1.2, inherit.aes = FALSE) +
  
  facet_grid(Model ~ True_Mig, labeller = labeller(True_Mig = label_both)) +
  
  scale_x_log10(
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1),
    labels = scales::scientific
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-5, base = 10),
    breaks = sort(c(-1,-0.1,-0.01,-0.001,-0.0001, 0, 1e-4, 1e-3, 1e-2, 1e-1)),
    labels = scales::scientific
  ) +
  
  scale_color_manual(values = c("TRON_ARES" = "#F8766D", "TRON_LEGACY" = "#00BFC4")) +
  labs(
    title = "Parameter Recovery: Selection Coefficient (s)",
    subtitle = "Red dashed line indicates perfect inference (y = x). Columns show True Migration Rate (D).",
    x = "Simulated Selection Coefficient (True s) [log10 scale]",
    # Aclaramos en el eje Y que el valor fue escalado
    y = "Estimated Selection Coefficient (Scaled s) [pseudo-log scale]" 
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("Publication_Recovery_Selection_Scaled.png", plot = plot_seleccion_recovery, width = 14, height = 8, dpi = 300)


# Filtramos la tabla resumen para enfocarnos solo en los escenarios bajo selección
df_error_seleccion <- summary_table %>% filter(Type == "Seleccion")

plot_error <- ggplot(df_error_seleccion, aes(x = True_Sel, y = MAE_s, color = Model, group = Model)) +
  
  # Líneas y puntos para ver la tendencia del error
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  
  # Dividimos por tasa de migración simulada en 4 columnas
  facet_wrap(~ True_Mig, labeller = labeller(True_Mig = label_both), ncol = 4) +
  
  # Mantenemos la escala logarítmica en X para coherencia con el plot anterior
  scale_x_log10(
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1),
    labels = scales::scientific
  ) +
  
  # Colores consistentes con tus gráficos previos
  scale_color_manual(values = c("TRON_ARES" = "#F8766D", "TRON_LEGACY" = "#00BFC4")) +
  
  # Etiquetas en inglés listas para el manuscrito
  labs(
    title = "Estimator Performance: Mean Absolute Error (MAE) of Selection",
    subtitle = "Lower values indicate better parameter recovery. Panels show True Migration Rate (D).",
    x = "Simulated Selection Coefficient (True s) [log10 scale]",
    y = "Mean Absolute Error of Scaled s"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave("Publication_MAE_Selection_Performance.png", plot = plot_error, width = 12, height = 5, dpi = 300)