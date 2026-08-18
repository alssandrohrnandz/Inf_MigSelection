library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(knitr) # Necesaria para la exportación a LaTeX

# Definimos la paleta de colores global y los nombres oficiales para publicación
mis_colores <- c("Classic" = "#00BFC4", "B-B" = "#F8766D", "BwS" = "#C77CFF")

# === 1. LECTURA Y ETIQUETADO DE ARCHIVOS ===
directorio <- "data/results_Discrete/outputs_LL"
rutas_archivos <- list.files(path = directorio, pattern = "\\.txt$", full.names = TRUE)

# Función para procesar un solo archivo
procesar_archivo <- function(ruta) {
  nombre_archivo <- basename(ruta)
  
  # Leer el archivo (usando read_tsv por las tabulaciones)
  df <- read_tsv(ruta, show_col_types = FALSE)
  
  # ¡LA SOLUCIÓN AQUÍ!: Forzar que la columna SNP sea siempre texto
  df <- df %>% mutate(SNP = as.character(SNP))
  
  # Extraer metadatos del nombre del archivo
  modelo <- case_when(
    str_detect(nombre_archivo, "B-B") | str_detect(nombre_archivo, "TRON_ARES") ~ "B-B",
    str_detect(nombre_archivo, "BwS") | str_detect(nombre_archivo, "TRON_SPIKES") ~ "BwS",
    TRUE ~ "Classic"
  )
  
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

# Forzamos a que Model sea un Factor ordenado para estandarizar todas las gráficas
df_raw <- df_raw %>%
  mutate(Model = factor(Model, levels = c("Classic", "B-B", "BwS")))

# === 2. TRANSFORMACIÓN, LIMPIEZA DE NAs Y COMPOSITE LIKELIHOOD ===

# 1. Transformamos a formato largo y limpiamos
df_long <- df_raw %>%
  # NUEVO: Filtro de seguridad explícito para aDNA
  filter(SNP != "INSUFFICIENT_DATA") %>% 
  pivot_longer(
    cols = starts_with("D_"), 
    names_to = "Grid_Param", 
    values_to = "Likelihood"
  ) %>%
  drop_na(Likelihood) %>% 
  mutate(
    Inferred_D = as.numeric(str_extract(Grid_Param, "(?<=D_)[0-9\\.]+")),
    Inferred_s = as.numeric(str_extract(Grid_Param, "(?<=s_)[-0-9\\.]+"))
  )

df_composite <- df_long %>%
  group_by(Model, Type, TaskID, Inferred_D, Inferred_s) %>%
  summarise(
    Composite_LL = sum(Likelihood),  
    N_SNPs_Usados = n(),             
    .groups = "drop"
  )

df_mle <- df_composite %>%
  group_by(Model, Type, TaskID) %>%
  slice_max(order_by = Composite_LL, n = 1, with_ties = FALSE) %>%
  ungroup()

mig_values <- c(0.0, 0.01, 0.05, 0.1)
sel_values <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
replicas_per_val <- 10

df_mle <- df_mle %>%
  mutate(
    idx = floor((TaskID - 1) / replicas_per_val),
    idx_mig = floor(idx / length(sel_values)) + 1,
    idx_sel = (idx %% length(sel_values)) + 1,
    True_Mig = mig_values[idx_mig],
    True_Sel = sel_values[idx_sel]
  ) %>%
  select(-idx, -idx_mig, -idx_sel)

df_mle <- df_mle %>%
  mutate(
    Inferred_s_scaled = Inferred_s * 2,
    Error_D = Inferred_D - True_Mig,
    Abs_Error_D = abs(Error_D),
    Error_s = Inferred_s_scaled - True_Sel,
    Abs_Error_s = abs(Error_s)
  )

# === 3. SEPARACIÓN EN LOS DF SOLICITADOS ===
df_neutros <- df_mle %>% filter(Type == "Neutro")
df_seleccion <- df_mle %>% filter(Type == "Seleccion")

# Generar una tabla resumen (MAE)
summary_table <- df_mle %>%
  group_by(Type, Model, True_Mig, True_Sel) %>%
  summarise(
    # NUEVO: Cuenta cuántas réplicas lograron sobrevivir a la máscara de aDNA
    Replicas_Exitosas = n(), 
    
    Mean_Inferred_D = mean(Inferred_D, na.rm = TRUE),
    MAE_D = mean(Abs_Error_D, na.rm = TRUE),
    Mean_Inferred_s = mean(Inferred_s, na.rm = TRUE),
    MAE_s = mean(Abs_Error_s, na.rm = TRUE),
    .groups = "drop"
  )

# Imprimir un resumen rápido en consola para que sepas cuántos datos perdiste
cat("\n=== REPORTE DE SUPERVIVENCIA aDNA ===\n")
print(summary_table %>% select(Type, Model, True_Sel, True_Mig, Replicas_Exitosas))
cat("=====================================\n\n")

write.csv(summary_table, "Estimator_Performance_Summary_3Models_aDNA.csv", row.names = FALSE)

print("¡Procesamiento finalizado! El Master DF de MLE está listo. Generando gráficos...")

# === 4. GRÁFICAS ===

# GRÁFICA 1: RECUPERACIÓN NEUTROS
plot_neutros_recovery <- ggplot(df_neutros, aes(x = True_Mig, y = Inferred_D, color = Model)) +
  geom_jitter(width = 0.003, height = 0, size = 2.5, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 1.2) +
  facet_wrap(~ Model) +
  scale_x_continuous(breaks = mig_values) +
  scale_color_manual(values = mis_colores) +
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

ggsave("Publication_Recovery_Neutros_3Models_aDNA.png", plot = plot_neutros_recovery, width = 14, height = 6, dpi = 300)

# GRÁFICA 2: RECUPERACIÓN SELECCIÓN
linea_perfecta <- data.frame(
  True_Sel = c(1e-4, 1e-1),
  Inferred_s_scaled = c(1e-4, 1e-1) 
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
  scale_color_manual(values = mis_colores) +
  labs(
    title = "Parameter Recovery: Selection Coefficient (s)",
    subtitle = "Red dashed line indicates perfect inference (y = x). Columns show True Migration Rate (D).",
    x = "Simulated Selection Coefficient (True s) [log10 scale]",
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

ggsave("Publication_Recovery_Selection_3Models_aDNA.png", plot = plot_seleccion_recovery, width = 14, height = 10, dpi = 300)

# GRÁFICA 3: DESEMPEÑO DE ERROR (MAE)
df_error_seleccion <- summary_table %>% filter(Type == "Seleccion") %>%
  mutate(
    True_Sel_Plot = case_when(
      Model == "B-B"     ~ True_Sel * 0.97,  
      Model == "Classic" ~ True_Sel,         
      Model == "BwS"     ~ True_Sel * 1.05   
    )
  )

plot_error <- ggplot(df_error_seleccion, aes(x = True_Sel_Plot, y = MAE_s, color = Model, group = Model, shape = Model)) +
  geom_line(linewidth = 1.2, alpha = 0.6) +
  geom_point(size = 3.5, alpha = 0.8) +
  facet_wrap(~ True_Mig, labeller = labeller(True_Mig = label_both), ncol = 4) +
  scale_x_log10(
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1),
    labels = scales::scientific
  ) +
  scale_color_manual(values = mis_colores) +
  scale_shape_manual(values = c("Classic" = 17, "B-B" = 16, "BwS" = 15)) +
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

ggsave("Publication_MAE_Selection_3Models_Fixed_aDNA.png", plot = plot_error, width = 14, height = 6, dpi = 300)

# GRÁFICA 4: DISTRIBUCIÓN GRID (CAJAS Y BIGOTES)
df_box_grid <- df_seleccion %>%
  mutate(
    True_Sel_Factor = factor(True_Sel, levels = sort(unique(True_Sel))),
    Model = factor(Model, levels = c("Classic", "B-B", "BwS"))
  )

plot_distribucion_grid <- ggplot(df_box_grid, aes(x = True_Sel_Factor, y = Inferred_s_scaled, fill = Model)) +
  geom_boxplot(
    width = 0.2, 
    alpha = 0.7,
    color = "black",
    outlier.size = 0.8,
    outlier.alpha = 0.4
  ) +
  geom_errorbar(aes(y = True_Sel, ymin = True_Sel, ymax = True_Sel), 
                color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_grid(Model ~ True_Mig, labeller = labeller(True_Mig = label_both)) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-5, base = 10),
    breaks = sort(c(-1,-0.1,-0.01,-0.001,-0.0001, 0, 1e-4, 1e-3, 1e-2, 1e-1)),
    labels = scales::scientific
  ) +
  scale_fill_manual(values = mis_colores) +
  labs(
    title = "Inference Precision and Bias across Models",
    subtitle = "Rows: Likelihood Models | Columns: Migration Rates. Red dashed segments indicate True Selection values.",
    x = "Simulated Selection Coefficient (True s)",
    y = "Estimated Selection Coefficient (Scaled s) [pseudo-log scale]"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none", 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave("Publication_Boxplot_Grid_3x4_aDNA.png", plot = plot_distribucion_grid, width = 15, height = 11, dpi = 300)

# === 5. GENERACIÓN DE TABLA LATEX PARA OVERLEAF ===
tabla_resumen <- df_raw %>%
  filter(Model == "Classic") %>% 
  select(Type, TaskID, SNP) %>%
  distinct() %>% 
  mutate(
    idx = floor((TaskID - 1) / replicas_per_val),
    idx_mig = floor(idx / length(sel_values)) + 1,
    idx_sel = (idx %% length(sel_values)) + 1,
    Migracion_D = mig_values[idx_mig],
    Seleccion_s = sel_values[idx_sel]
  ) %>%
  group_by(Migracion_D, Seleccion_s, Type) %>%
  summarise(Total_Mutaciones = n(), .groups = 'drop') %>%
  pivot_wider(
    names_from = Type, 
    values_from = Total_Mutaciones,
    values_fill = 0 
  ) %>%
  select(Migracion_D, Seleccion_s, Neutro, Seleccion) %>%
  arrange(Migracion_D, desc(Seleccion_s))

cat("\n=== COPIA EL SIGUIENTE CÓDIGO A OVERLEAF ===\n")
kable(
  tabla_resumen, 
  format = "latex", 
  booktabs = TRUE, 
  linesep = "",    
  col.names = c("Migración ($D$)", "Selección ($s$)", "Alelos Neutros", "Alelos Seleccionados"),
  caption = "Número total de mutaciones filtradas por escenario de simulación.",
  label = "resumen_mutaciones"
)

#======= Plot_aDNA_Timeline

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

# === 1. RUTAS DE LOS DATOS ===
dir_ll <- "data/results_Discrete/outputs_LL"
dir_masks <- "data/results_Discrete/aDNA_Masks"

# === 2. FILTRAR TASK_IDs EXITOSOS ===
archivos_ll <- list.files(dir_ll, pattern = "\\.txt$", full.names = TRUE)

obtener_taskid_valido <- function(ruta) {
  df_test <- suppressMessages(read_tsv(ruta, col_select = 1, n_max = 1, show_col_types = FALSE))
  if(nrow(df_test) > 0 && df_test$SNP[1] != "INSUFFICIENT_DATA") {
    task_id <- as.integer(str_extract(basename(ruta), "(?<=TaskID_)\\d+"))
    return(task_id)
  }
  return(NA)
}

cat("Analizando archivos de Likelihood para encontrar réplicas exitosas...\n")
task_ids_validos <- unique(na.omit(sapply(archivos_ll, obtener_taskid_valido)))
n_taskids <- length(task_ids_validos)
cat("Se encontraron", n_taskids, "TaskIDs con inferencias exitosas.\n")

# === 3. LEER MÁSCARAS ARQUEOLÓGICAS VALIDAS ===
leer_mascara <- function(tid) {
  ruta_mask <- file.path(dir_masks, paste0("aDNA_Global_Mask_", tid, ".csv"))
  if(file.exists(ruta_mask)) {
    df <- suppressMessages(read_csv(ruta_mask, show_col_types = FALSE))
    df$TaskID <- tid
    return(df)
  }
  return(NULL)
}

df_masks <- bind_rows(lapply(task_ids_validos, leer_mascara))

gen_max <- max(df_masks$Generation, na.rm = TRUE)
gen_min <- min(df_masks$Generation, na.rm = TRUE)
inflexion_neolitico <- gen_min + ((gen_max - gen_min) * 0.6)

# === 4. CÁLCULO DE PROMEDIOS POR VENTANA DE TIEMPO ===
bin_w <- 15 # Ventanas de 15 generaciones

df_promedios <- df_masks %>%
  # Asignamos cada registro a su "ventana" (bin) correspondiente
  mutate(Gen_Bin = floor(Generation / bin_w) * bin_w + (bin_w / 2)) %>%
  group_by(Gen_Bin) %>%
  summarise(
    Suma_Total = sum(Chr_Tot_Fijo),
    # Dividimos la suma total entre el número de réplicas para sacar el promedio esperado por corrida
    Promedio_por_Replica = sum(Chr_Tot_Fijo) / n_taskids, 
    .groups = "drop"
  )

# === 5. GENERAR LA LÍNEA DE TIEMPO (PLOT) ===
plot_timeline <- ggplot() +
  
  # FONDOS TEMÁTICOS (Paleolítico vs Neolítico)
  annotate("rect", xmin = gen_min, xmax = inflexion_neolitico, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "blue") +
  annotate("rect", xmin = inflexion_neolitico, xmax = gen_max, ymin = 0, ymax = Inf, 
           alpha = 0.1, fill = "green") +
  
  # Textos de los periodos
  annotate("text", x = gen_min + (inflexion_neolitico - gen_min)/2, y = Inf, 
           label = "Hunter-Gatherers\n(Sparse Data)", vjust = 1.5, fontface = "italic", color = "#2C3E50") +
  annotate("text", x = inflexion_neolitico + (gen_max - inflexion_neolitico)/2, y = Inf, 
           label = "Neolithic & Bronze Age\n(Data Boom)", vjust = 1.5, fontface = "italic", color = "#2C6E2F") +
  
  geom_vline(xintercept = inflexion_neolitico, linetype = "dashed", color = "black", linewidth = 0.8) +

  # BARRAS: Ahora muestran el Promedio en lugar de la Suma Bruta
  geom_col(data = df_promedios, aes(x = Gen_Bin, y = Promedio_por_Replica), 
           fill = "#2C3E50", color = "white", alpha = 0.7, width = bin_w) +
  
  # LÍNEA Y PUNTOS: Conecta los promedios directamente (reemplaza la línea roja de densidad)
  #geom_line(data = df_promedios, aes(x = Gen_Bin, y = Promedio_por_Replica), 
            #color = "#E74C3C", linewidth = 1.2) +
  #geom_point(data = df_promedios, aes(x = Gen_Bin, y = Promedio_por_Replica), 
             #color = "#C0392B", size = 2.5) +
  
  scale_x_continuous(breaks = seq(0, 500, by = 50), limits = c(0, 500)) +
  
  # ¡LA SOLUCIÓN!: Añade un 25% de espacio extra solo en la parte superior del eje Y
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  
  labs(
    title = "Ancient DNA Sampling Timeline: Expected Yield per Simulation",
    subtitle = "Bars and line represent the average number of successfully recovered chromosomes per simulation run.",
    x = "Generation (Past \u2192 Present)",
    y = "Average Chromosomes Recovered per Replica"
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 12)),
    axis.title.y = element_text(margin = margin(r = 12))
  )

ggsave("Publication_aDNA_Sampling_Timeline_Average.png", plot = plot_timeline, width = 12, height = 6, dpi = 300)

cat("Gráfica generada exitosamente: Publication_aDNA_Sampling_Timeline_Average.png\n")