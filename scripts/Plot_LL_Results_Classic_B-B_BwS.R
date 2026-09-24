library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

# Definimos la paleta de colores global y los nuevos nombres
mis_colores <- c("Classic" = "#004488", "B-B" = "#BB5566", "BwS" = "#DDAA33")
# === 1. LECTURA Y ETIQUETADO DE ARCHIVOS ===
directorio <- "data/results_Discrete/outputs_LL/independent_loci/neutros"
# * Este directorio contiene los archivos de salida de LL usando sitios independientes (loci) para los 3 modelos: Classic, B-B y BwS.
# * Descomentar para usar otro directorio si es necesario
# directorio <- "data/results_Discrete/outputs_LL/independent_loci/neutros"
rutas_archivos <- list.files(path = directorio, pattern = "\\.txt$", full.names = TRUE)

# Excluir los que empiezan con "P_TRON" o "P_"
rutas_archivos <- rutas_archivos[
  !grepl("^P_TRON|^P_", basename(rutas_archivos))
]

# Función para procesar un solo archivo
procesar_archivo <- function(ruta) {
  nombre_archivo <- basename(ruta)
  
  # Leer el archivo (usando read_tsv por las tabulaciones)
  df <- read_tsv(ruta, show_col_types = FALSE)
  
  # Extraer metadatos del nombre del archivo adaptado para 3 modelos
  modelo <- case_when(
    str_detect(nombre_archivo, "TRON_ARES") ~ "B-B",
    str_detect(nombre_archivo, "TRON_SPIKES") ~ "BwS",
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
# Retener solo SNPs que no tengan ningún NA en las columnas de la grilla
df_clean <- df_raw %>%
  filter(if_all(starts_with("D_"), ~ !is.na(.) & !is.infinite(.)))

df_clean <- df_clean %>%
  mutate(Model = factor(Model, levels = c("Classic", "B-B", "BwS")))
# === 2. TRANSFORMACIÓN, LIMPIEZA DE NAs Y COMPOSITE LIKELIHOOD ===

# 1. Transformamos a formato largo y limpiamos
df_long <- df_raw %>%
  pivot_longer(
    cols = starts_with("D_"), 
    names_to = "Grid_Param", 
    values_to = "Likelihood"
  ) %>%
  drop_na(Likelihood, SNP) %>% 
  separate_wider_delim(
    cols = Grid_Param,
    delim = "_",
    names = c(NA, "Inferred_D", NA, "Inferred_s"),
    cols_remove = FALSE
  ) %>% 
  mutate(
    Inferred_D = as.numeric(Inferred_D),
    Inferred_s = as.numeric(Inferred_s)
  )

# 2. CALCULAR EL COMPOSITE LIKELIHOOD (Sumar todos los SNPs por cada punto del Grid)
df_composite <- df_long %>%
  # Agrupamos por Archivo (TaskID/Model/Type) y por el Punto del Grid
  group_by(Model, Type, TaskID, Inferred_D, Inferred_s) %>%
  summarise(
    Composite_LL = sum(Likelihood),  
    N_SNPs_Usados = n(),             
    .groups = "drop"
  )

# 3. EXTRAER EL MLE GLOBAL POR ARCHIVO
df_mle <- df_composite %>%
  group_by(Model, Type, TaskID) %>%
  slice_max(order_by = Composite_LL, n = 1, with_ties = FALSE) %>%
  ungroup()

# * Estos valores de migracion y seleccion fueron definidos para sitios dependientes, para usar los valores
# * de simulación de los sitios independientes, se deben ajustar según el diseño experimental
mig_values <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125)
sel_values <- c(0.00)

replicas_per_val <- 50

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
    # Escalar el valor inferido de la EDP al modelo Clásico (Wright-Fisher)
    Inferred_s_scaled = Inferred_s * 2,
    
    # Calcular los errores
    Error_D = Inferred_D - True_Mig,
    Abs_Error_D = abs(Error_D),
    
    Error_s = Inferred_s_scaled - True_Sel,
    Abs_Error_s = abs(Error_s)
  )

# === 3. SEPARACIÓN EN LOS DF SOLICITADOS ===
# Separamos por Neutro y Selección
df_neutros <- df_mle %>% filter(Type == "Neutro")
df_seleccion <- df_mle %>% filter(Type == "Seleccion")

# Separación individual por si se necesita para inspección
df_neutros_ARES <- df_neutros %>% filter(Model == "TRON_ARES")
df_neutros_LEGACY <- df_neutros %>% filter(Model == "TRON_LEGACY")
df_neutros_SPIKES <- df_neutros %>% filter(Model == "TRON_SPIKES")

df_seleccion_ARES <- df_seleccion %>% filter(Model == "TRON_ARES")
df_seleccion_LEGACY <- df_seleccion %>% filter(Model == "TRON_LEGACY")
df_seleccion_SPIKES <- df_seleccion %>% filter(Model == "TRON_SPIKES")

# Generar una tabla resumen (MAE)
summary_table <- df_mle %>%
  group_by(Type, Model, True_Mig, True_Sel) %>%
  summarise(
    Mean_Inferred_D = mean(Inferred_D, na.rm = TRUE),
    MAE_D = mean(Abs_Error_D, na.rm = TRUE),
    Mean_Inferred_s = mean(Inferred_s, na.rm = TRUE),
    MAE_s = mean(Abs_Error_s, na.rm = TRUE),
    .groups = "drop"
  )

# Guardamos el resumen con un nombre actualizado
write.csv(summary_table, "Estimator_Performance_Summary_3Models_neutros.csv", row.names = FALSE)

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

# Anchura aumentada a 14 para acomodar 3 paneles
ggsave("Publication_Recovery_Neutros_p1.png", plot = plot_neutros_recovery, width = 14, height = 6, dpi = 300)

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

# Altura aumentada a 10 para acomodar las 3 filas
ggsave("Publication_Recovery_Selection_3Models.png", plot = plot_seleccion_recovery, width = 14, height = 10, dpi = 300)

# GRÁFICA 3: DESEMPEÑO DE ERROR (MAE)
# Filtramos la tabla resumen para enfocarnos solo en los escenarios bajo selección
df_error_seleccion <- summary_table %>% filter(Type == "Seleccion") %>%
  mutate(
    # Aplicamos el dodge usando los nuevos nombres
    True_Sel_Plot = case_when(
      Model == "B-B"     ~ True_Sel * 0.97,  # Desplaza ligeramente a la izquierda
      Model == "Classic" ~ True_Sel,         # Se mantiene exactamente en el centro
      Model == "BwS"     ~ True_Sel * 1.05   # Desplaza ligeramente a la derecha
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
  
  # Asignamos las formas a los nuevos nombres
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

# Guardamos con un nombre nuevo para que puedas comparar
ggsave("Publication_MAE_Selection_3Models_Fixed.png", plot = plot_error, width = 14, height = 6, dpi = 300)

# 1. Aseguramos que los datos de selección tengan los factores correctos
df_box_grid <- df_seleccion %>%
  mutate(
    # Convertimos True_Sel a factor para que las cajas no se encimen en el eje X
    True_Sel_Factor = factor(True_Sel, levels = sort(unique(True_Sel))),
    # Aseguramos que el orden de los modelos sea el deseado
    Model = factor(Model, levels = c("Classic", "B-B", "BwS"))
  )

# 2. Generamos el gráfico tipo Rejilla (Facet Grid)
plot_distribucion_grid <- ggplot(df_box_grid, aes(x = True_Sel_Factor, y = Inferred_s_scaled, fill = Model)) +
  
  # Violín de fondo (densidad de la inferencia)
  #geom_violin(alpha = 0.3, color = "transparent", scale = "width") +
  
  # Caja y bigotes (cuartiles y mediana)
  geom_boxplot(
    width = 0.2, 
    alpha = 0.7,
    color = "black",
    outlier.size = 0.8,
    outlier.alpha = 0.4
  ) +
  
  # LÍNEA DE REFERENCIA: El valor teórico (y = x)
  # Como X es un factor, esto requiere un truco: dibujar segmentos locales
  geom_errorbar(aes(y = True_Sel, ymin = True_Sel, ymax = True_Sel), 
                color = "red", linetype = "dashed", linewidth = 0.8) +

  # ESTRUCTURA DE REJILLA: Filas = Modelos, Columnas = Migración
  facet_grid(Model ~ True_Mig, labeller = labeller(True_Mig = label_both)) +
  
  # Escala Y en pseudo-log para capturar desde 0 hasta 0.1
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-5, base = 10),
    breaks = sort(c(-1,-0.1,-0.01,-0.001,-0.0001, 0, 1e-4, 1e-3, 1e-2, 1e-1)),
    labels = scales::scientific
  ) +
  
  # Colores consistentes
  scale_fill_manual(values = mis_colores) +
  
  labs(
    title = "Inference Precision and Bias across Models",
    subtitle = "Rows: Likelihood Models | Columns: Migration Rates. Red dashed segments indicate True Selection values.",
    x = "Simulated Selection Coefficient (True s)",
    y = "Estimated Selection Coefficient (Scaled s) [pseudo-log scale]"
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none", # Los nombres ya están en los títulos de las filas
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 11)
  )

# Guardamos con dimensiones optimizadas para el formato 3x4
ggsave("Publication_Boxplot_Grid_3x4.png", 
       plot = plot_distribucion_grid, 
       width = 15, 
       height = 11, 
       dpi = 300)


tabla_resumen <- df_raw %>%
  # Nos quedamos solo con un modelo (ej. TRON_LEGACY) para no contar los SNPs 3 veces
  filter(Model == "Classic") %>% 
  select(Type, TaskID, SNP) %>%
  distinct() %>% # Aseguramos no contar duplicados
  mutate(
    # Mapeamos el TaskID a los valores teóricos
    idx = floor((TaskID - 1) / replicas_per_val),
    idx_mig = floor(idx / length(sel_values)) + 1,
    idx_sel = (idx %% length(sel_values)) + 1,
    Migracion_D = mig_values[idx_mig],
    Seleccion_s = sel_values[idx_sel]
  ) %>%
  # Agrupamos por los parámetros de interés y contamos los SNPs
  group_by(Migracion_D, Seleccion_s, Type) %>%
  summarise(Total_Mutaciones = n(), .groups = 'drop') %>%
  
  # Pivotamos la tabla para que "Neutro" y "Seleccion" sean columnas separadas
  pivot_wider(
    names_from = Type, 
    values_from = Total_Mutaciones,
    values_fill = 0 # Rellena con 0 si en algún escenario no hubo mutaciones
  ) %>%
  
  # Ordenamos las columnas y las filas para que se vea limpio
  select(Migracion_D, Seleccion_s, Neutro, Seleccion) %>%
  arrange(Migracion_D, desc(Seleccion_s))

  kable(
  tabla_resumen, 
  format = "latex", 
  booktabs = TRUE, # Usa el formato formal de libros científicos
  linesep = "",    # Evita saltos de línea extraños
  col.names = c("Migración ($D$)", "Selección ($s$)", "Alelos Neutros", "Alelos Seleccionados"),
  caption = "Número total de mutaciones filtradas por escenario de simulación.",
  label = "resumen_mutaciones"
)

# 4. MAE D

df_error_dispersion <- summary_table %>% filter(Type == "Neutro") %>%
  mutate(
    # Aplicamos el dodge usando los nuevos nombres
    True_M_Plot = case_when(
      Model == "B-B"     ~ True_Mig * 0.97,  # Desplaza ligeramente a la izquierda
      Model == "Classic" ~ True_Mig,         # Se mantiene exactamente en el centro
      Model == "BwS"     ~ True_Mig * 1.05   # Desplaza ligeramente a la derecha
    )
  )

plot_error <- ggplot(df_error_dispersion, aes(x = True_M_Plot, y = MAE_D, color = Model, group = Model, shape = Model)) +
  
  geom_line(linewidth = 1.2, alpha = 0.6) +
  geom_point(size = 3.5, alpha = 0.8) +
  
  geom_vline(xintercept = 0.001, linetype = "dashed", color = "red", linewidth = 1.2) +
  

  scale_x_log10(
    breaks = mig_values,
    labels = scales::scientific
  ) +
  

  scale_color_manual(values = mis_colores) +
  
  # Asignamos las formas a los nuevos nombres
  scale_shape_manual(values = c("Classic" = 17, "B-B" = 16, "BwS" = 15)) +
  
  labs(
    title = "Estimator Performance: Mean Absolute Error (MAE) of Dispersion",
    subtitle = "Lower values indicate better parameter recovery. Red Dashed Line indicates Nm=1 individual per generation.",
    x = "Simulated Migration Rate (m)",
    y = "Mean Absolute Error of inferred D (in units of D)"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12)
  )

# Guardamos con un nombre nuevo para que puedas comparar
ggsave("Publication_MAE_Dispersion_3Models_Fixed.png", plot = plot_error, width = 14, height = 6, dpi = 300)

#########################################################
###  NUEVAS GRÁFICAS PARA PUBLICACIÓN ALELOS NEUTROS  ###
#########################################################

library(scales) # Necesario para pseudo_log_trans

# TODO: Esto se modifico para la NO INDEPENDENCIA de loci. Modificar en su defecto si se quiere usar para loci independientes.

# GRÁFICA 1: RECUPERACIÓN NEUTROS (MEJORADA)
plot_neutros_recovery <- ggplot(df_neutros, aes(x = True_Mig, y = Inferred_D, color = Model)) +
  # Reducimos un poco el jitter para no distorsionar en la escala logarítmica
  geom_jitter(width = 0, height = 0.05, size = 2.5, alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
  facet_wrap(~ Model) +
  # Aplicamos transformación pseudo-logarítmica para ver bien los valores diminutos y el cero
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1e-5, base = 10),
    breaks = c(0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125)
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-5, base = 10),
    breaks = c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125)
  ) +
  scale_color_manual(values = mis_colores) + # Asegúrate de tener 'mis_colores' definido
  labs(
    title = "Parameter Recovery: Diffusion (D) for Neutral Alleles (no independence)",
    subtitle = "Dashed line = perfect inference. Axes in pseudo-log scale.",
    x = "Simulated Diffusion Rate (True D)",
    y = "Estimated Diffusion Rate (Inferred D)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none", 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1), # Rota el texto para que no se encime
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold")
  )

# * Cuando usemos "Indepence_Loci" data, debemos usar la siguiente carpeta
ggsave("figures/independent_loci/Recovery_3Models_Neutral_Independent_Loci.png", plot = plot_neutros_recovery, width = 14, height = 6, dpi = 300)

#ggsave("figures/loci/Recovery_3Models_Neutral_No_Independence_Loci.png", plot = plot_neutros_recovery, width = 14, height = 6, dpi = 300)


###  Violin Neutros

# GRÁFICA 2: DISTRIBUCIÓN DE ERRORES (VIOLIN + BOXPLOT)
plot_neutros_violin <- ggplot(df_seleccion, aes(x = as.factor(True_Mig), y = Inferred_D, fill = Model)) +
  # Dibuja el violín para la densidad de los datos
  geom_violin(alpha = 0.5, color = "black", scale = "width", trim = FALSE) +
  # Agrega un boxplot pequeñito adentro para ver medianas y cuartiles
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
  # Dibuja un marcador rojo exactamente en el valor real de simulación
  geom_point(aes(y = True_Mig), color = "red", shape = 4, size = 3, stroke = 1.5) +
  # Si 0.001 es la categoría número 8 de tu eje X:
  geom_vline(xintercept = 7, linetype = "dashed", color = "red", linewidth = 1.2) +
  facet_wrap(~ Model) +
  # Escala pseudo-log en Y para no perder resolución en los errores pequeños
  # Escala pseudo-log en Y con etiquetas limpias
  # Reemplaza tu escala Y actual por esta:
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-6, base = 10),
    breaks = c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125), # Quitamos el 0.125
    labels = c("0", "0.0001", "0.0005", "0.001", "0.005", "0.01", "0.025", "0.05", "0.075", "0.1", "0.125")
  ) +
  scale_fill_manual(values = mis_colores) +

  # TODO: Esto se modifico para la NO INDEPENDENCIA de loci. Modificar en su defecto si se quiere usar para loci independientes.
  
  labs(
    title = "Dispersion of Estimated Diffusion Rates (MLE) for Neutral Alleles (no independence)",
    subtitle = "Red cross (X) indicates the True Simulated Diffusion (Target)",
    x = "Simulated Diffusion Rate (True D) - Categorical",
    y = "Estimated Diffusion Rate (Inferred D)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold")
  )


# ggsave("figures/independent_loci/Publication_Dispersion_Neutros_Violin.png", plot = plot_neutros_violin, width = 14, height = 6, dpi = 300)
ggsave("figures/loci/Publication_Dispersion_Neutros_Violin_No_Independence.png", plot = plot_neutros_violin, width = 14, height = 6, dpi = 300)


# Test comparación de medias con Wilcoxon y correción de Bonferri

# 1. Filtramos los datos a partir del umbral biológico Nm >= 1
df_filtrado <- df_neutros %>% filter(True_Mig >= 0.001)

# 2. Prueba Global de Kruskal-Wallis
cat("\n--- Prueba de Kruskal-Wallis Global ---\n")
kruskal_test <- kruskal.test(Abs_Error_D ~ Model, data = df_filtrado)
print(kruskal_test)

# Si el p-value es menor a 0.05, significa que un modelo es estadísticamente diferente a los otros.
if(kruskal_test$p.value < 0.05) {
  cat("\n--- Comparaciones Pareadas de Wilcoxon (Corrección Bonferroni) ---\n")
  # 3. Test post-hoc para ver quién le gana a quién
  pairwise_test <- pairwise.wilcox.test(
    df_filtrado$Abs_Error_D, 
    df_filtrado$Model, 
    p.adjust.method = "bonferroni", 
    exact = FALSE
  )
  print(pairwise_test)
}

# === 5. ANÁLISIS ESTADÍSTICO Y GRÁFICA DE RENDIMIENTO (Nm >= 1) ===

# 1. Filtramos los datos para la zona de buena señal biológica (Nm >= 1)
df_filtrado <- df_neutros %>% filter(True_Mig >= 0.001)

# 2. Calculamos el MAE promedio y el Error Estándar (SE) por modelo
df_summary_stats <- df_filtrado %>%
  group_by(Model) %>%
  summarise(
    Mean_MAE = mean(Abs_Error_D, na.rm = TRUE),
    SD_MAE = sd(Abs_Error_D, na.rm = TRUE),
    N = n(),
    SE_MAE = SD_MAE / sqrt(N),
    .groups = "drop"
  )

# 3. Generamos el gráfico de barras
plot_mae_stats <- ggplot(df_summary_stats, aes(x = Model, y = Mean_MAE, fill = Model)) +
  geom_bar(stat = "identity", color = "black", width = 0.6, alpha = 0.8) +
  # Agregamos las barras de error (Media +/- Error Estándar)
  geom_errorbar(aes(ymin = Mean_MAE - SE_MAE, ymax = Mean_MAE + SE_MAE), 
                width = 0.2, linewidth = 0.8) +
  scale_fill_manual(values = mis_colores) +
  labs(
    title = "Model Performance Comparison (True D \u2265 0.001)",
    subtitle = "Mean Absolute Error (MAE) \u00B1 Standard Error",
    x = "Inference Model",
    y = "Mean Absolute Error of Diffusion (MAE)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(), # Quitamos las líneas verticales para mayor limpieza
    axis.text.x = element_text(face = "bold", size = 12)
  )

# Guardamos la gráfica
ggsave("Supplementary_MAE_Comparison.png", plot = plot_mae_stats, width = 7, height = 6, dpi = 300)

####################################### LIKELIHOOD PROFILE ###################################################

library(ggplot2)
library(dplyr)
library(data.table)

# =========================================================================
# 1. CARGA Y ENSAMBLAJE DE DATOS (Resultados del Clúster)
# =========================================================================
# Supongamos que df_raw es la tabla donde ya juntaste todos los .txt de tu clúster
# (Si lo vas a hacer con datos empíricos, aquí cargas el txt de esos datos)
# df_raw <- fread("resultados_empiricos_completos.txt")

df_long

library(dplyr)

# =========================================================================
# INFERENCIA CONJUNTA (PROFILE LIKELIHOOD) A PARTIR DE DF_LONG
# =========================================================================

# Paso 1 y 2: Preparar y sumar la Verosimilitud Compuesta
df_composite <- df_long %>%
  # Limpieza: Evitar que un valor infinito arruine la suma de todo el genoma
  # (Veo un -1774161 en tu imagen para D=0.5, eso está bien, pero filtramos por si hay -Inf)
  filter(is.finite(Likelihood)) %>%
  
  # PERFILADO LOCAL: Agrupamos por SNP y valor D, y nos quedamos con el mejor 's'
  # (En tus neutros Inferred_s es 0, así que esto solo lo prepara para el futuro)
  group_by(TaskID, Model, SNP, Inferred_D) %>%
  slice_max(order_by = Likelihood, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  
  # VEROSIMILITUD CONJUNTA: Sumamos todos los SNPs para un mismo TaskID y D global
  group_by(TaskID, Model, Inferred_D) %>%
  summarise(
    Composite_LL = sum(Likelihood, na.rm = TRUE),
    N_SNPs_Usados = n_distinct(SNP),
    .groups = "drop"
  )

# Paso 3: Encontrar el MLE global (El pico de la montaña)
df_mle_global <- df_composite %>%
  group_by(TaskID, Model) %>%
  # Nos quedamos con la fila que tenga la Composite_LL más alta
  slice_max(order_by = Composite_LL, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(D_Estimado = Inferred_D, Max_LL = Composite_LL)

# Ver los resultados finales
cat("\n=== ESTIMACIONES GLOBALES POR TASK Y MODELO ===\n")
print(head(df_mle_global))

# =========================================================================
# 4. GRAFICAR LA CURVA DE PERFIL DE VEROSIMILITUD
# =========================================================================
plot_profile_likelihood <- ggplot(df_joint, aes(x = D, y = Composite_LL, color = Model, group = Model)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3) +
  
  # Marcar el pico (MLE) con una línea vertical para cada modelo
  geom_vline(data = df_mle_global, aes(xintercept = MLE_Global_D, color = Model), 
             linetype = "dashed", linewidth = 1, alpha = 0.7) +
  
  scale_x_log10(breaks = unique(df_joint$D), labels = scales::scientific) +
  
  theme_classic(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Joint Profile Likelihood for Diffusion Rate (D)",
    subtitle = "Aggregating spatial likelihoods across all independent SNPs",
    x = "Global Diffusion Rate (D)",
    y = "Joint Log-Likelihood",
    color = "Likelihood Framework"
  ) +
  theme(
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.major.x = element_line(color = "gray95", linetype = "dotted"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

print(plot_profile_likelihood)


# ============================================================
# 1. LISTADO Y PROCESAMIENTO INCREMENTAL (bajo uso de RAM)
# ============================================================

library(tidyverse)

directorio <- "data/results_Discrete/outputs_LL/independent_loci/seleccion"

rutas_archivos <- list.files(
  path       = directorio,
  pattern    = "_D_FULL_seleccion_m1\\.(txt|csv)$",
  full.names = TRUE
)
rutas_archivos <- rutas_archivos[!grepl("^P_TRON|^P_", basename(rutas_archivos))]

procesar_archivo <- function(ruta) {
  nombre_archivo <- basename(ruta)
  df <- read_tsv(ruta, col_select = starts_with("D_"), show_col_types = FALSE)

  cols_D <- names(df)[str_detect(names(df), "^D_")]
  if (length(cols_D) == 0) return(NULL)

  # Suma por punto de gradilla (ignora NA dentro del archivo)
  sumas <- colSums(as.matrix(df[, cols_D, drop = FALSE]), na.rm = TRUE)

  idx_mle        <- which.max(sumas)
  grid_param_mle <- cols_D[idx_mle]
  mle_valor      <- unname(sumas[idx_mle])

  # Parseo del punto MLE:  D_<val>_s_<val>
  D_val <- str_extract(grid_param_mle, "(?<=^D_)[^_]+")
  s_val <- str_extract(grid_param_mle, "(?<=_s_).*$")

  modelo <- case_when(
    str_detect(nombre_archivo, "TRON_ARES")   ~ "B-B",
    str_detect(nombre_archivo, "TRON_SPIKES") ~ "BwS",
    TRUE                                      ~ "Classic"
  )
  task_id <- as.integer(str_extract(nombre_archivo, "(?<=TaskID_)\\d+"))

  tibble(
    File           = nombre_archivo,
    Model          = modelo,
    TaskID         = task_id,
    Inferred_D     = as.numeric(D_val),
    Inferred_s     = as.numeric(s_val),
    Grid_Param_MLE = grid_param_mle,
    MLE_Sum        = mle_valor,
    n_SNPs         = sum(!is.na(df[[cols_D[1]]]))
  )
}

n          <- length(rutas_archivos)
resultados <- vector("list", n)

for (i in seq_len(n)) {
  resultados[[i]] <- tryCatch(
    procesar_archivo(rutas_archivos[i]),
    error = function(e) {
      message("Error en ", basename(rutas_archivos[i]), ": ", conditionMessage(e))
      NULL
    }
  )
  if (i %% 100 == 0) { gc(verbose = FALSE); message(sprintf("Procesados %d / %d", i, n)) }
}

tabla_final <- bind_rows(resultados)

# ============================================================
# 2. MAPEO TaskID -> (True_Mig, True_Sel)
#    Reproduce el ORDEN de los bucles en bash:
#       outer loop = MIG_VALUES
#       inner loop = SEL_VALUES  (rápido)
# ============================================================
# ORDEN EXACTO del bash (no reordenar):
# mig_values <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125)
# sel_values <- c(0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.0050, 0.0025, 0.001, 0.00075, 0.0005, 0.00025, 0.0001, 0.0)
replicas_per_val <- 50

mig_values <- c(0.0, 0.001, 0.1)
sel_values <- c(0.0, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1)

# Verificación rápida de consistencia con el bash
stopifnot(length(mig_values) * length(sel_values) * replicas_per_val ==
          max(tabla_final$TaskID, na.rm = TRUE))

tabla_final <- tabla_final %>%
  mutate(
    idx     = floor((TaskID - 1) / replicas_per_val),          # 0..(10*14-1)
    idx_mig = floor(idx / length(sel_values)) + 1,             # 1..10
    idx_sel = (idx %% length(sel_values)) + 1,                 # 1..14
    True_Mig = mig_values[idx_mig],
    True_Sel = sel_values[idx_sel],
    Replica  = (TaskID - 1) %% replicas_per_val + 1            # 1..50
  ) %>%
  select(-idx, -idx_mig, -idx_sel)

# ============================================================
# 3. ESCALADO Y ERRORES
#    (Wright-Fisher -> Classic: s se multiplica por 2)
# ============================================================
tabla_final <- tabla_final %>%
  mutate(
    Inferred_s_scaled = Inferred_s * 2,
    Error_D           = Inferred_D - True_Mig,
    Abs_Error_D       = abs(Error_D),
    Error_s           = Inferred_s_scaled - True_Sel,
    Abs_Error_s       = abs(Error_s)
  ) %>%
  relocate(Model, TaskID, Replica, True_Mig, True_Sel,
           Inferred_D, Inferred_s, Inferred_s_scaled,
           Error_D, Abs_Error_D, Error_s, Abs_Error_s,
           MLE_Sum, n_SNPs, Grid_Param_MLE, File)

# ============================================================
# 0. Preparación de tabla_final para el plot
# ============================================================
df_seleccion <- tabla_final %>%
  filter(!is.na(Inferred_s_scaled), !is.na(True_Sel)) %>%
  mutate(
    Model    = factor(Model, levels = c("Classic", "B-B", "BwS")),
    True_Mig = factor(True_Mig, levels = sort(unique(True_Mig)))
  )
)

# Línea de referencia y = x (cubre el rango simulado real: 0 → 0.1)
rango_max <- max(df_seleccion$True_Sel, na.rm = TRUE)
linea_perfecta <- data.frame(
  True_Sel          = c(0, rango_max),
  Inferred_s_scaled = c(0, rango_max)
)

# ============================================================
# 1. Plot: recuperación de s por modelo y por D (True_Mig)
# ============================================================
plot_seleccion_recovery <- ggplot(
    df_seleccion,
    aes(x = True_Sel, y = Inferred_s_scaled, color = Model)
  ) +
  geom_point(size = 1.8, alpha = 0.55) +
  # Línea identidad (perfect inference)
  geom_line(
    data = linea_perfecta,
    aes(x = True_Sel, y = Inferred_s_scaled),
    color = "red", linetype = "dashed", linewidth = 1,
    inherit.aes = FALSE
  ) +
  # Facetas: filas = Modelo, columnas = True_Mig (D)
  facet_grid(Model ~ True_Mig, labeller = labeller(True_Mig = label_both)) +
  # Eje X: pseudo-log para tolerar s = 0
  scale_x_continuous(
    trans  = pseudo_log_trans(sigma = 1e-5, base = 10),
    breaks = c(0, 1e-4, 1e-3, 1e-2, 1e-1),
    labels = scales::scientific
  ) +
  # Eje Y: pseudo-log para tolerar valores negativos y 0
  scale_y_continuous(
    trans  = pseudo_log_trans(sigma = 1e-5, base = 10),
    breaks = sort(c(-1e-2, -1e-3, -1e-4, 0,
                    1e-4, 1e-3, 1e-2, 1e-1)),
    labels = scales::scientific
  ) +
  scale_color_manual(values = mis_colores) +
  labs(
    title    = "Parameter Recovery: Selection Coefficient (s)",
    subtitle = "Red dashed line = perfect inference (y = x). Columns = True Migration Rate (D).",
    x        = "Simulated Selection Coefficient (True s)",
    y        = "Estimated Selection Coefficient (Scaled s)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "none",
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text       = element_text(face = "bold", size = 9)
  )

# ============================================================
# 2. Guardar
# ============================================================
ggsave(
  "figures/independent_loci/Exploracion/Selection_IndLoci_Test_Sel.png",
  plot   = plot_seleccion_recovery,
  width  = 16, height = 9, dpi = 300
)

# ============================================================
# EXPLORACION
#
# 0. MAPEO (mismo orden que el bash)
# ============================================================
mig_values <- c(0.0, 0.001, 0.1)
sel_values <- c(0.0, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1)
replicas_per_val <- 50

# Definimos los DOS grupos a explorar
grupos <- tibble(
  Grupo    = c("D0.001_s0.01", "D0.001_s0.025"),
  True_Mig = c(0.001,          0.001),
  True_Sel = c(0.01,          0.025),
  TaskID_ini = c(min(taskids_s005), min(taskids_s001)),
  TaskID_fin = c(max(taskids_s005), max(taskids_s001))
)

# D = 0.05  -> 7mo valor de mig_values
# s = 0.05  -> 3er valor de sel_values
# s = 0.01  -> 5to valor de sel_values
# (Verificado más abajo con stopifnot)

taskids_para <- function(D_objetivo, s_objetivo) {
  idx_mig <- which(mig_values == D_objetivo)
  idx_sel <- which(sel_values == s_objetivo)
  if (length(idx_mig) == 0 || length(idx_sel) == 0) stop("Valor no encontrado en el sweep")
  idx     <- (idx_mig - 1) * length(sel_values) + (idx_sel - 1)
  seq.int(idx * replicas_per_val + 1, (idx + 1) * replicas_per_val)
}

# Verificación del mapeo (por si cambian los vectores)
stopifnot(mig_values[2] == 0.001, sel_values[3] == 0.01, sel_values[4] == 0.025)

taskids_s005 <- taskids_para(0.001, 0.01)   # D=0.05, s=0.05
taskids_s001 <- taskids_para(0.001, 0.025)   # D=0.05, s=0.01

cat("D = 0.05 | s = 0.05  -> TaskIDs:", min(taskids_s005), "-", max(taskids_s005), "\n")
cat("D = 0.05 | s = 0.01  -> TaskIDs:", min(taskids_s001), "-", max(taskids_s001), "\n")
# Esperado:
#   D=0.05, s=0.05  -> 4301 - 4350
#   D=0.05, s=0.01  -> 4401 - 4450

info_rutas <- tibble(
  ruta   = rutas_todas,
  File   = basename(ruta),
  TaskID = as.integer(str_extract(File, "(?<=TaskID_)\\d+")),
  Model  = case_when(
    str_detect(File, "TRON_ARES")   ~ "B-B",
    str_detect(File, "TRON_SPIKES") ~ "BwS",
    TRUE                            ~ "Classic"
  )
)

# Filtrar solo los TaskIDs de interés
rutas_sel <- info_rutas %>%
  filter(TaskID %in% c(taskids_s005, taskids_s001)) %>%
  left_join(grupos, by = c("TaskID" = "TaskID_ini"))   # provisional
# Mejor unirlo explícitamente por rango:
rutas_sel <- info_rutas %>%
  rowwise() %>%
  mutate(
    Grupo = case_when(
      TaskID %in% taskids_s005 ~ "D_0.001_s0.01",
      TaskID %in% taskids_s001 ~ "D_0.001_s0.025",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(Grupo))

cat("\nArchivos encontrados por grupo y modelo:\n")
print(rutas_sel %>% count(Grupo, Model))


# ============================================================
# 2. LOCALIZAR ARCHIVOS DE ESOS TaskIDs (los 3 modelos)
# ============================================================
rutas_todas <- list.files(
  path       = directorio,
  pattern    = "_D_FULL_seleccion_m1\\.(txt|csv)$",
  full.names = TRUE
)
rutas_todas <- rutas_todas[!grepl("^P_TRON|^P_", basename(rutas_todas))]

info_rutas <- tibble(
  ruta   = rutas_todas,
  File   = basename(ruta),
  TaskID = as.integer(str_extract(File, "(?<=TaskID_)\\d+")),
  Model  = case_when(
    str_detect(File, "TRON_ARES")   ~ "B-B",
    str_detect(File, "TRON_SPIKES") ~ "BwS",
    TRUE                            ~ "Classic"
  )
)

# ============================================================
# 3. FUNCIÓN: PERFIL DE VEROSIMILITUD DE UN ARCHIVO
# ============================================================
perfilar_archivo <- function(ruta) {
  nombre <- basename(ruta)
  task   <- as.integer(str_extract(nombre, "(?<=TaskID_)\\d+"))
  model  <- case_when(
    str_detect(nombre, "TRON_ARES")   ~ "B-B",
    str_detect(nombre, "TRON_SPIKES") ~ "BwS",
    TRUE                              ~ "Classic"
  )

  df <- read_tsv(ruta, col_select = starts_with("D_"), show_col_types = FALSE)
  cols_D <- names(df)[str_detect(names(df), "^D_")]
  if (length(cols_D) == 0) return(NULL)

  mat        <- as.matrix(df[, cols_D, drop = FALSE])
  suma_total <- colSums(mat, na.rm = TRUE)
  n_validos  <- colSums(!is.na(mat))

  tibble(
    TaskID     = task,
    Model      = model,
    File       = nombre,
    Grid_Param = cols_D,
    Suma_LL    = as.numeric(suma_total),
    n_SNPs_val = as.numeric(n_validos)
  ) %>%
    mutate(
      D = as.numeric(str_extract(Grid_Param, "(?<=^D_)[^_]+")),
      s = as.numeric(str_extract(Grid_Param, "(?<=_s_).*$"))
    )
}

# ============================================================
# 4. LOOP SOBRE LOS ARCHIVOS SELECCIONADOS
# ============================================================
perfiles_lista <- vector("list", nrow(rutas_sel))

for (i in seq_len(nrow(rutas_sel))) {
  perfiles_lista[[i]] <- tryCatch(
    perfilar_archivo(rutas_sel$ruta[i]),
    error = function(e) { message("Error en ", rutas_sel$File[i]); NULL }
  )
  if (i %% 20 == 0) { gc(verbose = FALSE); message(sprintf("%d / %d", i, nrow(rutas_sel))) }
}

perfiles <- bind_rows(perfiles_lista) %>%
  mutate(
    Grupo = case_when(
      TaskID %in% taskids_s005 ~ "D_0.001_s0.01",
      TaskID %in% taskids_s001 ~ "D_0.001_s0.025",
      TRUE                     ~ NA_character_
    )
  ) %>%
  left_join(grupos %>% select(Grupo, True_Mig, True_Sel), by = "Grupo")

# Chequeo rápido
stopifnot(!any(is.na(perfiles$Grupo)),
          !any(is.na(perfiles$True_Mig)),
          !any(is.na(perfiles$True_Sel)))
cat("Grupos asignados correctamente:\n")

# ============================================================
# 5. MLE POR ARCHIVO
# ============================================================
mle_por_archivo <- perfiles %>%
  group_by(Grupo, Model, TaskID, File, True_Mig, True_Sel) %>%
  slice_max(Suma_LL, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    Grupo, Model, TaskID, File,
    True_Mig, True_Sel,
    Inferred_D        = D,
    Inferred_s        = s,
    Inferred_s_scaled = Inferred_s * FACTOR_EDP_TO_WF,
    MLE_Suma_LL       = Suma_LL,
    n_SNPs_val
  )

  # ============================================================
# 7. GRÁFICOS
# ============================================================

# --- (A) Perfiles individuales por grupo, facetado por modelo ---
p_A <- ggplot(perfiles,
              aes(x = s, y = Suma_LL, group = TaskID, color = factor(D))) +
  geom_line(linewidth = 0.3, alpha = 0.35) +
  facet_grid(Grupo ~ Model) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-5)) +
  labs(title    = "Perfiles log-LL — D = 0.05 con s ∈ {0.05, 0.01}",
       subtitle = "Cada línea es una réplica (n = 50 por celda). Colores = D del grid.",
       x = "s (grid EDP)", y = "Suma log-LL", color = "D") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# --- (B) Perfil promedio ± 1 SD por grupo y modelo ---
perfil_prom <- perfiles %>%
  group_by(Grupo, Model, True_Mig, True_Sel, D, s) %>%
  summarise(
    Suma_media = mean(Suma_LL, na.rm = TRUE),
    Suma_sd    = sd(Suma_LL,   na.rm = TRUE),
    .groups = "drop"
  )

p_B <- ggplot(perfil_prom,
              aes(x = s, y = Suma_media, color = factor(D), fill = factor(D))) +
  geom_ribbon(aes(ymin = Suma_media - Suma_sd, ymax = Suma_media + Suma_sd),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1) +
  geom_vline(data = grupos, aes(xintercept = True_Sel),
             linetype = "dotted", color = "black") +
  facet_grid(Grupo ~ Model) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-5)) +
  labs(title    = "Perfil promedio de log-LL (media ± 1 SD entre 50 réplicas)",
       subtitle = "Línea negra punteada = True_Sel del grupo",
       x = "s (grid EDP)", y = "Suma log-LL media", color = "D", fill = "D") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# --- (C) Distribución del Inferred_s por grupo y modelo ---
p_C <- ggplot(mle_por_archivo,
              aes(x = factor(Inferred_s_scaled), fill = Model)) +
  geom_bar(position = "dodge", alpha = 0.85) +
  geom_vline(data = grupos, aes(xintercept = which(sort(unique(mle_por_archivo$Inferred_s_scaled)) == True_Sel)),
             linetype = "dotted", color = "black") +
  facet_wrap(~ Grupo, ncol = 2, scales = "free_x") +
  scale_fill_manual(values = mis_colores) +
  labs(title    = "Distribución del MLE (Inferred_s_scaled) por modelo",
       subtitle = "Línea punteada = True_Sel del grupo",
       x = "Inferred_s_scaled (= 2 × Inferred_s)", y = "Nº réplicas") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# --- (D) Scatter: Inferred_s_scaled vs True_Sel, coloreado por modelo ---
p_D <- ggplot(mle_por_archivo,
              aes(x = True_Sel, y = Inferred_s_scaled, color = Model)) +
  geom_jitter(width = 0.002, height = 0, size = 2, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(~ Grupo) +
  scale_color_manual(values = mis_colores) +
  scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-5)) +
  scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5)) +
  labs(title    = "Recovery: Inferred_s_scaled vs True_Sel",
       subtitle = "Línea roja = y = x",
       x = "True_Sel", y = "Inferred_s_scaled") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_A); print(p_B); print(p_C); print(p_D)

ggsave("explorar_D005_s005_s001_A.png", p_A, width = 12, height = 7, dpi = 150)
ggsave("explorar_D005_s005_s001_B.png", p_B, width = 12, height = 7, dpi = 150)
ggsave("explorar_D005_s005_s001_C.png", p_C, width = 10, height = 5, dpi = 150)
ggsave("explorar_D005_s005_s001_D.png", p_D, width = 10, height = 5, dpi = 150)
