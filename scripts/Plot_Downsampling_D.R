library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(knitr)

# === 1. DEFINICIÓN DE TRUE VALUES (¡Ajusta estos valores!) ===
# Asumimos que para esta validación tienes un valor real de D (y/o s) 
# con el que simulaste los datos antes del downsampling.
mig_values <- c(0.0, 0.01, 0.05, 0.1) # Ejemplo: Los valores D que simulaste
sel_values <- c(0)                     # 0 porque son neutros
replicas_per_val <- 10

mis_colores <- c("ARES" = "#F8766D", "SPIKES" = "#00BFC4", "Classic" = "#C77CFF")

# === 2. LECTURA Y ETIQUETADO DE ARCHIVOS ===
directorio <- "results_validations" # Ajusta a tu directorio
rutas_archivos <- list.files(path = directorio, pattern = "\\.txt$", full.names = TRUE)

procesar_archivo_downsampling <- function(ruta) {
  nombre_archivo <- basename(ruta)
  df <- read_tsv(ruta, show_col_types = FALSE)
  
  # Forzar SNP a character por seguridad
  df <- df %>% mutate(SNP = as.character(SNP))
  
  # Extraer metadatos basados en la nueva nomenclatura:
  # Ej: SIM_OPTIM_ARES_Task_1_Prop_0.1_D_FULL_neutros.txt
  
  modelo <- str_extract(nombre_archivo, "(?<=SIM_OPTIM_)[A-Za-z]+") 
  task_id <- as.integer(str_extract(nombre_archivo, "(?<=Task_)\\d+"))
  proporcion <- as.numeric(str_extract(nombre_archivo, "(?<=Prop_)[0-9\\.]+"))
  tipo <- ifelse(str_detect(nombre_archivo, "neutros"), "Neutro", "Seleccion")
  
  df %>%
    mutate(
      Model = modelo,
      Type = tipo,
      TaskID = task_id,
      Proportion = proporcion,
      File = nombre_archivo
    )
}

df_raw <- map_dfr(rutas_archivos, procesar_archivo_downsampling)

df_plot <- df_raw %>%
  mutate(
    # Lógica de asignación basada en TaskID
    idx = floor((TaskID - 1) / replicas_per_val),
    idx_mig = floor(idx / length(sel_values)) + 1,
    idx_sel = (idx %% length(sel_values)) + 1,
    
    # Extraer los valores reales
    True_Mig = mig_values[idx_mig],
    True_Sel = sel_values[idx_sel],
    
    # Asegurar que SNP sea factor para agrupar las líneas correctamente
    SNP = as.factor(SNP)
  )

df_mae_clean <- df_mae %>%
  # Eliminar cualquier fila donde True_Mig no haya sido asignado correctamente
  filter(!is.na(True_Mig))
nombres_viejos <- c("LEGACY", "ARES", "SPIKES")
nombres_nuevos <- c("Binomial", "B-B", "BwS")
df_mae_clean <- df_mae_clean %>%
  mutate(Model = factor(Model, 
                        levels = nombres_viejos, 
                        labels = nombres_nuevos))
# =========================================================================
# 2. GRÁFICO UNIFICADO DE SENSIBILIDAD (Cuadrícula 2x2)
# =========================================================================
plot_mae_sensitivity <- ggplot(df_mae_clean, aes(x = Prop_Data, y = MAE, color = Model, group = Model)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 3, alpha = 0.9) +
  
  scale_x_continuous(breaks = unique(df_mae_clean$Prop_Data)) +
  
  # Forzar el layout de 2x2 paneles
  facet_wrap(~ paste("True Migration (m) =", True_Mig), 
             scales = "free_y", 
             ncol = 2, 
             nrow = 2) +
  
  theme_classic(base_size = 14) +
  scale_color_brewer(palette = "Set1") + 
  labs(
    title = "Model Error Sensitivity Across Data Downsampling",
    subtitle = "Comparison of MAE as a function of the proportion of retained archaeogenomic data",
    x = "Proportion of Retained Data (prop)",
    y = "Mean Absolute Error (MAE)",
    color = "Likelihood Framework"
  ) +
  theme(
    panel.grid.major.y = element_line(color = "gray93", linetype = "dashed"),
    panel.grid.major.x = element_line(color = "gray95", linetype = "dotted"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 16),
    strip.background = element_rect(fill = "gray95", color = NA),
    strip.text = element_text(face = "bold", size = 12)
  )

# Desplegar el gráfico en el visor de R (opcional)
print(plot_mae_sensitivity)

# =========================================================================
# 3. EXPORTACIÓN A PDF PARA PUBLICACIÓN
# =========================================================================
# Guardar en alta calidad. Las dimensiones 10x8 pulgadas suelen ser ideales 
# para encajar dos columnas en el formato estándar de revistas científicas.
ggsave("MAE_Sensitivity_Plot_2x2.pdf", 
       plot = plot_mae_sensitivity, 
       device = "pdf", 
       width = 10, 
       height = 8, 
       units = "in", 
       dpi = 300)

cat("Gráfico exportado exitosamente como 'MAE_Sensitivity_Plot_2x2.pdf'\n")
