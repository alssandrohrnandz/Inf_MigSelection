library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(knitr)

args <- commandArgs(trailingOnly = TRUE)

input_dir   <- "/mnt/data/dortega/hlopezh/Inf_MigSelection/data/results_Discrete/outputs_LL/independent_loci"
output_dir  <- "../../data/results_Discrete/figures/Global_Sel_D"
prefix <- "P_TRON_LEGACY_Grid_TaskID"
# === 1. DEFINICIÓN DE TRUE VALUES ===
mig_values <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125) 
sel_values <- c(0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.0050, 0.0025, 0.001, 0.00075, 0.0005, 0.00025, 0.0001, 0.0)                   
replicas_per_val <- 50

prefix <- "P_TRON_LEGACY_Grid_TaskID"
input_dir <- "/mnt/data/dortega/hlopezh/Inf_MigSelection/data/results_Discrete/outputs_LL/independent_loci"

# Opción 1: Empieza con el prefijo y termina en .txt (Recomendada)
ruta <- list.files(
  path = input_dir, 
  pattern = paste0("^", prefix, ".*\\.txt$"), 
  full.names = TRUE
)

# 1. Función de lectura adaptada a tu nomenclatura actual
procesar_archivo <- function(ruta) {
  nombre_archivo <- basename(ruta)
  df <- read_tsv(ruta, show_col_types = FALSE)
  
  # Adaptar regex a: P_TRON_LEGACY_Grid_TaskID_4535_D_FULL_seleccion_m1.txt
  modelo <- str_extract(nombre_archivo, "(?<=P_TRON_)[A-Za-z]+") # Extraerá "LEGACY"
  task_id <- as.integer(str_extract(nombre_archivo, "(?<=TaskID_)\\d+")) # Extraerá 4535
  tipo <- ifelse(str_detect(nombre_archivo, "neutros"), "Neutro", "Seleccion")
  
  df %>%
    mutate(
      Model = modelo,
      Type = tipo,
      TaskID = task_id,
      File = nombre_archivo
    )
}

# 2. Leer todos los archivos
df_raw <- map_dfr(ruta, procesar_archivo)
saveRDS(df_raw, file = "df_raw.rds")

# 3. Asignar los verdaderos valores de Migración y Selección según el TaskID
df_plot <- df_raw %>%
  mutate(
    # La matemática de índices está perfecta
    idx = floor((TaskID - 1) / replicas_per_val),
    idx_mig = floor(idx / length(sel_values)) + 1,
    idx_sel = (idx %% length(sel_values)) + 1,
    
    True_Mig = mig_values[idx_mig],
    True_Sel = sel_values[idx_sel]
  )

# 4. Calcular Composite Likelihood y obtener el MLE por archivo (TaskID)
df_mle_per_file <- df_plot %>%
  # Paso A: Composite Likelihood (Sumamos el LL de todos los SNPs para cada combinación de D y S evaluada en el grid)
  # IMPORTANTE: Asumo que en tu TXT tienes las columnas D, S y LL
  group_by(TaskID, File, Model, Type, True_Mig, True_Sel) %>%
  summarise(across(starts_with("D_"), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  
  # Paso B: Transformar de formato ancho a largo para poder buscar el máximo
  # Esto convertirá todas las columnas "D_X_s_Y" en una columna llamada "Grid_Point" 
  # y sus valores en "Composite_LL"
  pivot_longer(
    cols = starts_with("D_"),
    names_to = "Grid_Point",
    values_to = "Composite_LL"
  ) %>%
  
  # Paso C: Encontrar el MLE (La combinación de D y S que da el Max Composite_LL por archivo)
  group_by(TaskID, File, Model, Type, True_Mig, True_Sel) %>%
  filter(Composite_LL == max(Composite_LL)) %>%
  slice(1) %>% # Desempate si dos puntos del grid tienen exactamente el mismo LL
  ungroup() %>%
  
  # Paso D: Extraer los valores numéricos de D y S del nombre de la columna
  # Usamos expresiones regulares para sacar los números (incluyendo notación científica como 1e-04)
  mutate(
    D_EST = as.numeric(str_extract(Grid_Point, "(?<=D_)[0-9e\\.-]+")),
    S_EST = as.numeric(str_extract(Grid_Point, "(?<=s_)[0-9e\\.-]+"))
  ) %>%
  
  # Limpieza: quitamos la columna de texto original del grid point para dejarlo limpio
  select(-Grid_Point)

# Aquí puedes seguir con tu bloque de renombramiento de modelos:
nombres_viejos <- c("LEGACY", "ARES", "SPIKES")
nombres_nuevos <- c("Binomial", "B-B", "BwS")

df_final <- df_mle_per_file %>%
  mutate(Model = factor(Model, 
                        levels = nombres_viejos, 
                        labels = nombres_nuevos))

saveRDS(df_mle_per_file, file = "df_mle_grid_results.rds")

## ANALISIS

summary_tbl <- df_mle_per_file %>%
  group_by(True_Mig, True_Sel, Model, Type) %>%
  summarise(
    bias_D = mean(D_EST - True_Mig, na.rm = TRUE),
    bias_S = mean(S_EST - True_Sel, na.rm = TRUE),
    rmse_D = sqrt(mean((D_EST - True_Mig)^2, na.rm = TRUE)),
    rmse_S = sqrt(mean((S_EST - True_Sel)^2, na.rm = TRUE)),
    n_rep  = n(),
    .groups = "drop"
  )

df_mle_per_file_corregido <- df_mle_per_file %>%
  mutate(
    error_D    = D_EST - True_Mig,
    error_S    = S_EST - True_Sel,
    True_Mig_f = factor(True_Mig, levels = sort(unique(True_Mig))),
    True_Sel_f = factor(True_Sel, levels = sort(unique(True_Sel)))
  )

# Heatmap de sesgo en D
p1 <- ggplot(summary_tbl, aes(x = factor(True_Mig), y = factor(True_Sel), fill = bias_D)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  facet_wrap(~ Model) +
  labs(x = "D real", y = "s real", fill = "Sesgo (D_est - D_real)") +
  theme_minimal()

# ggsave ajustado con el nombre correcto de variable y archivo representativo
ggsave(
  filename = "Heatmap_Bias_D_SpatialGrid.pdf", 
  plot = p1, 
  width = 10, 
  height = 8, 
  device = "pdf"
)

library(ggplot2)

# 4A. Error de D, por valor real de D, faceteado por modelo
p2 <- ggplot(df_mle_per_file_corregido, aes(x = True_Mig_f, y = error_D)) +
    geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Model) +
    labs(
        x = "D real", 
        y = "Error (D_est - D_real)",
        title = "Distribución del error de D por valor real"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Guardar p2
ggsave(
    filename = "Boxplot_Error_D_Distribution.pdf",
    plot = p2,
    width = 9,
    height = 6,
    device = "pdf"
)

# 4B. Error de s, por valor real de s, faceteado por modelo
p3 <- ggplot(df_mle_per_file_corregido, aes(x = True_Sel_f, y = error_S)) +
    geom_boxplot(outlier.alpha = 0.3, outlier.size = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~ Model) +
    labs(
        x = "s real", 
        y = "Error (s_est - s_real)",
        title = "Distribución del error de s por valor real"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Guardar p3
ggsave(
    filename = "Boxplot_Error_S_Distribution.pdf",
    plot = p3,
    width = 9,
    height = 6,
    device = "pdf"
)