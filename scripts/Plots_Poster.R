library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(knitr)

# === 1. DEFINICIÓN DE TRUE VALUES ===
mig_values <- c(0.0, 0.01, 0.05, 0.1) 
sel_values <- c(0)                     
replicas_per_val <- 10

# === PALETA DE COLORES NEÓN (TRON) ===
# Asignamos colores neón específicos a los nuevos nombres de los modelos
tron_colors <- c("Binomial" = "#ff0055",   # Rojo/Rosa Neón
                 "B-B"      = "#00BFC4",   # Cian Brillante (o puedes usar #00ff66 para verde)
                 "BwS"      = "#ffaa00")   # Naranja Encendido

# === 2. LECTURA Y ETIQUETADO DE ARCHIVOS ===
directorio <- "results_validations" 
rutas_archivos <- list.files(path = directorio, pattern = "\\.txt$", full.names = TRUE)

procesar_archivo_downsampling <- function(ruta) {
  nombre_archivo <- basename(ruta)
  df <- read_tsv(ruta, show_col_types = FALSE)
  
  df <- df %>% mutate(SNP = as.character(SNP))
  
  modelo <- str_extract(nombre_archivo, "(?<=SIM_OPTIM_)[A-Za-z]+") 
  task_id <- as.integer(str_extract(nombre_archivo, "(?<=Task_)\\d+"))
  proporcion <- as.numeric(str_extract(nombre_archivo, "(?<=Prop_)[0-9\\.]+"))
  tipo <- ifelse(str_detect(nombre_archivo, "neutros"), "Neutro", "Seleccion")
  
  df %>%
    mutate(
      Model = modelo,
      Type = tipo,
      TaskID = task_id,
      Prop_Data = proporcion, # Asegurando que la columna se llame Prop_Data para el gráfico
      File = nombre_archivo
    )
}

df_raw <- map_dfr(rutas_archivos, procesar_archivo_downsampling)

df_plot <- df_raw %>%
  mutate(
    idx = floor((TaskID - 1) / replicas_per_val),
    idx_mig = floor(idx / length(sel_values)) + 1,
    idx_sel = (idx %% length(sel_values)) + 1,
    
    True_Mig = mig_values[idx_mig],
    True_Sel = sel_values[idx_sel],
    
    SNP = as.factor(SNP)
  )
df_plot_clean<-df_plot[df_plot[, 12] == "Neutro", ]
df_mae <- df_plot_clean %>%
  # Paso 1: Composite Likelihood (Suma de LL de todos los SNPs por Job/Réplica)
  group_by(TaskID, True_Mig, Model, Prop_Data, D) %>%
  summarise(Composite_LL = sum(LL, na.rm = TRUE), .groups = "drop") %>%
  
  # Paso 2: MLE (Encontrar el valor D que maximiza la Composite Likelihood)
  group_by(TaskID, True_Mig, Model, Prop_Data) %>%
  filter(Composite_LL == max(Composite_LL)) %>%
  slice(1) %>% # Desempate en caso de valores idénticos de LL
  ungroup() %>%
  rename(D_EST = D) %>%
  
  # Paso 3: Error Absoluto por réplica
  mutate(Absolute_Error = abs(D_EST - True_Mig)) %>%
  
  # Paso 4: Mean Absolute Error (MAE) por grupo de simulación (True_Mig)
  group_by(True_Mig, Model, Prop_Data) %>%
  summarise(MAE = mean(Absolute_Error, na.rm = TRUE), .groups = "drop")

df_mae_clean <- df_mae %>%
  filter(!is.na(True_Mig))

# APLICAR NUEVA NOMENCLATURA (Regla del Historial)
nombres_viejos <- c("LEGACY", "ARES", "SPIKES")
nombres_nuevos <- c("Binomial", "B-B", "BwS")

df_mae_clean <- df_mae_clean %>%
  mutate(Model = factor(Model, 
                        levels = nombres_viejos, 
                        labels = nombres_nuevos))

# =========================================================================
# 2. GRÁFICO UNIFICADO DE SENSIBILIDAD (ESTÉTICA TRON)
# =========================================================================
plot_mae_sensitivity <- ggplot(df_mae_clean, aes(x = Prop_Data, y = MAE, color = Model, group = Model)) +
  
  # Efecto Neón: Línea base más gruesa (Glow)
  geom_line(linewidth = 1.8, alpha = 0.5) +
  # Efecto Neón: Núcleo de la línea brillante
  geom_line(linewidth = 0.6, color = "#ffffff") +
  
  # Puntos brillantes
  geom_point(size = 3.5, alpha = 0.8) +
  geom_point(size = 1.5, color = "#ffffff") +
  
  scale_x_continuous(breaks = unique(df_mae_clean$Prop_Data)) +
  
  facet_wrap(~ paste("True Dispersal (\u03C3) =", True_Mig), # Cambiado a sigma para coherencia con tus fórmulas
             scales = "free_y", 
             ncol = 2, 
             nrow = 2) +
  
  # Asignar la paleta Neón
  scale_color_manual(values = tron_colors) + 
  
  labs(
    title = "Model Error Sensitivity Across Data Downsampling",
    subtitle = "Comparison of MAE as a function of the proportion of retained archaeogenomic data",
    x = "Proportion of Retained Data",
    y = "Mean Absolute Error (MAE)",
    color = "Likelihood Framework"
  ) +
  
  # === TEMA OSCURO CIBERNÉTICO ===
  theme_minimal(base_size = 14) +
  theme(
    # Fondos
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    
    # Cuadrícula (Rejilla de simulación)
    panel.grid.major.y = element_line(color = "#00ffff15", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "#00ffff15", linewidth = 0.3, linetype = "dashed"),
    panel.grid.minor   = element_blank(),
    
    # Ejes y Textos
    text             = element_text(color = "#00ffff"),
    axis.title       = element_text(color = "#00ffffcc", face = "bold"),
    axis.text        = element_text(color = "#00ffff66"),
    axis.line        = element_line(color = "#00ffff33"),
    
    # Leyenda
    legend.position   = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.text       = element_text(color = "#00ffffcc"),
    legend.title      = element_text(color = "#00ffff", face = "bold"),
    
    # Títulos y Facetas
    plot.title       = element_text(face = "bold", size = 16, color = "#00ffff"),
    plot.subtitle    = element_text(size = 12, color = "#00ffff99"),
    strip.background = element_rect(fill = "#00ffff15", color = "#00ffff33", linewidth = 0.5),
    strip.text       = element_text(face = "bold", size = 12, color = "#00ffff")
  )

# Desplegar
print(plot_mae_sensitivity)

# =========================================================================
# 3. EXPORTACIÓN A PDF PARA INKSCAPE
# =========================================================================
ggsave("Frame_D_MAE_Sensitivity_TRON.pdf", 
       plot = plot_mae_sensitivity, 
       device = "pdf", 
       width = 10, 
       height = 8, 
       units = "in", 
       dpi = 300,
       useDingbats = FALSE) # Vital para exportar vectores limpios a Inkscape

cat("Gráfico TRON exportado exitosamente como 'Frame_D_MAE_Sensitivity_TRON.pdf'\n")

# ========================================================================
# 4. GRÁFICO DE CURVAS DE VEROSIMILITUD COMPUESTAS (Ejemplo Conceptual)
# =======================================================================
library(ggplot2)
library(dplyr)

# Definir el eje X de parámetros (sigma o D)
x_vals <- seq(0, 1, length.out = 300)
# Valor real inferido
true_d <- 0.6 

# Generar varias distribuciones (locis)
set.seed(42)
n_loci <- 15
df_curves <- data.frame()
for(i in 1:n_loci) {
  mu <- true_d + rnorm(1, 0, 0.1) # Cada locus "piensa" que está en un lugar ligeramente distinto
  y <- exp(-(x_vals - mu)^2 / (2 * 0.05^2))
  df_curves <- rbind(df_curves, data.frame(x = x_vals, y = y, Locus = as.factor(i)))
}

# Calcular la suma (Composite Likelihood)
df_sum <- df_curves %>% group_by(x) %>% summarise(y = sum(y))

# Graficar
ggplot() +
  geom_line(data = df_curves, aes(x = x, y = y, group = Locus), 
            color = "#00ffff", alpha = 0.15, linewidth = 0.5) +
  geom_line(data = df_sum, aes(x = x, y = y), 
            color = "#ffaa00", linewidth = 2) +
  geom_vline(xintercept = true_d, linetype = "dashed", color = "#ff0055") +
  annotate("text", x = true_d, y = max(df_sum$y), label = expression(hat(D)), 
           color = "white", vjust = -1) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "#02050a"),
        plot.background = element_rect(fill = "#02050a"),
        panel.grid = element_blank(),
        axis.text = element_text(color = "#00ffff"),
        axis.title = element_text(color = "#00ffff"))
# Exportacion a PDF
ggsave("Poster_Composite_Likelihood_Curves_TRON.pdf", 
       width = 8, height = 6, units = "in", dpi = 300, useDingbats = FALSE)
cat("Gráfico de curvas de verosimilitud compuesto exportado como 'Composite_Likelihood_Curves_TRON.pdf'\n")

# ========================================================================
# 5. VIOLIN NEUTROS
# ========================================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
# Filtrar solo los datos neutros

# === 1. LECTURA Y ETIQUETADO DE ARCHIVOS ===
directorio <- "data/results_Discrete/outputs_LL/independent_loci"
rutas_archivos <- list.files(
  path = directorio,
  pattern = "^TRON_LEGACY_.*\\.txt$",
  full.names = TRUE
)

rutas_neutros <- grep("neutros", rutas_archivos, value = TRUE)

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

df_raw <- map_dfr(rutas_neutros, procesar_archivo)
df_raw <- df_raw %>%
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
  mutate(
    Inferred_D = as.numeric(str_extract(Grid_Param, "(?<=D_)[0-9\\.]+")),
    Inferred_s = as.numeric(str_extract(Grid_Param, "(?<=s_)[-0-9\\.]+"))
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

mig_values <- c(0.0, 0.000005, 0.00001, 0.00005,  0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05 ,0.075, 0.1, 0.125)
sel_values <- c(0.0)#c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
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

# Filtrar solo los datos neutros
df_neutros <- df_mle %>% filter(Type == "Neutro")

# Filtrar solo modelo "Classis"

N_pob <- 1000

df_classic <- df_neutros %>% 
  filter(Model == "Classic") %>%
  mutate(
    # Convertimos la proporción a individuos absolutos
    True_Humans     = True_Mig * N_pob,
    Inferred_Humans = Inferred_D * N_pob
  )
# 2. GRÁFICA DE VIOLÍN TRON (ESCALADA)
plot_neutros_box <- ggplot(df_classic, aes(x = as.factor(True_Humans), y = Inferred_Humans)) +
  
  geom_boxplot(
  width = 0.3, 
  fill = "#00BFC4",     # Relleno cian
  alpha = 0.2,          # Transparencia para que el rojo de la X resalte
  color = "#00BFC4",    # Borde brillante
  outlier.shape = NA, 
  linewidth = 0.8
) +
  
  geom_segment(aes(x = as.numeric(as.factor(True_Humans)) - 0.2, 
                 xend = as.numeric(as.factor(True_Humans)) + 0.2, 
                 y = True_Humans, 
                 yend = True_Humans), 
             color = "#ff0055", linewidth = 0.9) +

  geom_vline(xintercept = 7, linetype = "dashed", color = "#ff0055", linewidth = 0.8, alpha = 0.8) +

  facet_wrap(~ Model) +
  
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-3, base = 10),
    breaks = c(0, 0.01, 0.1, 1, 10, 100),
    labels = c("0", "0.01", "0.1", "1", "10", "100")
  ) +
  
  labs(
    title = "Dispersion Scaled to Human Movement",
    subtitle = "Boxplots show distribution of inferred movement",
    x = "Simulated Movement (Individuals per generation)",
    y = "Estimated Movement (Individuals per generation)"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    panel.grid.major.y = element_line(color = "#00ffff15", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "#00ffff15", linewidth = 0.3, linetype = "dashed"),
    panel.grid.minor   = element_blank(),
    text             = element_text(color = "#00ffff"),
    axis.title       = element_text(color = "#00ffffcc", face = "bold"),
    axis.text.y      = element_text(color = "#00ffff66"),
    axis.text.x      = element_text(angle = 45, hjust = 1, color = "#00ffff66"),
    axis.line        = element_line(color = "#00ffff33"),
    plot.title       = element_text(face = "bold", size = 16, color = "#00ffff"),
    plot.subtitle    = element_text(size = 12, color = "#00ffff99", face = "italic"),
    strip.background = element_rect(fill = "#00ffff10", color = "#00ffff33", linewidth = 0.5),
    strip.text       = element_text(face = "bold", size = 14, color = "#00BFC4")
  )


print(plot_neutros_box)

ggsave("Dispersion_Classic_TRON.pdf", plot = plot_neutros_box, width = 8, height = 6, device = "pdf", useDingbats = FALSE)
# ========================================================================
# 6. ERROR DIRECCIONAL: SOBRE vs SUBESTIMACIÓN (ESTILO TRON)
# ========================================================================

library(ggplot2)
library(dplyr)
library(scales)

# 1. CALCULAR EL ERROR RELATIVO (%)
df_relativo <- df_classic %>%
  # Excluimos el 0 para evitar la división por cero (Inf)
  filter(True_Humans > 0) %>%
  mutate(
    # Calculamos el porcentaje de error
    Relative_Error_Pct = ((Inferred_Humans - True_Humans) / True_Humans) * 100
  )

# 2. GENERAR LA GRÁFICA DE ERROR RELATIVO
plot_error_relativo <- ggplot(df_relativo, aes(x = as.factor(True_Humans), y = Relative_Error_Pct)) +
  
  # Línea en 0% (Estimación perfecta)
  geom_hline(yintercept = 0, linetype = "dashed", color = "#ff0055", linewidth = 1.2) +
  
  # Boxplots con relleno cian translúcido
  geom_boxplot(
    width = 0.5, 
    fill = "#00BFC4", 
    alpha = 0.2, 
    color = "#00BFC4", 
    outlier.color = "#00ffff", 
    outlier.alpha = 0.4,
    linewidth = 0.8
  ) +
  
  # Escala pseudo-log para manejar porcentajes positivos gigantes y negativos hasta -100%
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 10, base = 10),
    breaks = c(-100, -50, 0, 50, 100, 500, 1000),
    labels = function(x) paste0(x, "%") # Añade el símbolo de % automáticamente
  ) +
  
  labs(
    title = "Relative Estimation Error (%)",
    subtitle = "Percentage of deviation from the true simulated movement | 0% = Perfect Match",
    x = "Simulated Movement (Individuals per generation)",
    y = "Relative Error (%)"
  ) +
  
  # === TEMA OSCURO CIBERNÉTICO ===
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    
    # Cuadrícula
    panel.grid.major.y = element_line(color = "#00ffff15", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "#00ffff15", linewidth = 0.3, linetype = "dashed"),
    panel.grid.minor   = element_blank(),
    
    # Textos y Ejes
    text             = element_text(color = "#00ffff"),
    axis.title       = element_text(color = "#00ffffcc", face = "bold"),
    axis.text.y      = element_text(color = "#00ffff66", face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1, color = "#00ffff66"),
    
    # Títulos
    plot.title       = element_text(face = "bold", size = 16, color = "#00ffff"),
    plot.subtitle    = element_text(size = 12, color = "#ff0055", face = "italic")
  )

# 3. EXPORTAR A PDF
ggsave("Relative_Error_TRON.pdf", plot = plot_error_relativo, width = 8, height = 6, device = "pdf", useDingbats = FALSE)



###
# Seleccion
##

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

# Directorio donde están los txt (Ajusta esta ruta)
directorio <- "data/results_Discrete/outputs_LL/0505"

# 1. OBTENER LISTA DINÁMICA DE ARCHIVOS
rutas_seleccion <- list.files(
  path = directorio,
  pattern = "^TRON_LEGACY_.*\\_m2.txt$",
  full.names = TRUE
)

todos_los_ids <- as.numeric(gsub(".*TaskID_([0-9]+).*", "\\1", basename(rutas_seleccion)))

rutas_filtradas <- rutas_seleccion[todos_los_ids >= 1 & todos_los_ids <= 14000]

# Función para procesar un solo archivo
procesar_archivo <- function(ruta) {
  # Leer el archivo (usando read_tsv por las tabulaciones)
  df <- read_tsv(ruta, show_col_types = FALSE)
  nombre_archivo <- basename(ruta)
  # Extraer metadatos del nombre del archivo adaptado para 3 modelos
  task_id <- as.numeric(gsub(".*TaskID_([0-9]+).*", "\\1", nombre_base))
  
  tipo <- ifelse(str_detect(nombre_archivo, "neutros"), "Neutro", "Seleccion")
  task_id <- str_extract(nombre_archivo, "(?<=TaskID_)\\d+")
  
  # Agregar metadatos como columnas nuevas
  df %>%
    mutate(
      Type = tipo,
      TaskID = as.integer(task_id),
      File = nombre_archivo
    )
}

df_raw <- map_dfr(rutas_filtradas,procesar_archivo)

# Guardar
#saveRDS(df_raw, file = "df_raw.rds")

# Cargar
#df_raw <- readRDS("df_raw.rds")


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

# Definir tus vectores originales y N (Idéntico a tu Bash)
# mig_values <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125)
# sel_values <- c(0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.0050, 0.0025, 0.001, 0.00075, 0.0005, 0.00025, 0.0001, 0.0)
# replicas_per_val <- 10
# num_sel <- length(sel_values)
# N_pob <- 1000 

mig_values <- c(0.0, 0.01, 0.05, 0.1)
sel_values <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
replicas_per_val <- 10
num_sel <- length(sel_values)
N_pob <- 1000 

# 2. CALCULAR EL COMPOSITE LIKELIHOOD (Sumar todos los SNPs por cada punto del Grid)
df_composite <- df_long %>%
  # Agrupamos por Archivo (TaskID/Model/Type) y por el Punto del Grid
  group_by(Type, TaskID, Inferred_D, Inferred_s) %>%
  summarise(
    Composite_LL = sum(Likelihood),  
    N_SNPs_Usados = n(),             
    .groups = "drop"
  )

df_mle <- df_composite %>%
  group_by(Type, TaskID) %>%
  slice_max(order_by = Composite_LL, n = 1, with_ties = FALSE) %>%
  ungroup()

df_mle <- df_mle %>%
  mutate(
    idx = floor((TaskID - 1) / replicas_per_val),
    idx_mig = floor(idx / length(sel_values)) + 1,
    idx_sel = (idx %% length(sel_values)) + 1,
    True_Humans = mig_values[idx_mig] * N_pob,
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

df_seleccion <- df_mle %>% filter(Type == "Seleccion")



# Función para traducir el Task ID a parámetros reales
get_true_params <- function(task_id) {
  idx <- floor((task_id - 1) / replicas_per_val)
  idx_mig <- floor(idx / num_sel) + 1
  idx_sel <- (idx %% num_sel) + 1
  replica <- ((task_id - 1) %% replicas_per_val) + 1
  
  return(list(
    True_Mig = mig_values[idx_mig],
    True_Sel = sel_values[idx_sel],
    True_Humans = mig_values[idx_mig] * N_pob,
    Replica = replica
  ))
}

# 2. BUCLE SOBRE LOS ARCHIVOS REALES
resultados_lista <- list()

for (archivo in rutas_filtradas) {
  
  # Identificación del Job (ya sabemos que está entre 21 y 30)
  nombre_base <- basename(archivo)
  task_id <- as.numeric(gsub(".*TaskID_([0-9]+).*", "\\1", nombre_base))
  
  # Extraer parámetros reales simulados por el Proceso Cortana
  params <- get_true_params(task_id)
  
  # Lectura de la matriz de Verosimilitud
  datos_crudos <- read.table(archivo, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  # PROCESO DURANDAL: CÁLCULO DEL COMPOSITE LIKELIHOOD
  matriz_ll <- datos_crudos[, -1]
  composite_likelihoods <- colSums(matriz_ll, na.rm = TRUE)
  
  # Encontrar el argmax para las estimaciones conjuntas
  columna_mle <- names(composite_likelihoods)[which.max(composite_likelihoods)]
  
  # Extracción de parámetros MLE con Regex
  partes_nombre <- strsplit(columna_mle, "_s_")[[1]]
  s_hat <- as.numeric(partes_nombre[2])
  d_hat <- as.numeric(gsub("D_", "", partes_nombre[1]))
  
  # GUARDAR LOS RESULTADOS
  resultados_lista[[task_id]] <- data.frame(
    Task_ID = task_id,
    Replica = params$Replica,
    True_Mig = params$True_Mig,
    True_Humans = params$True_Humans,
    True_Sel = params$True_Sel,
    Inferred_Sel = s_hat,
    Inferred_Mig = d_hat
  )
}

# Unificar todo en el Data Frame final
df_seleccion <- bind_rows(resultados_lista)

# === SCRIPT DEL GRÁFICO FINAL MODIFICADO (SELECCIÓN CONTINUA) ===

plot_seleccion_grid <- ggplot(df_seleccion, aes(x = as.factor(True_Sel), y = Inferred_s)) +
  
  # 1. Boxplots restaurados: Al usar as.factor(x), el width = 0.6 funciona perfecto
  geom_boxplot(
    width = 0.6, 
    fill = "#00BFC4", 
    alpha = 0.2, 
    color = "#00BFC4", 
    outlier.color = "#00ffff", 
    outlier.alpha = 0.4,
    linewidth = 0.4
  ) +
  
  # 2. Diana segura: Usamos el guion (shape = 95) en lugar del segmento para evitar explosiones logarítmicas
  geom_point(aes(y = True_Sel), color = "#ff0055", shape = 95, size = 6, stroke = 1.5) +
  
  # Facetado por régimen de dispersión
  facet_wrap(
    ~ True_Humans, 
    ncol = 5, 
    labeller = labeller(True_Humans = function(x) paste(x, "migrants/gen"))
  ) +
  
  # 3. Escala Pseudo-Log SOLO en el eje Y (El eje X ahora se ordena uniformemente)
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1),
    labels = c("0", "1e-4", "0.001", "0.01", "0.05", "0.1")
  ) +
  
  labs(
    title = "Accuracy of Selection Inference Across Dispersion Regimes",
    subtitle = "Magenta dashes indicate the true simulated selection coefficient (s)",
    x = "True Simulated Selection Coefficient (s)",
    y = "Inferred Selection Coefficient (MLE)"
  ) +
  
  # === TEMA TRON ===
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    
    panel.grid.major.y = element_line(color = "#00ffff15", linewidth = 0.3),
    panel.grid.major.x = element_line(color = "#00ffff15", linewidth = 0.3, linetype = "dashed"),
    panel.grid.minor   = element_blank(),
    
    text             = element_text(color = "#00ffff"),
    axis.title       = element_text(color = "#00ffffcc", face = "bold", size = 14),
    axis.text.y      = element_text(color = "#00ffff66"),
    axis.text.x      = element_text(angle = 45, hjust = 1, color = "#00ffff66", size = 8),
    
    plot.title       = element_text(face = "bold", size = 18, color = "#00ffff"),
    plot.subtitle    = element_text(size = 14, color = "#ff0055", face = "italic"),
    
    strip.background = element_rect(fill = "#00ffff10", color = "#00ffff33", linewidth = 0.5),
    strip.text       = element_text(face = "bold", size = 11, color = "#00BFC4")
  )

ggsave("Selection_Inference_Grid_TRON_Fixed.pdf", plot = plot_seleccion_grid, width = 16, height = 8, device = "pdf", useDingbats = FALSE)

########################################################################################################################
########################################################################################################################
## GRAFICO PARA DATOS EMPIRICOS (NEUTROS)
########################################################################################################################
########################################################################################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

directorio_empirico <- "results/Empirical_data/"
datos_alelos <- "AADNA_data/CADD_Data/Stats_CADD_Ancestrales_0_5.csv"
rutas_empiricas <- list.files(path=directorio_empirico, pattern="Analysis_LL_.*\\.txt$", full.names = TRUE)

procesar_archivo <- function(ruta) {
  nombre_archivo <- basename(ruta)
  
  # Leer el archivo (usando read_tsv por las tabulaciones)
  df <- read_tsv(ruta, show_col_types = FALSE)
  
  tipo <- ifelse(str_detect(nombre_archivo, "neutrales"), "Neutro", "Seleccion")
  task_id <- as.numeric(str_extract(nombre_archivo, "(?<=_LL_)\\d+"))
  
  # Agregar metadatos como columnas nuevas
  df %>%
    mutate(
      Type = tipo,
      TaskID = as.integer(task_id),
      File = nombre_archivo
    )
}

df_raw_empiricos <- map_dfr(rutas_empiricas, procesar_archivo)
datos_alelos_df <- read_tsv(datos_alelos, show_col_types = FALSE) 

texto <- datos_alelos_df[[1]]

df <- read.table(text = texto, sep = ",", header = FALSE, stringsAsFactors = FALSE)
colnames(df) <- c("CHR", "POS", "A1", "A2", "NI_IDEAD","CADD", "SNP", "CHR_x", "POS_y", "Ref", "Anc")
df_sub <- df[, c("SNP", "POS")]

data_merged_ll_pos<- df_raw_empiricos %>%
  left_join(df_sub, by = "SNP") %>%
  filter(!is.na(POS)) %>%
  select(-File) %>%
  rename(Position = POS)
head(data_merged_ll_pos)

# Ordenar el data.frame por la columna Position
data_merged_ll_pos <- data_merged_ll_pos[order(data_merged_ll_pos$Position), ]

df_long <- data_merged_ll_pos %>%
  pivot_longer(
    cols = starts_with("D_"), 
    names_to = "Grid_Param", 
    values_to = "Likelihood"
  ) %>%
  drop_na(Likelihood, SNP) 

df_composite <- df_long %>%
  # Agrupamos por Archivo (TaskID/Model/Type) y por el Punto del Grid
  group_by(Grid_Param) %>%
  summarise(
    Composite_LL = sum(Likelihood),  
    N_SNPs_Usados = n(),             
    .groups = "drop"
  )

df_mle <- df_composite %>%
  slice_max(order_by = Grid_Param, n = 1, with_ties = FALSE) %>%
  ungroup()

# 1. CALCULAR MLE POR CADA SNP INDIVIDUAL (Para el Plot 1)
df_mle_individual <- df_long %>%
  group_by(SNP, Position, Type) %>%
  # Encontramos la Verosimilitud máxima para CADA SNP
  slice_max(order_by = Likelihood, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    # Extraer el valor numérico de D para poder graficarlo en el eje Y
    D_MLE = as.numeric(gsub("D_", "", Grid_Param))
  )

# 2. PREPARAR DATOS PARA LA CURVA COMPOSITE (Con Filtro de Integridad)
df_composite_clean <- df_composite %>%
  mutate(
    D_val = as.numeric(gsub("D_", "", Grid_Param))
  ) %>%
  # FILTRO INTELIGENTE: Conservamos estrictamente los puntos del grid 
  # que utilizaron el máximo número de SNPs disponibles.
  filter(N_SNPs_Usados == max(N_SNPs_Usados)) %>% 
  # Ordenamos numéricamente para la curva
  arrange(D_val)

# 3. OBTENER EL MLE GLOBAL CORREGIDO
# Esto te devolverá un data.frame con DOS filas en lugar de una
mle_globales <- df_composite_clean %>%
  slice_max(order_by = Composite_LL, n = 1, with_ties = TRUE)

print(mle_globales)

# Valor exacto del MLE para usarlo en las gráficas
best_D <- mle_globales$D_val 
max_LL <- mle_globales$Composite_LL

plot_mle_cromosoma <- ggplot(df_mle_individual, aes(x = Position, y = D_MLE)) +
  
  # Puntos para cada SNP
  geom_point(color = "#00FF66", alpha = 0.6, size = 2) +
  
  # Línea magenta horizontal que indica cuál fue el MLE GLOBAL de todo el cromosoma
  geom_hline(yintercept = best_D, color = "#ff0055", linetype = "dashed", linewidth = 1) +
  
  # Usamos pseudo_log porque D puede ser 0
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1),
    labels = c("0", "1e-4", "0.001", "0.01", "0.05", "0.1")
  ) +
  
  # Formatear el eje X (Posición genómica) para que use comas (ej. 150,000,000)
  scale_x_continuous(labels = comma) +
  
  labs(
    title = "Local Dispersal Inference across Chromosome 1",
    subtitle = "Dashed magenta line represents the Global Composite MLE",
    x = "Genomic Position (bp)",
    y = "Individual MLE for Dispersion (D)"
  ) +
  
  # TEMA VERDE NEÓN
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    panel.grid.major.y = element_line(color = "#00FF6633", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "#00FF6633", linewidth = 0.5, linetype = "dashed"),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "#00FF6655", linewidth = 0.6),
    text             = element_text(color = "#00FF66"),
    axis.title       = element_text(color = "#00FF66CC", face = "bold"),
    axis.text        = element_text(color = "#00FF6688", face = "bold"),
    plot.title       = element_text(face = "bold", size = 18, color = "#00FF66"),
    plot.subtitle    = element_text(size = 14, color = "#ff0055", face = "italic")
  )

ggsave("Empirical_CHR1_MLE_Distribution.pdf", plot = plot_mle_cromosoma, width = 12, height = 6)

plot_composite_curve <- ggplot(df_composite_clean, aes(x = D_val, y = Composite_LL)) +
  
  # Dibujar la curva de verosimilitud
  geom_line(color = "#00FF66", linewidth = 1.2) +
  geom_point(color = "#00FF66", size = 3) +
  
  # Marcar el pico máximo (MLE) con una línea vertical magenta
  geom_vline(xintercept = best_D, color = "#ff0055", linetype = "dashed", linewidth = 1.2) +
  
  # Marcar el punto exacto del MLE
  geom_point(data = mle_globales, aes(x = D_val, y = Composite_LL), 
             color = "#ff0055", size = 5, shape = 18) + # Rombo magenta
  
  # Escala logarítmica en X para que los valores pequeños (1e-04) no se aplasten contra el 0
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1),
    labels = c("0", "1e-4", "0.001", "0.01", "0.05", "0.1")
  ) +
  
  labs(
    title = "Composite Likelihood Surface for Empirical Data",
    subtitle = paste0("Maximum Likelihood Estimate (MLE): D = ", best_D),
    x = "Evaluated Dispersion Parameter (D)",
    y = "Log-Composite Likelihood"
  ) +
  
  # TEMA VERDE NEÓN
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    panel.grid.major.y = element_line(color = "#00FF6633", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "#00FF6633", linewidth = 0.5, linetype = "dashed"),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "#00FF6655", linewidth = 0.6),
    text             = element_text(color = "#00FF66"),
    axis.title       = element_text(color = "#00FF66CC", face = "bold"),
    axis.text        = element_text(color = "#00FF6688", face = "bold"),
    plot.title       = element_text(face = "bold", size = 18, color = "#00FF66"),
    plot.subtitle    = element_text(size = 14, color = "#ff0055", face = "italic")
  )

ggsave("Empirical_Composite_Likelihood_Curve.pdf", plot = plot_composite_curve, width = 10, height = 6)

library(ggplot2)
library(dplyr)
library(scales)

# Asegurar que manejamos los empates en el MLE global
mle_globales <- df_composite_clean %>%
  slice_max(order_by = Composite_LL, n = 1, with_ties = TRUE)

# Extraer los valores para las líneas guía (pueden ser uno o dos puntos debido al empate)
best_D_values <- mle_globales$D_val


# ==============================================================================
# OPCON A PARA EL CROMOSOMA: Puntos ultra-finos con Jittering (Evita solapamiento)
# ==============================================================================
plot_mle_cromosoma_puntos <- ggplot(df_mle_individual, aes(x = Position, y = D_MLE)) +
  # Usamos position_jitter para despegar ligeramente los puntos encimados en el eje Y
  geom_jitter(color = "#00BFC4", alpha = 0.15, size = 0.4, height = 0.05, width = 0) +
  
  # Líneas de los MLE globales empatados en magenta
  geom_hline(yintercept = best_D_values, color = "#ff0055", linetype = "dashed", linewidth = 0.8) +
  
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1),
    labels = c("0", "1e-4", "0.001", "0.01", "0.05", "0.1")
  ) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Local Dispersal Inference across Chromosome 1 (Optimized Points)",
    subtitle = "Dashed magenta lines represent the tied Global Composite MLEs",
    x = "Genomic Position (bp)", y = "Individual MLE for Dispersion (D)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    panel.grid.major.y = element_line(color = "#00BFC433", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "#00BFC433", linewidth = 0.5, linetype = "dashed"),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "#00BFC455", linewidth = 0.6),
    text             = element_text(color = "#00BFC4"),
    axis.title       = element_text(color = "#00BFC4CC", face = "bold"),
    axis.text        = element_text(color = "#00BFC488", face = "bold"),
    plot.title       = element_text(face = "bold", size = 18, color = "#00BFC4"),
    plot.subtitle    = element_text(size = 14, color = "#ff0055", face = "italic")
  )


# ==============================================================================
# OPCIÓN B PARA EL CROMOSOMA: Mapa de Densidad Hexagonal (RECOMENDADO para >20k SNPs)
# ==============================================================================
# Nota: Requiere instalar el paquete 'hexbin' (install.packages("hexbin"))
plot_mle_cromosoma_hex <- ggplot(df_mle_individual, aes(x = Position, y = D_MLE)) +
  geom_hex(bins = 60) +
  
  # Paleta de degradado que va desde el fondo negro pasando por el azul hasta el verde neón en el pico de densidad
  scale_fill_gradientn(
    colors = c("#02050a", "#00BFC4", "#00FF66"),
    values = c(0, 0.3, 1),
    name = "SNP Count"
  ) +
  
  geom_hline(yintercept = best_D_values, color = "#ff0055", linetype = "dashed", linewidth = 0.9) +
  
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1),
    labels = c("0", "1e-4", "0.001", "0.01", "0.05", "0.1")
  ) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Genomic Density of Dispersal Inference (Chromosome 1)",
    subtitle = "Dashed magenta lines represent the tied Global Composite MLEs",
    x = "Genomic Position (bp)", y = "Individual MLE for Dispersion (D)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    panel.grid.major = element_blank(), # Quitamos rejillas para resaltar los hexágonos
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "#00BFC455", linewidth = 0.6),
    text             = element_text(color = "#00BFC4"),
    axis.title       = element_text(color = "#00BFC4CC", face = "bold"),
    axis.text        = element_text(color = "#00BFC488", face = "bold"),
    legend.text      = element_text(color = "#00BFC4"),
    plot.title       = element_text(face = "bold", size = 18, color = "#00BFC4"),
    plot.subtitle    = element_text(size = 14, color = "#ff0055", face = "italic")
  )


# ==============================================================================
# GRÁFICA DE CURVA COMPOSITE RESTAURADA A AZUL (`#00BFC4`)
# ==============================================================================
plot_composite_curve <- ggplot(df_composite_clean, aes(x = D_val, y = Composite_LL)) +
  geom_line(color = "#00BFC4", linewidth = 1.2) +
  geom_point(color = "#00BFC4", size = 3) +
  
  # Dibujar líneas para TODOS los valores del MLE global empatados
  geom_vline(xintercept = best_D_values, color = "#ff0055", linetype = "dashed", linewidth = 1.2) +
  geom_point(data = mle_globales, aes(x = D_val, y = Composite_LL), 
             color = "#ff0055", size = 5, shape = 18) +
  
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1),
    labels = c("0", "1e-4", "0.001", "0.01", "0.05", "0.1")
  ) +
  labs(
    title = "Composite Likelihood Surface for Empirical Data",
    subtitle = paste0("Tied Maximum Likelihood Estimates (MLE): D = ", paste(best_D_values, collapse = " and ")),
    x = "Evaluated Dispersion Parameter (D)", y = "Log-Composite Likelihood"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    panel.grid.major.y = element_line(color = "#00BFC433", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "#00BFC433", linewidth = 0.5, linetype = "dashed"),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "#00BFC455", linewidth = 0.6),
    text             = element_text(color = "#00BFC4"),
    axis.title       = element_text(color = "#00BFC4CC", face = "bold"),
    axis.text        = element_text(color = "#00BFC488", face = "bold"),
    plot.title       = element_text(face = "bold", size = 18, color = "#00BFC4"),
    plot.subtitle    = element_text(size = 14, color = "#ff0055", face = "italic")
  )

# Guardar los archivos resultantes
ggsave("Empirical_Composite_Likelihood_Curve.pdf", plot = plot_composite_curve, width = 10, height = 6)
ggsave("Empirical_CHR1_MLE_Points_Optimized.pdf", plot = plot_mle_cromosoma_puntos, width = 12, height = 6)
ggsave("Empirical_CHR1_MLE_HexDensity.pdf", plot = plot_mle_cromosoma_hex, width = 12, height = 6)

# Instala ggrepel si no lo tienes: install.packages("ggrepel")
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel)

# 1. LEER Y PROCESAR LOS DATOS DE LA MATRIZ DE VENTAJOSOS
df_ventajosos <- read.table("Analysis_LL_Ventajosos.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Convertir a formato largo y extraer el MLE para cada SNP
df_mle_vent <- df_ventajosos %>%
  pivot_longer(
    cols = -SNP, 
    names_to = "Param_Grid", 
    values_to = "Likelihood"
  ) %>%
  group_by(SNP) %>%
  slice_max(order_by = Likelihood, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Magia Regex para separar D y s de los encabezados (ej. "D_1e-04_s_0.05")
  mutate(
    D_hat = as.numeric(sub("D_([0-9e.-]+)_s_.*", "\\1", Param_Grid)),
    s_hat = as.numeric(sub(".*_s_([0-9e.-]+)", "\\1", Param_Grid))
  )

# 2. DEFINIR UNA PALETA DE COLORES NEÓN DISTINTIVA
# Generamos suficientes colores de alto contraste para destacar cada alelo en el fondo negro
neon_palette <- c(
  "#FF0055", # Magenta TRON
  "#39FF14", # Verde Radioactivo
  "#FFD700", # Amarillo Oro
  "#00FFFF", # Cian Puro
  "#FF8C00", # Naranja Neón
  "#B026FF", # Púrpura Eléctrico
  "#FE4164", # Rosa Neón
  "#FDFEFE", # Blanco Brillante
  "#00FF66", # Verde Ventajoso
  "#0088FF"  # Azul Intenso
)

# 3. GENERAR EL GRÁFICO 2D ESTILO TRON
plot_ventajosos <- ggplot(df_mle_vent, aes(x = s_hat, y = D_hat, color = SNP)) +
  
  # Puntos grandes y brillantes
  geom_point(size = 6, alpha = 0.9) +
  
  # Etiquetas inteligentes que repelen a otras para no encimarse
  geom_text_repel(
    aes(label = SNP),
    size = 5,
    fontface = "bold",
    color = "#FDFEFE",         # Texto en blanco para máxima legibilidad
    bg.color = "#02050a",      # Borde negro alrededor del texto para que contraste
    bg.r = 0.15,
    box.padding = 0.8,
    point.padding = 0.5,
    segment.color = "#00BFC488", # Líneas guía en azul cian translúcido
    segment.size = 0.6
  ) +
  
  # Aplicar la paleta neón
  scale_color_manual(values = neon_palette) +
  
  # Escalas Pseudo-Log para ambos ejes (maneja ceros y notación científica)
  scale_x_continuous(
    trans = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5),
    labels = c("0", "1e-4", "0.001", "0.01", "0.05", "0.1", "0.5")
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = c(0, 1e-4, 1e-3, 0.01, 0.05, 0.1),
    labels = c("0", "1e-4", "0.001", "0.01", "0.05", "0.1")
  ) +
  
  labs(
    title = "Joint Inference of Selection and Dispersion (Beneficial Alleles)",
    subtitle = "Spatial mapping of Maximum Likelihood Estimates (\U0001D412\U0302 and \U0001D403\U0302)",
    x = "Inferred Selection Coefficient (s)",
    y = "Inferred Dispersion Rate (D)"
  ) +
  
  # TEMA TRON (Fondo oscuro, ejes azules)
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none", # Quitamos la leyenda porque las etiquetas ya dicen qué es cada punto
    panel.background = element_rect(fill = "#02050a", color = NA),
    plot.background  = element_rect(fill = "#02050a", color = NA),
    
    # Cuadrícula y ejes en tu azul cian
    panel.grid.major = element_line(color = "#00BFC433", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "#00BFC4", linewidth = 0.8),
    
    text             = element_text(color = "#00BFC4"),
    axis.title       = element_text(color = "#00BFC4", face = "bold", size = 15),
    axis.text        = element_text(color = "#00BFC4AA", face = "bold"),
    
    plot.title       = element_text(face = "bold", size = 18, color = "#00BFC4"),
    plot.subtitle    = element_text(size = 14, color = "#ff0055", face = "italic")
  )

ggsave("Empirical_Beneficial_Alleles_2DSpace.pdf", plot = plot_ventajosos, width = 10, height = 8, device = "pdf")