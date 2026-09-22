library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

# === 1. LECTURA Y PROCESAMIENTO DE ARCHIVOS SLIM ===
directorio <- "data/results_Discrete/outputs_slim/independent_loci/seleccion"
# Cambiamos el patrón para buscar los .csv de SLiM
rutas_archivos <- list.files(path = directorio, pattern = "\\.csv$", full.names = TRUE)

# Función adaptada para extraer la duración (Lifespan) de los alelos
procesar_archivo_slim <- function(ruta) {
  nombre_archivo <- basename(ruta)
  
  # Usamos read.csv porque la imagen indica delimitación por comas
  df <- read.csv(ruta, header = TRUE, sep = ",")
  
  # CORRECCIÓN 1: Extraer el TaskID (números justo antes de .csv)
  tipo <- ifelse(str_detect(nombre_archivo, "neutros"), "Neutro", "Seleccion")
  task_id <- str_extract(nombre_archivo, "\\d+(?=\\.csv$)")
  
  df_resumen <- df %>%
    # CORRECCIÓN 2: Homogeneizar los tipos de datos antes de operar
    mutate(
      Generation = as.numeric(Generation),
      MutationID = as.character(MutationID)
    ) %>%
    # Eliminar posibles NAs en Generation para evitar warnings en min()/max()
    filter(!is.na(Generation)) %>%
    group_by(MutationID) %>%
    summarise(
      # max - min + 1 nos da la duración total en generaciones
      Lifespan = max(Generation, na.rm = TRUE) - min(Generation, na.rm = TRUE) + 1, 
      .groups = "drop"
    ) %>%
    summarise(
      Mean_Lifespan = mean(Lifespan, na.rm = TRUE),
      Total_Mutations = n(), 
      se = sd(Lifespan, na.rm = TRUE) / sqrt(n())
    ) %>%
    mutate(
      Type = tipo,
      TaskID = as.integer(task_id),
      File = nombre_archivo
    )
  
  return(df_resumen)
}

# Iterar sobre todos los .csv y colapsarlos en un Master DF ligero
df_raw <- map_dfr(rutas_archivos, procesar_archivo_slim)


# === 2. ASIGNACIÓN DEL GRID DE PARÁMETROS (TRUE_MIG, TRUE_SEL) ===
mig_values <- c(0.0, 0.01, 0.05, 0.1)
sel_values <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001)
replicas_per_val <- 10

df_lifespans <- df_raw %>%
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


# === 3. VISUALIZACIÓN ===
# Para ver el efecto de la selección, filtramos los datos no neutrales
df_plot_seleccion <- df_lifespans %>%
  filter(!is.na(True_Mig))

plot_lifespan <- ggplot(df_plot_seleccion, aes(x = True_Sel, y = Mean_Lifespan)) +
  geom_boxplot(
    aes(group = True_Sel),      # Agrupa las cajas por cada valor del coeficiente de selección
    alpha = 0.3,              # Hace las cajas translúcidas para que sigas viendo los puntos
    width = 0.2,              # Ajusta el ancho de las cajas
    outlier.shape = NA        # Oculta los outliers por defecto para no duplicar los puntos translúcidos que ya tienes
  ) +
  # Puntos individuales de cada réplica de SLiM con algo de transparencia
  geom_jitter(width = 0, height = 0, size = 1.5, alpha = 0.3) +
  
  # Líneas y puntos gruesos que marcan el promedio general de las 10 réplicas
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # Facet_wrap por tasa de migración simulada
  facet_wrap(~ True_Mig, labeller = labeller(True_Mig = label_both), ncol = 2) +
  
  # Escala logarítmica en el eje X para s
  scale_x_log10(
    breaks = sel_values,
    labels = scales::scientific
  ) +
  
  labs(
    title = "Average Allele Lifespan by Selection Coefficient",
    subtitle = "Faded points are individual replicates; solid lines are means. Panels show True Migration Rate (D).",
    x = "Simulated Selection Coefficient (True s) [log10 scale]",
    y = "Mean Lifespan (Generations)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12)
  )

print(plot_lifespan)
ggsave("Publication_Allele_Lifespan.png", plot = plot_lifespan, width = 12, height = 8, dpi = 300)


library(tidyverse)
library(data.table)   # fread + fwrite: más rápido y menor RAM que read.csv

# ============================================================
# 0. CONFIGURACIÓN
# ============================================================
directorio <- "data/results_Discrete/outputs_slim/independent_loci/seleccion"
archivo_salida  <- "results/slim_lifespan_summary.csv"

rutas_archivos <- list.files(directorio, pattern = "\\.csv$", full.names = TRUE)
rutas_archivos <- rutas_archivos[!grepl("^P_TRON|^P_", basename(rutas_archivos))]

cat("Archivos a procesar:", length(rutas_archivos), "\n")

# ============================================================
# 1. FUNCIÓN: leer UN archivo -> 1 fila de resumen
# ============================================================
# Columnas asumidas: Generation, MutationID, Frequency
# Si tu SLiM output usa otros nombres, ajusta el `select =`.
procesar_archivo_slim <- function(ruta) {
  nombre <- basename(ruta)

  # Lectura eficiente: solo las columnas necesarias
  dt <- tryCatch(
    fread(ruta, select = c("Generation", "MutationID", "Frequency"),
          showProgress = FALSE, data.table = TRUE),
    error = function(e) { message("Fallo lectura ", nombre, ": ", e$message); NULL }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  # Homogeneizar tipos
  dt[, Generation := as.numeric(Generation)]
  dt[, MutationID := as.character(MutationID)]
  dt <- dt[!is.na(Generation)]

  if (nrow(dt) == 0) return(NULL)

  gen_max <- max(dt$Generation)

  # ---- Resumen por mutación ----
  per_mut <- dt[, .(
    first_gen = min(Generation),
    last_gen  = max(Generation),
    lifespan  = max(Generation) - min(Generation) + 1,
    max_freq  = max(Frequency, na.rm = TRUE),
    last_freq = Frequency[which.max(Generation)]
  ), by = MutationID]

  # ---- Clasificación fijado / perdido / segregando ----
  per_mut[, status := fifelse(
    last_freq >= 0.99,                        "Fixed",
    fifelse(last_gen < gen_max | last_freq == 0, "Lost",
                                                 "Segregating")
  )]

  # ---- Metadatos ----
  tipo    <- ifelse(str_detect(nombre, "neutros"), "Neutro", "Seleccion")
  task_id <- as.integer(str_extract(nombre, "\\d+(?=\\.csv$)"))

  # ---- Resumen a nivel de archivo (UNA fila) ----
  tibble(
    File                = nombre,
    TaskID              = task_id,
    Type                = tipo,
    Total_Mutations     = nrow(per_mut),
    Mean_Lifespan       = mean(per_mut$lifespan,  na.rm = TRUE),
    Median_Lifespan     = median(per_mut$lifespan, na.rm = TRUE),
    SE_Lifespan         = sd(per_mut$lifespan, na.rm = TRUE) / sqrt(nrow(per_mut)),
    N_Fixed             = sum(per_mut$status == "Fixed"),
    N_Lost              = sum(per_mut$status == "Lost"),
    N_Segregating       = sum(per_mut$status == "Segregating"),
    Pct_Fixed           = 100 * mean(per_mut$status == "Fixed"),
    Pct_Lost            = 100 * mean(per_mut$status == "Lost"),
    Mean_Lifespan_Fixed = mean(per_mut$lifespan[per_mut$status == "Fixed"], na.rm = TRUE),
    Mean_Lifespan_Lost  = mean(per_mut$lifespan[per_mut$status == "Lost"],  na.rm = TRUE)
  )
}

# ============================================================
# 2. LOOP INCREMENTAL (con reanudación)
# ============================================================
# Si el CSV ya existe, saltamos los archivos ya procesados
if (file.exists(archivo_salida)) {
  ya_procesados <- fread(archivo_salida, select = "File")$File
  rutas_archivos <- rutas_archivos[!basename(rutas_archivos) %in% ya_procesados]
  cat("Reanudando. Faltan:", length(rutas_archivos), "archivos.\n")
}

for (i in seq_along(rutas_archivos)) {
  res <- tryCatch(
    procesar_archivo_slim(rutas_archivos[i]),
    error = function(e) {
      message("Error en ", basename(rutas_archivos[i]), ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(res)) {
    # Escribir INMEDIATAMENTE en disco (append)
    existe <- file.exists(archivo_salida)
    fwrite(res, archivo_salida, append = existe, col.names = !existe)
  }

  # Liberar memoria cada 50 archivos
  if (i %% 50 == 0) {
    gc(verbose = FALSE)
    message(sprintf("Procesados %d / %d", i, length(rutas_archivos)))
  }
}

# ============================================================
# 3. CARGAR EL RESULTADO FINAL (ya ligero en RAM)
# ============================================================
df_raw <- fread(archivo_salida) %>% as_tibble()
cat("Filas totales:", nrow(df_raw), "\n")
print(head(df_raw))