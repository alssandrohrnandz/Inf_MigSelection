# ============================================================
# 0. CONFIGURACIÓN GLOBAL
# ============================================================
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)

base_dir      <- "data/results_Discrete/outputs_LL/independent_loci"
dir_neutros   <- file.path("data/results_Discrete/outputs_LL/independent_loci/neutros/")
dir_seleccion <- file.path("data/results_Discrete/outputs_LL/independent_loci/seleccion")

replicas_per_val <- 50
FACTOR_EDP_TO_WF <- 2   # Wright-Fisher -> Classic: s se multiplica por 2

mis_colores <- c("Classic" = "#004488", "B-B" = "#BB5566", "BwS" = "#DDAA33")

# Crear directorios de salida si no existen
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# >>> REVISAR: vectores del sweep (ORDEN EXACTO DEL BASH)
# ------------------------------------------------------------
mig_values_neutros <- c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.0)
sel_values_neutros <- c(0.0)

mig_values_seleccion <- c(0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)
sel_values_seleccion <- c(0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.0)

# ============================================================
# 1. FUNCIONES
# ============================================================

# --- 1.1 Leer un archivo y devolver MLE ---------------------
procesar_archivo <- function(ruta) {
  nombre_archivo <- basename(ruta)
  
  df <- read_tsv(ruta, show_col_types = FALSE)
  cols_D <- names(df)[str_detect(names(df), "^D_")]
  if (length(cols_D) == 0) return(NULL)
  
  # Descartar SNPs con algún NA en las columnas D_
  n_SNPs_original <- nrow(df)
  df <- df %>% filter(if_all(all_of(cols_D), ~ !is.na(.x)))
  n_SNPs        <- nrow(df)
  n_descartados <- n_SNPs_original - n_SNPs
  
  if (n_SNPs == 0) {
    warning("Todos los SNPs de ", nombre_archivo, " tienen NA. Se omite.")
    return(NULL)
  }
  
  sumas <- colSums(as.matrix(df[, cols_D, drop = FALSE]))
  
  idx_mle        <- which.max(sumas)
  grid_param_mle <- cols_D[idx_mle]
  mle_valor      <- unname(sumas[idx_mle])
  
  D_val <- str_extract(grid_param_mle, "(?<=^D_)[^_]+")
  s_val <- str_extract(grid_param_mle, "(?<=_s_).*$")
  
  modelo <- case_when(
    str_detect(nombre_archivo, "TRON_ARES")   ~ "B-B",
    str_detect(nombre_archivo, "TRON_SPIKES") ~ "BwS",
    TRUE                                      ~ "Classic"
  )
  task_id <- as.integer(str_extract(nombre_archivo, "(?<=TaskID_)\\d+"))
  
  # Extraer el tipo de filtrado del nombre del archivo
  filtro <- case_when(
    str_detect(nombre_archivo, "_all\\.(txt|csv)$")          ~ "all",
    str_detect(nombre_archivo, "_nuevos_400\\.(txt|csv)$")   ~ "nuevos_400",
    str_detect(nombre_archivo, "_vigentes_400\\.(txt|csv)$") ~ "vigentes_400",
    TRUE                                                     ~ "unknown"
  )
  
  tibble(
    File               = nombre_archivo,
    Model              = modelo,
    Filtro             = filtro,
    TaskID             = task_id,
    Inferred_D         = as.numeric(D_val),
    Inferred_s         = as.numeric(s_val),
    Grid_Param_MLE     = grid_param_mle,
    MLE_Sum            = mle_valor,
    n_SNPs             = n_SNPs,
    n_SNPs_descartados = n_descartados
  )
}

# --- 1.2 Mapear TaskID -> (True_Mig, True_Sel) -------------
mapear_taskid <- function(tabla, mig_values, sel_values, replicas = replicas_per_val, etiqueta = "") {
  n_esperado <- length(mig_values) * length(sel_values) * replicas
  max_task   <- max(tabla$TaskID, na.rm = TRUE)
  if (n_esperado != max_task) {
    warning(sprintf("[%s] n esperado (%d) != max TaskID (%d). Revisa vectores de sweep.", etiqueta, n_esperado, max_task))
  }
  
  tabla %>%
    mutate(
      idx      = floor((TaskID - 1) / replicas),
      idx_mig  = floor(idx / length(sel_values)) + 1,
      idx_sel  = (idx %% length(sel_values)) + 1,
      True_Mig = mig_values[idx_mig],
      True_Sel = sel_values[idx_sel],
      Replica  = (TaskID - 1) %% replicas + 1
    ) %>%
    select(-idx, -idx_mig, -idx_sel)
}

# --- 1.3 Pipeline completo para un directorio --------------
procesar_directorio <- function(directorio, patron, mig_values, sel_values, etiqueta = "") {
  rutas <- list.files(directorio, pattern = patron, full.names = TRUE)
  rutas <- rutas[!grepl("^P_TRON|^P_", basename(rutas))]
  
  n <- length(rutas)
  message(sprintf("[%s] %d archivos encontrados en %s", etiqueta, n, directorio))
  
  resultados <- vector("list", n)
  for (i in seq_len(n)) {
    resultados[[i]] <- tryCatch(
      procesar_archivo(rutas[i]),
      error = function(e) { message("Error en ", basename(rutas[i]), ": ", conditionMessage(e)); NULL }
    )
    if (i %% 100 == 0) { gc(verbose = FALSE); message(sprintf("  %d / %d", i, n)) }
  }
  
  tabla <- bind_rows(resultados)
  tabla <- mapear_taskid(tabla, mig_values, sel_values, etiqueta = etiqueta)
  
  tabla %>%
    mutate(
      Inferred_s_scaled = Inferred_s * FACTOR_EDP_TO_WF,
      Error_D           = Inferred_D - True_Mig,
      Abs_Error_D       = abs(Error_D),
      Error_s           = Inferred_s_scaled - True_Sel,
      Abs_Error_s       = abs(Error_s),
      Model             = factor(Model, levels = c("Classic", "B-B", "BwS"))
    ) %>%
    relocate(Model, Filtro, TaskID, Replica, True_Mig, True_Sel,
             Inferred_D, Inferred_s, Inferred_s_scaled,
             Error_D, Abs_Error_D, Error_s, Abs_Error_s,
             MLE_Sum, n_SNPs, Grid_Param_MLE, File)
}

# ============================================================
# 2. EJECUCIÓN: NEUTROS Y SELECCIÓN (Todos los filtros)
# ============================================================

message("=== Procesando NEUTROS ===")
# Se actualizó el patrón Regex para atrapar all, nuevos_400 y vigentes_400
tabla_neutros <- procesar_directorio(
  directorio = dir_neutros,
  patron     = "D_FULL_neutros_m1_(all|nuevos_400|vigentes_400)\\.(txt|csv)$",
  mig_values = mig_values_neutros,
  sel_values = sel_values_neutros,
  etiqueta   = "neutros"
)

message("=== Procesando SELECCION ===")
tabla_seleccion <- procesar_directorio(
  directorio = dir_seleccion,
  patron     = "D_FULL_seleccion_m1_(all|nuevos_400|vigentes_400)\\.(txt|csv)$",
  mig_values = mig_values_seleccion,
  sel_values = sel_values_seleccion,
  etiqueta   = "seleccion"
)

# Guardar bases de datos globales procesadas
saveRDS(tabla_neutros,   file = "data/processed/tabla_neutros_all_filters.rds")
saveRDS(tabla_seleccion, file = "data/processed/tabla_seleccion_all_filters.rds")
message("Tablas guardadas exitosamente en formato .rds en data/processed/")

tabla_neutros <- readRDS("data/processed/tabla_neutros_all_filters.rds")
tabla_seleccion <- readRDS("data/processed/tabla_seleccion_all_filters.rds")


# ============================================================
# 3. BUCLE DE PLOTS POR TIPO DE FILTRO
# ============================================================

filtros_disponibles <- unique(c(tabla_neutros$Filtro, tabla_seleccion$Filtro))
filtros_disponibles <- filtros_disponibles[!is.na(filtros_disponibles) & filtros_disponibles != "unknown"]

mis_breaks <- c(0.0, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1)

for (f in filtros_disponibles) {
  message(sprintf("=== Generando gráficas para el filtro: %s ===", f))
  
  # ----------------------------------------------------------
  # A. PLOTS: NEUTROS
  # ----------------------------------------------------------
  df_neutros <- tabla_neutros %>%
    filter(Filtro == f, !is.na(Inferred_D), !is.na(True_Mig)) %>%
    filter(!True_Mig %in% c(0.025, 0.075, 0.125))
  
  if(nrow(df_neutros) > 0) {
    plot_neutros_recovery <- ggplot(df_neutros, aes(x = True_Mig, y = Inferred_D, color = Model)) +
      geom_jitter(width = 0, height = 0.05, size = 2.5, alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) +
      facet_wrap(~ Model) +
      scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-5, base = 10), breaks = mis_breaks) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5, base = 10), breaks = mis_breaks) +
      scale_color_manual(values = mis_colores) +
      labs(
        title = sprintf("Parameter Recovery: Diffusion (D) - Neutral [%s]", f),
        subtitle = "Dashed line = perfect inference. Axes in pseudo-log scale.",
        x = "Simulated Diffusion Rate (True D)",
        y = "Estimated Diffusion Rate (Inferred D)"
      ) +
      theme_bw(base_size = 14) +
      theme(legend.position = "none", panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "gray90", color = "black"),
            strip.text = element_text(face = "bold"))
    
    ggsave(sprintf("figures/Recovery_Neutral_IndLoci_%s.png", f), plot = plot_neutros_recovery, width = 14, height = 6, dpi = 300)
    
    plot_neutros_recovery_box <- ggplot(df_neutros, aes(x = factor(True_Mig), y = Inferred_D, fill = Model)) +
      geom_boxplot(alpha = 0.7, outlier.size = 1.5, outlier.alpha = 0.5, width = 0.5, color = "black") +
      facet_wrap(~ Model) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5, base = 10), breaks = mis_breaks) +
      scale_fill_manual(values = mis_colores) +
      labs(
        title = sprintf("Parameter Recovery: Diffusion (D) - Neutral [%s]", f),
        subtitle = "Distribution of inferred D. Y-axis in pseudo-log scale.",
        x = "Simulated Diffusion Rate (True D)",
        y = "Estimated Diffusion Rate (Inferred D)"
      ) +
      theme_bw(base_size = 14) +
      theme(legend.position = "none", panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_rect(fill = "gray90", color = "black"),
            strip.text = element_text(face = "bold"))
    
    ggsave(sprintf("figures/Recovery_Neutral_IndLoci_Boxplot_%s.png", f), plot = plot_neutros_recovery_box, width = 14, height = 6, dpi = 300)
  }

  df_error_D_summary <- df_neutros %>%
      group_by(Model, True_Mig) %>%
      summarise(
        N_reps        = n(),
        MAE_D         = mean(Abs_Error_D, na.rm = TRUE),
        Min_Abs_Error = min(Abs_Error_D, na.rm = TRUE),
        .groups       = "drop"
      )
    
    plot_error_rate_D <- ggplot(df_error_D_summary, aes(x = True_Mig, y = MAE_D, color = Model, group = Model)) +
      geom_line(linewidth = 1) + 
      geom_point(size = 2.5) +
      # Cinta sombreada para el límite del error mínimo
      geom_ribbon(aes(ymin = Min_Abs_Error, ymax = MAE_D, fill = Model), alpha = 0.15, color = NA) +
      scale_x_continuous(
        trans = pseudo_log_trans(sigma = 1e-5, base = 10), 
        breaks = mis_breaks
      ) +
      scale_y_continuous(
        trans = pseudo_log_trans(sigma = 1e-5, base = 10),
        labels = scales::scientific
      ) +
      scale_color_manual(values = mis_colores) +
      scale_fill_manual(values = mis_colores) +
      labs(
        title = sprintf("Mean Absolute Error (MAE) for Diffusion (D) - Neutral [%s]", f),
        subtitle = "Lines = Mean Absolute Error (MAE). Shaded area = Minimum Absolute Error.",
        x = "Simulated Diffusion Rate (True D)",
        y = "Absolute Error (D_inferred - D_true)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position  = "bottom", 
        axis.text.x      = element_text(angle = 45, hjust = 1), 
        panel.grid.minor = element_blank()
      )
    
    # Guardar la gráfica de error de D
    ggsave(
      sprintf("figures/Error_Rate_D_Neutral_%s.png", f), 
      plot = plot_error_rate_D, 
      width = 10, 
      height = 6, 
      dpi = 300
    )
  
  # ----------------------------------------------------------
  # B. PLOTS: SELECCIÓN
  # ----------------------------------------------------------
  df_seleccion <- tabla_seleccion %>%
    filter(Filtro == f, !is.na(Inferred_s_scaled), !is.na(True_Sel))
  
  if(nrow(df_seleccion) > 0) {
    rango_max      <- max(df_seleccion$True_Sel, na.rm = TRUE)
    linea_perfecta <- data.frame(True_Sel = c(0, rango_max), Inferred_s_scaled = c(0, rango_max))
    
    plot_seleccion_recovery <- ggplot(df_seleccion, aes(x = True_Sel, y = Inferred_s_scaled, color = Model)) +
      geom_point(size = 1.8, alpha = 0.55) +
      geom_line(data = linea_perfecta, aes(x = True_Sel, y = Inferred_s_scaled),
                color = "red", linetype = "dashed", linewidth = 1, inherit.aes = FALSE) +
      facet_grid(Model ~ True_Mig, labeller = labeller(True_Mig = label_both)) +
      scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-5, base = 10), breaks = c(0, 1e-4, 1e-3, 1e-2, 1e-1), labels = scales::scientific) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5, base = 10), breaks = sort(c(-1e-2, -1e-3, -1e-4, 0, 1e-4, 1e-3, 1e-2, 1e-1)), labels = scales::scientific) +
      scale_color_manual(values = mis_colores) +
      labs(
        title = sprintf("Parameter Recovery: Selection Coefficient (s) [%s]", f),
        subtitle = "Red dashed line = perfect inference (y = x). Columns = True D.",
        x = "Simulated Selection Coefficient (True s)",
        y = "Estimated Selection Coefficient (Scaled s)"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8),
            panel.grid.minor = element_blank(), strip.background = element_rect(fill = "gray90", color = "black"), strip.text = element_text(face = "bold", size = 9))
    
    ggsave(sprintf("figures/Recovery_Seleccion_IndLoci_%s.png", f), plot = plot_seleccion_recovery, width = 16, height = 9, dpi = 300)
    
    plot_seleccion_recovery_boxplot <- ggplot(df_seleccion, aes(x = factor(True_Sel), y = Inferred_s_scaled, fill = Model, color = Model)) +
      geom_boxplot(alpha = 0.6, outlier.size = 1, outlier.alpha = 0.5, width = 0.6) +
      facet_grid(Model ~ True_Mig, labeller = labeller(True_Mig = label_both)) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5, base = 10), breaks = sort(c(-1e-2, -1e-3, -1e-4, 0, 1e-4, 1e-3, 1e-2, 1e-1)), labels = scales::scientific) +
      scale_color_manual(values = mis_colores) +
      scale_fill_manual(values = mis_colores) +
      labs(
        title = sprintf("Parameter Recovery: Selection Coefficient (s) Boxplot [%s]", f),
        x = "Simulated Selection Coefficient (True s)",
        y = "Estimated Selection Coefficient (Scaled s)"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.text.y = element_text(size = 8),
            panel.grid.minor = element_blank(), strip.background = element_rect(fill = "gray90", color = "black"), strip.text = element_text(face = "bold", size = 9))
    
    ggsave(sprintf("figures/Recovery_Seleccion_IndLoci_Boxplot_%s.png", f), plot = plot_seleccion_recovery_boxplot, width = 16, height = 9, dpi = 300)
    
    # ----------------------------------------------------------
    # C. CÁLCULO Y GRÁFICO DE TASA DE ERROR
    # ----------------------------------------------------------
    df_error_summary <- df_seleccion %>%
      group_by(Model, True_Mig, True_Sel) %>%
      summarise(
        N_reps        = n(),
        MAE_s         = mean(Abs_Error_s, na.rm = TRUE),
        Min_Abs_Error = min(Abs_Error_s, na.rm = TRUE),
        RMSE_s        = sqrt(mean(Error_s^2, na.rm = TRUE)),
        Misclass_Rate = mean((True_Sel > 0 & Inferred_s_scaled <= 1e-4) | (True_Sel == 0 & Inferred_s_scaled > 1e-4), na.rm = TRUE),
        .groups       = "drop"
      )
    
    plot_error_rate_s <- ggplot(df_error_summary, aes(x = True_Sel, y = MAE_s, color = Model, group = Model)) +
      geom_line(linewidth = 1) + geom_point(size = 2.5) +
      geom_ribbon(aes(ymin = Min_Abs_Error, ymax = MAE_s, fill = Model), alpha = 0.15, color = NA) +
      facet_wrap(~ True_Mig, labeller = labeller(True_Mig = label_both), ncol = 4) +
      scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-5, base = 10), breaks = c(0, 1e-4, 1e-3, 1e-2, 1e-1), labels = scales::scientific) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5, base = 10), labels = scales::scientific) +
      scale_color_manual(values = mis_colores) +
      scale_fill_manual(values = mis_colores) +
      labs(
        title = sprintf("Minimum & Mean Absolute Error Rate for Selection (s) [%s]", f),
        subtitle = "Lines = MAE. Lower shaded bound = Minimum Absolute Error.",
        x = "Simulated Selection Coefficient (True s)",
        y = "Absolute Error |s_inferred - s_true|"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), panel.grid.minor = element_blank(), strip.background = element_rect(fill = "gray90", color = "black"), strip.text = element_text(face = "bold", size = 9))
    
    ggsave(sprintf("figures/Error_Rate_Selection_IndLoci_%s.png", f), plot = plot_error_rate_s, width = 14, height = 8, dpi = 300)
    ### ESTE ULTIMO HEATMAP NO VALE LA PENA ###
    plot_heatmap_min_error <- ggplot(df_error_summary, aes(x = factor(True_Sel), y = factor(True_Mig), fill = Min_Abs_Error)) +
      geom_tile(color = "white", linewidth = 0.5) +
      facet_wrap(~ Model) +
      scale_fill_viridis_c(option = "magma", direction = -1, name = "Min Abs Error") +
      labs(
        title = sprintf("Minimum Absolute Error Grid Across Sweep Parameters [%s]", f),
        x = "True Selection (s)",
        y = "True Migration/Diffusion (D)"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), strip.background = element_rect(fill = "gray90", color = "black"), strip.text = element_text(face = "bold"))
    
    ggsave(sprintf("figures/Heatmap_Min_Error_Selection_%s.png", f), plot = plot_heatmap_min_error, width = 14, height = 5, dpi = 300)
  }
}

# ============================================================
# 6. EXPLORACIÓN DE PERFILES DE LL (Todos los filtros)
# ============================================================

mig_values_expl <- c(0.0, 0.001, 0.1)
sel_values_expl <- c(0.0, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1)

taskids_para <- function(D_objetivo, s_objetivo, mig_values = mig_values_expl, sel_values = sel_values_expl, replicas = replicas_per_val) {
  idx_mig <- which(mig_values == D_objetivo)
  idx_sel <- which(sel_values == s_objetivo)
  if (length(idx_mig) == 0 || length(idx_sel) == 0) stop("Valor no encontrado en el sweep")
  idx <- (idx_mig - 1) * length(sel_values) + (idx_sel - 1)
  seq.int(idx * replicas + 1, (idx + 1) * replicas)
}

taskids_s001  <- taskids_para(0.001, 0.01)
taskids_s0025 <- taskids_para(0.001, 0.025)

grupos <- tibble(
  Grupo    = c("D0.001_s0.01", "D0.001_s0.025"),
  True_Mig = c(0.001,           0.001),
  True_Sel = c(0.01,            0.025)
)

rutas_todas <- list.files(
  path       = dir_seleccion,
  pattern    = "D_FULL_seleccion_m1_(all|nuevos_400|vigentes_400)\\.(txt|csv)$",
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
  ),
  Filtro = case_when(
    str_detect(File, "_all\\.(txt|csv)$")          ~ "all",
    str_detect(File, "_nuevos_400\\.(txt|csv)$")   ~ "nuevos_400",
    str_detect(File, "_vigentes_400\\.(txt|csv)$") ~ "vigentes_400",
    TRUE                                           ~ "unknown"
  )
)

rutas_sel <- info_rutas %>%
  filter(TaskID %in% c(taskids_s001, taskids_s0025)) %>%
  mutate(
    Grupo = case_when(
      TaskID %in% taskids_s001  ~ "D0.001_s0.01",
      TaskID %in% taskids_s0025 ~ "D0.001_s0.025",
      TRUE                      ~ NA_character_
    )
  )

perfilar_archivo <- function(ruta) {
  nombre <- basename(ruta)
  task   <- as.integer(str_extract(nombre, "(?<=TaskID_)\\d+"))
  model  <- case_when(str_detect(nombre, "TRON_ARES") ~ "B-B", str_detect(nombre, "TRON_SPIKES") ~ "BwS", TRUE ~ "Classic")
  
  df <- read_tsv(ruta, col_select = starts_with("D_"), show_col_types = FALSE)
  cols_D <- names(df)[str_detect(names(df), "^D_")]
  if (length(cols_D) == 0) return(NULL)
  
  df <- df %>% filter(if_all(all_of(cols_D), ~ !is.na(.x)))
  if (nrow(df) == 0) return(NULL)
  
  mat        <- as.matrix(df[, cols_D, drop = FALSE])
  suma_total <- colSums(mat)
  n_validos  <- rep(nrow(df), ncol(mat))
  
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

perfiles_lista <- vector("list", nrow(rutas_sel))
for (i in seq_len(nrow(rutas_sel))) {
  perfiles_lista[[i]] <- tryCatch(
    perfilar_archivo(rutas_sel$ruta[i]),
    error = function(e) { NULL }
  )
}

perfiles_global <- bind_rows(perfiles_lista) %>%
  left_join(info_rutas %>% select(File, Filtro), by = "File") %>%
  mutate(
    Grupo = case_when(
      TaskID %in% taskids_s001  ~ "D0.001_s0.01",
      TaskID %in% taskids_s0025 ~ "D0.001_s0.025",
      TRUE                      ~ NA_character_
    )
  ) %>%
  left_join(grupos, by = "Grupo")

for (f in filtros_disponibles) {
  perfiles <- perfiles_global %>% filter(Filtro == f)
  
  if (nrow(perfiles) > 0) {
    mle_por_archivo <- perfiles %>%
      group_by(Grupo, Model, TaskID, File, True_Mig, True_Sel) %>%
      slice_max(Suma_LL, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(Grupo, Model, TaskID, File, True_Mig, True_Sel, Inferred_D = D, Inferred_s = s, Inferred_s_scaled = Inferred_s * FACTOR_EDP_TO_WF, MLE_Suma_LL = Suma_LL, n_SNPs_val)
    
    p_A <- ggplot(perfiles, aes(x = s, y = Suma_LL, group = TaskID, color = factor(D))) +
      geom_line(linewidth = 0.3, alpha = 0.35) + facet_grid(Grupo ~ Model) +
      scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-5)) +
      labs(title = sprintf("Perfiles log-LL por grupo [%s]", f), subtitle = "Cada línea = una réplica.", x = "s (grid EDP)", y = "Suma log-LL", color = "D") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
    
    perfil_prom <- perfiles %>%
      group_by(Grupo, Model, True_Mig, True_Sel, D, s) %>%
      summarise(Suma_media = mean(Suma_LL, na.rm = TRUE), Suma_sd = sd(Suma_LL, na.rm = TRUE), .groups = "drop")
    
    p_B <- ggplot(perfil_prom, aes(x = s, y = Suma_media, color = factor(D), fill = factor(D))) +
      geom_ribbon(aes(ymin = Suma_media - Suma_sd, ymax = Suma_media + Suma_sd), alpha = 0.15, color = NA) +
      geom_line(linewidth = 0.6) + geom_point(size = 1) +
      geom_vline(data = grupos, aes(xintercept = True_Sel), linetype = "dotted", color = "black") +
      facet_grid(Grupo ~ Model) + scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-5)) +
      labs(title = sprintf("Perfil promedio log-LL (media ± 1 SD) [%s]", f), subtitle = "Línea punteada = True_Sel", x = "s", y = "Suma log-LL media") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
    
    p_C <- ggplot(mle_por_archivo, aes(x = factor(Inferred_s_scaled), fill = Model)) +
      geom_bar(position = "dodge", alpha = 0.85) + facet_wrap(~ Grupo, ncol = 2, scales = "free_x") +
      scale_fill_manual(values = mis_colores) + labs(title = sprintf("Distribución del MLE [%s]", f), x = "Inferred_s_scaled", y = "Nº réplicas") + theme_bw(base_size = 11) + theme(legend.position = "bottom")
    
    p_D <- ggplot(mle_por_archivo, aes(x = True_Sel, y = Inferred_s_scaled, color = Model)) +
      geom_jitter(width = 0.002, height = 0, size = 2, alpha = 0.6) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      facet_wrap(~ Grupo) + scale_color_manual(values = mis_colores) + scale_x_continuous(trans = pseudo_log_trans(sigma = 1e-5)) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1e-5)) + labs(title = sprintf("Recovery: Inferred_s_scaled vs True_Sel [%s]", f), x = "True_Sel", y = "Inferred_s_scaled") + theme_bw(base_size = 12) + theme(legend.position = "bottom")
    
    ggsave(sprintf("figures/explorar_A_%s.png", f), p_A, width = 12, height = 7, dpi = 150)
    ggsave(sprintf("figures/explorar_B_%s.png", f), p_B, width = 12, height = 7, dpi = 150)
    ggsave(sprintf("figures/explorar_C_%s.png", f), p_C, width = 10, height = 5, dpi = 150)
    ggsave(sprintf("figures/explorar_D_%s.png", f), p_D, width = 10, height = 5, dpi = 150)
  }
}