library(tidyverse)
library(data.table)
library(scales)

# ============================================================
# 0. CARGAR TODOS LOS CHUNKS Y COMBINAR
# ============================================================
archivos_resumen <- Sys.glob("slim_lifespan_summary_chunk*.csv")
if (length(archivos_resumen) == 0) archivos_resumen <- "slim_lifespan_summary.csv"

df_raw <- rbindlist(lapply(archivos_resumen, fread), fill = TRUE) %>% as_tibble()
cat("Filas cargadas:", nrow(df_raw), "\n")

# ============================================================
# 1. MAPEO TaskID -> (True_Mig, True_Sel)
# ============================================================
mig_values      <- c(0.0001, 0.0005, 0.001, 0.005, 0.01,
                     0.025, 0.05, 0.075, 0.1, 0.125)
sel_values_sel  <- c(0.1, 0.075, 0.05, 0.025, 0.01, 0.0075, 0.0050, 0.0025,
                     0.001, 0.00075, 0.0005, 0.00025, 0.0001, 0.0)
sel_values_neu  <- c(0.0)
replicas_per_val <- 50

mapear_params <- function(task_id, tipo) {
  sel_values <- if (tipo == "Neutro") sel_values_neu else sel_values_sel
  n_sel <- length(sel_values)
  idx      <- floor((task_id - 1) / replicas_per_val)
  idx_mig  <- floor(idx / n_sel) + 1
  idx_sel  <- (idx %% n_sel) + 1
  tibble(
    True_Mig = mig_values[idx_mig],
    True_Sel = sel_values[idx_sel],
    Replica  = (task_id - 1) %% replicas_per_val + 1
  )
}

df_params <- df_raw %>%
  rowwise() %>%
  mutate(params = list(mapear_params(TaskID, Type))) %>%
  unnest(params) %>%
  ungroup()

# --- Verificación rápida ---
cat("\n=== Chequeo de mapeo ===\n")
cat("TaskID=1    -> Mig=0.0001, Sel=0.1 (Seleccion)\n")
cat("TaskID=51   -> Mig=0.0001, Sel=0.075\n")
cat("TaskID=701  -> Mig=0.0005, Sel=0.1\n")
cat("TaskID=1 (Neutro) -> Mig=0.0001, Sel=0.0\n\n")

print(df_params %>% filter(Type == "Seleccion") %>%
        slice(c(1, 51, 701)) %>%
        select(TaskID, Type, True_Mig, True_Sel, Replica))

# ============================================================
# 2. FILTRAR Y PREPARAR DATOS PARA PLOTS
# ============================================================
# Ratio Lifespan / (1/s) -- escala de referencia "tiempo de fijación esperado"
df_plot <- df_params %>%
  filter(!is.na(True_Sel), !is.na(True_Mig), !is.na(Mean_Lifespan)) %>%
  mutate(
    Type     = factor(Type, levels = c("Seleccion", "Neutro")),
    True_Mig = factor(True_Mig, levels = sort(unique(True_Mig))),
    Sel_Label = ifelse(Type == "Neutro", "Neutro",
                       sprintf("%.4g", True_Sel)),
    # Escala relativa: lifespan en unidades de 1/s (tiempo de fijación teórico)
    Lifespan_over_1s = ifelse(True_Sel > 0, Mean_Lifespan * True_Sel, NA_real_)
  )

# ============================================================
# 3. RESUMEN: media y SE agregada por combinación
# ============================================================
resumen <- df_plot %>%
  group_by(Type, True_Mig, True_Sel, Sel_Label) %>%
  summarise(
    n_replicas          = n(),
    Mean_Lifespan_agg   = mean(Mean_Lifespan, na.rm = TRUE),
    SD_Lifespan         = sd(Mean_Lifespan,   na.rm = TRUE),
    SE_Lifespan         = SD_Lifespan / sqrt(n_replicas),
    Mean_Lifespan_Fixed = mean(Mean_Lifespan_Fixed, na.rm = TRUE),
    Mean_Lifespan_Lost  = mean(Mean_Lifespan_Lost,  na.rm = TRUE),
    Mean_Pct_Fixed      = mean(Pct_Fixed, na.rm = TRUE),
    Mean_Pct_Lost       = mean(Pct_Lost,  na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 4. GRÁFICA 1: Lifespan promedio en función de s simulada
#               (coloreado por D, facetado por Type)
# ============================================================
p1 <- ggplot(
    resumen %>% filter(Type == "Seleccion"),
    aes(x = True_Sel, y = Mean_Lifespan_agg,
        color = True_Mig, group = True_Mig)
  ) +
  geom_ribbon(aes(ymin = Mean_Lifespan_agg - SE_Lifespan,
                  ymax = Mean_Lifespan_agg + SE_Lifespan,
                  fill = True_Mig),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_x_continuous(
    trans  = pseudo_log_trans(sigma = 1e-5, base = 10),
    breaks = sort(unique(df_plot$True_Sel)),
    labels = scales::label_number(accuracy = 1e-5)
  ) +
  scale_y_log10() +
  scale_color_viridis_d(option = "plasma", name = "True_Mig") +
  scale_fill_viridis_d(option  = "plasma", name = "True_Mig") +
  labs(
    title    = "Lifespan promedio de alelos por coeficiente de selección",
    subtitle = "Media ± SE sobre 50 réplicas. Bandas = error estándar. Eje Y en escala log.",
    x        = "Coeficiente de selección simulado (True_Sel)",
    y        = "Lifespan promedio (generaciones)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

print(p1)
ggsave("lifespan_vs_seleccion.png", p1, width = 12, height = 6, dpi = 300)

# ============================================================
# 5. GRÁFICA 2: Distribución del lifespan (boxplot) por s
#               (facetado por D)
# ============================================================
p2 <- ggplot(
    df_plot %>% filter(Type == "Seleccion"),
    aes(x = factor(True_Sel), y = Mean_Lifespan,
        fill = factor(True_Sel))
  ) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4,
               linewidth = 0.3, alpha = 0.85) +
  facet_grid(~ True_Mig, labeller = labeller(True_Mig = label_both)) +
  scale_y_log10() +
  scale_fill_viridis_d(guide = "none") +
  labs(
    title    = "Distribución del lifespan por True_Sel y True_Mig",
    subtitle = "Cada caja resume 50 réplicas. Eje Y en escala log.",
    x        = "Coeficiente de selección simulado (True_Sel)",
    y        = "Lifespan promedio por archivo (generaciones)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 9)
  )

print(p2)
ggsave("lifespan_boxplot_vs_seleccion.png", p2, width = 16, height = 5, dpi = 300)

# ============================================================
# 6. (OPCIONAL) Añadir comparación con neutros
# ============================================================
# Línea base: lifespan promedio de los neutros para cada D
neutros_base <- resumen %>%
  filter(Type == "Neutro") %>%
  select(True_Mig, Lifespan_Neutro = Mean_Lifespan_agg)

p3 <- ggplot(
    resumen %>% filter(Type == "Seleccion"),
    aes(x = True_Sel, y = Mean_Lifespan_agg,
        color = True_Mig, group = True_Mig)
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_hline(data = neutros_base,
             aes(yintercept = Lifespan_Neutro, color = True_Mig),
             linetype = "dashed", alpha = 0.5, linewidth = 0.4) +
  scale_x_continuous(
    trans  = pseudo_log_trans(sigma = 1e-5, base = 10),
    breaks = sort(unique(df_plot$True_Sel)),
    labels = scales::label_number(accuracy = 1e-5)
  ) +
  scale_y_log10() +
  scale_color_viridis_d(option = "plasma", name = "True_Mig") +
  labs(
    title    = "Lifespan bajo selección vs línea base neutra",
    subtitle = "Líneas punteadas = lifespan esperado bajo neutralidad (Type = Neutro) para cada D.",
    x        = "Coeficiente de selección simulado (True_Sel)",
    y        = "Lifespan promedio (generaciones)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(p3)
ggsave("lifespan_vs_neutros.png", p3, width = 12, height = 6, dpi = 300)

# ============================================================
# 7. Guardar tabla resumen
# ============================================================
write_csv(resumen, "lifespan_resumen_por_parametros.csv")
cat("\nResumen guardado en 'lifespan_resumen_por_parametros.csv'\n")

# Chequeo rápido
cat("\n=== Filas por combinación (debe ser 50) ===\n")
print(resumen %>% count(Type, True_Mig, True_Sel) %>%
        filter(n != 50) %>% head())