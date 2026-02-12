#!/usr/bin/env Rscript

# Carga silenciosa de librerías
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(viridis) # Para escalas de color científicas
})

# ==========================================
# 1. Lectura de Argumentos
# ==========================================
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
  stop("Uso: Rscript Plot_Joint_Likelihood.R <Input_Dir> <Output_Dir> <Regex_IDs> <Valor_Migracion> [Prefijo]")
}

input_dir   <- args[1]
output_dir  <- args[2]
task_id_reg <- args[3]
m_value     <- as.numeric(args[4])
prefix      <- if(!is.na(args[5])) args[5] else "Analysis"
s_value     <- as.numeric(args[6])

print(paste("Analizando Superficie Conjunta (D vs s). Migración teórica:", m_value))
print(paste("Analizando Superficie Conjunta (D vs s). Selección teórica:", s_value))

# ==========================================
# 2. Carga y Composite Likelihood
# ==========================================

pattern_to_search <- paste0("^", prefix, ".*_", task_id_reg, "_SNP_.*\\.txt$")
file_list <- list.files(path = input_dir, pattern = pattern_to_search, full.names = TRUE)

if(length(file_list) == 0) stop("No se encontraron archivos.")

print(paste("Procesando", length(file_list), "archivos para generar la superficie..."))

# Función de lectura robusta
read_and_process <- function(f) {
  tryCatch({
    d <- read.table(f, header = TRUE)
    return(d) # Debe tener columnas D, s, LL
  }, error = function(e) return(NULL))
}

# Leemos y combinamos todo
raw_data_list <- lapply(file_list, read_and_process)
full_data <- bind_rows(raw_data_list)

if(nrow(full_data) == 0) stop("Datos vacíos.")

# --- CÁLCULO DEL COMPOSITE LIKELIHOOD ---
# Sumamos el LL para cada combinación única de D y s
composite_surface <- full_data %>%
  group_by(D, s) %>%
  summarise(LL_sum = sum(LL, na.rm = TRUE), .groups = "drop")

# ==========================================
# 3. Encontrar el Máximo (MLE)
# ==========================================

# El punto con el LL más alto en toda la superficie
mle_point <- composite_surface %>%
  filter(LL_sum == max(LL_sum)) %>%
  slice(1)

best_D <- mle_point$D[1]
best_s <- mle_point$s[1]
max_LL <- mle_point$LL_sum[1]

inferred_D <- best_D
inferred_s <- best_s
max_LL_val <- max_LL

print(paste("MLE Encontrado -> D:", best_D, "| s:", best_s))

# Creamos una columna de Delta Likelihood para mejorar el color del plot
# (Hacemos que el máximo sea 0 y el resto negativo)
composite_surface <- composite_surface %>%
  mutate(Delta_LL = LL_sum - max_LL)

# Filtramos valores extremadamente bajos para que el gráfico no se vea plano
# Solo mostramos valores dentro de un rango razonable del pico (ej. los top 100 log-units)
limit_LL <- -100 
plot_data <- composite_surface %>%
  filter(Delta_LL > limit_LL)

# ==========================================
# 4. Generación de Gráficos
# ==========================================

# --- GRÁFICO 1: HEATMAP CONJUNTO (D vs s) ---
# Este es el gráfico más importante para ver la correlación

p1 <- ggplot(plot_data, aes(x = D, y = s, fill = Delta_LL)) +
  geom_tile() + # Dibuja la cuadrícula
  scale_x_log10(
    breaks = unique(composite_surface$D), # Usamos los cortes de tu grid
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_fill_viridis(option = "magma", name = expression(Delta*LL)) +
  
  # Marcamos el punto máximo con una estrella o punto blanco
  geom_point(data = mle_point, aes(x = D, y = s), shape = 21, fill = "white", color = "black", size = 3, stroke = 1.5) +
  
  # Líneas teóricas (si quieres ver dónde cae la verdad)
  geom_vline(xintercept = m_value, linetype = "dashed", color = "cyan", alpha=0.5) +
  # Si conoces el s teórico, podrías agregar geom_hline aquí también
  
  labs(
    title = paste0("Joint Likelihood Surface (D= ",m_value," vs s= ",s_value,")"),
    subtitle = paste0("MLE: D=", best_D, ", s=", best_s),
    x = expression(paste("Diffusion Coefficient ", italic(D))),
    y = expression(paste("Selection Coefficient ", italic(s)))
  ) +
  theme_dark() + # Tema oscuro resalta mejor el heatmap
  theme(panel.grid = element_blank())

# --- GRÁFICO 2: PERFIL DE SELECCIÓN (s) ---
# Marginalizamos D (tomamos el mejor D para cada s)
profile_s <- composite_surface %>%
  group_by(s) %>%
  summarise(Profile_LL = max(LL_sum), .groups = "drop")

p2 <- ggplot(profile_s, aes(x = s, y = Profile_LL)) +
  geom_vline(aes(xintercept = inferred_s, color = "Inferred MLE", linetype = "Inferred MLE"), linewidth = 1) +
  geom_vline(aes(xintercept = s_value, color = "Theoretical", linetype = "Theoretical"), linewidth = 1) +
  geom_line(color = "#e74c3c", size = 1) +
  geom_point(data = filter(profile_s, Profile_LL == max(Profile_LL)), aes(x=s, y=Profile_LL), size=3) +
  labs(
    title = paste0("Profile Likelihood for Selection (s= ",s_value,")"),
    subtitle = paste0("Marginalized over Diffusion (D= ",m_value ,")"),
    x = "Selection Coefficient (s)",
    y = "Log-Likelihood"
  ) +
  theme_bw()

# --- GRÁFICO 3: PERFIL DE DIFUSIÓN (D) ---
# Marginalizamos s (tomamos el mejor s para cada D)
profile_D <- composite_surface %>%
  group_by(D) %>%
  summarise(Profile_LL = max(LL_sum)) 

p3 <- ggplot(profile_D, aes(x = D, y = Profile_LL)) +
  geom_line(color = "#3498db", size = 1) +
  scale_x_log10() +
  geom_point(data = filter(profile_s, Profile_LL == max(Profile_LL)), aes(x=D, y=Profile_LL), size=3) +
  labs(
    title = paste0("Profile Likelihood for Diffusion (D= ",m_value,")"),
    subtitle = paste0("Marginalized over Selection (s= ",s_value,")"),
    x = "Diffusion Coefficient (D)",
    y = "Log-Likelihood"
  ) +
  theme_bw()

# ==========================================
# 5. Guardado
# ==========================================
filename_base <- paste0(prefix, "_Joint_Analysis_Mig_", m_value,"_Sel",s_value)

ggsave(filename = file.path(output_dir, paste0(filename_base, "_Heatmap.png")), plot = p1, width = 8, height = 6)
ggsave(filename = file.path(output_dir, paste0(filename_base, "_Profile_s.png")), plot = p2, width = 6, height = 4)
ggsave(filename = file.path(output_dir, paste0(filename_base, "_Profile_D.png")), plot = p3, width = 6, height = 4)

print(paste("Gráficos guardados en:", output_dir))