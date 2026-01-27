#!/usr/bin/env Rscript

# Load libraries silently
suppressPackageStartupMessages({
  library(scales)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
})

# ==========================================
# 1. Argument Parsing
# ==========================================
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
  stop("Usage: Rscript Plot_Likelihood.R <Input_Dir> <Output_Dir> <Task_ID_Regex> <Migration_Value> [Prefix]")
}

input_dir   <- args[1]
output_dir  <- args[2]
task_id_reg <- args[3]
# CORRECCIÓN: Convertir a numérico inmediatamente para evitar errores matemáticos posteriores
m_value     <- as.numeric(args[4]) 
prefix      <- if(!is.na(args[5])) args[5] else "Analysis" 

print(paste("Migration Rate (m):", m_value))
print("Regex Pattern loaded for file search.")
print(paste("Guardar en: ", output_dir))

# ==========================================
# 2. File Discovery & Grouping
# ==========================================

pattern_to_search <- paste0("^", prefix, ".*_", task_id_reg, "_SNP_.*\\.txt$")

# substr solo para el print, para no saturar la consola
print(paste("Searching pattern:", substr(pattern_to_search, 1, 50), "..."))

file_list <- list.files(
  path = input_dir,
  pattern = pattern_to_search,
  full.names = TRUE
)

num_files <- length(file_list)
print(paste("Found", num_files, "files matching the Task Group."))

# --- LÓGICA DE GRUPOS ---
TARGET_GROUPS <- 50

if (num_files < TARGET_GROUPS) {
    print(paste("ADVERTENCIA: Solo hay", num_files, "archivos. Se reduce el número de grupos."))
    actual_groups <- num_files
} else {
    actual_groups <- TARGET_GROUPS
}

if (actual_groups == 0) {
    stop("No se encontraron archivos para procesar.")
}

# 3. Barajar y Dividir
set.seed(12345) 

shuffled_files <- sample(file_list)
group_indices <- cut(seq_along(shuffled_files), breaks = actual_groups, labels = FALSE)
file_groups <- split(shuffled_files, group_indices)

print(paste("Generados", length(file_groups), "grupos Composite."))
print(paste("Promedio archivos por grupo:", round(mean(sapply(file_groups, length)), 1)))

# ==========================================
# 3. Processing
# ==========================================

process_group <- function(files, group_idx) {
  group_name <- paste("Group", group_idx)
  
  # Read all files
  raw_data_list <- lapply(files, function(f) {
    tryCatch({
      d <- read.table(f, header = TRUE)
      return(d)
    }, error = function(e) {
      warning(paste("Error reading:", f))
      return(NULL)
    })
  })
  
  full_data <- bind_rows(raw_data_list)
  if(nrow(full_data) == 0) return(NULL)
  
  # COMPOSITE LIKELIHOOD CALCULATION
  composite_surface <- full_data %>%
    group_by(D, s) %>%
    summarise(LL_sum = sum(LL, na.rm = TRUE), .groups = "drop") %>%
    mutate(Group = group_name)
  
  return(composite_surface)
}

print("Calculating Composite Likelihoods...")
group_summaries <- lapply(seq_along(file_groups), function(i) process_group(file_groups[[i]], i))
combined_data <- bind_rows(group_summaries)

if(nrow(combined_data) == 0) stop("No valid data could be loaded.")

# ==========================================
# 4. Profile Likelihood Extraction
# ==========================================

profile_LL <- combined_data %>%
  group_by(Group, D) %>%
  summarise(
    Profile_LL = max(LL_sum),
    # OPTIMIZACIÓN: which.max es ligeramente más seguro/rápido
    Best_s = s[which.max(LL_sum)], 
    .groups = "drop"
  )

profile_LL <- profile_LL %>%
  group_by(Group) %>%
  mutate(delta_LL = Profile_LL - max(Profile_LL)) %>%
  ungroup()

# Find global maximum points (MLE) per group
max_points <- profile_LL %>%
  group_by(Group) %>%
  slice_max(Profile_LL, n = 1, with_ties = FALSE) %>%
  ungroup()

median_D <- median(max_points$D, na.rm = TRUE)

# Debug prints
head(profile_LL, 15)
head(max_points, 20)

# Statistics for Caption


bias_median_log <- median(max_points$D, na.rm = TRUE) - m_value
rmse_log  <- sqrt(mean((max_points$D - m_value)^2))

bias_rmse_text <- paste0(
  "Bias = ", round(bias_median_log, 3),
  " | RMSE = ", round(rmse_log, 3)
)



# A. Profile Likelihood Curve for D (Line Plot)
p1 <- ggplot(profile_LL, aes(x = D, y = Profile_LL, color = Group)) +
  geom_line(linewidth = 0.8, alpha = 0.8) + # linewidth actualizado
  geom_point(
    data = max_points,
    aes(x = D, y = Profile_LL),
    size = 3, shape = 21, fill = "white", stroke = 1.5
  ) +
  scale_x_log10() +          
  labs(
    title = "Composite Likelihood Profile (Diffusion)",
    subtitle = paste0("Migration Rate (m): ", m_value, " | Groups: ", actual_groups),
    x = expression(paste("Diffusion Coefficient D (", m^2, "/gen)")),
    y = "Sum of Log-Likelihood",
    color = "Subset Group"
  ) +
  theme_minimal() +
  theme(legend.position = "none") 

# B. Distribution of Estimated D (Violin/Boxplot)
p2 <- ggplot(max_points, aes(x = "MLE Estimate", y = D)) + 
  geom_violin(fill = "grey95", color = "grey60") +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.05, height = 0, color = "#2c7bb6", size = 2, alpha = 0.7) +
  
  # Líneas horizontales con mapping para leyenda
  geom_hline(
    aes(yintercept = median_D, color = "Median value", linetype = "Median value"),
    linewidth = 1
  ) +
  geom_hline(
    aes(yintercept = m_value, color = "Theoric value", linetype = "Theoric value"),
    linewidth = 1
  ) +
  
  scale_y_log10(
    breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 0.5, 1),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  
  annotation_logticks(sides = "l") +
  
  # --- CORRECCIÓN DE LA LEYENDA ---
  # Definimos 'breaks' para asegurar el orden exacto de las etiquetas
  scale_color_manual(
    name = "",
    breaks = c("Theoric value", "Median value"), # Orden forzado
    values = c(
      "Theoric value" = "#d7191c",
      "Median value"  = "blue"
    ),
    labels = c(
      paste0("Theoric value = ", m_value),
      paste0("Median value = ", format(median_D, scientific = FALSE, digits = 3))
    )
  ) +
  
  scale_linetype_manual(
    name = "",
    breaks = c("Theoric value", "Median value"), # Orden forzado
    values = c(
      "Theoric value" = "dashed",
      "Median value"  = "dotted"
    ),
    labels = c(
      paste0("Theoric value = ", m_value),
      paste0("Median value = ", format(median_D, scientific = FALSE, digits = 3))
    )
  ) +
  
  labs(
    title = "Variance of D Estimation",
    subtitle = paste0("Distribution of LL across ", actual_groups, " random Composite groups (2Ne=1000)"), #TODO: Modificar 500 o 1000
    x = "",
    y = expression(paste("Estimated D (", m^2, "/gen) - Log Scale")),
    caption = bias_rmse_text
  ) +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "bottom" # Mover leyenda abajo para mejor lectura
  )

# ==========================================
# 6. Saving
# ==========================================

clean_filename <- paste0(prefix, "_Composite_N1000_Mig_", m_value) #TODO:Editar a 500 o 100

ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Profile.png")), 
       plot = p1, width = 8, height = 6)

ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Distribution.png")), 
       plot = p2, width = 6, height = 6)

print(paste("Plots saved successfully in:", output_dir))