#!/usr/bin/env Rscript

# Load libraries silently
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
})

# ==========================================
# 1. Argument Parsing & Configuration
# ==========================================
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 4) {
  stop("Usage: Rscript Plot_Likelihood.R <Input_Dir> <Output_Dir> <Task_ID_Regex> <Migration_Value> [Prefix]")
}

input_dir   <- args[1]
output_dir  <- args[2]
task_id_reg <- args[3] # Ahora esto es un REGEX: "(1|2|...|50)"
m_value     <- args[4] 
prefix      <- if(!is.na(args[5])) args[5] else "Analysis" 

print(paste("Migration Rate (m):", m_value))
print("Regex Pattern loaded for file search.")

# ==========================================
# 2. File Discovery & Loading
# ==========================================

# Pattern: Matches strictly the current prefix and the group of task_ids
pattern_to_search <- paste0("^", prefix, ".*_", task_id_reg, "_SNP_.*\\.txt$")

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

# 1. Definimos cuántos grupos haremos realmente
if (num_files < TARGET_GROUPS) {
    print(paste("ADVERTENCIA: Solo hay", num_files, "archivos. Se reduce el número de grupos."))
    actual_groups <- num_files
} else {
    actual_groups <- TARGET_GROUPS
}

# 2. Validación de seguridad
if (actual_groups == 0) {
    stop("No se encontraron archivos para procesar.")
}

# 3. Barajar y Dividir
# CORRECCIÓN: No podemos usar task_id_reg como entero para la seed.
# Usamos el valor de migración (convertido a factor numérico) o una seed fija.
set.seed(12345) 

shuffled_files <- sample(file_list)

group_indices <- cut(seq_along(shuffled_files), breaks = actual_groups, labels = FALSE)
file_groups <- split(shuffled_files, group_indices)

print(paste("Generados", length(file_groups), "grupos Composite."))
print(paste("Promedio archivos por grupo:", round(mean(sapply(file_groups, length)), 1)))

# ==========================================
# 3. Data Processing (Composite Likelihood)
# ==========================================

process_group <- function(files, group_idx) {
  group_name <- paste("Group", group_idx)
  
  # Read all files in the group
  raw_data_list <- lapply(files, function(f) {
    tryCatch({
      d <- read.table(f, header = TRUE)
      return(d)
    }, error = function(e) {
      warning(paste("Error reading:", f))
      return(NULL)
    })
  })
  
  # Combine
  full_data <- bind_rows(raw_data_list)
  if(nrow(full_data) == 0) return(NULL)
  
  # COMPOSITE LIKELIHOOD CALCULATION
  # Sum Log-Likelihoods (LL) for each combination of parameters (D and s)
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

profile_D <- combined_data %>%
  group_by(Group, D) %>%
  summarise(
    Profile_LL = max(LL_sum),      # Best LL for this D (optimizing s)
    Best_s = s[which.max(LL_sum)], # The s value associated with that max
    .groups = "drop"
  )

# Find global maximum points (MLE) per group
max_points <- profile_D %>%
  group_by(Group) %>%
  filter(Profile_LL == max(Profile_LL)) %>%
  ungroup()

# ==========================================
# 5. Plotting (English Titles)
# ==========================================

# A. Profile Likelihood Curve for D (Line Plot)
p1 <- ggplot(profile_D, aes(x = D, y = Profile_LL, color = Group)) +
  geom_line(size = 0.8, alpha = 0.8) +
  geom_point(data = max_points, aes(x = D, y = Profile_LL), 
             size = 3, shape = 21, fill = "white", stroke = 1.5) +
  labs(
    title = "Composite Likelihood Profile (Diffusion)",
    # CORRECCIÓN: Usamos m_value en el título, no el regex feo
    subtitle = paste0("Migration Rate (m): ", m_value, " | Groups: ", actual_groups),
    x = expression(paste("Diffusion Coefficient D (", m^2, "/gen)")),
    y = "Sum of Log-Likelihood",
    color = "Subset Group"
  ) +
  theme_minimal() +
  theme(legend.position = "none") # Ocultamos leyenda si son 50 grupos (mucho ruido)

# B. Distribution of Estimated D (Violin/Boxplot)
teoric_value <- as.numeric(m_value)
p2 <- ggplot(max_points, aes(x = "MLE Estimate", y = D)) +
  geom_violin(fill = "grey95", color = "grey60", draw_quantiles = 0.5) +
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.05, height = 0, color = "#2c7bb6", size = 2, alpha = 0.7) +
  geom_hline(yintercept = teoric_value, linetype = "dashed", color = "#d7191c", linewidth = 1) +
  labs(
    title = "Variance of D Estimation",
    # CORRECCIÓN: Usamos actual_groups en el subtítulo
    subtitle = paste0("Distribution of Maxima across ", actual_groups, " random Composite groups"),
    x = "",
    y = "Estimated D"
  ) +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

# ==========================================
# 6. Saving Outputs
# ==========================================

# CORRECCIÓN: Limpieza del nombre de archivo.
# No podemos usar el regex con barras '|' en el nombre del archivo.
# Usamos el valor de migración para nombrar el archivo limpiamente.
clean_filename <- paste0(prefix, "_Composite_Mig_", m_value)

ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Profile.png")), 
       plot = p1, width = 8, height = 6)

ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Distribution.png")), 
       plot = p2, width = 6, height = 6)

print(paste("Plots saved successfully in:", output_dir))