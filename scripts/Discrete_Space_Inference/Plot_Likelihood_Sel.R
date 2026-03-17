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

print(paste("Searching pattern:", pattern_to_search))

file_list <- list.files(
  path = input_dir,
  pattern = pattern_to_search,
  full.names = TRUE
)

num_files <- length(file_list)
if (num_files == 0) stop("No se encontraron archivos para procesar.")

print(paste("Found", num_files, "files matching the Task Group."))

replica_ids <- str_extract(basename(file_list), "_[0-9]+_SNP_") %>% 
               str_extract("[0-9]+") %>% 
               as.numeric()
print(paste("Replica IDs found:", paste(unique(replica_ids), collapse = ", ")))

file_groups <- split(file_list, replica_ids)

actual_groups <- length(file_groups)
print(paste("Se identificaron", actual_groups, "simulaciones/réplicas distintas."))
print(paste("Promedio de SNPs por réplica:", round(mean(sapply(file_groups, length)), 1)))

# ==========================================
# 3. Data Processing (Composite & Single-SNP MLE)
# ==========================================

process_group <- function(files, rep_id) {
  group_name <- paste("Replica", rep_id)
  
  # 1. Leer archivos y ETIQUETAR de qué SNP/archivo viene cada dato
  raw_data_list <- lapply(files, function(f) {
    tryCatch({
      d <- read.table(f, header = TRUE)
      # Extraemos el número de SNP del nombre del archivo (ej. "..._SNP_125.txt" -> "125")
      snp_num <- stringr::str_extract(basename(f), "(?<=_SNP_)[0-9]+")
      d$SNP_ID <- ifelse(is.na(snp_num), basename(f), snp_num) # Fallback al nombre completo si falla
      return(d)
    }, error = function(e) {
      warning(paste("Error reading:", f))
      return(NULL)
    })
  })
  
  # Combinar todo
  full_data <- bind_rows(raw_data_list)
  if(nrow(full_data) == 0) return(NULL)
  
  # A) CÁLCULO COMPOSITE LIKELIHOOD (Para toda la réplica junta)
  composite_surface <- full_data %>%
    group_by(D, s) %>%
    summarise(LL_sum = sum(LL, na.rm = TRUE), .groups = "drop") %>%
    mutate(Group = group_name, Replica_ID = rep_id)
  
  # B) EXTRACCIÓN SINGLE-SNP MLE (El máximo absoluto de cada archivo individual)
  single_snp_mle <- full_data %>%
    group_by(SNP_ID) %>%
    # Nos quedamos con la fila que tenga el LL más alto para este SNP
    slice_max(order_by = LL, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Group = group_name, Replica_ID = rep_id)
  
  # Devolvemos ambos en forma de lista
  return(list(composite = composite_surface, singles = single_snp_mle))
}

print("Calculating Composite Likelihoods & Extracting Single-SNP MLEs...")

# Iteramos sobre las réplicas
group_summaries <- lapply(names(file_groups), function(id) {
  process_group(file_groups[[id]], id)
})

# Como ahora es una lista doble, separamos y unimos los dataframes correspondientes
combined_data <- bind_rows(lapply(group_summaries, `[[`, "composite"))
all_singles_mle <- bind_rows(lapply(group_summaries, `[[`, "singles"))

if(nrow(combined_data) == 0) stop("No valid data could be loaded.")

# ==========================================
# 4. Extracción de Máximos (MLE por Réplica)
# ==========================================
print("Extrayendo MLE (Maximum Likelihood Estimate) por cada réplica...")

# Para cada réplica, buscamos la fila exacta que tiene el LL_sum más alto
mle_per_replica <- combined_data %>%
  group_by(Replica_ID) %>%
  slice_max(order_by = LL_sum, n = 1, with_ties = FALSE) %>%
  ungroup()

print(paste("Se extrajeron", nrow(mle_per_replica), "puntos máximos (uno por simulación)."))

# ==========================================
# 5. Generación de Gráficos de Distribución
# ==========================================
print("Generando gráficos...")

# Definimos una paleta coherente
color_theoric <- "#d7191c"  # Rojo para el valor real
color_puntos_D <- "#3498db" # Azul para D
color_puntos_s <- "#9b59b6" # Morado para s

# --- GRÁFICO A: Distribución de D ---
p_dist_D <- ggplot(mle_per_replica, aes(x = "MLE Estimates", y = D)) +
  geom_violin(fill = "grey95", color = "grey60", alpha = 0.5) +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, color = color_puntos_D, size = 3, alpha = 0.7) +
  geom_hline(aes(yintercept = m_value, color = "Theoretical"), linetype = "dashed", linewidth = 1) +
  scale_y_log10(
    breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 0.5, 1),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks(sides = "l") +
  scale_color_manual(name = "", values = c("Theoretical" = color_theoric), 
                     labels = paste0("Theoric D = ", m_value)) +
  labs(
    title = "Variance of Diffusion (D) Estimation",
    subtitle = paste0("Distribution of MLEs across ", actual_groups, " independent simulations\nFixed s parameter = ", s_value),
    x = "",
    y = expression(paste("Estimated D (", m^2, "/gen) - Log Scale"))
  ) +
  theme_bw() + theme(legend.position = "bottom", axis.ticks.x = element_blank())

# --- GRÁFICO B: Distribución de s ---
p_dist_s <- ggplot(mle_per_replica, aes(x = "MLE Estimates", y = s)) +
  geom_violin(fill = "grey95", color = "grey60", alpha = 0.5) +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0, color = color_puntos_s, size = 3, alpha = 0.7) +
  geom_hline(aes(yintercept = s_value, color = "Theoretical"), linetype = "dashed", linewidth = 1) +
  scale_color_manual(name = "", values = c("Theoretical" = color_theoric), 
                     labels = paste0("Theoric s = ", s_value)) +
  labs(
    title = "Variance of Selection (s) Estimation",
    subtitle = paste0("Distribution of MLEs across ", actual_groups, " independent simulations\nFixed D parameter = ", m_value),
    x = "",
    y = "Estimated Selection Coefficient (s)"
  ) +
  theme_bw() + theme(legend.position = "bottom", axis.ticks.x = element_blank())

# --- GRÁFICO C: Diagnóstico de SNPs Individuales ---
total_snps <- nrow(all_singles_mle)

p_diag_snps <- ggplot(all_singles_mle, aes(x = s, y = LL)) +
  geom_jitter(aes(fill = as.factor(Replica_ID)), width = 0.002, size = 2.5, shape = 21, color = "black", alpha = 0.6) +
  geom_vline(xintercept = s_value, linetype = "dashed", color = "red", linewidth = 1) +
  labs(
    title = "Single-Locus MLE Diagnostics: LL vs Inferred Selection (s)",
    subtitle = paste0("Showing max Likelihood points for ", total_snps, " individual SNPs across all replicas.\nDashed line = Theoretical s (", s_value, ")"),
    x = "Inferred Selection Coefficient (s)",
    y = "Maximum Log-Likelihood (LL)",
    fill = "Replica ID"
  ) +
  theme_bw() +
  theme(legend.position = "right")

# ==========================================
# 6. Guardado de Resultados
# ==========================================

# Definimos el nombre PRIMERO
clean_filename <- paste0(prefix, "_Performance_Mig_", m_value, "_Sel_", s_value)

# Guardamos los 3 gráficos
ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Dist_D.png")), 
       plot = p_dist_D, width = 6, height = 6)

ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Dist_s.png")), 
       plot = p_dist_s, width = 6, height = 6)

ggsave(filename = file.path(output_dir, paste0(clean_filename, "_SingleSNP_Diagnostics.png")), 
       plot = p_diag_snps, width = 8, height = 6)

print(paste("Gráficos guardados exitosamente en:", output_dir))