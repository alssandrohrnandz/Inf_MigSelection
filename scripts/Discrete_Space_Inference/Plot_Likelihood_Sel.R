library(dplyr)
library(stringr)
library(ggplot2)

# ==========================================
# 1. Argumentos y Setup (Tu código original)
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
s_value    <- as.numeric(args[6]) 
s_value <- s_value/2
print(paste("Analizando Superficie Conjunta (D vs s). Migración teórica:", m_value))
print(paste("Analizando Superficie Conjunta (D vs s). Selección teórica:", s_value , "(Novembre*2)"))

# ==========================================
# 2. Carga de Archivos
# ==========================================
pattern_to_search <- paste0("^", prefix, "_TaskID_", task_id_reg, "_All_SNPs.txt")
print(pattern_to_search)
file_list <- list.files(path = input_dir, pattern = pattern_to_search, full.names = TRUE)

num_files <- length(file_list)
if (num_files == 0) stop("No se encontraron archivos para procesar.")

print(paste("Found", num_files, "files matching the Task Group."))

replica_ids <- str_extract(basename(file_list), "(?<=TaskID_)[0-9]+") %>% as.numeric()
print(paste("Replica IDs found:", paste(unique(replica_ids), collapse = ", ")))

# ==========================================
# 3. Cálculo MLE (Máxima Verosimilitud por Réplica)
# ==========================================
print("Calculando sumas de verosimilitud y extrayendo el punto máximo por archivo...")

process_file <- function(f, rep_id) {
  # check.names = FALSE evita que R reemplace caracteres raros en los nombres de columnas
  d <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  # Separar la columna SNP de los valores numéricos de Log-Likelihood
  ll_matrix <- d[, colnames(d) != "SNP", drop = FALSE]
  
  # A) CÁLCULO COMPOSITE LIKELIHOOD (Para toda la réplica)
  # Sumamos todas las filas (SNPs) para cada columna (combinación de D y s)
  ll_sums <- colSums(ll_matrix, na.rm = TRUE)
  
  # Encontrar el nombre de la columna que dio la suma más alta
  best_col <- names(which.max(ll_sums))
  best_LL_sum <- max(ll_sums, na.rm = TRUE)
  
  # Extraer D y s del nombre de la columna (ej. "D_0.01_s_0.05")
  best_D <- as.numeric(str_extract(best_col, "(?<=D_)[0-9.eE+-]+"))
  best_s <- as.numeric(str_extract(best_col, "(?<=s_)[0-9.eE+-]+"))
  
  composite_mle <- data.frame(
    Replica_ID = rep_id,
    D = best_D,
    s = best_s,
    LL_sum = best_LL_sum
  )
  
  # B) EXTRACCIÓN SINGLE-SNP MLE (Opcional para el Gráfico C)
  # Identificamos la mejor columna para cada fila individualmente
  max_idx <- max.col(ll_matrix, ties.method = "first")
  best_cols_snps <- colnames(ll_matrix)[max_idx]
  
  single_snp_mle <- data.frame(
    SNP_ID = d$SNP,
    Replica_ID = rep_id,
    D = as.numeric(str_extract(best_cols_snps, "(?<=D_)[0-9.eE+-]+")),
    s = as.numeric(str_extract(best_cols_snps, "(?<=s_)[0-9.eE+-]+")),
    LL = apply(ll_matrix, 1, max, na.rm = TRUE)
  )
  
  return(list(composite = composite_mle, singles = single_snp_mle))
}

# Procesar todos los archivos de la lista
results_list <- lapply(seq_along(file_list), function(i) {
  process_file(file_list[i], replica_ids[i])
})

# Unir resultados en DataFrames
mle_per_replica <- bind_rows(lapply(results_list, `[[`, "composite"))
all_singles_mle <- bind_rows(lapply(results_list, `[[`, "singles"))

actual_groups <- nrow(mle_per_replica)
print(paste("Se extrajeron", actual_groups, "puntos máximos (uno por archivo)."))

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
margen_s <- 0.05

p_diag_snps <- ggplot(all_singles_mle, aes(x = s, y = LL)) +
  geom_jitter(aes(fill = as.factor(Replica_ID)), width = 0.002, size = 2.5, shape = 21, color = "black", alpha = 0.6) +
  geom_vline(xintercept = s_value, linetype = "dashed", color = "red", linewidth = 1) +
  # Agregamos el zoom aquí:
  coord_cartesian(xlim = c(s_value - margen_s, s_value + margen_s)) +
  labs(
    title = "Single-Locus MLE Diagnostics: LL vs Inferred Selection (s)",
    subtitle = paste0("Zoom cerca del valor teórico. Límite: +/- ", margen_s),
    x = "Inferred Selection Coefficient (s)",
    y = "Maximum Log-Likelihood (LL)",
    fill = "Replica ID"
  ) +
  theme_bw() +
  theme(legend.position = "right")

# ==========================================
# 6. Guardado de Resultados
# ==========================================
clean_filename <- paste0(prefix, "_Performance_Mig_", m_value, "_Sel_", s_value)

#ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Dist_D.png")), 
#       plot = p_dist_D, width = 6, height = 6)

ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Dist_s.png")), 
       plot = p_dist_s, width = 6, height = 6)

ggsave(filename = file.path(output_dir, paste0(clean_filename, "_SingleSNP_Diagnostics.png")), 
       plot = p_diag_snps, width = 8, height = 6)

print(paste("Gráficos guardados exitosamente en:", output_dir))