library(dplyr)
library(stringr)
library(ggplot2)
library(scales)
# ==========================================
# 1. Argumentos y Setup
# ==========================================
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 5) {
  stop("Uso: Rscript Plot_Joint_Likelihood.R <Input_Dir> <Output_Dir> <Regex_IDs> <Valor_Migracion> <Valor_Seleccion>")
}

input_dir   <- args[1]
output_dir  <- args[2]
task_id_reg <- args[3]
m_value     <- as.numeric(args[4])
s_value     <- as.numeric(args[5]) 
s_value     <-  s_value/2
print(paste("Analizando Superficie Conjunta (D vs s). Migración teórica:", m_value))
print(paste("Analizando Superficie Conjunta (D vs s). Selección teórica:", s_value))

# ==========================================
# 2. Función de Extracción de MLE
# ==========================================
process_file <- function(f, rep_id, method_name) {
  d <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  
  ll_matrix <- d[, colnames(d) != "SNP", drop = FALSE]
  
  # A) CÁLCULO COMPOSITE LIKELIHOOD 
  ll_sums <- colSums(ll_matrix, na.rm = TRUE)
  best_col <- names(which.max(ll_sums))
  best_LL_sum <- max(ll_sums, na.rm = TRUE)
  
  best_D <- as.numeric(str_extract(best_col, "(?<=D_)[0-9.eE+-]+"))
  best_s <- as.numeric(str_extract(best_col, "(?<=s_)[0-9.eE+-]+"))
  
  composite_mle <- data.frame(
    Method = method_name,
    Replica_ID = rep_id,
    D = best_D,
    s = best_s,
    LL_sum = best_LL_sum
  )
  
  # B) EXTRACCIÓN SINGLE-SNP MLE
  max_idx <- max.col(ll_matrix, ties.method = "first")
  best_cols_snps <- colnames(ll_matrix)[max_idx]
  
  single_snp_mle <- data.frame(
    Method = method_name,
    SNP_ID = d$SNP,
    Replica_ID = rep_id,
    D = as.numeric(str_extract(best_cols_snps, "(?<=D_)[0-9.eE+-]+")),
    s = as.numeric(str_extract(best_cols_snps, "(?<=s_)[0-9.eE+-]+")),
    LL = apply(ll_matrix, 1, max, na.rm = TRUE)
  )
  
  return(list(composite = composite_mle, singles = single_snp_mle))
}

# ==========================================
# 3. Carga y Procesamiento de Archivos
# ==========================================
print("Buscando y procesando archivos TRON...")

methods_to_find <- c("TRON_LEGACY_Grid", "TRON_ARES_Grid", "TRON_SPIKES_Grid")
all_composite <- list()
all_singles <- list()

for (method in methods_to_find) {
  # MODIFICACIÓN: Usamos .*\\.txt$ para que capture tanto archivos 'neutros' como 'seleccion'
  pattern_to_search <- paste0("^", method, "_TaskID_", task_id_reg, "_.*\\.txt$")
  
  file_list <- list.files(path = input_dir, pattern = pattern_to_search, full.names = TRUE)
  print(file_list)
  
  if (length(file_list) > 0) {
    # El resto se mantiene exactamente igual
    replica_ids <- str_extract(basename(file_list), "(?<=TaskID_)[0-9A-Za-z_]+(?=_)") 
    
    for (i in seq_along(file_list)) {
      res <- process_file(file_list[i], replica_ids[i], method)
      all_composite[[length(all_composite) + 1]] <- res$composite
      all_singles[[length(all_singles) + 1]] <- res$singles
    }
  } else {
    print(paste("Advertencia: No se encontraron archivos para", method, "con el patrón:", pattern_to_search))
  }
}

if (length(all_composite) == 0) stop("No se procesó ningún archivo. Revisa las rutas y Regex.")

mle_per_replica_raw <- bind_rows(all_composite)
all_singles_mle_raw <- bind_rows(all_singles)

# Formatear el nombre del método para los gráficos
mle_per_replica_raw$Method <- gsub("_Grid", "", mle_per_replica_raw$Method)
all_singles_mle_raw$Method <- gsub("_Grid", "", all_singles_mle_raw$Method)

# ==========================================
# 4. Agrupación por Lotes (Tu Lógica)
# ==========================================
GROUP_SIZE <- 10

mle_per_replica <- mle_per_replica_raw %>%
  mutate(
    TaskID_num   = as.integer(str_extract(Replica_ID, "[0-9]+")),
    TaskID_Group = paste0(
      "Group_", sprintf("%02d", (TaskID_num - 1) %/% GROUP_SIZE + 1),
      "_IDs_", ((TaskID_num - 1) %/% GROUP_SIZE) * GROUP_SIZE + 1,
      "-", ((TaskID_num - 1) %/% GROUP_SIZE + 1) * GROUP_SIZE
    )
  ) %>% filter(!is.na(TaskID_num))

all_singles_mle <- all_singles_mle_raw %>%
  mutate(
    TaskID_num   = as.integer(str_extract(Replica_ID, "[0-9]+")),
    TaskID_Group = paste0(
      "Group_", sprintf("%02d", (TaskID_num - 1) %/% GROUP_SIZE + 1),
      "_IDs_", ((TaskID_num - 1) %/% GROUP_SIZE) * GROUP_SIZE + 1,
      "-", ((TaskID_num - 1) %/% GROUP_SIZE + 1) * GROUP_SIZE
    )
  ) %>% filter(!is.na(TaskID_num))

task_groups <- sort(unique(mle_per_replica$TaskID_Group))
cat("\nGrupos detectados:", length(task_groups), "\n")
print(task_groups)

# ==========================================
# 5. Generación de Gráficos POR GRUPO
# ==========================================
color_theoric <- "#d7191c"
color_legacy <- "#3498db" 
color_ares <- "#2ecc71"
color_spikes <- "#C77CFF"   
margen_s <- 0.05

for (current_group in task_groups) {
  
  cat("Generando gráficos para:", current_group, "...\n")
  
  # Filtrar datos solo para este grupo
  group_mle <- mle_per_replica %>% filter(TaskID_Group == current_group)
  group_singles <- all_singles_mle %>% filter(TaskID_Group == current_group)
  actual_replicas <- length(unique(group_mle$Replica_ID))
  
  # --- GRÁFICO A: Distribución de D ---
  p_dist_D <- ggplot(group_mle, aes(x = Method, y = D, fill = Method)) +
    geom_violin(color = "black", alpha = 0.5) +
    geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, size = 3, alpha = 0.7, color = "black", shape = 21) +
    geom_hline(aes(yintercept = m_value, color = "Theoretical"), linetype = "dashed", linewidth = 1) +
    scale_y_continuous(
    trans = pseudo_log_trans(sigma = 1e-10, base = 10),
    breaks = c(1e-10, 1e-08, 1e-06, 1e-04, 0.01, 1)) +
    annotation_logticks(sides = "l") +
    scale_fill_manual(values = c("TRON_LEGACY" = color_legacy, "TRON_ARES" = color_ares, "TRON_SPIKES"= color_spikes)) +
    scale_color_manual(name = "", values = c("Theoretical" = color_theoric), 
                       labels = paste0("Theoric D = ", m_value)) +
    labs(
      title = paste("Variance of Diffusion (D) Estimation -", gsub("_", " ", current_group)),
      subtitle = paste0("Distribution of MLEs across ", actual_replicas, " replicas. Fixed s = ", s_value),
      x = "Likelihood Model",
      y = expression(paste("Estimated D (", m^2, "/gen) - Log Scale"))
    ) +
    theme_bw() + theme(legend.position = "bottom")

  # --- GRÁFICO B: Distribución de s ---
  p_dist_s <- ggplot(group_mle, aes(x = Method, y = s, fill = Method)) +
    geom_violin(color = "black", alpha = 0.5) +
    geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, size = 3, alpha = 0.7, color = "black", shape = 21) +
    geom_hline(aes(yintercept = s_value, color = "Theoretical"), linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values = c("TRON_LEGACY" = color_legacy, "TRON_ARES" = color_ares, "TRON_SPIKES"= color_spikes)) +
    scale_color_manual(name = "", values = c("Theoretical" = color_theoric), 
                       labels = paste0("Theoric s = ", s_value)) +
    labs(
      title = paste("Variance of Selection (s) Estimation -", gsub("_", " ", current_group)),
      subtitle = paste0("Distribution of MLEs across ", actual_replicas, " replicas. Fixed D = ", m_value),
      x = "Likelihood Model",
      y = "Estimated Selection Coefficient (s)"
    ) +
    theme_bw() + theme(legend.position = "bottom")

  # --- GRÁFICO C: Diagnóstico de SNPs Individuales ---
  p_diag_snps <- ggplot(group_singles, aes(x = s, y = LL)) +
    geom_jitter(aes(fill = as.factor(Replica_ID)), width = 0.002, size = 2.5, shape = 21, color = "black", alpha = 0.6) +
    geom_vline(xintercept = s_value, linetype = "dashed", color = "red", linewidth = 1) +
    coord_cartesian(xlim = c(s_value - margen_s, s_value + margen_s)) +
    facet_wrap(~ Method, scales = "free_y") + 
    labs(
      title = paste("Single-Locus MLE Diagnostics -", gsub("_", " ", current_group)),
      subtitle = paste0("Comparación directa. Línea roja = s teórico (", s_value, ")"),
      x = "Inferred Selection Coefficient (s)",
      y = "Maximum Log-Likelihood (LL)",
      fill = "Replica ID"
    ) +
    theme_bw() +
    theme(legend.position = "right")

  # ==========================================
  # 6. Guardado de Resultados Limpio
  # ==========================================
  # El nombre ahora está basado en el grupo (ej. Group_01_IDs_1-10) y no en el regex
  clean_filename <- paste0("Comparison_m1_", current_group, "_Mig_", m_value, "_Sel_", s_value)
  
  ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Dist_D.png")), 
         plot = p_dist_D, width = 8, height = 6)
  
  ggsave(filename = file.path(output_dir, paste0(clean_filename, "_Dist_s.png")), 
         plot = p_dist_s, width = 8, height = 6)
  
  ggsave(filename = file.path(output_dir, paste0(clean_filename, "_SingleSNP_Diagnostics.png")), 
         plot = p_diag_snps, width = 10, height = 6)
}

print("Todos los gráficos comparativos agrupados han sido guardados exitosamente.")