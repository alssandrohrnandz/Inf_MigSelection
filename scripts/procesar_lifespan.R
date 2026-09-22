#!/usr/bin/env Rscript
# ============================================================
# procesar_lifespan.R
# Procesa outputs de SLiM y extrae:
#   - Resumen por archivo (1 fila)
#   - Lifespan por mutación (1 fila por MutID) -- incluye MutID
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
  library(tibble)
})

# ---------------- Configuración ----------------
directorio       <- "data/results_Discrete/outputs_slim/independent_loci/seleccion"   # <-- AJUSTA
archivo_resumen  <- "slim_lifespan_summary.csv"
archivo_permut   <- "slim_lifespan_per_mutation.csv"
umbral_fix       <- 0.99

# ---------------- Listado de archivos ----------------
rutas_archivos <- list.files(directorio, pattern = "\\.csv$", full.names = TRUE)
rutas_archivos <- rutas_archivos[!grepl("^P_TRON|^P_", basename(rutas_archivos))]

# --- Soporte para SLURM job array ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
  chunk_id  <- as.integer(args[1])
  n_chunks  <- as.integer(args[2])
  idx <- seq(chunk_id, length(rutas_archivos), by = n_chunks)
  rutas_archivos <- rutas_archivos[idx]
  # Sufijo para evitar colisiones entre chunks
  archivo_resumen <- sprintf("slim_lifespan_summary_chunk%03d.csv", chunk_id)
  archivo_permut  <- sprintf("slim_lifespan_per_mutation_chunk%03d.csv", chunk_id)
  cat(sprintf("Chunk %d/%d -> %d archivos\n", chunk_id, n_chunks, length(rutas_archivos)))
}

cat(sprintf("[%s] Archivos encontrados: %d\n", Sys.time(), length(rutas_archivos)))

# ---------------- Reanudación ----------------
if (file.exists(archivo_resumen)) {
  ya <- fread(archivo_resumen, select = "File")$File
  rutas_archivos <- rutas_archivos[!basename(rutas_archivos) %in% ya]
  cat(sprintf("Reanudando. Faltan: %d\n", length(rutas_archivos)))
}

# ---------------- Función: procesa 1 archivo ----------------
procesar_archivo_slim <- function(ruta) {
  nombre <- basename(ruta)

  dt <- tryCatch(
    fread(ruta,
          select = c("Cycle", "MutID", "Freq"),
          showProgress = FALSE, data.table = TRUE),
    error = function(e) { message("Fallo lectura ", nombre, ": ", e$message); NULL }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  # Tipos homogéneos
  dt[, Cycle := as.numeric(Cycle)]
  dt[, MutID := as.character(MutID)]
  dt[, Freq  := as.numeric(Freq)]
  dt <- dt[!is.na(Cycle)]

  # Regla: MutID == 0 corresponde a la raíz / no mutante -> excluir
  dt <- dt[MutID != "0" & MutID != "0.0"]
  if (nrow(dt) == 0) return(NULL)

  gen_max <- max(dt$Cycle)

  # --- Resumen por mutación ---
  per_mut <- dt[, .(
    First_Gen  = min(Cycle),
    Last_Gen   = max(Cycle),
    Lifespan   = max(Cycle) - min(Cycle) + 1,
    Max_Freq   = max(Freq, na.rm = TRUE),
    Last_Freq  = Freq[which.max(Cycle)]
  ), by = MutID]

  # Clasificación
  per_mut[, Status := fifelse(
    Last_Freq >= umbral_fix,                  "Fixed",
    fifelse(Last_Gen < gen_max | Last_Freq == 0, "Lost",
                                               "Segregating")
  )]

  # Metadatos
  tipo    <- ifelse(str_detect(nombre, "neutros"), "Neutro", "Seleccion")
  task_id <- as.integer(str_extract(nombre, "\\d+(?=\\.csv$)"))

  per_mut[, `:=`(File = nombre, TaskID = task_id, Type = tipo)]

  # --- Resumen de archivo (UNA fila) ---
  resumen_archivo <- tibble(
    File                = nombre,
    TaskID              = task_id,
    Type                = tipo,
    Total_Mutations     = nrow(per_mut),
    Mean_Lifespan       = mean(per_mut$Lifespan,   na.rm = TRUE),
    Median_Lifespan     = median(per_mut$Lifespan, na.rm = TRUE),
    SE_Lifespan         = sd(per_mut$Lifespan, na.rm = TRUE) / sqrt(nrow(per_mut)),
    N_Fixed             = sum(per_mut$Status == "Fixed"),
    N_Lost              = sum(per_mut$Status == "Lost"),
    N_Segregating       = sum(per_mut$Status == "Segregating"),
    Pct_Fixed           = 100 * mean(per_mut$Status == "Fixed"),
    Pct_Lost            = 100 * mean(per_mut$Status == "Lost"),
    Mean_Lifespan_Fixed = mean(per_mut$Lifespan[per_mut$Status == "Fixed"], na.rm = TRUE),
    Mean_Lifespan_Lost  = mean(per_mut$Lifespan[per_mut$Status == "Lost"],  na.rm = TRUE)
  )

  # Ordenar columnas del per-mut para CSV
  setcolorder(per_mut, c("File", "TaskID", "Type", "MutID",
                         "First_Gen", "Last_Gen", "Lifespan",
                         "Max_Freq", "Last_Freq", "Status"))

  list(resumen = resumen_archivo, per_mut = as_tibble(per_mut))
}

# ---------------- Loop incremental ----------------
n <- length(rutas_archivos)
t0 <- Sys.time()

for (i in seq_len(n)) {
  out <- tryCatch(
    procesar_archivo_slim(rutas_archivos[i]),
    error = function(e) {
      message("Error en ", basename(rutas_archivos[i]), ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(out)) {
    # Append al resumen
    existe_res <- file.exists(archivo_resumen)
    fwrite(out$resumen, archivo_resumen,
           append = existe_res, col.names = !existe_res)

    # Append al per-mut (puede ser grande, se escribe incremental)
    existe_pm <- file.exists(archivo_permut)
    fwrite(out$per_mut, archivo_permut,
           append = existe_pm, col.names = !existe_pm)
  }

  if (i %% 50 == 0) {
    gc(verbose = FALSE)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("[%s] %d / %d  (%.1f min, %.2f s/archivo)\n",
                Sys.time(), i, n, elapsed, 60 * elapsed / i))
  }
}

cat(sprintf("[%s] LISTO. Resumen: %s | Per-mut: %s\n",
            Sys.time(), archivo_resumen, archivo_permut))