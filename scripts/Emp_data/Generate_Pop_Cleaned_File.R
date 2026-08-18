# Cargar librerías necesarias
library(dplyr)
library(stringr)

# Definir las rutas de tus archivos
fam_file   <- "AADNA_data/1240k.pedind"
ind_file   <- "AADNA_data/v66.1240K.aadr.PUB.ind"
csv_file   <- "AADNA_data/Eurasia_Populations_Periods_Grid.csv"
output_pop <- "AADNA_data/Populations_Cleaned.pop"

# 1. Leer el archivo .fam (para jalar el Family ID numérico correcto)
fam_data <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE) %>%
  select(FID = V1, IID = V2)

# 2. Leer el archivo .ind limpiando los espacios en blanco
# Los archivos .ind a veces usan múltiples espacios como separadores
ind_lines <- readLines(ind_file)
ind_clean <- lapply(ind_lines, function(line) {
  parts <- str_split(str_trim(line), "\\s+")[[1]]
  # Asegurar que tenga al menos las 3 columnas (IID, Sex, Population)
  if(length(parts) >= 3) {
    return(data.frame(IID = parts[1], Population = parts[3], stringsAsFactors = FALSE))
  }
  return(NULL)
})
ind_data <- bind_rows(ind_clean)

# 3. Leer tu CSV de poblaciones objetivo
pop_grid <- read.csv(csv_file, stringsAsFactors = FALSE)

# 4. Cruzar la información y filtrar
pop_cleaned <- fam_data %>%
  # Cruza por el identificador del individuo (IID) para obtener su población
  inner_join(ind_data, by = "IID") %>%
  # Filtra para mantener solo las poblaciones que están en la columna CLST del CSV
  filter(Population %in% pop_grid$CLST) %>%
  # Estructura exacta solicitada por PLINK: 1. FID, 2. IID, 3. Cluster
  select(FID, IID, Population)

# 5. Guardar el archivo sin encabezados ni comillas (formato plano)
write.table(pop_cleaned, file = output_pop,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("¡Listo! El archivo", output_pop, "se generó con", nrow(pop_cleaned), "muestras coordinadas.\n")
