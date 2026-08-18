# === Archivo: merge_and_filter_by_chr.R === #

# 1. Recibir argumentos de la línea de comandos desde Bash
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Error. Uso: Rscript merge_and_filter_by_chr.R <CHR> <pop_info> <freq_pattern> <advantageous> <output_prefix>")
}

chr_num           <- args[1]
pop_info_file     <- args[2]
freq_pattern      <- args[3]
advantageous_file <- args[4]
output_prefix     <- args[5]
cadd_file         <- args[6]  # <-- Nuevo argumento

# Formatear el nombre del archivo de entrada cambiando {CHR} por el número actual
input_freq_file <- gsub("\\{CHR\\}", chr_num, freq_pattern)

cat("=========================================\n")
cat("Procesando Cromosoma:", chr_num, "\n")
cat("=========================================\n")

# 2. Cargar data.table
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("El paquete 'data.table' no está instalado. Ejecuta install.packages('data.table') en R primero.")
}
library(data.table)

# 3. Cargar información de poblaciones
cat("Cargando info de población...\n")
pop_info <- fread(pop_info_file)
# Renombrar columna para que coincida con el archivo de frecuencias
setnames(pop_info, "Group_Master_ID", "CLST") 

# 4. Cargar lista de alelos ventajosos
cat("Cargando lista de alelos bajo selección...\n")
# Según tu imagen, ventajosos.txt no parece tener encabezado, el rsID está en la columna 2
adv_data <- fread(advantageous_file, header = FALSE) 
adv_alleles <- unique(adv_data$V2) # Extraer los rsID únicos

# 5. Cargar el archivo de frecuencias estratificadas del cromosoma
cat("Leyendo datos genéticos:", input_freq_file, "...\n")
# fread es rapidísimo y solo cargará en memoria la fracción del cromosoma (unos cientos de MB)
freq_data <- fread(input_freq_file)

# 6. Unir (Merge) los datos
cat("Realizando la unión (merge) de frecuencias con metadata...\n")
# Hacemos un inner join eficiente. Solo traemos YBP y Total_Muestras
merged_data <- merge(
  freq_data, 
  pop_info[, .(CLST, Media_YBP, Media_Lat, Media_Long)], 
  by = "CLST", 
  all.x = FALSE # Ignora poblaciones en el .frq que no estén en tu txt
)

# Reordenar las columnas para dejarlo presentable
setcolorder(merged_data, c("CHR", "SNP", "CLST", "A1", "A2", "MAF", "MAC", "NCHROBS", "Media_YBP", "Media_Lat", "Media_Long"))

# 7. Filtrar alelos bajo selección vs neutrales
cat("Separando alelos...\n")
data_sel <- merged_data[SNP %in% adv_alleles]
data_neu <- merged_data[!SNP %in% adv_alleles]

# 8. Integrar puntajes CADD como filtro para los neutrales
cat("Cargando puntajes CADD y filtrando datos neutrales...\n")
# Usamos 'select' para leer solo las dos columnas que nos importan y ahorrar RAM
cadd_data <- fread(cadd_file, select = c("ID", "CADD"))
setnames(cadd_data, "ID", "SNP") # Renombramos para que el merge funcione automáticamente

# Inner join: all = FALSE elimina cualquier SNP en data_neu que no esté en cadd_data
data_neu_filtrado <- merge(data_neu, cadd_data, by = "SNP", all = FALSE)

# 9. Escribir resultados
output_neu_file <- paste0(output_prefix, "_CHR", chr_num, "_neutral_CADD.csv")
output_sel_file <- paste0(output_prefix, "_CHR", chr_num, "_selection.csv")

cat("Guardando:", output_neu_file, "\n")
fwrite(data_neu_filtrado, output_neu_file, quote = FALSE)

cat("Guardando:", output_sel_file, "\n")
fwrite(data_sel, output_sel_file, quote = FALSE)

cat("¡Cromosoma", chr_num, "completado con éxito!\n")