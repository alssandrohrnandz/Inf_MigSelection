# ==============================================================================
# SCRIPT PARA GENERAR ARCHIVO MAESTRO DE FRECUENCIAS (SÓLO CHR 1)
# ==============================================================================
library(data.table)
library(readr)

args <- commandArgs(trailingOnly = TRUE)

#Mensaje para que entre mis argumentos
if (length(args)==0){
    cat("Inserte cromosoma(s) a estudiar: ")
    CHR<-readLines("stdin",n=1)
} else {
   CHR <- args[1]
}

#si metemos por error espacios:
CHR<-trimws(CHR)

# RUTAS DE ARCHIVOS
CADD_FILE <- "/mnt/data/dortega/hlopezh/Inf_MigSelection/AADNA_data/CADD_Data/OrderedCADDDatasets1.txt"
BIM_FILE  <- "/mnt/data/dortega/hlopezh/Inf_MigSelection/AADNA_data/1240k_bin.bim"
FRQ_FILE  <- paste0("/mnt/data/dortega/hlopezh/Inf_MigSelection/chr_frequencies/1240k_freq_chr_",CHR,".frq.strat")
if (!file.exists(FRQ_FILE)) {
   stop("Archivo no encontrado. Hay que llamar a BATMAN. Ejecución suspendida", call. = TRUE)
}
cat(paste0("Analizando sitios del Cromosoma ", CHR, " en archivo \"", FRQ_FILE, "\"\n" ))
#OUT_FILE  <- "data/Master_Chr1_Frequencies.tsv"
OUTPUT_DIR <- "AADNA_data/frq_chunks"
#stop("Fin del test")
print("Leyendo OrderedCADDDatasets1.txt (sitios neutros)")

cadd <- fread(CADD_FILE, header = TRUE)
colnames(cadd) <- c("RowID", "Chrom", "Pos", "Ref", "Alt", "RawScore", "PHRED") #los nombres deben estar en este orden en el cadd

cadd_chr1 <- unique(cadd[Chrom == 1, .(Chrom, Pos)]) #filtrao

print(paste("-> Se encontraron", nrow(cadd_chr1), "posiciones únicas en CADD para el Chr 1."))

print("Filtrando y mapeando con data.bim ")

bim <- fread(BIM_FILE, select = c(1, 2, 4), col.names = c("Chrom", "ID", "Pos"))
bim_chr1 <- bim[Chrom == 1]
SNPS_no<-nrow(bim_chr1)

print(paste0("Se encontraron ",SNPS_no, " SNPs únicos en el archivo BIM"))

stop("Fin del test")

# TODO: verificar de calculos para alelos ancestrales

cadd_with_rsid <- merge(cadd_chr1, bim_chr1, by = c("Chrom", "Pos"), all = FALSE) #merge del archivo cadd y archivo bim

print("=== 3. Cruzando con el archivo gigante frq_v4.txt ===")
frq <- fread(FRQ_FILE, header = FALSE)
colnames(frq) <- c('Chrom_frq', 'ID', 'Population', 'A1', 'A2', 'AlleleFrequency', 'DerivedAlleles', 'AlleleCount')

frq <- frq[ID != "SNP" & Population != "CLST"]

master_final <- merge(frq, cadd_with_rsid, by = "ID", all = FALSE)

master_final <- master_final[, .(Chrom, Pos, ID, Population, A1, A2, AlleleFrequency, DerivedAlleles, AlleleCount)] # reordenamiento

print(paste("-> Matriz maestra unificada creada con", nrow(master_final), "filas totales."))

print("=== 4. Subdividiendo en bloques de 500 SNPs únicos ===") #extraemos los id unicos
unique_alleles <- unique(master_final$ID)
total_alleles  <- length(unique_alleles)
chunk_size     <- 500 # cada archivo tiene 500 SNP
num_chunks     <- ceiling(total_alleles / chunk_size)

print(paste("Total de SNPs únicos detectados:", total_alleles))

for (i in 1:num_chunks) {
    start_idx <- ((i - 1) * chunk_size) + 1
    end_idx   <- min(i * chunk_size, total_alleles)
    
    current_alleles <- unique_alleles[start_idx:end_idx]
    
    chunk_data <- master_final[ID %in% current_alleles]
    
    chunk_file <- file.path(OUTPUT_DIR, paste0("Master_Chunk_", i, ".rds"))
    
    saveRDS(chunk_data, file = chunk_file)
    
    if (i %% 10 == 0 || i == num_chunks) {
        print(paste("Progreso:", i, "/", num_chunks, "bloques guardados."))
    }
}

print("=== ¡PROCESO DE SUBDIVISIÓN FINALIZADO CON ÉXITO! ===")
#los archivos se guardaron en AADNA_data/frq_chunks
# continua en RUN_PDE_Legacy.sh > infLilelihood_muitations_binomial.R