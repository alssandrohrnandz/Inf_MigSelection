library(dplyr) #######CON JOB_ID

output_dir <- "/Users/miel_/Documents/Doctorado/output/"

for (i in 1:10){

group_results <- list()
pattern_to_search <- paste0("^SelectionDiffusion_Analysis_.*_Job_",i,"\\.txt$") 
pattern_to_search

file_list <- list.files(
  path = "/Users/miel_/Documents/Doctorado/LL/m_05",   #empirical_data",
  pattern = pattern_to_search,
  full.names = TRUE
)
if (length(file_list) == 0) {
  print(paste("ADVERTENCIA: No se encontraron archivos para el Job/Grupo", i))
  next # Pasar a la siguiente iteración si no hay archivos
}

group_data <- do.call(rbind, lapply(file_list, function(f) {
  # Usar tryCatch es una buena práctica para archivos corruptos/vacíos
  df <- tryCatch(read.table(f, header = TRUE), error = function(e) {
    warning(paste("Error leyendo archivo:", f, e$message))
    return(NULL)
  })
  return(df)
}))

if (is.null(group_data) || nrow(group_data) == 0) {
  print(paste("ADVERTENCIA: Datos vacíos o fallidos para el Job/Grupo", i))
  next
}

summary <- group_data %>%
  group_by(D) %>%
  summarise(LL_sum = sum(LL, na.rm = TRUE), .groups = "drop") # Usar na.rm=TRUE

new_summary<-summary[2:62,]

plot(new_summary$D,new_summary$LL_sum)

max_ll_index <- which.max(new_summary$LL_sum)
max_D <- new_summary$D[max_ll_index]
max_LL <- new_summary$LL_sum[max_ll_index]

file_name <- paste0(output_dir,"ll_empirical_data_job_",i,".png")

png(filename = file_name, width = 730, height = 400)

plot(new_summary$D, new_summary$LL_sum,
     xlab = "Parámetro D",
     ylab = "Verosimilitud Compuesta (300 alelos)",
     main = paste("Verosimilitud Compuesta vs. Parámetro D\n (300 alelos aleatorios). Grupo=", i),
     pch = 19,
     col = "blue")
lines(new_summary$D, new_summary$LL_sum, col = "red", lwd = 2)
points(max_D, max_LL, col = "green", pch = 17, cex = 1.5) # Triángulo grande
text(max_D, max_LL,
     labels = paste0("D=", max_D),
     pos = 1, # 4 = a la derecha
     col = "green")
file_name <- paste(output_dir, "/ll_job_", i,"_v2", ".png")
dev.off()
}

# Inicializar una lista para guardar los resultados máximos de cada Job/Grupo
max_results <- list() 
output_dir <- "/Users/miel_/Documents/Doctorado/output/"

for (i in 1:10){
  
  group_results <- list()
  pattern_to_search <- paste0("^SelectionDiffusion_Analysis_.*_Job_",i,"\\.txt$")
  
  file_list <- list.files(
    path = "/Users/miel_/Documents/Doctorado/LL/m_05/",
    pattern = pattern_to_search,
    full.names = TRUE
  )
  
  if (length(file_list) == 0) {
    print(paste("ADVERTENCIA: No se encontraron archivos para el Job/Grupo", i))
    next # Pasar a la siguiente iteración si no hay archivos
  }
  
  group_data <- do.call(rbind, lapply(file_list, function(f) {
    df <- tryCatch(read.table(f, header = TRUE), error = function(e) {
      warning(paste("Error leyendo archivo:", f, e$message))
      return(NULL)
    })
    return(df)
  }))
  
  if (is.null(group_data) || nrow(group_data) == 0) {
    print(paste("ADVERTENCIA: Datos vacíos o fallidos para el Job/Grupo", i))
    next
  }
  
  summary <- group_data %>%
    group_by(D) %>%
    summarise(LL_sum = sum(LL, na.rm = TRUE), .groups = "drop")
  
  new_summary<-summary[2:62,] # Filtrar los datos como lo haces
  
  # ----------------------------------------------------
  # PASO CLAVE: Encontrar y Guardar el Máximo
  # ----------------------------------------------------
  max_ll_index <- which.max(new_summary$LL_sum)
  max_D <- new_summary$D[max_ll_index]
  max_LL <- new_summary$LL_sum[max_ll_index]
  
  # Almacenar el resultado en la lista
  max_results[[i]] <- data.frame(
    Job = i,
    D_max = max_D,
    LL_max = max_LL
  )
  
  # --- El código para graficar el LL vs D individual (que ya tenías)
  # file_name <- paste0(output_dir,"ll_job_",i,".png")
  # png(filename = file_name, width = 730, height = 400)
  # ... (resto del código de plot individual)
  # dev.off()
}

# ----------------------------------------------------
# PASO FINAL: Combinar los resultados en un único data frame
# ----------------------------------------------------
final_max_df <- do.call(rbind, max_results)
print(head(final_max_df)) # Muestra los primeros resultados

library(ggplot2)

# Crear el gráfico de densidad de los 50 puntos (D_max y LL_max)
p <- ggplot(final_max_df, aes(x = D_max, y = LL_max)) +
  
  # Añadir puntos de dispersión
  geom_point(alpha = 0.7, # Transparencia para ver superposición
             size = 3,
             color = "black") +
  
  # Opcional pero recomendado para ver dónde caen los puntos: Añadir una capa de densidad 2D
  stat_density_2d(aes(fill = after_stat(level)), 
                  geom = "polygon",
                  alpha = 0.3,
                  show.legend = FALSE) +
  scale_fill_gradientn(colors = c("lightblue", "yellow", "red")) + # Escala de color para la densidad
  
  # Etiquetas y Títulos
  labs(title = "Distribución de los Parámetros Óptimos (D vs. LL) en 50 Grupos",
       subtitle = paste0("Máxima Verosimilitud Compuesta (LL) en función del Parámetro D."),
       x = "Parámetro D de Máxima Verosimilitud (D_max)",
       y = "Máxima Verosimilitud Compuesta (LL_max)") +
  
  # Tema limpio
  theme_minimal() +
  
  # Ajustar límites del eje X si quieres que vaya de 0 a 30 (como mencionaste)
  xlim(0, 30)

# Mostrar el gráfico
print(p)

# Opcional: Guardar el gráfico final
ggsave(paste0(output_dir, "distribucion_maximos_50_grupos.png"), 
       plot = p, width = 10, height = 6)
# ----------------------------------------------------
# SIN JOB ID
# ----------------------------------------------------
library(dplyr) 
library(tidyverse)

output_dir <- "/Users/alessandrohernandez/Documents/Posgrado/Doctorado/Discrete/output"

# Parámetros
group_size <- 1
N <- 7
pattern_to_search <- "^SelectionDiffusion_Analysis_.*\\.txt$"
path_to_files <- "/Users/alessandrohernandez/Documents/Posgrado/Doctorado/Discrete/LL/Selection"

# 1. Listar y barajar archivos
file_list <- list.files(
  path = path_to_files,
  pattern = pattern_to_search,
  full.names = TRUE
)

if (length(file_list) < group_size * N) {
  stop("No hay suficientes archivos para formar ", N, " grupos de ", group_size, " archivos.")
}

set.seed(123)
shuffled_files <- sample(file_list)

# 2. Dividir en N grupos de tamaño fijo
file_groups <- split(shuffled_files[1:(group_size * N)], rep(1:N, each = group_size))

# 3. Procesar cada grupo
group_summaries <- lapply(seq_along(file_groups), function(i) {
  files <- file_groups[[i]]
  group_name <- paste0("Grupo ", i)
  
  data <- do.call(rbind, lapply(files, function(f) {
    tryCatch(read.table(f, header = TRUE), error = function(e) {
      warning(paste("Error leyendo archivo:", f, e$message))
      return(NULL)
    })
  }))
  
  if (is.null(data)) return(NULL)
  
  summary <- data %>%
    group_by(s) %>%
    summarise(LL_sum = sum(LL, na.rm = TRUE), .groups = "drop") %>%
    mutate(Group = group_name)
  
  return(summary)
})

# 4. Combinar y graficar

combined_data <- bind_rows(group_summaries)
max_points <- combined_data %>%
  group_by(Group) %>%
  filter(LL_sum == max(LL_sum)) %>%
  ungroup()


ggplot(combined_data, aes(x = D, y = LL_sum, fill = Group)) +
  geom_line() +
  labs(title = "Comparación de LL por grupo", x = "D", y = "Suma de LL") +
  theme_minimal()

# --- DEFINE VALOR ---
valor_teorico <- 0.01 
ggplot(max_points, aes(x = "Inferencia", y = s)) +
  geom_violin(fill = "grey90", color = "grey40", draw_quantiles = c(0.5)) + # Línea central es mediana
  # 2. Puntos (Jitter)
  geom_jitter(width = 0.15, alpha = 0.3, color = "#2c7bb6", size = 1.0) +
  # 3. Boxplot Angosto (Resumen estadístico robusto)
  geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA, alpha = 0.5) +
  # 4. Línea Roja Teórica
  geom_hline(yintercept = valor_teorico, linetype = "dashed", color = "#d7191c", size = 1) +
  labs(
    title = "Inferencia del Coeficiente de Difusión",
    subtitle = bquote(paste("Valor esperado de ", italic(s), " = ", .(valor_teorico))),
    y = expression(paste("D estimado (", m^2, "/gen)")) # Notación matemática en el eje
  ) +
  theme_bw() + # Tema blanco y negro, muy limpio
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    text = element_text(family = "sans", size = 12)
  )

ggplot(combined_data, aes(x = D, y = LL_sum, color = Group, group = Group)) +
 # geom_line(size = 1) +
  geom_point(data = max_points, aes(x = D, y = LL_sum), size = 3, shape = 21, fill = "white", stroke = 1.5) +
  labs(title = "LL por grupo con máximos destacados", x = "D", y = "Suma de LL") +
  theme_minimal()

ggplot(max_points, aes(x = "", y = D)) +
  geom_violin(fill = "skyblue", color = "darkblue") +           # Violin
  geom_jitter(width = 0.1, alpha = 0.5, color = "red") +       # Puntos individuales
  stat_summary(fun = median, geom = "point", color = "black", size = 3) +  # Mediana
  stat_summary(fun.data = function(x) {
    data.frame(
      y = mean(x),
      ymin = quantile(x, 0.25),
      ymax = quantile(x, 0.75)
    )
  }, geom = "errorbar", width = 0.1, color = "darkgreen") +   # Cuartiles y media
  labs(
    title = "Distribución de D con media y cuartiles (m=0.0005",
    x = "",
    y = "Valor de D"
  ) +
  theme_minimal()

ggplot(combined_data, aes(x = D)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  labs(title = "Densidad de valores de D", x = "D", y = "Densidad") +
  theme_minimal()

#############
###SELECCION
##########
# -------------------------------------------------------------------------
# 1. CARGA DE LIBRERÍAS Y CONFIGURACIÓN
# -------------------------------------------------------------------------
library(tidyverse) # Necesario para map_dfr, ggplot, pipes, etc.

# Valores verdaderos a inferir
true_D <- 0.01
true_s <- 0.01

# -------------------------------------------------------------------------
# 2. LECTURA DE ARCHIVOS (Tu código integrado)
# -------------------------------------------------------------------------
pattern_to_search <- "^SelectionDiffusion_Analysis_.*\\.txt$"
path_to_files <- "/Users/alessandrohernandez/Documents/Posgrado/Doctorado/Discrete/LL/Selection"

file_list <- list.files(
  path = path_to_files,
  pattern = pattern_to_search,
  full.names = TRUE
)

# Verificación de seguridad: ¿Encontró archivos?
if(length(file_list) == 0) {
  stop("No se encontraron archivos. Verifica la ruta o el patrón de búsqueda.")
} else {
  message(paste("Se encontraron", length(file_list), "archivos para procesar."))
}

# -------------------------------------------------------------------------
# 3. EXTRACCIÓN DEL MEJOR VALOR (MLE) DE CADA ARCHIVO
# -------------------------------------------------------------------------

# Esta función lee un archivo y devuelve solo la fila con el mejor Log-Likelihood
process_file <- function(filepath) {
  # Leer el archivo (asumiendo separador por espacio o tabulador)
  df <- read.table(filepath, header = TRUE)
  
  # Filtrar la fila con el máximo LL
  best_row <- df %>% 
    filter(LL == max(LL)) %>% 
    # Si hay empates en el máximo, tomamos el primero para evitar duplicados
    slice(1) 
  
  return(best_row)
}

# Aplicamos la función a todos los archivos y combinamos en un solo dataframe
# map_dfr itera sobre la lista y une los resultados en un data.frame automáticamente
results_df <- map_dfr(file_list, process_file, .id = "File_ID")

# Opcional: ver las primeras filas de los resultados
head(results_df)

# -------------------------------------------------------------------------
# 4. VISUALIZACIÓN (BOXPLOT DE DISPERSIÓN)
# -------------------------------------------------------------------------

# Transformamos los datos a formato largo para graficar D y s juntos
plot_data <- results_df %>%
  select(D, s) %>% # Seleccionamos solo las columnas de interés
  pivot_longer(cols = c("D", "s"), names_to = "Parametro", values_to = "Valor_Inferido")

# Generar el gráfico
p <- ggplot(plot_data, aes(x = Parametro, y = Valor_Inferido, fill = Parametro)) +
  # 1. Boxplot para ver la distribución (mediana y cuartiles)
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  
  # 2. Jitter para ver cada archivo individualmente (puntos dispersos)
  geom_jitter(width = 0.2, alpha = 0.4, size = 2, color = "#2c3e50") +
  
  # 3. Línea roja indicando el VALOR VERDADERO (0.01)
  geom_hline(yintercept = true_D, color = "red", linetype = "dashed", linewidth = 1) +
  
  # 4. Facetas: Separa los gráficos si las escalas son muy diferentes
  # "scales = free" permite que el eje Y se ajuste independientemente para D y s
  facet_wrap(~Parametro, scales = "free") + 
  
  # Estética
  scale_fill_manual(values = c("D" = "#3498db", "s" = "#e67e22")) +
  theme_bw() +
  labs(
    title = paste("Dispersión de la inferencia en", length(file_list), "archivos simulados"),
    subtitle = "La línea roja discontinua indica el valor verdadero (0.01)",
    y = "Valor Estimado (MLE)",
    x = "Parámetro"
  ) +
  theme(legend.position = "none") # Ocultar leyenda redundante

print(p)


dev.off()

# Inicializar una lista para guardar los resultados máximos de cada Job/Grupo
max_results <- list() 
output_dir <- "/Users/miel_/Documents/Doctorado/output/"

for (i in 1:50){
  
  group_results <- list()
  pattern_to_search <- paste0("^SelectionDiffusion_Analysis_.*_Job_",i,"\\.txt$")
  
  file_list <- list.files(
    path = "/Users/miel_/Documents/Doctorado/LL",
    pattern = pattern_to_search,
    full.names = TRUE
  )
  
  if (length(file_list) == 0) {
    print(paste("ADVERTENCIA: No se encontraron archivos para el Job/Grupo", i))
    next # Pasar a la siguiente iteración si no hay archivos
  }
  
  group_data <- do.call(rbind, lapply(file_list, function(f) {
    df <- tryCatch(read.table(f, header = TRUE), error = function(e) {
      warning(paste("Error leyendo archivo:", f, e$message))
      return(NULL)
    })
    return(df)
  }))
  
  if (is.null(group_data) || nrow(group_data) == 0) {
    print(paste("ADVERTENCIA: Datos vacíos o fallidos para el Job/Grupo", i))
    next
  }
  
  summary <- group_data %>%
    group_by(D) %>%
    summarise(LL_sum = sum(LL, na.rm = TRUE), .groups = "drop")
  
  new_summary<-summary[2:62,] # Filtrar los datos como lo haces
  
  # ----------------------------------------------------
  # PASO CLAVE: Encontrar y Guardar el Máximo
  # ----------------------------------------------------
  max_ll_index <- which.max(new_summary$LL_sum)
  max_D <- new_summary$D[max_ll_index]
  max_LL <- new_summary$LL_sum[max_ll_index]
  
  # Almacenar el resultado en la lista
  max_results[[i]] <- data.frame(
    Job = i,
    D_max = max_D,
    LL_max = max_LL
  )
  
  # --- El código para graficar el LL vs D individual (que ya tenías)
  # file_name <- paste0(output_dir,"ll_job_",i,".png")
  # png(filename = file_name, width = 730, height = 400)
  # ... (resto del código de plot individual)
  # dev.off()
}

# ----------------------------------------------------
# PASO FINAL: Combinar los resultados en un único data frame
# ----------------------------------------------------
final_max_df <- do.call(rbind, max_results)
print(head(final_max_df)) # Muestra los primeros resultados

library(ggplot2)

# Crear el gráfico de densidad de los 50 puntos (D_max y LL_max)
p <- ggplot(final_max_df, aes(x = D_max, y = LL_max)) +
  
  # Añadir puntos de dispersión
  geom_point(alpha = 0.7, # Transparencia para ver superposición
             size = 3,
             color = "black") +
  
  # Opcional pero recomendado para ver dónde caen los puntos: Añadir una capa de densidad 2D
  stat_density_2d(aes(fill = after_stat(level)), 
                  geom = "polygon",
                  alpha = 0.3,
                  show.legend = FALSE) +
  scale_fill_gradientn(colors = c("lightblue", "yellow", "red")) + # Escala de color para la densidad
  
  # Etiquetas y Títulos
  labs(title = "Distribución de los Parámetros Óptimos (D vs. LL) en 50 Grupos",
       subtitle = paste0("Máxima Verosimilitud Compuesta (LL) en función del Parámetro D."),
       x = "Parámetro D de Máxima Verosimilitud (D_max)",
       y = "Máxima Verosimilitud Compuesta (LL_max)") +
  
  # Tema limpio
  theme_minimal() +
  
  # Ajustar límites del eje X si quieres que vaya de 0 a 30 (como mencionaste)
  xlim(0, 30)

# Mostrar el gráfico
print(p)

# Opcional: Guardar el gráfico final
ggsave(paste0(output_dir, "distribucion_maximos_50_grupos.png"), 
       plot = p, width = 10, height = 6)


