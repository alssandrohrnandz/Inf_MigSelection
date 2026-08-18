library(data.table)
library(tidyverse)
# Definir la ruta del archivo
archivo_anno <- "AADNA_data/v66.1240K.aadr.PUB.anno"

# Leer solo el encabezado y las primeras filas
preview <- fread(archivo_anno, nrows = 5, sep = "\t")
# print(preview)
# header <- names(fread("AADNA_data/v66.1240K.aadr.PUB.anno", nrows = 0))
# Mostrar la estructura básica
# data.frame(Index = 1:length(header), Column_Name = header)
cols_to_keep <- c(1, 11, 15, 18,19)
anno_subset <- fread("AADNA_data/v66.1240K.aadr.PUB.anno", 
                     select = cols_to_keep, 
                     sep = "\t")

setnames(anno_subset, 
         old = names(anno_subset), 
         new = c("ID", "Date_BP", "CLST", "Lat", "Lon"))

cols_to_fix <- c("Date_BP", "Lat", "Lon")
anno_subset[, (cols_to_fix) := lapply(.SD, function(x) as.numeric(as.character(x))), .SDcols = cols_to_fix]
anno_eurasia <- anno_subset[Lat >= 30 & Lat <= 75 & Lon >= -15 & Lon <= 55 & Date_BP <= 11700]
# Filtrar por las coordenadas de tu cuadrícula de Eurasia
cat("Muestras originales:", nrow(anno_subset), "\n")
cat("Muestras en Eurasia:", nrow(anno_eurasia), "\n")
print(anno_eurasia)



# 2. Agrupar por CLST y calcular el promedio (Mean)
# Usamos na.rm = TRUE para ignorar los valores faltantes en el cálculo
resumen_clst <- anno_eurasia[, .(
  Mean_Date_BP = mean(Date_BP, na.rm = TRUE),
  Mean_Lat     = mean(Lat, na.rm = TRUE),
  Mean_Lon     = mean(Lon, na.rm = TRUE),
  N_Samples    = .N  # Contamos cuántas muestras hay por cada grupo
), by = CLST]

resumen_clst <- resumen_clst[order(-N_Samples)]
print(resumen_clst)
print(nrow(resumen_clst))
head(resumen_clst)

library(ggplot2)
library(scales)

# 1. Asumimos que ya corriste tu filtro del Holoceno:
# anno_eurasia <- anno_subset[Lat >= 30 & Lat <= 75 & Lon >= -15 & Lon <= 140 & Date_BP <= 11700]

# 2. Calcular la Generación (Forward in time: Gen 1 = 11700 BP)
resumen_clst[, Generation := floor((11700 - Mean_Date_BP) / 29) + 1]

# 3. Calcular cromosomas totales por generación
# Nota metodológica: Aunque los datos aDNA (1240k) suelen ser "pseudohaploides" 
# (se selecciona un solo alelo al azar), la muestra proviene de un individuo diploide. 
# Si tu PDE modela cromosomas totales de la población, multiplicamos por 2.
timeline_dt <- resumen_clst[, .(
  N_Individuals = .N,
  N_Chromosomes = .N * 2
), by = Generation][order(Generation)]

# Completar generaciones sin datos (rellenar con ceros)
all_gens <- data.table(Generation = 1:max(timeline_dt$Generation))
timeline_dt <- merge(all_gens, timeline_dt, by = "Generation", all.x = TRUE)
timeline_dt[is.na(N_Chromosomes), N_Chromosomes := 0]

# 1. Definir los periodos arqueológicos (Fechas muy generales para Eurasia)
# Puedes ajustar los valores Inicio_BP y Fin_BP si tu tesis se centra en una zona específica
periodos <- data.table(
  Periodo = c("Mesolithic", "Neolithic", "Bronze Age", "Iron Age", "Historical"),
  Inicio_BP = c(11700, 8000, 5000, 3200, 2000),
  Fin_BP = c(8000, 5000, 3200, 2000, 0),
  # Colores pastel tenues (Azul, Verde, Naranja, Gris, Morado)
  Color = c("#e0f7fa", "#f1f8e9", "#fff3e0", "#eceff1", "#f3e5f5") 
)

# Convertir los tiempos BP a tu escala de Generaciones (eje X)
periodos[, xmin := (11700 - Inicio_BP) / 29]
periodos[, xmax := (11700 - Fin_BP) / 29]

# 4. Generación del gráfico (Estilo Nature / MBE)
plot_timeline <- ggplot(timeline_dt, aes(x = Generation, y = N_Chromosomes)) +
  
  # A. CAPA DE FONDO: Rectángulos de los periodos arqueológicos
  geom_rect(data = periodos,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Periodo),
            inherit.aes = FALSE, alpha = 0.5) +
  # Asignar los colores exactos definidos en la tabla y ocultar la leyenda lateral
  scale_fill_manual(values = setNames(periodos$Color, periodos$Periodo), guide = "none") +
  
  # B. TEXTOS: Nombres de los periodos en la parte superior
  # Nota: y = 170 es porque nuestro techo visual será 180
  geom_text(data = periodos,
            aes(x = (xmin + xmax) / 2, y = 110, label = Periodo),
            inherit.aes = FALSE, color = "gray40", size = 4.5, fontface = "italic") +
  
  # C. DATOS CRUDOS: Barras 
  geom_col(fill = "gray60", width = 1, alpha = 0.8) +
  
  # D. CURVA LOESS: Tendencia
  geom_smooth(method = "loess", span = 0.1, color = "#2c3e50", 
              fill = "#3498db", alpha = 0.3, linewidth = 1.2) +
  
  # E. EJES (Con la sintaxis actualizada para evitar warnings)
  scale_x_continuous(
    name = "Time (Generations since 11,700 BP)",
    expand = c(0, 0),
    sec.axis = sec_axis(
      transform = ~ 11700 - (. - 1) * 29, 
      name = "Years Before Present (BP)",
      labels = comma 
    )
  ) +
  
  scale_y_continuous(
    name = "Total Chromosomes Sampled",
    expand = expansion(mult = c(0, 0)),
    labels = comma
  ) +
  
  # Recorte visual para lidiar con el sesgo moderno (techo en 450)
  coord_cartesian(ylim = c(0, 120)) +
  
  # F. ESTÉTICA GENERAL
  theme_classic(base_size = 14) +
  labs(
    title = "Spatiotemporal Distribution of Ancient Eurasian Genomes",
    subtitle = "Holocene Epoch (11,700 BP - Present) | Assumes 29-year generation time (Fenner et al., 2005)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.x.top = element_text(color = "#c0392b", face = "bold", margin = margin(b = 10)),
    axis.text.x.top = element_text(color = "#c0392b"),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    axis.line = element_line(color = "black"),
    # Ajustar márgenes para que los nombres de los periodos no se corten
    plot.margin = margin(t = 20, r = 20, b = 10, l = 10) 
  )

# 5. Visualizar y Guardar
#print(plot_timeline)
ggsave("FigS1_Holocene_Sampling_Timeline.pdf", plot_timeline, width = 12, height = 6)

library(data.table)

# 1. Definir los periodos arqueológicos utilizando fcase (más rápido y limpio que ifelse)
resumen_clst[, Period := fcase(
  Mean_Date_BP <= 11700 & Mean_Date_BP > 8000, "Mesolithic",
  Mean_Date_BP <= 8000 & Mean_Date_BP > 5000, "Neolithic",
  Mean_Date_BP <= 5000 & Mean_Date_BP > 3200, "Bronze Age",
  Mean_Date_BP <= 3200 & Mean_Date_BP > 2000, "Iron Age",
  Mean_Date_BP <= 2000 & Mean_Date_BP >= 0,    "Historical",
  default = NA_character_
)]

# 2. Calcular las coordenadas en la cuadrícula de la PDE (nx=62, ny=18). Cada gradilla mide 2.5°
resumen_clst[, X := as.integer(cut(Mean_Lon, breaks = seq(-15, 55, length.out = 62 + 1), include.lowest = TRUE))]
resumen_clst[, Y := as.integer(cut(Mean_Lat, breaks = seq(30, 75, length.out = 18 + 1), include.lowest = TRUE))]
cat("Muestras en Eurasia:", nrow(resumen_clst), "\n")
# 3. Resumir la información a nivel de población (CLST)
poblaciones_dt <- resumen_clst

# 4. Eliminar poblaciones que hayan quedado sin periodo (ej. fuera de rango de fechas)
poblaciones_dt <- poblaciones_dt[!is.na(Period)]
cat("Muestras en Eurasia:", nrow(poblaciones_dt), "\n")
# 5. Exportar el CSV para tus registros y material suplementario
fwrite(poblaciones_dt, "Eurasia_Populations_Periods_Grid.csv")
cat("Archivo CSV exportado exitosamente con", nrow(poblaciones_dt), "poblaciones.\n")

library(ggplot2)
library(maps)
library(viridis)

# 1. Obtener el mapa base del mundo
world_map <- map_data("world")

# 2. Generar las líneas divisorias exactas de tu PDE
lon_breaks <- seq(-15, 55, length.out = 62)#31 + 1)
lat_breaks <- seq(30, 75, length.out = 18)#9 + 1)

# 3. El truco de la Facetación: Crear el panel "All Holocene"
poblaciones_total <- copy(poblaciones_dt)
poblaciones_total[, Period := "All Holocene"]

# Unir los datos originales (separados) con la copia (agrupada)
plot_data <- rbind(poblaciones_total, poblaciones_dt)

# Ordenar los niveles cronológicamente para que "All Holocene" salga primero,
# seguido del Mesolítico hasta el periodo Histórico.
plot_data[, Period := factor(Period, levels = c(
  "All Holocene", "Mesolithic", "Neolithic", "Bronze Age", "Iron Age", "Historical"
))]

# 4. Construcción del Mapa Principal
plot_maps <- ggplot() +
  
  # A. Dibujar la masa continental
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
               fill = "gray85", color = "white", linewidth = 0.2) +
  
  # B. Dibujar la Cuadrícula Espacial (Grid) de tu PDE
  geom_vline(xintercept = lon_breaks, color = "gray60", linewidth = 0.3, linetype = "dotted") +
  geom_hline(yintercept = lat_breaks, color = "gray60", linewidth = 0.3, linetype = "dotted") +
  
  # C. Colocar las poblaciones (CLST)
  geom_point(data = plot_data, 
             aes(x = Mean_Lon, y = Mean_Lat, color = Period, size = N_Samples), 
             alpha = 0.75, shape = 16) +
  
  # D. Restringir el enfoque a tu Bounding Box de Eurasia
  coord_quickmap(xlim = c(-15, 55), ylim = c(30, 75)) +
  
  # E. Dividir en los 6 paneles (2 filas x 3 columnas)
  facet_wrap(~ Period, ncol = 3) +
  
  # F. Estética: Escala de colores y tamaños
  scale_color_viridis_d(option = "turbo", guide = "none") + # turbo es genial para separar épocas
  scale_size_continuous(range = c(1, 5), name = "Sample Size (N)") +
  
  # G. Limpieza general para formato revista científica
  theme_minimal(base_size = 14) +
  labs(
    title = "Spatiotemporal Discretization of Ancient Eurasian Populations",
    subtitle = "Dotted lines represent the 31x9 spatial grid applied in the Partial Differential Equation (PDE) model.",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    strip.background = element_rect(fill = "gray30", color = NA),
    strip.text = element_text(color = "white", face = "bold", size = 12),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "aliceblue", color = NA), # Color azul agua muy tenue para el mar
    legend.position = "bottom"
  )

# 5. Visualizar y Guardar en alta resolución
print(plot_maps)
ggsave("Fig_2_Spatiotemporal_Grid_Maps.pdf", plot_maps, width = 14, height = 8, dpi = 300)