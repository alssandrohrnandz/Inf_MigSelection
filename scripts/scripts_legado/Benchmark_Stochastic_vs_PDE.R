# ==============================================================================
# COMPARACIÓN DE EFICIENCIA: ESTOCÁSTICO (MONTE CARLO) VS PDE + BETA-BINOMIAL
# ==============================================================================

# 1. Parámetros Globales
grid_size <- 10        # Grilla NxN
generations <- 100      # Generaciones a simular
m <- 0.05              # Tasa de migración por vecino
s <- 0.01              # Coeficiente de selección
Ne <- 1000             # Tamaño de subpoblación
replicas_mc <- 100     # Réplicas para el método Monte Carlo

# Crear matriz inicial (Mutación surge en el centro)
P_init <- matrix(0, nrow = grid_size, ncol = grid_size)
P_init[5, 5] <- 0.1  

# Datos simulados "observados" al final de la generación 50 para evaluar verosimilitud
# (Asumimos conteos hipotéticos para el ejemplo)
conteos_obs <- matrix(rbinom(grid_size^2, 2*Ne, P_init), grid_size, grid_size)
cromosomas_totales <- matrix(2*Ne, grid_size, grid_size)

# ==============================================================================
# FUNCIÓN AUXILIAR: Difusión Espacial 2D rápida (Laplaciano Discreto)
# ==============================================================================
difusion_2D <- function(P, m) {
  n <- nrow(P)
  # Bordes reflectantes (simplificados) para mantener masa
  P_up    <- rbind(P[1, ], P[1:(n-1), ])
  P_down  <- rbind(P[2:n, ], P[n, ])
  P_left  <- cbind(P[, 1], P[, 1:(n-1)])
  P_right <- cbind(P[, 2:n], P[, n])
  
  # Matriz actualizada por migración
  P_new <- P * (1 - 4*m) + m * (P_up + P_down + P_left + P_right)
  return(P_new)
}

# ==============================================================================
# ENFOQUE 1: SIMULACIÓN ESTOCÁSTICA EXPLÍCITA (MONTE CARLO)
# Simula la deriva (Fokker-Planck implícito) en cada paso espacial y temporal.
# ==============================================================================
run_stochastic_approach <- function() {
  operaciones <- 0
  
  # Para calcular la verosimilitud, necesitamos correr muchas réplicas
  # y promediar los resultados espaciales.
  P_sum <- matrix(0, grid_size, grid_size)
  
  for (rep in 1:replicas_mc) {
    P_curr <- P_init
    for (t in 1:generations) {
      # 1. Migración y Selección (Determinista temporal)
      P_mig <- difusion_2D(P_curr, m)
      P_sel <- P_mig + (s * P_mig * (1 - P_mig))
      P_sel <- pmax(pmin(P_sel, 1), 0) # Limitar entre 0 y 1
      
      # 2. Deriva Génica (Estocástica - Muestreo Binomial por deme)
      # Aquí está el cuello de botella computacional
      P_curr <- matrix(rbinom(grid_size^2, 2*Ne, P_sel) / (2*Ne), grid_size, grid_size)
      operaciones <- operaciones + grid_size^2
    }
    P_sum <- P_sum + P_curr
  }
  
  # Promedio de frecuencias esperadas tras R réplicas
  P_expected <- P_sum / replicas_mc
  
  # Verosimilitud Binomial estándar asumiendo el promedio
  ll <- sum(dbinom(as.vector(conteos_obs), as.vector(cromosomas_totales), 
                   as.vector(P_expected), log = TRUE), na.rm = TRUE)
  
  return(list(logLik = ll, ops_deriva = operaciones))
}

# ==============================================================================
# ENFOQUE 2: Fisher-KPP (PDE) + BETA-BINOMIAL (VEROSIMILITUD ANALÍTICA)
# Tu propuesta: Determinar esperanza matemática 1 vez, usar matemática para el ruido.
# ==============================================================================
# Función auxiliar Beta-Binomial
dbbinom_custom <- function(x, size, alpha, beta) {
  # Log-likelihood de la Beta-Binomial usando lgamma para evitar desbordamientos
  log_num <- lgamma(size + 1) + lgamma(x + alpha) + lgamma(size - x + beta) + lgamma(alpha + beta)
  log_den <- lgamma(x + 1) + lgamma(size - x + 1) + lgamma(size + alpha + beta) + lgamma(alpha) + lgamma(beta)
  return(log_num - log_den)
}

run_analytical_approach <- function() {
  operaciones <- 0
  P_curr <- P_init
  
  # Se corre UNA sola vez
  for (t in 1:generations) {
    # Migración y Selección
    P_mig <- difusion_2D(P_curr, m)
    P_curr <- P_mig + (s * P_mig * (1 - P_mig))
    P_curr <- pmax(pmin(P_curr, 1), 0)
    # 0 llamadas a rbinom (cero simulaciones de deriva)
  }
  
  # Absorber ruido en la Verosimilitud Beta-Binomial
  theta_drift <- 4 * Ne # Ajuste teórico de varianza
  alpha_param <- P_curr * theta_drift
  beta_param  <- (1 - P_curr) * theta_drift
  
  # Prevenir ceros absolutos en parámetros beta
  alpha_param <- pmax(alpha_param, 1e-6)
  beta_param  <- pmax(beta_param, 1e-6)
  
  ll <- sum(dbbinom_custom(as.vector(conteos_obs), as.vector(cromosomas_totales), 
                           as.vector(alpha_param), as.vector(beta_param)))
  
  operaciones <- 1 # Se resolvió la EDP 1 vez
  
  return(list(logLik = ll, ops_deriva = operaciones))
}

# ==============================================================================
# EJECUCIÓN Y COMPARACIÓN DE TIEMPOS
# ==============================================================================
cat("\n--- INICIANDO BENCHMARK ---\n\n")

# Medir Enfoque 1 (Monte Carlo Estocástico)
start_time <- Sys.time()
res_stochastic <- run_stochastic_approach()
end_time <- Sys.time()
time_stochastic <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Medir Enfoque 2 (EDP + Beta Binomial)
start_time <- Sys.time()
res_analytical <- run_analytical_approach()
end_time <- Sys.time()
time_analytical <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Resultados
cat(sprintf("1. MÉTODOS ESTOCÁSTICOS MÚLTIPLES (Monte Carlo con %d réplicas)\n", replicas_mc))
cat(sprintf("   - Tiempo de ejecución : %.4f segundos\n", time_stochastic))
cat(sprintf("   - Evaluaciones de red : %d veces (Cuello de botella)\n", res_stochastic$ops_deriva))
cat(sprintf("   - Log-Verosimilitud   : %.2f\n\n", res_stochastic$logLik))

cat("2. MÉTODO FISHER-KPP + BETA-BINOMIAL (Tu enfoque)\n")
cat(sprintf("   - Tiempo de ejecución : %.4f segundos\n", time_analytical))
cat(sprintf("   - Evaluaciones de red : 0 simulaciones estocásticas (1 resolución analítica)\n"))
cat(sprintf("   - Log-Verosimilitud   : %.2f\n\n", res_analytical$logLik))

aceleracion <- time_stochastic / time_analytical
cat(sprintf(">> CONCLUSIÓN: Tu enfoque PDE+Beta-Binomial es aprox. %.1f VECES MÁS RÁPIDO.\n", aceleracion))
cat("   Imagina multiplicar esto por los 1071 puntos de tu Grid Search.\n")

# ==============================================================================
# VISUALIZACIÓN ESPACIAL: Estocástico vs Analítico (Fisher-KPP)
# ==============================================================================

# 1. Parámetros Globales
grid_size <- 10
generations <- 100
m <- 0.05
s <- 0.01
Ne <- 1000
replicas_mc <- 100

# Matriz inicial (mutación al centro)
P_init <- matrix(0, nrow = grid_size, ncol = grid_size)
P_init[5, 5] <- 0.1  

# Función de difusión
difusion_2D <- function(P, m) {
  n <- nrow(P)
  P_up    <- rbind(P[1, ], P[1:(n-1), ])
  P_down  <- rbind(P[2:n, ], P[n, ])
  P_left  <- cbind(P[, 1], P[, 1:(n-1)])
  P_right <- cbind(P[, 2:n], P[, n])
  return(P * (1 - 4*m) + m * (P_up + P_down + P_left + P_right))
}

# 2. Obtener matriz final - ESTOCÁSTICO (Promedio de Monte Carlo)
P_sum <- matrix(0, grid_size, grid_size)
for (rep in 1:replicas_mc) {
  P_curr <- P_init
  for (t in 1:generations) {
    P_mig <- difusion_2D(P_curr, m)
    P_sel <- P_mig + (s * P_mig * (1 - P_mig))
    P_sel <- pmax(pmin(P_sel, 1), 0)
    P_curr <- matrix(rbinom(grid_size^2, 2*Ne, P_sel) / (2*Ne), grid_size, grid_size)
  }
  P_sum <- P_sum + P_curr
}
P_estocastico_final <- P_sum / replicas_mc

# 3. Obtener matriz final - ANALÍTICO (EDP Fisher-KPP)
P_analitico_final <- P_init
for (t in 1:generations) {
  P_mig <- difusion_2D(P_analitico_final, m)
  P_analitico_final <- P_mig + (s * P_mig * (1 - P_mig))
  P_analitico_final <- pmax(pmin(P_analitico_final, 1), 0)
}

# ==============================================================================
# GRAFICACIÓN (Mapas de calor)
# ==============================================================================

# Rota la matriz 90 grados para que la función image() la dibuje 
# tal cual se orienta matemáticamente (eje X e Y)
rotar_matriz <- function(m) { t(m)[, nrow(m):1] }

# Definir colores (del amarillo claro al rojo oscuro)
paleta_colores <- hcl.colors(20, "YlOrRd", rev = TRUE)

# Escala máxima compartida para que los colores signifiquen lo mismo en ambas gráficas
max_freq <- max(P_estocastico_final, P_analitico_final)

# Configurar la ventana gráfica (1 fila, 2 columnas)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

# Gráfico 1: Estocástico
image(1:grid_size, 1:grid_size, rotar_matriz(P_estocastico_final), 
      col = paleta_colores, zlim = c(0, max_freq),
      main = "Estocástico (Monte Carlo)",
      xlab = "Coordenada X", ylab = "Coordenada Y")
contour(1:grid_size, 1:grid_size, rotar_matriz(P_estocastico_final), 
        add = TRUE, col = "black", alpha = 0.5)

# Gráfico 2: Analítico
image(1:grid_size, 1:grid_size, rotar_matriz(P_analitico_final), 
      col = paleta_colores, zlim = c(0, max_freq),
      main = "Analítico (Fisher-KPP)",
      xlab = "Coordenada X", ylab = "Coordenada Y")
contour(1:grid_size, 1:grid_size, rotar_matriz(P_analitico_final), 
        add = TRUE, col = "black", alpha = 0.5)

# Restaurar la configuración de la ventana gráfica
par(mfrow = c(1, 1))
