# ==============================================================================
# PERIODOS DE RETORNO DE PRECIPITACIÓN
# ==============================================================================

# Carga paquetes necesarion

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("evd")) install.packages("evd")
library(ggplot2)
library(dplyr)
library(evd)

# 1. Carga de datos diarios
# ------------------------------------------------------------------------------
df_diario <- read.csv("VALORES ESTACION METEOROLOGICA REUS AEROPUERTO.csv", 
                  check.names = FALSE, 
                  encoding = "UTF-8")


# Extraer y convertir las fechas al formato Date de R (día-mes-año)
fechas <- as.Date(df_diario$FECHA, format = "%d-%m-%Y")

n_dias <- length(fechas)

# Recorta datos si es necesario
df_diario <- subset(df_diario, 
       fechas >= as.Date("1940-01-01") & FECHA <= as.Date("2026-12-31"))
fechas <- as.Date(df_diario$FECHA, format = "%d-%m-%Y")

# Extraer las precipitaciones asegurándote de que sean numéricas
lluvia_diaria <- as.numeric(df_diario[["Precipitación: l/m2"]])

df_diario <- data.frame(Fecha = fechas, Precip_Diaria = round(lluvia_diaria, 1))

# 2. Extracción de Máximos Anuales
# ------------------------------------------------------------------------------
df_diario$Ano <- format(df_diario$Fecha, "%Y")

df_anual <- df_diario %>%
  group_by(Ano) %>%
  summarise(Prec_Max = max(Precip_Diaria, na.rm = TRUE)) %>%
  ungroup()

P <- sort(df_anual$Prec_Max) # Datos ordenados
n <- length(P)
T_target <- c(10, 100, 500)  # Periodos de retorno objetivo

media_P <- mean(P)
sd_P    <- sd(P)
cv_P    <- sd_P / media_P

# 3. Posición de ploteo empírica (Weibull)
# ------------------------------------------------------------------------------
df_observado <- data.frame(
  T_obs = 1 / (1 - (1:n) / (n + 1)),
  P_obs = P
)

# 4. Ajustes de Distribuciones
# ------------------------------------------------------------------------------

# --- MÉTODO 1: Gumbel ---
alpha_gumbel <- (sqrt(6) / pi) * sd_P
beta_gumbel  <- media_P - 0.5772156649 * alpha_gumbel
calc_gumbel  <- function(T_val) beta_gumbel - alpha_gumbel * log(-log(1 - 1 / T_val))

# --- MÉTODO 2: Gauss (Normal) ---
calc_gauss   <- function(T_val) qnorm(1 - (1 / T_val), mean = media_P, sd = sd_P)

# --- MÉTODO 3: SQRT-ET max (Etoh et al., 1986) ---
# Funciones para el cálculo numérico de momentos de SQRT-ET max
momento_u <- function(u, k, pow) {
  # Sustitución u = sqrt(alpha * x)
  k * (u^pow) * exp(-u) * exp(-k * (1 + u) * exp(-u))
}

# Coeficiente de variación teórico en función de k (independiente de alpha)
cv_sqrt_et <- function(k) {
  # Usamos sub-intervalos para mejorar la tolerancia de integración en k elevados
  m1 <- integrate(momento_u, 0, Inf, k = k, pow = 3, rel.tol = 1e-7)$value
  m2 <- integrate(momento_u, 0, Inf, k = k, pow = 5, rel.tol = 1e-7)$value
  var <- m2 - (m1^2)
  if (var < 0) return(NA) # Evitar imprecisiones numéricas en límites extremos
  return(sqrt(var) / m1)
}

# A. Hallar k tal que Cv_teórico(k) == Cv_observado (Intervalo ampliado + extensión automática)
obj_k <- function(k) cv_sqrt_et(k) - cv_P

# Buscamos k en un rango amplio de 1e-4 a 1e7
sol_k <- uniroot(obj_k, interval = c(1e-4, 1e7), extendInt = "yes")
k_opt <- sol_k$root

# B. Hallar alpha tal que Media_teórica(k, alpha) == Media_observada
m1_ref <- integrate(momento_u, 0, Inf, k = k_opt, pow = 3, rel.tol = 1e-7)$value
alpha_opt <- m1_ref / media_P

# C. Función para resolver el cuantil x_T para SQRT-ET max
calc_sqrt_et <- function(T_val) {
  target <- -log(1 - (1 / T_val))
  # Resolver para u = sqrt(alpha * x) en: k * (1 + u) * exp(-u) = target
  eq_u <- function(u) k_opt * (1 + u) * exp(-u) - target
  
  # Límite superior amplio para u en función de k
  u_max <- max(50, log(k_opt * 100))
  u_val <- uniroot(eq_u, interval = c(0, u_max), extendInt = "yes")$root
  x_val <- (u_val^2) / alpha_opt
  return(x_val)
}

# 5. Cálculo para T = 10, 100 y 500 años
# ------------------------------------------------------------------------------
P_gumbel  <- sapply(T_target, calc_gumbel)
P_gauss   <- sapply(T_target, calc_gauss)
P_sqrt_et <- sapply(T_target, calc_sqrt_et)

resultados <- data.frame(
  Periodo_Retorno_Anos = T_target,
  Gumbel_mm      = round(P_gumbel, 2),
  Gauss_mm       = round(P_gauss, 2),
  SQRT_ET_max_mm = round(P_sqrt_et, 2)
)

print("=== PRECIPITACIÓN MÁXIMA ESTIMADA (mm) ===")
print(resultados)

# 6. Representación Gráfica
# ------------------------------------------------------------------------------
T_seq <- exp(seq(log(1.05), log(500), length.out = 400))

df_curvas <- data.frame(
  T = rep(T_seq, 3),
  P = c(sapply(T_seq, calc_gumbel), 
        sapply(T_seq, calc_gauss), 
        sapply(T_seq, calc_sqrt_et)),
  Metodo = factor(rep(c("Gumbel", "Gauss", "SQRT-ET max"), each = length(T_seq)))
)

df_puntos_destacados <- data.frame(
  T = rep(T_target, 3),
  P = c(P_gumbel, P_gauss, P_sqrt_et),
  Metodo = factor(rep(c("Gumbel", "Gauss", "SQRT-ET max"), each = length(T_target)))
)

g <- ggplot() +
  geom_point(data = df_observado, aes(x = T_obs, y = P_obs), 
             color = "black", size = 2, alpha = 0.6) +
  geom_line(data = df_curvas, aes(x = T, y = P, color = Metodo), linewidth = 1) +
  geom_point(data = df_puntos_destacados, aes(x = T, y = P, color = Metodo), 
             size = 3.5, shape = 18) +
  scale_x_log10(breaks = c(1/365,1,2, 5, 10, 25, 50, 100, 250, 500), limits = c(1, 550)) +
  scale_color_manual(values = c("Gumbel" = "#E41A1C", "Gauss" = "#377EB8", "SQRT-ET max" = "#4DAF4A")) +
  labs(
    title = "Precipitaciones Máximas Anuales (Reus Aeropuerto)",
    subtitle = "Ajuste SQRT-ET max (Etoh et al., 1986), Gumbel y Gauss",
    x = "Periodo de Retorno T (Años - Escala Logarítmica)",
    y = "Precipitación Máxima Diaria (mm)",
    color = "Método de Ajuste",
    caption = "Puntos negros: Datos observados."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_line(linewidth = 0.2, color = "gray90"),
    panel.grid.major = element_line(linewidth = 0.4, color = "gray80")
  )

print(g)



