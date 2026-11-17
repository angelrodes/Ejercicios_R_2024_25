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

# Carga de datos diarios
# ------------------------------------------------------------------------------
df_diario <- read.csv("VALORES ESTACION METEOROLOGICA REUS AEROPUERTO.csv", 
                  check.names = FALSE, 
                  encoding = "UTF-8")


# Extraer y convertir las fechas al formato Date de R (día-mes-año)
fechas <- as.Date(df_diario$FECHA, format = "%d-%m-%Y")

# Recorta datos si es necesario
df_diario <- subset(df_diario, 
       fechas >= as.Date("1940-01-01") & FECHA <= as.Date("2026-12-31"))
fechas <- as.Date(df_diario$FECHA, format = "%d-%m-%Y")

# Extraer las precipitaciones asegurándote de que sean numéricas
lluvia_diaria <- as.numeric(df_diario[["Precipitación: l/m2"]])

# Un vistazo a los datos
# ------------------------------------------------------------------------------

# datos "en bruto"
plot(fechas,lluvia_diaria)

# precipitación media anual
anios <- as.numeric(format(fechas, "%Y"))
lluvia_anual <- tapply(lluvia_diaria, anios, sum, na.rm = TRUE)
lista_anios <- as.numeric(names(lluvia_anual))
plot(lista_anios, lluvia_anual,      type = "b")

# Días con precipitaciones intensas
dias_mas_50 <- tapply(lluvia_diaria > 50, anios, sum, na.rm = TRUE)
plot(lista_anios, dias_mas_50,
     type = "h", lwd = 3, col = "darkred",
     main = "Días al año con precipitación > 50 l/m²",
     xlab = "Año",
     ylab = "Número de días",
     las = 1)

# Periodos de retorno observados
# ------------------------------------------------------------------------------
P_obs <- sort(lluvia_diaria) # Datos ordenados
n <- length(P_obs)
total_anios <- (n /365)

Tr_obs=total_anios/ (1+n-(1:n))

plot(Tr_obs,P_obs, 
     log = "x")

# Filtrar solo periodos de retorno de más de un año
Tr_mas_de_un_anno <- Tr_obs>1
P_obs <- P_obs[Tr_mas_de_un_anno]
n <- length(P_obs)
Tr_obs <- Tr_obs[Tr_mas_de_un_anno]

plot(Tr_obs,P_obs, 
     log = "x")

# Estimación lineal (con x logarítmico)
# Extracción de Máximos Anuales
# ------------------------------------------------------------------------------
# Ajustar la regresión lineal usando el logaritmo de Tr_obs
modelo <- lm(P_obs ~ log(Tr_obs))


# Dibujar el gráfico base ampliando xlim hasta 1000 para que quepan todos los puntos de interés
plot(Tr_obs, P_obs,
     log = "x",
     xlim = c(min(Tr_obs), 1000),
     ylim = c(min(P_obs), predict(modelo, newdata = data.frame(Tr_obs = 1000))),
     pch = 19,
     cex = 1.5,       # Tamaño mayor para resaltarlos
     col = "steelblue",
     main = "Curva de Periodo de Retorno",
     xlab = "Periodo de retorno Tr (años, escala log)",
     ylab = "Precipitación (l/m²)",
     las = 1)
grid()

# Dibujar la línea de tendencia recta sobre el eje logarítmico
x_seq <- seq(min(Tr_obs), 500, length.out = 200)
y_pred <- predict(modelo, newdata = data.frame(Tr_obs = x_seq))
lines(x_seq, y_pred, col = "red", lwd = 2)

# Calcular y resaltar los puntos para Tr = 10, 100 y 500 años
tr_destacados <- c(10, 100, 500)
p_destacados <- predict(modelo, newdata = data.frame(Tr_obs = tr_destacados))

points(tr_destacados, p_destacados, 
       pch = 19,        # Círculo relleno
       col = "darkred")

# Añadir etiquetas con el valor numérico junto a cada punto resaltado
text(tr_destacados, p_destacados,
     labels = paste0("Tr=", tr_destacados, ": ", round(p_destacados, 1), " l/m²"),
     pos = 3,           # Posición: 3 = encima del punto (usa 2 para izquierda, 4 para derecha)
     cex = 0.85,
     font = 2,
     col = "darkred")



