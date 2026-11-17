# Instalar  y carga paquetes si es necesario
# ------------------------------------------------------------------------------
if (!require("terra")) install.packages("terra")
if (!require("sf")) install.packages("sf")
if (!require("giscoR")) install.packages("giscoR")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("tidyterra")) install.packages("tidyterra")

# Cargar paquetes
library(terra)
library(sf)
library(giscoR)
library(tidyverse)
library(tidyterra)

# Opción A: Datos de copernicus (Si no puedes descargar los datos de Copernicus)
# ------------------------------------------------------------------------------
# Visitar https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php
# Descargar datos 0.1 deg. regular grid, Ensemble mean: RR
# Unos 3 GB (unos 5 minutos en los ordenadores de la sala)

# Importar datos de precipitación de Europa (OJO con el nombre del archivo. Puede ser distinto si descargais otra versión)
precipitation <- rast("rr_ens_mean_0.1deg_reg_v34.0e.nc")

# Recortar mapa
limites_cataluna <- gisco_get_nuts(
  country = "Spain",
  nuts_level = 2,
  resolution = "03",
  nuts_id = "ES51"
)
plot(limites_cataluna)
precipitacion_cat <- crop(precipitation, limites_cataluna) #este proceso podría tardar algunos minutos
writeRaster(precipitacion_cat, filename = "precipitacion_cataluna.tif", overwrite = TRUE) # Guardar archivo por si lo necesitamos
plot(precipitacion_cat)

# Opción B: datos recortados
# ------------------------------------------------------------------------------
# Si no podeis descargar datos de Copernicus, podeis descargar los datos recortados de Catalunya desde aquí:
# https://github.com/angelrodes/Ejercicios_R_2024_25/raw/refs/heads/main/precipitacion_cataluna.tif

# Recuperar datos de precipitación de Cataluña
precipitacion_cat <- rast("precipitacion_cataluna.tif")
tiempo_en_segundos <- precipitacion_cat@pntr[["time"]] # ¿Qué formato de tiempo es?
fechas <- as.Date(as.POSIXct(tiempo_en_segundos, origin="1970-01-01"))
tabla_precip <- as.data.frame(precipitacion_cat, xy = TRUE)
colnames(tabla_precip)[3:ncol(tabla_precip)] <-   format(fechas, "%Y-%m-%d")

# Limitar a datos por fecha
# ------------------------------------------------------------------------------
# Identificar las columnas que NO son 'x' ni 'y'
columnas_fechas <- setdiff(names(tabla_precip), c("x", "y"))
# Extraer el año de esas columnas (primeros 4 caracteres) y pasarlo a número
anios <- as.numeric(substr(columnas_fechas, 1, 4))
# Filtrar cuáles de esas fechas son del año 2030 o anteriores (es decir: todas)
filtrado_fecha <- anios <= 2030
fechas_validas <- columnas_fechas[filtrado_fecha]
# Reconstruir el data frame manteniendo 'x', 'y' y las fechas filtradas
fechas <- fechas[filtrado_fecha]
anios <- anios[filtrado_fecha]
tabla_precip <- tabla_precip[, c("x", "y", fechas_validas)]

# Extensión temporal
intervalo <- interval(min(fechas, na.rm = TRUE), max(fechas, na.rm = TRUE))
intervalo_anios <- time_length(intervalo, unit = "years") # numero de años que tenemos datos

# Precipitación media anual
# ------------------------------------------------------------------------------
prec_mean<-app(precipitacion_cat,mean)*365.25
plot(prec_mean, main = paste("Precipitación media anual - ", intervalo))
isolineas <- as.contour(prec_mean, levels = seq(0, 3000, by = 50))
lines(isolineas, col = "black", lwd = 1)

# Precipitación media anual por décadas
# ------------------------------------------------------------------------------
decadas <- unique(floor(as.numeric(format(fechas, "%Y")) / 10) * 10)
rango_comun <- range(prec_mean, na.rm = TRUE)
for(dec in decadas) {
  precip_decada <- mean(precipitacion_cat[[which(format(fechas, "%Y") >= dec & format(fechas, "%Y") <= dec+9)]]) * 365.25
  plot(precip_decada, main = paste("Precipitación Media Anual - Década", dec))
}

# Calculamos los periodos de retorno para todos los puntos del mapa
# ------------------------------------------------------------------------------

# Identificamos cuáles son las columnas de fechas (todas excepto 'x' e 'y')
columnas_fechas <- setdiff(names(tabla_precip), c("x", "y"))
numero_filas <- nrow(tabla_precip)

# Generamos variables para rellenar
resultados_Tr10 <- numeric(nrow(tabla_precip))
resultados_Tr100 <- numeric(nrow(tabla_precip))
resultados_Tr500 <- numeric(nrow(tabla_precip))


# Iniciamos los cálculos de Periodos de retorno para CADA PUNTO del mapa
for (i in 1:numero_filas) {
porcentaje <- (i / numero_filas) * 100
cat(sprintf("\rProgreso: %.1f%% (%d de %d)", porcentaje, i, numero_filas))
precipitaciones_punto <- as.numeric(tabla_precip[i, columnas_fechas]) # Datos de precipitación de un punto
P_obs <- sort(precipitaciones_punto) # Datos ordenados
n <- length(P_obs)
total_anios <- (n /365)
Tr_obs=total_anios/ (1+n-(1:n))
Tr_mas_de_un_anno <- Tr_obs>1 # Filtrar solo periodos de retorno de más de un año
P_obs <- P_obs[Tr_mas_de_un_anno]
n <- length(P_obs)
Tr_obs <- Tr_obs[Tr_mas_de_un_anno]
# Ajustar la regresión lineal usando el logaritmo de Tr_obs
modelo <- lm(P_obs ~ log(Tr_obs))
x_seq <- c(10, 100, 500) # Periodos de retorno de interés
y_pred <- predict(modelo, newdata = data.frame(Tr_obs = x_seq))

# Gurdamos los datos que queremos representar
resultados_Tr10[i] <- y_pred[1] 
resultados_Tr100[i] <- y_pred[2] 
resultados_Tr500[i] <- y_pred[3] 
}

# Generamos las tablas de resultados con x,yTr
mapa_prediccion_Tr10 <- data.frame(
  x = tabla_precip$x,
  y = tabla_precip$y,
  resultado = resultados_Tr10
)

mapa_prediccion_Tr100 <- data.frame(
  x = tabla_precip$x,
  y = tabla_precip$y,
  resultado = resultados_Tr100
)

mapa_prediccion_Tr500 <- data.frame(
  x = tabla_precip$x,
  y = tabla_precip$y,
  resultado = resultados_Tr500
)

# Los transformamos a mapas SpatRaster
raster_resultados10 <- rast(mapa_prediccion_Tr10, type = "xyz")
raster_resultados100 <- rast(mapa_prediccion_Tr100, type = "xyz")
raster_resultados500 <- rast(mapa_prediccion_Tr500, type = "xyz")

# Representamos los mapas

# Nos aseguramos que los mismos valores corresponden a los mismos colores para poder comparar mapas
# Definir los cortes exactos desde 0 hasta 300
cortes <- seq(0, 300, length.out = 31)
paleta_precip <- colorRampPalette(c("forestgreen", "yellow", "orange", "red"))(100)

plot(raster_resultados10, main = "Precipitación máxima diaria T10",  col = paleta_precip, breaks = cortes)
isolineas <- as.contour(raster_resultados10, levels = seq(0, 3000, by = 10))
lines(isolineas, col = "black", lwd = 1)

plot(raster_resultados100, main = "Precipitación máxima diaria T100",  col = paleta_precip, breaks = cortes)
isolineas <- as.contour(raster_resultados100, levels = seq(0, 3000, by = 10))
lines(isolineas, col = "black", lwd = 1)

plot(raster_resultados500, main = "Precipitación máxima diaria T500",  col = paleta_precip, breaks = cortes)
isolineas <- as.contour(raster_resultados500, levels = seq(0, 3000, by = 10))
lines(isolineas, col = "black", lwd = 1)

# Ejercicios:
# - Comparar mapa de precipitación media anual y periodos de retorno de precipitación máxima diaria
# - Comparar mpas de peridos de retorno con los datos de INUNCAT
# - Repetir unsando datos de distinto rango de fechas y observar si se puede deducir efectos del cambio climático


