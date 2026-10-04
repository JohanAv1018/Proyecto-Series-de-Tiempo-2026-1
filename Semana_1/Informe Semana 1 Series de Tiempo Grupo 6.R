# ==========================================================
# INFORME SEMANA 1 SERIES DE TIEMPO GRUPO 6
#
# Ana Sofía Salazar Álvarez
# Winston Obeymar Lucano Villota
# Julian Camilo Tobo Tenen
# Johan Steven Avilan Peñaloza
#
# Universidad Nacional de Colombia - Departamento de Estadística
#
# Serie: MSCI COLCAP (cierre diario, en puntos)
# Fecha de corte de los datos: viernes 2 de octubre de 2026
# Pronóstico: lunes 5 al viernes 9 de octubre de 2026 (h = 1,...,5)
#
# Modelos: ETS, Naive, sNaive (lag 5) y Drift
#
# Validación: holdout repetido sin traslape (bloques de h = 5),
#             ventana expansible y fija,
#             con y sin recalibración
#
# Métrica: RMSE por paso y RMSE multi-paso
#
# Datos: Banco de la República
# ==========================================================


# ==========================================================
# PREELIMINARES
# ==========================================================

install.packages("pacman")  #Solo una vez
library(pacman)
p_load(fpp3, dplyr, tidyr, purrr, readxl, lubridate, readr, ggplot2)

# Mostrar decimales en las tablas (por defecto tibble muestra solo 3 cifras significativas)
options(pillar.sigfigs = 7)


# ==========================================================
# 1. DATOS
# ==========================================================

# Descarga del archivo de datos
url <- "https://raw.githubusercontent.com/JohanAv1018/Proyecto-Series-de-Tiempo-2026-1/main/Semana_1/colcap_semana1.xlsx"
tmp <- tempfile(fileext = ".xlsx")
download.file(url, tmp, mode = "wb")


# Carga y transformación de los datos
datos <- read_excel(tmp) %>%
  rename(fecha = 1, colcap = 2) %>%
  mutate(
    fecha  = dmy(fecha),
    colcap = parse_number(as.character(colcap),
                          locale = locale(decimal_mark = ",", grouping_mark = "."))
  ) %>%
  arrange(fecha) %>%
  filter(fecha <= as.Date("2026-10-02"))


# ==========================================================
# 2. PARÁMETROS GENERALES Y PREPARACION DE LA SERIE
# ==========================================================

inicio_ventana <- as.Date("2021-05-28")
h <- 5      
n_bloques_eval <- 100    
ancho_fijo <- 500 

serie <- datos %>%
  filter(fecha >= inicio_ventana) %>%
  mutate(t = row_number()) %>%
  as_tsibble(index = t)

autoplot(serie, colcap)

inicio_train <- nrow(serie) - h * n_bloques_eval

# ==========================================================
# 3. FUNCIÓN GENERAL
# ==========================================================

evaluar_modelo <- function(serie,
                           tipo_ventana = "expansible",
                           recalibrar = TRUE,
                           modelo = "ETS") {
  
  y <- serie$colcap
  n <- length(y)
  
  errores <- list()
  
  i <- 1
  t0 <- inicio_train
  
  while((t0 + h) <= n) {
    
    # -------------------------
    # Definir ventana
    # -------------------------
    
    if(tipo_ventana == "expansible") {
      
      inicio <- 1
      fin <- t0
      
    }
    
    if(tipo_ventana == "fija") {
      
      inicio <- max(1, t0 - ancho_fijo + 1)
      fin <- t0
      
    }
    
    train <- serie[inicio:fin, ]
    
    # -------------------------
    # Ajustar modelo
    # -------------------------
    
    if(recalibrar || i == 1) {
      
      if(modelo == "ETS") {
        
        fit <- train |>
          model(ETS(colcap))
        
      }
      
      if(modelo == "Naive") {
        
        fit <- train |>
          model(NAIVE(colcap))
        
      }
      
      if(modelo == "sNaive") {
        
        fit <- train |>
          model(SNAIVE(colcap ~ lag(5)))
        
      }
      
      if(modelo == "Drift") {
        
        fit <- train |>
          model(RW(colcap ~ drift()))
        
      }
      
    } else {
      
      fit <- refit(fit, train, reestimate = FALSE)
      
    }
    
    # -------------------------
    # Pronóstico
    # -------------------------
    
    fc <- forecast(fit, h = h)
    
    y_real <- y[(t0+1):(t0+h)]
    
    y_pred <- fc$.mean
    
    # -------------------------
    # Error cuadrático
    # -------------------------
    
    err2 <- (y_real - y_pred)^2
    
    errores[[i]] <- err2
    
    t0 <- t0 + h
    
    i <- i + 1
    
  }
  
  errores_mat <- do.call(rbind, errores)
  
  rmse_horizonte <- sqrt(colMeans(errores_mat))
  
  rmse_global <- sqrt(mean(errores_mat))
  
  return(list(
    rmse_horizonte = rmse_horizonte,
    rmse_global = rmse_global
  ))
  
}

# ==========================================================
# 4. ESCENARIOS
# ==========================================================

escenarios <- tribble(
  ~nombre, ~ventana, ~recalibrar, ~modelo,
  "ETS_expansible_recalibra", "expansible", TRUE, "ETS",
  
  "ETS_expansible_fijo", "expansible", FALSE, "ETS",
  
  "ETS_fija_recalibra", "fija", TRUE, "ETS",
  
  "ETS_fija_fijo", "fija", FALSE, "ETS",
  
  "Naive_expansible", "expansible", TRUE, "Naive",
  
  "sNaive_expansible", "expansible", TRUE, "sNaive",
  
  "sNaive_fija", "fija", TRUE, "sNaive",
  
  "Drift_expansible", "expansible", TRUE, "Drift",
  
  "Drift_fija", "fija", TRUE, "Drift"
)



# ==========================================================
# 5. EJECUCIÓN
# ==========================================================

resultados <- list()

for(i in 1:nrow(escenarios)) {
  
  esc <- escenarios[i, ]
  
  cat("Ejecutando:", esc$nombre, "\n")
  
  res <- evaluar_modelo(
    serie,
    tipo_ventana = esc$ventana,
    recalibrar = esc$recalibrar,
    modelo = esc$modelo
  )
  
  resultados[[esc$nombre]] <- res
  
}


# ==========================================================
# 6. TABLA FINAL
# ==========================================================

tabla_final <- map_dfr(
  names(resultados),
  function(nombre) {
    
    res <- resultados[[nombre]]
    
    tibble(
      Modelo = nombre,
      
      RMSE_h1 = res$rmse_horizonte[1],
      RMSE_h2 = res$rmse_horizonte[2],
      RMSE_h3 = res$rmse_horizonte[3],
      RMSE_h4 = res$rmse_horizonte[4],
      RMSE_h5 = res$rmse_horizonte[5],
      
      RMSE_global = res$rmse_global
    )
    
  }
)

print(tabla_final, n = Inf)


# ==========================================================
# 7. MODELO GANADOR
# ==========================================================

mejor_modelo <- tabla_final |>
  arrange(RMSE_global) |>
  slice(1)

print(mejor_modelo)

nombre_ganador <- mejor_modelo$Modelo
nombre_ganador

# ==========================================================
# 8. AJUSTE FINAL
# ==========================================================

if(grepl("fija", nombre_ganador)) {
  
  serie_final <- tail(serie, ancho_fijo)
  
} else {
  
  serie_final <- serie
  
}

if(grepl("ETS", nombre_ganador)) {
  
  modelo_final <- serie_final |>
    model(ETS(colcap))
  
}

if(grepl("Naive", nombre_ganador)) {
  
  modelo_final <- serie_final |>
    model(NAIVE(colcap))
  
}

if(grepl("sNaive", nombre_ganador)) {
  
  modelo_final <- serie_final |>
    model(SNAIVE(colcap ~ lag(5)))
  
}

if(grepl("Drift", nombre_ganador)) {
  
  modelo_final <- serie_final |>
    model(RW(colcap ~ drift()))
  
}

forecast_final <- forecast(modelo_final, h = h)

# Pronósticos con fecha (lunes 5 a viernes 9 de octubre de 2026)
pronosticos <- tibble(
  fecha = seq(as.Date("2026-10-05"), as.Date("2026-10-09"), by = "day"),
  pronostico = round(forecast_final$.mean, 2)
)

# Se imprime como texto con 2 decimales fijos (el tibble recorta decimales al mostrar numericos)
print(pronosticos |> mutate(pronostico = sprintf("%.2f", pronostico)), n = Inf)

# ==========================================================
# 9. GRÁFICO FINAL
# ==========================================================

autoplot(forecast_final, serie) +
  labs(
    title = "Pronóstico final (5 pasos adelante)",
    y = "Colcap"
  )



