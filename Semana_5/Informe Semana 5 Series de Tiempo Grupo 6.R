# ============================================================
# COLCAP (MSCI, desde 2021-05-28): suavizamiento/ETS + ARMA/ARIMA/SARIMA
# ============================================================

library(pacman)
p_load(dplyr, readxl, forecast, lubridate, ggplot2, tseries)

FECHA_CAMBIO <- as.Date("2021-05-28")
W_FINAL      <- NULL   

# ---- 1. Datos ----

datos <- read_xlsx("Semana_5/colcap_semana5.xlsx", sheet = 1, col_types = "text")
names(datos) <- c("Fecha", "COLCAP")
datos$Fecha  <- suppressWarnings(dmy(datos$Fecha))
datos        <- datos |> filter(!is.na(Fecha))

limpiar_colcap <- function(x) {
  ifelse(grepl(",", x),
         as.numeric(gsub(",", ".", gsub("\\.", "", x))),
         as.numeric(x))
}

datos$COLCAP <- limpiar_colcap(datos$COLCAP)
datos <- datos |> arrange(Fecha) |> filter(Fecha >= FECHA_CAMBIO)

stopifnot(sum(is.na(datos$COLCAP)) == 0)
stopifnot(all(datos$COLCAP > 100 & datos$COLCAP < 5000))

tail(datos, 5)
nrow(datos)

# ---- 2. Identificacion: estacionariedad y ACF/PACF ----

y_log <- log(datos$COLCAP)
ret   <- diff(y_log)

adf.test(y_log)
kpss.test(y_log)
ndiffs(y_log, test = "kpss")

adf.test(ret)
kpss.test(ret)

ggtsdisplay(ret, lag.max = 30)

nsdiffs(ts(y_log, frequency = 5))

# ---- 3. Candidatos ----

f_naive  <- function(x, h) naive(as.numeric(x), h = h)
f_drift  <- function(x, h) rwf(as.numeric(x), h = h, drift = TRUE)
f_ets    <- function(x, h) forecast(ets(as.numeric(x)), h = h)

f_arima <- function(x, h) {
  fit <- auto.arima(as.numeric(x), d = 1, seasonal = FALSE, lambda = 0)
  forecast(fit, h = h)
}

f_arima210 <- function(x, h) {
  fit <- Arima(as.numeric(x), order = c(2, 1, 0), include.drift = TRUE, lambda = 0)
  forecast(fit, h = h)
}

modelos <- list(
  naive = f_naive, drift = f_drift, ets = f_ets,
  arima = f_arima, arima210 = f_arima210
)

# ---- 4. Validacion ----

generar_origenes <- function(n, h, estrategia, gap = 0, min_train = 250) {
  ultimos <- n - h - gap
  if (estrategia == "rep_no_overlap")     return(rev(seq(ultimos, min_train, by = -h)))
  if (estrategia == "rep_no_overlap_gap") return(rev(seq(ultimos, min_train, by = -(h + gap))))
}

evaluar_modelo <- function(y, f_modelo, h, origenes, ventana = NULL) {
  errores <- t(sapply(origenes, function(t) {
    inicio <- if (is.null(ventana)) 1 else max(1, t - ventana + 1)
    train  <- y[inicio:t]
    real   <- y[(t + 1):(t + h)]
    pred   <- tryCatch(as.numeric(f_modelo(ts(train), h)$mean),
                       error = function(e) rep(NA, h))
    real - pred
  }))
  colnames(errores) <- paste0("h", 1:h)
  errores
}

dm_p <- function(e_mod, e_ref, k) {
  ok <- complete.cases(e_mod[, k], e_ref[, k])
  tryCatch(dm.test(e_ref[ok, k], e_mod[ok, k], h = 1, power = 2)$p.value,
           error = function(e) NA_real_)
}

y <- datos$COLCAP
n <- length(y)
h <- 4
gap_dias <- 4
hs <- c(1, 2, 4)

disenos <- list(
  list(muestreo = "rep_no_overlap_gap", ventana = "expansiva", w = NULL),
  list(muestreo = "rep_no_overlap_gap", ventana = "fija",      w = 250),
  list(muestreo = "rep_no_overlap_gap", ventana = "fija",      w = 500),
  list(muestreo = "rep_no_overlap",     ventana = "expansiva", w = NULL)
)

resultados <- list()

for (d in disenos) {
  origenes <- generar_origenes(n, h, d$muestreo, gap = gap_dias)
  errs  <- lapply(modelos, function(f) evaluar_modelo(y, f, h, origenes, ventana = d$w))
  e_ref <- errs$naive
  rmse_ref <- sqrt(colMeans(e_ref[, hs]^2, na.rm = TRUE))
  
  for (m in names(errs)) {
    e <- errs[[m]]
    rmse <- sqrt(colMeans(e[, hs]^2, na.rm = TRUE))
    resultados[[length(resultados) + 1]] <- data.frame(
      muestreo   = d$muestreo,
      ventana    = d$ventana,
      w          = ifelse(is.null(d$w), NA, d$w),
      modelo     = m,
      horizonte  = hs,
      rmse       = unname(rmse),
      rmse_rel   = unname(rmse / rmse_ref),
      p_dm       = sapply(hs, function(k) dm_p(e, e_ref, k)),
      n_origenes = length(origenes)
    )
  }
}

tabla_val <- do.call(rbind, resultados)

tabla_val |>
  filter(horizonte == 4) |>
  arrange(muestreo, ventana, w, rmse)

tabla_val |>
  filter(horizonte == 1) |>
  arrange(muestreo, ventana, w, rmse)

# ---- 5. Ajuste y diagnostico del ARIMA ----

fit_arima <- auto.arima(y_log, d = 1, seasonal = FALSE,
                        stepwise = FALSE, approximation = FALSE)
summary(fit_arima)
checkresiduals(fit_arima, lag = 10)

# Residuales de ARIMA(2,1,0) (Modelo elegido)
fit_210 <- Arima(y_log, order = c(2, 1, 0), include.drift = TRUE)
checkresiduals(fit_210, lag = 10)

# ---- 6. Pronostico de la semana ----

ultima <- max(datos$Fecha)
fechas_pron <- as.Date(c("2026-11-03", "2026-11-04", "2026-11-05", "2026-11-06"))  # lunes 2-nov festivo
fechas_pron

y_final <- if (is.null(W_FINAL)) y else tail(y, W_FINAL)

pron <- sapply(modelos, function(f) as.numeric(f(ts(y_final), h)$mean))
tabla_pron <- data.frame(fecha = fechas_pron, round(pron, 2))
print(tabla_pron)



