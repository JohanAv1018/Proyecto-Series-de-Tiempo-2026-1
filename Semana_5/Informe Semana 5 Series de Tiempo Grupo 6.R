# ============================================================
# COLCAP (MSCI, desde 2021-05-28): suavizamiento/ETS + ARMA/ARIMA/SARIMA
# ============================================================

library(pacman)
p_load(dplyr, readxl, forecast, lubridate, ggplot2, tseries)

FECHA_CAMBIO <- as.Date("2021-05-28")
W_FINAL      <- NULL   

# ---- 1. Datos ----

datos <- read_xlsx("COLCAP_datos_usados_Semana4_corte_2026-10-23.xlsx", sheet = 1, col_types = "text")
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
f_snaive <- function(x, h) snaive(ts(as.numeric(x), frequency = 5), h = h)
f_drift  <- function(x, h) rwf(as.numeric(x), h = h, drift = TRUE)
f_ses    <- function(x, h) ses(as.numeric(x), h = h)
f_holt   <- function(x, h) holt(as.numeric(x), h = h)
f_damp   <- function(x, h) holt(as.numeric(x), h = h, damped = TRUE)
f_ets    <- function(x, h) forecast(ets(as.numeric(x)), h = h)

f_arma <- function(x, h) {
  x  <- as.numeric(x)
  r  <- diff(log(x))
  fit <- auto.arima(r, d = 0, seasonal = FALSE, stationary = TRUE)
  pr <- as.numeric(forecast(fit, h = h)$mean)
  list(mean = x[length(x)] * exp(cumsum(pr)))
}

f_arima <- function(x, h) {
  fit <- auto.arima(as.numeric(x), d = 1, seasonal = FALSE, lambda = 0)
  forecast(fit, h = h)
}

f_arima210 <- function(x, h) {
  fit <- Arima(as.numeric(x), order = c(2, 1, 0), include.drift = TRUE, lambda = 0)
  forecast(fit, h = h)
}

f_sarima <- function(x, h) {
  fit <- auto.arima(ts(as.numeric(x), frequency = 5), d = 1, seasonal = TRUE, lambda = 0)
  forecast(fit, h = h)
}

# ---- 3.1 Combinaciones (modelos en paralelo, se promedian pronosticos) ----

pron_vec <- function(f, x, h) as.numeric(f(x, h)$mean)

f_combo_media <- function(x, h) {
  P <- cbind(pron_vec(f_ets, x, h), pron_vec(f_arima, x, h), pron_vec(f_naive, x, h))
  list(mean = rowMeans(P))
}

f_combo_mediana <- function(x, h) {
  P <- cbind(pron_vec(f_ets, x, h), pron_vec(f_arima, x, h), pron_vec(f_naive, x, h),
             pron_vec(f_drift, x, h), pron_vec(f_damp, x, h))
  list(mean = apply(P, 1, median))
}

# ---- 3.2 Hibridos por etapas (cada modelo trabaja sobre lo que dejo el anterior) ----

# etapa 1: ETS sobre log del nivel; etapa 2: ARMA sobre sus residuos
f_hib_ets_arma <- function(x, h) {
  lx  <- log(as.numeric(x))
  fit1 <- ets(lx, additive.only = TRUE)
  p1   <- as.numeric(forecast(fit1, h = h)$mean)
  fit2 <- auto.arima(as.numeric(residuals(fit1)), d = 0, seasonal = FALSE,
                     stationary = TRUE)
  p2   <- as.numeric(forecast(fit2, h = h)$mean)
  list(mean = exp(p1 + p2))
}

# etapa 1: STL con ciclo semanal; etapa 2: ETS sobre la serie desestacionalizada
f_hib_stl_ets <- function(x, h) {
  stlf(ts(as.numeric(x), frequency = 5), h = h, method = "ets", lambda = 0)
}

modelos <- list(
  naive = f_naive, snaive = f_snaive, drift = f_drift, ses = f_ses, holt = f_holt,
  damped = f_damp, ets = f_ets,
  arma = f_arma, arima = f_arima, arima210 = f_arima210, sarima = f_sarima,
  combo_media = f_combo_media, combo_mediana = f_combo_mediana,
  hib_ets_arma = f_hib_ets_arma, hib_stl_ets = f_hib_stl_ets
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
h <- 5
gap_dias <- 5
hs <- c(1, 3, 5)

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
  filter(horizonte == 5) |>
  arrange(muestreo, ventana, w, rmse)

tabla_val |>
  filter(horizonte == 1) |>
  arrange(muestreo, ventana, w, rmse)

# ---- 5. Ajuste y diagnostico del ARIMA ----

fit_arima <- auto.arima(y_log, d = 1, seasonal = FALSE,
                        stepwise = FALSE, approximation = FALSE)
summary(fit_arima)
checkresiduals(fit_arima, lag = 10)

fit_arma_ret <- auto.arima(ret, d = 0, seasonal = FALSE, stationary = TRUE,
                           stepwise = FALSE, approximation = FALSE)
summary(fit_arma_ret)

# ---- 6. Pronostico de la semana ----

ultima <- max(datos$Fecha)
cand   <- seq(ultima + 1, by = "day", length.out = 12)
fechas_pron <- cand[wday(cand, week_start = 1) <= 5][1:5]
fechas_pron

y_final <- if (is.null(W_FINAL)) y else tail(y, W_FINAL)

pron <- sapply(modelos, function(f) as.numeric(f(ts(y_final), 5)$mean))
tabla_pron <- data.frame(fecha = fechas_pron, round(pron, 2))
print(tabla_pron)

dir.create("resultados", showWarnings = FALSE)
saveRDS(tabla_val,  "resultados/tabla_val_arima.rds")
saveRDS(tabla_pron, "resultados/tabla_pron_arima.rds")

# ---- 7. Paseo aleatorio con drift y transformaciones ----

lam <- BoxCox.lambda(y, method = "guerrero")
lam

f_rw_drift <- function(x, h) {
  fit <- Arima(as.numeric(x), order = c(0, 1, 0), include.drift = TRUE, lambda = 0)
  forecast(fit, h = h)
}

f_rw_drift_bc <- function(x, h) {
  fit <- Arima(as.numeric(x), order = c(0, 1, 0), include.drift = TRUE, lambda = lam)
  forecast(fit, h = h)
}

f_rw_drift_adj <- function(x, h) {
  fit <- Arima(as.numeric(x), order = c(0, 1, 0), include.drift = TRUE, lambda = 0)
  forecast(fit, h = h, biasadj = TRUE)
}

modelos_cmp <- list(naive = f_naive, arima = f_arima, drift = f_drift,
                    rw_drift = f_rw_drift, rw_drift_bc = f_rw_drift_bc,
                    rw_drift_adj = f_rw_drift_adj)

cmp <- list()
for (d in disenos) {
  origenes <- generar_origenes(n, h, d$muestreo, gap = gap_dias)
  errs <- lapply(modelos_cmp, function(f) evaluar_modelo(y, f, h, origenes, ventana = d$w))
  rmse_ref <- sqrt(colMeans(errs$naive[, hs]^2, na.rm = TRUE))
  for (m in names(errs)) {
    rmse <- sqrt(colMeans(errs[[m]][, hs]^2, na.rm = TRUE))
    cmp[[length(cmp) + 1]] <- data.frame(
      muestreo = d$muestreo, w = ifelse(is.null(d$w), NA, d$w), modelo = m,
      horizonte = hs, rmse_rel = round(rmse / rmse_ref, 4),
      p_dm = round(sapply(hs, function(k) dm_p(errs[[m]], errs$naive, k)), 3)
    )
  }
}
tabla_cmp <- do.call(rbind, cmp)

tabla_cmp |> filter(horizonte == 1) |> arrange(muestreo, w, rmse_rel)
tabla_cmp |> filter(horizonte == 5) |> arrange(muestreo, w, rmse_rel)

sapply(modelos_cmp, function(f) round(as.numeric(f(ts(y), 5)$mean), 2))

cambio5 <- diff(log(y), lag = 5)
mean(cambio5 > 0)