# ============================================================
# COLCAP semana 6: variables exogenas rezagadas + estrategia directa
# Periodo MSCI COLCAP (desde 2021-05-28)
# ============================================================

library(pacman)
p_load(dplyr, readxl, forecast, lubridate, quantmod, zoo, lmtest, ranger)

FECHA_CAMBIO <- as.Date("2021-05-28")
H <- 5

# ---- 1. COLCAP ----
# NOTA: ajuste el nombre del archivo al de su descarga (p.ej. "COLCAP_1.xlsx")
datos <- read_xlsx("COLCAP_1.xlsx", sheet = 1, col_types = "text")
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

CORTE <- max(datos$Fecha)
CORTE
tail(datos, 3)

# ---- 2. Exogenas (Yahoo Finance) ----

simbolos <- c(sp500 = "^GSPC", eem = "EEM", brent = "BZ=F", usdcop = "COP=X", vix = "^VIX")

bajar <- function(sym) {
  x <- getSymbols(sym, src = "yahoo", from = FECHA_CAMBIO - 30, to = CORTE + 1,
                  auto.assign = FALSE)
  tibble(Fecha = as.Date(index(x)), valor = as.numeric(Ad(x))) |>
    filter(!is.na(valor), Fecha <= CORTE)
}

# a cada fecha del COLCAP: ultimo valor disponible en o antes de esa fecha, nunca posterior
alinear <- function(fechas, df) {
  todo <- full_join(tibble(Fecha = fechas), df, by = "Fecha") |> arrange(Fecha)
  todo$valor <- na.locf(todo$valor, na.rm = FALSE)
  todo$valor[match(fechas, todo$Fecha)]
}

# Copia point-in-time: si ya existe, se reutiliza (reproducible sin depender de Yahoo)
if (file.exists("exogenas_semana6.rds")) {
  base <- readRDS("exogenas_semana6.rds")
} else {
  base <- datos
  for (nm in names(simbolos)) base[[nm]] <- alinear(base$Fecha, bajar(simbolos[[nm]]))
  saveRDS(base, "exogenas_semana6.rds")
}

summary(base)
tail(base, 3)

# ---- 3. Rasgos del dia t y objetivos t -> t+k ----

lr <- function(v) c(NA, diff(log(v)))

ly <- log(base$COLCAP)
n  <- nrow(base)

feat <- data.frame(
  colcap1 = lr(base$COLCAP),
  colcap5 = c(rep(NA, 5), diff(ly, lag = 5)),
  sp500   = lr(base$sp500),
  eem     = lr(base$eem),
  brent   = lr(base$brent),
  usdcop  = lr(base$usdcop),
  vix     = lr(base$vix)
)
feat$colcap2 <- c(NA, head(feat$colcap1, -1))
X <- as.matrix(feat)

Y <- sapply(1:H, function(k) c(ly[(1 + k):n] - ly[1:(n - k)], rep(NA, k)))

# ---- 4. Seleccion de exogenas ----

ok  <- complete.cases(X)
exo <- c("sp500", "eem", "brent", "usdcop", "vix")

granger_p <- sapply(exo, function(v) {
  grangertest(x = X[ok, v], y = X[ok, "colcap1"], order = 2)$`Pr(>F)`[2]
})
round(granger_p, 4)

cor_pred <- sapply(colnames(X), function(v) {
  c(h1 = cor(X[, v], Y[, 1], use = "complete.obs"),
    h5 = cor(X[, v], Y[, 5], use = "complete.obs"))
})
round(cor_pred, 3)
2 / sqrt(sum(ok))

# ---- 5. Modelos: cada uno devuelve el log-retorno acumulado de t a t+1..t+H ----

filas_train <- function(t, k) {
  idx <- seq_len(t - k)
  idx[complete.cases(X[idx, , drop = FALSE]) & !is.na(Y[idx, k])]
}

m_naive <- function(t) rep(0, H)

m_rw_drift <- function(t) (1:H) * mean(diff(ly[1:t]))

m_arima <- function(t) {
  fit <- auto.arima(ly[1:t], d = 1, seasonal = FALSE)
  as.numeric(forecast(fit, h = H)$mean) - ly[t]
}

hacer_lm <- function(vars) function(t) {
  sapply(1:H, function(k) {
    idx <- filas_train(t, k)
    b <- lm.fit(cbind(1, X[idx, vars, drop = FALSE]), Y[idx, k])$coefficients
    b[is.na(b)] <- 0
    sum(c(1, X[t, vars]) * b)
  })
}

m_lm_colcap <- hacer_lm(c("colcap1", "colcap2"))
m_lm_exo    <- hacer_lm(c("colcap1", "colcap2", exo))

m_rf <- function(t) {
  sapply(1:H, function(k) {
    idx <- filas_train(t, k)
    d   <- data.frame(y = Y[idx, k], X[idx, , drop = FALSE])
    fit <- ranger(y ~ ., data = d, num.trees = 300, min.node.size = 20, seed = 123)
    predict(fit, data.frame(X[t, , drop = FALSE]))$predictions
  })
}

modelos <- list(naive = m_naive, rw_drift = m_rw_drift, arima = m_arima,
                lm_colcap = m_lm_colcap, lm_exo = m_lm_exo, rf = m_rf)

# ---- 6. Validacion: bloques de 5 dias sin solapamiento (puede tardar unos minutos) ----

origenes <- rev(seq(n - H, 250, by = -H))
length(origenes)

pred <- lapply(modelos, function(f) {
  t(sapply(origenes, function(t) tryCatch(f(t), error = function(e) rep(NA, H))))
})
pred$combo <- (pred$arima + pred$lm_exo + pred$rf) / 3

real <- t(sapply(origenes, function(t) base$COLCAP[t + (1:H)]))
ult  <- base$COLCAP[origenes]
err  <- lapply(pred, function(p) real - ult * exp(p))

dm_p <- function(e, e0, k) {
  ok <- complete.cases(e[, k], e0[, k])
  tryCatch(dm.test(e0[ok, k], e[ok, k], h = 1, power = 2)$p.value,
           error = function(x) NA_real_)
}

rmse0 <- sqrt(colMeans(err$naive^2, na.rm = TRUE))

tabla_val <- do.call(rbind, lapply(names(err), function(m) {
  e    <- err[[m]]
  rmse <- sqrt(colMeans(e^2, na.rm = TRUE))
  data.frame(modelo = m, horizonte = 1:H, rmse = round(rmse, 2),
             rmse_rel = round(rmse / rmse0, 4),
             p_dm = round(sapply(1:H, function(k) dm_p(e, err$naive, k)), 3))
}))

tabla_val |> filter(horizonte %in% c(1, 5)) |> arrange(horizonte, rmse_rel)

# coeficientes de la regresion directa a 1 dia con toda la muestra
idx1 <- filas_train(n, 1)
summary(lm(Y[idx1, 1] ~ X[idx1, c("colcap1", "colcap2", exo)]))

# ---- 7. Pronostico de la semana ----

cand <- seq(CORTE + 1, by = "day", length.out = 12)
fechas_pron <- cand[wday(cand, week_start = 1) <= 5][1:H]
fechas_pron

p_fin <- sapply(modelos, function(f) f(n))
p_fin <- cbind(p_fin, combo = (p_fin[, "arima"] + p_fin[, "lm_exo"] + p_fin[, "rf"]) / 3)

tabla_pron <- data.frame(fecha = fechas_pron, round(base$COLCAP[n] * exp(p_fin), 2))
tabla_pron

# Modelo entregado esta semana: paseo aleatorio con deriva (rw_drift)
dir.create("resultados", showWarnings = FALSE)
saveRDS(tabla_val,  "resultados/tabla_val_semana6.rds")
saveRDS(tabla_pron, "resultados/tabla_pron_semana6.rds")
