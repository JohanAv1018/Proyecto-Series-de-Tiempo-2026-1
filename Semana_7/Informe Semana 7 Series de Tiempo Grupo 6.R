# ============================================================
# COLCAP semana 7: pronostico puntual + intervalos al 95%
# Periodo MSCI COLCAP (desde 2021-05-28)
# Pronostico: martes 17 a viernes 20 de noviembre (lunes 16 festivo)
# ============================================================

library(pacman)
p_load("dplyr", "tidyr", "readxl", "forecast", "lubridate", "ggplot2", "TSA")

set.seed(6)

fecha_cambio <- as.Date("2021-05-28")
corte        <- as.Date("2026-11-13")
alpha        <- 0.05
min_train    <- 250
n_boot       <- 200
bloque_l     <- 6         # longitud de bloque ~ horizonte (antes 10)

# ---- 1. COLCAP ----
datos <- read_xlsx("Semana_7/colcap_semana7.xlsx", sheet = 1, col_types = "text")
datos <- datos[, 1:2]
names(datos) <- c("Fecha", "COLCAP")
datos$Fecha  <- suppressWarnings(dmy(datos$Fecha))
datos        <- datos |> filter(!is.na(Fecha))

limpiar_colcap <- function(x) {
  ifelse(grepl(",", x),
         as.numeric(gsub(",", ".", gsub("\\.", "", x))),
         as.numeric(x))
}

datos$COLCAP <- limpiar_colcap(datos$COLCAP)
datos <- datos |> arrange(Fecha) |> filter(Fecha >= fecha_cambio, Fecha <= corte)
stopifnot(sum(is.na(datos$COLCAP)) == 0,
          max(datos$Fecha) == corte,
          !any(duplicated(datos$Fecha)),
          all(wday(datos$Fecha, week_start = 1) <= 5))


x  <- datos$COLCAP
ly <- log(x)
r  <- diff(ly)
n  <- length(x)

# fechas a pronosticar: martes 17 a viernes 20 (lunes 16 es festivo)
fechas_pron <- seq(as.Date("2026-11-17"), by = "day", length.out = 4)
fechas_pron
h <- length(fechas_pron)

# agrupamiento de volatilidad: p-valor Ljung-Box (lag 10) de r y de r^2
c(r  = Box.test(r,   10, "Ljung-Box")$p.value,
  r2 = Box.test(r^2, 10, "Ljung-Box")$p.value)

# ---- 2. Funciones auxiliares ----

winkler <- function(y, l, u) {
  (u - l) + (2 / alpha) * pmax(l - y, 0) + (2 / alpha) * pmax(y - u, 0)
}

# desviacion estandar diaria de los log-retornos
sd_full <- function(rr) sd(rr)
sd_w    <- function(w) function(rr) sd(tail(rr, w))

# ---- 3. Metodos de intervalo: cada uno devuelve centro, li, ls (h = 1..h) ----
# t = indice de la ultima observacion usada (solo se usa informacion hasta t)

# 3.1 Normal sobre log-retornos: centro con deriva, escala sqrt(h)
int_sigma <- function(t, sigma_fun, z = qnorm(1 - alpha / 2)) {
  rr     <- r[1:(t - 1)]
  s      <- sigma_fun(rr)
  ks     <- 1:h
  centro <- x[t] * exp(ks * mean(rr))
  list(centro = centro,
       li = centro * exp(-z * s * sqrt(ks)),
       ls = centro * exp( z * s * sqrt(ks)))
}

# 3.2 Cuantiles empiricos de los retornos acumulados a k pasos (ultimos 500 dias)
int_emp <- function(t, ventana = 500) {
  ks <- 1:h
  q <- sapply(ks, function(k) {
    rk <- tail(ly[(k + 1):t] - ly[1:(t - k)], ventana)
    quantile(rk, c(alpha / 2, 1 - alpha / 2), names = FALSE)
  })
  list(centro = x[t] * exp(ks * mean(r[1:(t - 1)])),
       li = x[t] * exp(q[1, ]), ls = x[t] * exp(q[2, ]))
}

# 3.3 ARIMA sobre log(COLCAP), d = 1
#     fijo: orden y coeficientes de los primeros min_train datos (como M3_4b_uso.R)
#     recal: auto.arima reestimado en cada origen
fit_ini <- auto.arima(ly[1:min_train], d = 1, seasonal = FALSE,
                      stepwise = FALSE, approximation = FALSE)
fit_ini

# ---- 3.3b Diagnostico de residuales del ARIMA (estilo M3_3): ARCH / McLeod-Li ----
# Justifica usar intervalos NO normales (bootstrap/empirico) cuando hay
# agrupamiento de volatilidad en los residuales.
fit_diag <- auto.arima(ly, d = 1, seasonal = FALSE,
                       stepwise = FALSE, approximation = FALSE)
res_diag <- as.numeric(residuals(fit_diag))

# Ljung-Box sobre residuos (correlacion lineal) y sobre residuos^2 (heterocedasticidad)
cat("Ljung-Box residuos    (lag 10): p =",
    round(Box.test(res_diag,   10, "Ljung-Box")$p.value, 4), "\n")
cat("Ljung-Box residuos^2  (lag 10): p =",
    round(Box.test(res_diag^2, 10, "Ljung-Box")$p.value, 4), "\n")

# Test de McLeod-Li (ARCH): p-valores bajos (< 0.05) => efecto ARCH / volatilidad agrupada
McLeod.Li.test(y = res_diag)

ajustar_arima <- function(t, recal) {
  yt <- ly[1:t]
  if (recal) auto.arima(yt, d = 1, seasonal = FALSE) else Arima(yt, model = fit_ini)
}

int_arima <- function(t, recal = FALSE) {
  fc <- forecast(ajustar_arima(t, recal), h = h, level = 100 * (1 - alpha))
  list(centro = exp(as.numeric(fc$mean)),
       li = exp(as.numeric(fc$lower[, 1])),
       ls = exp(as.numeric(fc$upper[, 1])))
} 

# 3.4 ARIMA con bootstrap de residuos (estilo M3_4a_uso.R).
#     Trayectoria futura = media + sum_j psi_{k-j} e_j, con e_j remuestreados:
#     iid = bootstrap de innovaciones; mbb = moving block; cbb = circular;
#     sieve = AR ajustado a los residuos y remuestreo de sus innovaciones
conv <- function(a, b) {
  out <- rep(0, length(a) + length(b) - 1)
  for (i in seq_along(a)) out[i:(i + length(b) - 1)] <- out[i:(i + length(b) - 1)] + a[i] * b
  out
}

psi_arima <- function(fit) {                      # pesos psi del proceso integrado (d = 1)
  a_full <- conv(c(1, -fit$model$phi), c(1, -1))
  c(1, ARMAtoMA(ar = -a_full[-1], ma = fit$model$theta, lag.max = h - 1))
}

innov_boot <- function(res, B, tipo) {
  m <- length(res); l <- bloque_l
  if (tipo == "sieve") {
    a <- ar(res, aic = TRUE, order.max = 10, method = "yule-walker")
    e <- as.numeric(na.omit(a$resid)); phi <- as.numeric(a$ar)
  }
  sapply(seq_len(B), function(b) {
    switch(tipo,
           iid   = sample(res, h, replace = TRUE),
           mbb   = { st <- sample(1:(m - l + 1), ceiling(h / l), replace = TRUE)
           unlist(lapply(st, function(s) res[s:(s + l - 1)]))[1:h] },
           cbb   = { st <- sample(1:m, ceiling(h / l), replace = TRUE)
           unlist(lapply(st, function(s) res[((s:(s + l - 1)) - 1) %% m + 1]))[1:h] },
           sieve = as.numeric(arima.sim(model = if (length(phi)) list(ar = phi) else list(),
                                        n = h, innov = sample(e, h, TRUE),
                                        n.start = 50, start.innov = sample(e, 50, TRUE)))
    )
  })
}

int_arima_boot <- function(t, tipo, recal = FALSE) {
  fit  <- ajustar_arima(t, recal)
  res  <- as.numeric(residuals(fit))[-(1:10)]     # se omiten residuos de inicializacion
  mu   <- as.numeric(forecast(fit, h = h)$mean)
  psi  <- psi_arima(fit)
  Psi  <- outer(1:h, 1:h, function(k, j) ifelse(k >= j, psi[pmax(k - j + 1, 1)], 0))
  sims <- exp(mu + Psi %*% innov_boot(res, n_boot, tipo))   # h x n_boot, en puntos
  list(centro = exp(mu),
       li = apply(sims, 1, quantile, alpha / 2, names = FALSE),
       ls = apply(sims, 1, quantile, 1 - alpha / 2, names = FALSE))
}

# 3.5 ETS sobre log(COLCAP) (sin estacionalidad), normal y bootstrap
int_ets <- function(t) {
  fc <- forecast(ets(ly[1:t], model = "ZZN"), h = h, level = 100 * (1 - alpha))
  list(centro = exp(as.numeric(fc$mean)),
       li = exp(as.numeric(fc$lower[, 1])),
       ls = exp(as.numeric(fc$upper[, 1])))
}

int_ets_boot <- function(t) {
  fit  <- ets(ly[1:t], model = "ZZN")
  sims <- exp(sapply(seq_len(n_boot), function(b)
    as.numeric(simulate(fit, nsim = h, future = TRUE, bootstrap = TRUE))))
  list(centro = exp(as.numeric(forecast(fit, h = h)$mean)),
       li = apply(sims, 1, quantile, alpha / 2, names = FALSE),
       ls = apply(sims, 1, quantile, 1 - alpha / 2, names = FALSE))
}

modelos <- list(
  normal_full      = function(t) int_sigma(t, sd_full),
  rw_drift_IC_normal60       = function(t) int_sigma(t, sd_w(60)),
  empirico500      = function(t) int_emp(t),
  arima_norm       = function(t) int_arima(t),
  arima_recal_norm = function(t) int_arima(t, recal = TRUE),
  arima_boot_iid   = function(t) int_arima_boot(t, "iid"),
  arima_boot_mbb   = function(t) int_arima_boot(t, "mbb"),
  arima_boot_cbb   = function(t) int_arima_boot(t, "cbb"),
  arima_boot_sieve = function(t) int_arima_boot(t, "sieve"),
  ets_norm         = function(t) int_ets(t),
  ets_boot         = function(t) int_ets_boot(t)
)

# ---- 4. Validacion: bloques de 4 dias sin solapamiento, ventana expansiva ----
# (puede tardar unos minutos por el bootstrap)

origenes <- rev(seq(n - h, min_train, by = -h))
length(origenes)

ev <- do.call(rbind, lapply(names(modelos), function(m) {
  do.call(rbind, lapply(origenes, function(t) {
    iv <- tryCatch(modelos[[m]](t), error = function(e) NULL)
    if (is.null(iv)) return(NULL)
    data.frame(modelo = m, origen = t, paso = 1:h, obs = x[t + (1:h)],
               li = iv$li, ls = iv$ls)
  }))
}))
ev$w   <- winkler(ev$obs, ev$li, ev$ls)
ev$cub <- ev$obs >= ev$li & ev$obs <= ev$ls
stopifnot(all(table(ev$modelo) == length(origenes) * h))   # ningun metodo fallo

tabla_val <- ev |>
  group_by(modelo, horizonte = paso) |>
  summarise(winkler = round(mean(w), 2), cobertura = round(mean(cub), 3),
            ancho = round(mean(ls - li), 2), .groups = "drop")

ranking <- ev |>
  group_by(modelo) |>
  summarise(winkler = round(mean(w), 2), cobertura = round(mean(cub), 3),
            ancho = round(mean(ls - li), 2), .groups = "drop") |>
  arrange(winkler)

# comparacion pareada contra el mejor metodo (orientativa: hay dependencia entre origenes)
mejor <- ranking$modelo[1]
w_or <- ev |> group_by(modelo, origen) |> summarise(w = mean(w), .groups = "drop") |>
  pivot_wider(names_from = modelo, values_from = w)
ranking$p_vs_mejor <- round(sapply(ranking$modelo, function(m)
  if (m == mejor) NA_real_ else t.test(w_or[[m]], w_or[[mejor]], paired = TRUE)$p.value), 3)

as.data.frame(ranking)

# ---- 5. Pronostico de la semana ----

ganador <- ranking$modelo[1]      # menor Winkler
ganador

set.seed(6)                       # reproducibilidad del pronostico final
n_boot <- 2000                    # mas trayectorias para el pronostico final
iv <- modelos[[ganador]](n)

tabla_pron <- data.frame(fecha = fechas_pron, paso = 1:h, modelo = ganador,
                         puntual = round(iv$centro, 2),
                         LI95 = round(iv$li, 2), LS95 = round(iv$ls, 2))
tabla_pron

dir.create("Semana_7/resultados", showWarnings = FALSE, recursive = TRUE)
saveRDS(tabla_pron, "Semana_7/resultados/tabla_pron_semana7.rds")
