# ============================================================
# COLCAP semana 8: pronostico puntual + intervalos al 95%
# (EWMA / GARCH) en lugar de varianza constante
# ============================================================

library(pacman)
p_load("dplyr", "tidyr", "readxl", "lubridate", "forecast", "fGarch", "ggplot2")

set.seed(6)

fecha_cambio <- as.Date("2021-05-28")
corte        <- as.Date("2026-11-20")
alpha        <- 0.05
z            <- qnorm(1 - alpha / 2)
min_train    <- 250
refit_cada   <- 1         # reestimacion de los GARCH en cada origen (con 20 el resultado casi no cambia y corre ~15 veces mas rapido)
lambda_ewma  <- 0.94      # RiskMetrics

# ---- 1. COLCAP ----
datos <- read_xlsx("colcap_semana8.xlsx", sheet = 1, col_types = "text")
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

x <- datos$COLCAP
r <- 100 * diff(log(x))   # log-retornos en %
n <- length(x)
n

fechas_pron <- seq(as.Date("2026-11-23"), by = "day", length.out = 5)
h  <- length(fechas_pron)
ks <- 1:h

# ---- 2. Diagnostico: agrupamiento de volatilidad y colas pesadas ----
# Ljung-Box sobre r^2 = estadistico de McLeod-Li
data.frame(rezago   = c(5, 10, 20),
           p_r      = sapply(c(5, 10, 20), function(l) Box.test(r,   l, "Ljung-Box")$p.value),
           p_r2     = sapply(c(5, 10, 20), function(l) Box.test(r^2, l, "Ljung-Box")$p.value))

c(curtosis_exceso = mean((r - mean(r))^4) / var(r)^2 - 3,
  sd_20d = sd(tail(r, 20)), sd_60d = sd(tail(r, 60)), sd_total = sd(r))

# ---- 3. Modelos de varianza ----
# Todos usan el mismo centro: paseo aleatorio con deriva, c_k = k * media(r).
# Cambia solo la varianza del retorno acumulado a k dias, V_k = sum_{j<=k} sigma^2_{t+j}.
# Intervalo: x_t * exp((c_k -/+ q * sqrt(V_k)) / 100)
# q = 1.96 para errores normales; para errores t, cuantil de la t estandarizada a varianza 1
# (sigma_t del GARCH es la desviacion, no la escala de la t)

# 3.1 GARCH(1,1) y GJR-GARCH(1,1) con fGarch.
#     fGarch parametriza el GJR como APARCH con delta = 2:
#     sigma2_t = omega + alpha (|e| - gamma e)^2 + beta sigma2_{t-1}
ajustar_garch <- function(rr, tipo) {
  f <- switch(tipo,
              garch_n = garchFit(~ garch(1, 1), data = rr, cond.dist = "norm", trace = FALSE),
              garch_t = garchFit(~ garch(1, 1), data = rr, cond.dist = "std",  trace = FALSE),
              gjr_t   = garchFit(~ aparch(1, 1), data = rr, cond.dist = "std",
                                 include.delta = FALSE, delta = 2, trace = FALSE))
  cf <- coef(f)
  list(mu = cf[["mu"]], omega = cf[["omega"]], alpha = cf[["alpha1"]],
       beta = cf[["beta1"]], gamma = if ("gamma1" %in% names(cf)) cf[["gamma1"]] else 0,
       nu = if ("shape" %in% names(cf)) cf[["shape"]] else Inf)
}

# Filtro de la varianza con coeficientes dados y pronostico a h pasos.
# E[(|e| - gamma e)^2] = sigma^2 (1 + gamma^2), de ahi la persistencia del GJR.
var_garch <- function(rr, p) {
  e  <- rr - p$mu
  s2 <- var(rr)
  for (i in seq_along(e)) s2 <- p$omega + p$alpha * (abs(e[i]) - p$gamma * e[i])^2 + p$beta * s2
  pers <- p$alpha * (1 + p$gamma^2) + p$beta
  v <- numeric(h); v[1] <- s2
  for (k in 2:h) v[k] <- p$omega + pers * v[k - 1]
  v
}

# 3.2 EWMA (RiskMetrics): sigma2_{t+1} = lambda sigma2_t + (1 - lambda) r_t^2, plana en k
var_ewma <- function(rr) {
  w <- lambda_ewma^(rev(seq_along(rr)) - 1)
  rep(sum(w * rr^2) / sum(w), h)
}

# 3.3 Varianzas constantes (semana 7)
var_sd60  <- function(rr) rep(var(tail(rr, 60)), h)

q_std <- function(nu) if (is.finite(nu)) qt(1 - alpha / 2, nu) * sqrt((nu - 2) / nu) else z
var_total <- function(rr) rep(var(rr), h)

# ---- 4. Validacion: origen movil diario, ventana expansiva ----
# Origen t = indice del ultimo cierre usado; se pronostican x[t+1], ..., x[t+5].
# Los GARCH se reestiman cada 'refit_cada' origenes; entre reestimaciones
# los coeficientes quedan fijos y solo se actualiza el filtro con los datos nuevos.

origenes <- min_train:(n - h)
length(origenes)

modelos_var <- c("sd_total", "sd60", "ewma", "garch_n", "garch_t", "gjr_t")
coefs <- list()
V <- lapply(modelos_var, function(m) matrix(NA, length(origenes), h))
names(V) <- modelos_var
Q <- lapply(modelos_var, function(m) rep(z, length(origenes)))
names(Q) <- modelos_var
centro <- matrix(NA, length(origenes), h)

for (i in seq_along(origenes)) {
  t  <- origenes[i]
  rr <- r[1:(t - 1)]
  if ((i - 1) %% refit_cada == 0)
    for (m in c("garch_n", "garch_t", "gjr_t")) coefs[[m]] <- ajustar_garch(rr, m)
  centro[i, ] <- ks * mean(rr)
  V$sd_total[i, ] <- cumsum(var_total(rr))
  V$sd60[i, ]     <- cumsum(var_sd60(rr))
  V$ewma[i, ]     <- cumsum(var_ewma(rr))
  for (m in c("garch_n", "garch_t", "gjr_t")) {
    V[[m]][i, ] <- cumsum(var_garch(rr, coefs[[m]]))
    Q[[m]][i]   <- q_std(coefs[[m]]$nu)
  }
}

winkler <- function(y, l, u) (u - l) + (2 / alpha) * pmax(l - y, 0) + (2 / alpha) * pmax(y - u, 0)

obs <- t(sapply(origenes, function(t) x[t + ks]))
xt  <- x[origenes]

bandas <- lapply(modelos_var, function(m) list(li = xt * exp((centro - Q[[m]] * sqrt(V[[m]])) / 100),
                                               ls = xt * exp((centro + Q[[m]] * sqrt(V[[m]])) / 100)))
names(bandas) <- modelos_var

# Combinacion: promedio de los limites de GARCH-t y EWMA
bandas$comb_garcht_ewma <- list(li = (bandas$garch_t$li + bandas$ewma$li) / 2,
                                ls = (bandas$garch_t$ls + bandas$ewma$ls) / 2)

ev <- do.call(rbind, lapply(names(bandas), function(m) {
  b <- bandas[[m]]
  data.frame(modelo = m, origen = rep(origenes, h), paso = rep(ks, each = length(origenes)),
             obs = as.vector(obs), li = as.vector(b$li), ls = as.vector(b$ls))
}))
ev$w   <- winkler(ev$obs, ev$li, ev$ls)
ev$cub <- ev$obs >= ev$li & ev$obs <= ev$ls

# 4.1 Todos los origenes (con solape)
ranking <- ev |>
  group_by(modelo) |>
  summarise(winkler = mean(w), cobertura = mean(cub), ancho = mean(ls - li), .groups = "drop") |>
  arrange(winkler)

# 4.2 Bloques de 5 dias sin solape (comparable con semanas anteriores)
orig_bloque <- rev(seq(n - h, min_train, by = -h))
ranking_bloques <- ev |>
  filter(origen %in% orig_bloque) |>
  group_by(modelo) |>
  summarise(winkler_bloques = mean(w), cobertura_bloques = mean(cub), .groups = "drop")

ranking <- ranking |> left_join(ranking_bloques, by = "modelo")

# 4.3 Diebold-Mariano sobre la perdida de Winkler (HAC con h-1 rezagos por el solape)
#     dm.test con power = 1 sobre perdidas no negativas compara directamente el Winkler
w_origen <- ev |> group_by(modelo, origen) |> summarise(w = mean(w), .groups = "drop") |>
  pivot_wider(names_from = modelo, values_from = w)
mejor <- ranking$modelo[1]
ranking$p_DM_vs_mejor <- sapply(ranking$modelo, function(m)
  if (m == mejor) NA_real_ else
    dm.test(w_origen[[m]], w_origen[[mejor]], h = h, power = 1, alternative = "greater")$p.value)

ranking |> mutate(across(where(is.numeric), ~ round(.x, 3))) |> as.data.frame()

tabla_horizonte <- ev |>
  group_by(modelo, paso) |>
  summarise(winkler = round(mean(w), 1), cobertura = round(mean(cub), 3), .groups = "drop") |>
  pivot_wider(names_from = paso, values_from = c(winkler, cobertura))
as.data.frame(tabla_horizonte)

# 4.4 Calibracion condicional: cobertura en dias de volatilidad alta vs baja
#     (un intervalo de varianza constante sub-cubre justo cuando la volatilidad sube)
vol_ref <- data.frame(origen = origenes, s20 = sapply(origenes, function(t) sd(r[(t - 20):(t - 1)])))
vol_ref$regimen <- ifelse(vol_ref$s20 > median(vol_ref$s20), "alta", "baja")
cobertura_regimen <- ev |>
  left_join(vol_ref, by = "origen") |>
  group_by(modelo, regimen) |>
  summarise(cobertura = round(mean(cub), 3), winkler = round(mean(w), 1), .groups = "drop") |>
  pivot_wider(names_from = regimen, values_from = c(cobertura, winkler))
as.data.frame(cobertura_regimen)

# 4.5 Contrafactual de la semana 7 (origen 13-nov, 4 dias): Winkler realizado
t7    <- which(datos$Fecha == as.Date("2026-11-13"))
rr7   <- r[1:(t7 - 1)]
obs7  <- x[t7 + 1:4]
c7    <- (1:4) * mean(rr7)
coef7 <- lapply(c(garch_n = "garch_n", garch_t = "garch_t", gjr_t = "gjr_t"), function(m) ajustar_garch(rr7, m))
V7 <- list(sd_total = cumsum(var_total(rr7)), sd60 = cumsum(var_sd60(rr7)), ewma = cumsum(var_ewma(rr7)),
           garch_n = cumsum(var_garch(rr7, coef7$garch_n)), garch_t = cumsum(var_garch(rr7, coef7$garch_t)),
           gjr_t = cumsum(var_garch(rr7, coef7$gjr_t)))
Q7 <- c(sd_total = z, sd60 = z, ewma = z, sapply(coef7, function(p) q_std(p$nu)))
b7 <- lapply(names(V7), function(m) list(li = x[t7] * exp((c7 - Q7[[m]] * sqrt(V7[[m]][1:4])) / 100),
                                         ls = x[t7] * exp((c7 + Q7[[m]] * sqrt(V7[[m]][1:4])) / 100)))
names(b7) <- names(V7)
b7$comb_garcht_ewma <- list(li = (b7$garch_t$li + b7$ewma$li) / 2, ls = (b7$garch_t$ls + b7$ewma$ls) / 2)
contraf <- data.frame(modelo = names(b7),
                      winkler_sem7 = sapply(b7, function(b) round(mean(winkler(obs7, b$li, b$ls)), 1)),
                      cubiertos = sapply(b7, function(b) sum(obs7 >= b$li & obs7 <= b$ls))) |>
  arrange(winkler_sem7)
contraf

# ---- 5. Pronostico de la semana ----
ganador <- ranking$modelo[1]
ganador

rr_fin <- r
c_fin  <- ks * mean(rr_fin)
coef_fin <- lapply(c(garch_n = "garch_n", garch_t = "garch_t", gjr_t = "gjr_t"),
                   function(m) ajustar_garch(rr_fin, m))
coef_fin$garch_t

V_fin <- list(sd_total = cumsum(var_total(rr_fin)),
              sd60     = cumsum(var_sd60(rr_fin)),
              ewma     = cumsum(var_ewma(rr_fin)),
              garch_n  = cumsum(var_garch(rr_fin, coef_fin$garch_n)),
              garch_t  = cumsum(var_garch(rr_fin, coef_fin$garch_t)),
              gjr_t    = cumsum(var_garch(rr_fin, coef_fin$gjr_t)))

Q_fin <- c(sd_total = z, sd60 = z, ewma = z, sapply(coef_fin, function(p) q_std(p$nu)))
Q_fin
b_fin <- lapply(names(V_fin), function(m) list(li = x[n] * exp((c_fin - Q_fin[[m]] * sqrt(V_fin[[m]])) / 100),
                                               ls = x[n] * exp((c_fin + Q_fin[[m]] * sqrt(V_fin[[m]])) / 100)))
names(b_fin) <- names(V_fin)
b_fin$comb_garcht_ewma <- list(li = (b_fin$garch_t$li + b_fin$ewma$li) / 2,
                               ls = (b_fin$garch_t$ls + b_fin$ewma$ls) / 2)

# volatilidad diaria pronosticada (%) por modelo
sapply(V_fin, function(v) round(sqrt(diff(c(0, v))), 3))

tabla_pron <- data.frame(fecha = fechas_pron, paso = ks, modelo = ganador,
                         puntual = round(x[n] * exp(c_fin / 100), 2),
                         LI95 = round(b_fin[[ganador]]$li, 2),
                         LS95 = round(b_fin[[ganador]]$ls, 2))
tabla_pron

# El puntual exp(c_k) es la mediana bajo log-normalidad; la media agrega V_k/2.
# Diferencia en puntos (referencia; en la validacion la mediana tuvo RMSE igual o menor):
V_comb <- (V_fin$garch_t + V_fin$ewma) / 2
round(x[n] * exp(c_fin / 100 + V_comb / (2 * 100^2)) - tabla_pron$puntual, 2)

dir.create("Semana_8/resultados", showWarnings = FALSE, recursive = TRUE)
saveRDS(tabla_pron, "Semana_8/resultados/tabla_pron_semana8.rds")
saveRDS(ranking,    "Semana_8/resultados/ranking_semana8.rds")

# ---- 6. Grafico: ultimos 120 dias y pronostico ----
hist_df <- tail(datos, 120)
ggplot() +
  geom_line(data = hist_df, aes(Fecha, COLCAP)) +
  geom_ribbon(data = tabla_pron, aes(fecha, ymin = LI95, ymax = LS95), fill = "steelblue", alpha = 0.3) +
  geom_line(data = tabla_pron, aes(fecha, puntual), colour = "steelblue") +
  labs(x = NULL, y = "COLCAP", title = "Pronostico 23-27 nov 2026 con IC 95% (GARCH-t + EWMA)") +
  theme_minimal()
 