############################################################################
##  PROYECTO 1 - SERIES DE TIEMPO UNIVARIADAS
##  Universidad Nacional de Colombia - Depto. de Estadistica
##
##  SEMANA 3 : pronostico del MSCI COLCAP
##  Corte de datos : viernes 16 de octubre de 2026
##  A pronosticar  : lunes 19 a viernes 23 de octubre de 2026
##
##  Restriccion de la semana: se permiten suavizamiento, ETS, benchmarks y
##  ahora tambien ARMA/ARIMA/SARIMA. En ESTE grupo decidimos NO usar ARIMA
##  y aprovechar la opcion de Transformaciones
##  (log, Box-Cox) y Medias moviles.
##  SIN variables exogenas.
##
##  Autor/grupo : Grupo 6 
##  Fecha        : 2026-10-18
############################################################################


library(fpp3)        
library(readxl)    
library(purrr)      

set.seed(123)

## ==========================================================================
## 1. LECTURA Y LIMPIEZA DE LOS DATOS
## ==========================================================================

archivo <- file.choose()      
datos <- read.csv(archivo) |>
  mutate(Fecha = as.Date(Fecha)) |>
  arrange(Fecha) |>
  distinct(Fecha, .keep_all = TRUE)

#chequeo

chequeo <- datos |> filter(Fecha %in% as.Date(c("2026-09-17","2026-09-18")))
print(chequeo)
stopifnot(
  nrow(chequeo) == 2,
  abs(chequeo$COLCAP[chequeo$Fecha == as.Date("2026-09-17")] - 2521.85) < 1e-6,
  abs(chequeo$COLCAP[chequeo$Fecha == as.Date("2026-09-18")] - 2548.22) < 1e-6
)

corte <- as.Date("2026-10-16")
datos <- datos |> filter(Fecha <= corte)
stopifnot(max(datos$Fecha) == corte)  

cat("\nObservaciones (toda la historia):", nrow(datos),
    "| del", format(min(datos$Fecha)), "al", format(max(datos$Fecha)),
    "| ultimo valor =", tail(datos$COLCAP, 1), "\n")

## ---- Copia identificable de los datos utilizados (regla 5e) ---------------
write.csv(datos, "COLCAP_datos_usados_Semana3_corte_2026-10-16.csv",
          row.names = FALSE)

## ==========================================================================
## 2. VENTANA DE ENTRENAMIENTO Y SERIE IGUALMENTE ESPACIADA
## ==========================================================================
## Quiebre metodologico: el MSCI COLCAP reemplazo al bvc COLCAP el 28/05/2021.
## Usamos como ventana base la historia del MSCI COLCAP (>= 2021-05-28) para
## trabajar sobre una sola metodologia.
inicio_ventana <- as.Date("2021-05-28")

serie <- datos |>
  filter(Fecha >= inicio_ventana) |>
  arrange(Fecha) |>
  mutate(t = row_number()) |>
  as_tsibble(index = t)

cat("Ventana de entrenamiento:", format(min(serie$Fecha)), "a",
    format(max(serie$Fecha)), "|", nrow(serie), "observaciones\n\n")

## Parametro Box-Cox
lambda_bc <- serie |>
  features(COLCAP, features = guerrero) |>
  pull(lambda_guerrero)
cat("Lambda Box-Cox (Guerrero):", round(lambda_bc, 4), "\n")

autoplot(serie, COLCAP) +
  labs(title = "MSCI COLCAP - ventana de entrenamiento (dias de negociacion)",
       x = "Observacion de negociacion (t)", y = "Indice (puntos)")

## ==========================================================================
## 3. HORIZONTE Y FECHAS A PRONOSTICAR
## ==========================================================================
h <- 5L
fechas_obj <- as.Date(c("2026-10-19","2026-10-20","2026-10-21",
                        "2026-10-22","2026-10-23")) 

## ==========================================================================
## 4. DEFINICION DE MODELOS (Semana 3: suavizamiento/ETS/benchmarks +
##    transformaciones y medias moviles; SIN ARIMA, SIN exogenas)
## ==========================================================================
ajustar <- list(
  ## --- benchmarks de clase ---
  "Naive"        = function(d) model(d, m = NAIVE(COLCAP)),
  "Drift"        = function(d) model(d, m = RW(COLCAP ~ drift())),
  "Media"        = function(d) model(d, m = MEAN(COLCAP)),
  ## --- medias moviles (transformacion por promedio local) ---
  "MM20"         = function(d) model(d, m = MEAN(COLCAP ~ window(size = 20))),
  "MM60"         = function(d) model(d, m = MEAN(COLCAP ~ window(size = 60))),
  ## --- suavizamiento / ETS ---
  "SES_ANN"      = function(d) model(d, m = ETS(COLCAP ~ error("A")+trend("N")+season("N"))),
  "Holt_AAN"     = function(d) model(d, m = ETS(COLCAP ~ error("A")+trend("A")+season("N"))),
  "Damped_AAdN"  = function(d) model(d, m = ETS(COLCAP ~ error("A")+trend("Ad")+season("N"))),
  "ETS_auto"     = function(d) model(d, m = ETS(COLCAP)),
  ## --- ETS sobre transformaciones ---
  "ETS_log"      = function(d) model(d, m = ETS(log(COLCAP))),
  "ETS_boxcox"   = function(d) model(d, m = ETS(box_cox(COLCAP, lambda_bc)))
)

## ==========================================================================
## 5. FUNCION DE EVALUACION (rolling origin / origen movil)  -- h = 5, RMSE
## ==========================================================================
evaluar_modelo <- function(serie, modelo,
                           tipo_ventana = c("expansible","fija"),
                           recalibrar   = TRUE,
                           ancho_fijo   = 500L,
                           inicio_train = NULL,
                           paso         = 1L,
                           h            = 5L) {

  tipo_ventana <- match.arg(tipo_ventana)
  fun_fit <- ajustar[[modelo]]
  n <- nrow(serie)
  if (is.null(inicio_train)) inicio_train <- floor(0.6 * n)

  errores <- list(); i <- 1; t0 <- inicio_train; fit_fijo <- NULL

  while ((t0 + h) <= n) {
    if (tipo_ventana == "expansible") { ini <- 1L; fin <- t0 }
    else { ini <- max(1L, t0 - ancho_fijo + 1L); fin <- t0 }
    train <- serie[ini:fin, ]

    pred <- tryCatch({
      if (recalibrar) {
        fit <- fun_fit(train)
      } else {
        if (is.null(fit_fijo)) fit_fijo <- fun_fit(train)
        fit <- refit(fit_fijo, train, reestimate = FALSE)
      }
      forecast(fit, h = h)$.mean
    }, error = function(e) rep(NA_real_, h))

    y_real <- serie$COLCAP[(t0 + 1):(t0 + h)]
    errores[[i]] <- (y_real - pred)^2
    t0 <- t0 + paso; i <- i + 1
  }

  E <- do.call(rbind, errores)
  list(rmse_paso  = sqrt(colMeans(E, na.rm = TRUE)),
       rmse_multi = sqrt(mean(E,     na.rm = TRUE)),
       n_origenes = nrow(E))
}

## ==========================================================================
## 6. ESCENARIOS DE VALIDACION
## ==========================================================================
escenarios <- tribble(
  ~nombre,                   ~modelo,        ~ventana,      ~recalibrar,
  "Naive",                   "Naive",        "expansible",  TRUE,
  "Drift",                   "Drift",        "expansible",  TRUE,
  "Media",                   "Media",        "expansible",  TRUE,
  "Media movil 20",          "MM20",         "expansible",  TRUE,
  "Media movil 60",          "MM60",         "expansible",  TRUE,
  "SES (A,N,N)",             "SES_ANN",      "expansible",  TRUE,
  "Holt (A,A,N)",            "Holt_AAN",     "expansible",  TRUE,
  "Holt amortig. (A,Ad,N)",  "Damped_AAdN",  "expansible",  TRUE,
  "ETS auto (AICc)",         "ETS_auto",     "expansible",  TRUE,
  "ETS sobre log",           "ETS_log",      "expansible",  TRUE,
  "ETS sobre Box-Cox",       "ETS_boxcox",   "expansible",  TRUE,
  "ETS auto - ventana fija", "ETS_auto",     "fija",        TRUE,
  "Drift - ventana fija",    "Drift",        "fija",        TRUE
)

## ==========================================================================
## 7. EJECUCION DE LA VALIDACION
## ==========================================================================

paso_cv <- 1L

resultados <- list()
for (k in seq_len(nrow(escenarios))) {
  e <- escenarios[k, ]
  cat("Evaluando:", e$nombre, "...\n")
  resultados[[e$nombre]] <- evaluar_modelo(
    serie, modelo = e$modelo, tipo_ventana = e$ventana,
    recalibrar = e$recalibrar, paso = paso_cv, h = h
  )
}

## ==========================================================================
## 8. TABLA RESUMEN (RMSE por paso y multi-paso)
## ==========================================================================
tabla_cv <- imap_dfr(resultados, function(res, nombre) {
  tibble(Modelo = nombre,
         RMSE_h1 = res$rmse_paso[1], RMSE_h2 = res$rmse_paso[2],
         RMSE_h3 = res$rmse_paso[3], RMSE_h4 = res$rmse_paso[4],
         RMSE_h5 = res$rmse_paso[5], RMSE_multi = res$rmse_multi,
         n_orig = res$n_origenes)
}) |> arrange(RMSE_multi)

cat("RESULTADOS DE VALIDACION")
print(as.data.frame(tabla_cv), digits = 5)
write.csv(tabla_cv, "COLCAP_validacion_CV_Semana3.csv", row.names = FALSE)

## ==========================================================================
## 9. MODELO GANADOR
## ==========================================================================
ganador <- tabla_cv$Modelo[1]
modelo_ganador <- escenarios$modelo[escenarios$nombre == ganador]
cat("\n>> Modelo seleccionado por RMSE multi-paso:", ganador,
    "(familia:", modelo_ganador, ")\n")

## ==========================================================================
## 10. AJUSTE FINAL SOBRE TODA LA VENTANA Y PRONOSTICO h = 5
## ==========================================================================
modelos_finales <- serie |>
  model(
    Naive       = NAIVE(COLCAP),
    Drift       = RW(COLCAP ~ drift()),
    MM20        = MEAN(COLCAP ~ window(size = 20)),
    SES_ANN     = ETS(COLCAP ~ error("A")+trend("N")+season("N")),
    Holt_AAN    = ETS(COLCAP ~ error("A")+trend("A")+season("N")),
    Damped_AAdN = ETS(COLCAP ~ error("A")+trend("Ad")+season("N")),
    ETS_auto    = ETS(COLCAP),
    ETS_log     = ETS(log(COLCAP)),
    ETS_boxcox  = ETS(box_cox(COLCAP, lambda_bc))
  )

cat("Especificacion del ETS automatico elegido")
print(report(select(modelos_finales, ETS_auto)))

fc_todos <- forecast(modelos_finales, h = h)

tabla_pron <- fc_todos |>
  as_tibble() |>
  mutate(paso = t - max(serie$t)) |>
  select(.model, paso, .mean) |>
  pivot_wider(names_from = paso, values_from = .mean, names_prefix = "h") |>
  rename(Modelo = .model) |>
  mutate(across(starts_with("h"), ~round(.x, 2)))

cat("PRONOSTICOS PUNTUALES (lun19 - vie23 de octubre)")
cat("h1=19  h2=20  h3=21  h4=22  h5=23 (octubre 2026)\n")
print(as.data.frame(tabla_pron))

## Pronostico FINAL a reportar = modelo ganador en CV
entrega <- tibble(
  Fecha      = fechas_obj,
  Paso_h     = 1:5,
  Modelo     = ganador,
  Prediccion = round(forecast(select(modelos_finales, all_of(modelo_ganador)),
                              h = h)$.mean, 2)
)
cat("PRONOSTICO (Semana 3)")
print(as.data.frame(entrega))
write.csv(entrega, "COLCAP_pronostico_Semana3.csv", row.names = FALSE)


## 11. GRAFICO DEL PRONOSTICO FINAL

autoplot(filter(fc_todos, .model == modelo_ganador), tail(serie, 60)) +
  labs(title = paste0("MSCI COLCAP - pronostico Semana 3 (", ganador, ")"),
       subtitle = "Ultimas 60 observaciones + 5 pasos (19-23 oct 2026)",
       x = "Observacion de negociacion (t)", y = "Indice (puntos)")


## 12. Reproducibilidad

cat("sessionInfo (conservar para reproducibilidad)")
print(sessionInfo())

