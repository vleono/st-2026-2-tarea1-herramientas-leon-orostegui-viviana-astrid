# =========================================================================
# ejemplos/ejemplos.R
# Corre los ocho ejemplos del literal 4 y el contraejemplo del 4(e).
# Se ejecuta desde la raiz del repositorio con:
#   source("ejemplos/ejemplos.R")
# =========================================================================

# ---------- 0. Librerias ------------------------------------------------
library(ggplot2)
library(patchwork)
library(tibble)
library(dplyr)

# ---------- 0c. Deteccion de la raiz del proyecto ----------------------
# Permite correr el script tanto desde la raiz como desde informe/
if (file.exists("R/00-lectura.R")) {
  raiz <- "."
} else if (file.exists("../R/00-lectura.R")) {
  raiz <- ".."
} else {
  stop("No se encuentra la carpeta R/. Ejecuta desde la raiz del repositorio.")
}

# ---------- 0b. Funciones del paquete ----------------------------------
source(file.path(raiz, "R/00-lectura.R"))
source(file.path(raiz, "R/01-graficos.R"))
source(file.path(raiz, "R/02-metodos.R"))
source(file.path(raiz, "R/03-evaluacion.R"))

# =========================================================================
# Protocolo del literal 4(a)
# =========================================================================
ejecutar_ejemplo <- function(serie_ts,
                             nombre_serie,
                             unidad,
                             metodo,
                             rejilla    = NULL,
                             par_tend   = NULL,
                             m_ing      = 1L) {
  
  # ---------- Verificaciones defensivas --------------------------------
  if (!inherits(serie_ts, "ts")) {
    stop("serie_ts debe ser un objeto ts.")
  }
  if (!is.null(dim(serie_ts)) && ncol(serie_ts) > 1) {
    stop("serie_ts debe ser univariado. Extrae una columna antes, p. ej.: ",
         "serie <- Seatbelts[, \"DriversKilled\"]")
  }
  metodo <- match.arg(metodo,
                      c("media", "mm", "ses", "dmm",
                        "tendencia", "holt"))
  
  # ---------- 1. Descripcion de la serie ------------------------------
  y_full <- as.numeric(serie_ts)
  T_full <- length(y_full)
  
  info <- list(
    nombre     = nombre_serie,
    fuente     = paste(utils::help(nombre_serie), collapse = ""),
    unidad     = unidad,
    frecuencia = frequency(serie_ts),
    inicio     = start(serie_ts),
    fin        = end(serie_ts),
    n          = T_full
  )
  
  # ---------- 2. Grafico de la serie y correlograma -------------------
  datos   <- leer_serie(serie_ts, fuente = info$fuente, unidad = unidad)
  g_serie <- graficar_serie(datos, titulo = paste("Serie", nombre_serie))
  
  m_acf    <- min(floor(T_full / 4), 24)
  g_acf    <- correlograma(y_full, m = m_acf)
  r_serie  <- acf_muestral(y_full, lag.max = m_acf)
  lb_serie <- ljung_box(r_serie, T = T_full, m = m_acf, p = 0)
  
  # ---------- 3. Particion estimacion / validacion --------------------
  h <- min(12, floor(0.2 * T_full))
  s <- frequency(serie_ts)
  if (s > 1) h <- max(h, s)
  
  y_est <- y_full[1:(T_full - h)]
  y_val <- y_full[(T_full - h + 1):T_full]
  
  # ---------- 4. Optimizacion y ajuste --------------------------------
  if (!is.null(rejilla)) {
    opt <- optimizar(y_est, metodo, rejilla)
    if (metodo %in% c("mm", "dmm")) {
      par_opt <- list(k = opt$optimo$k)
    } else if (metodo == "ses") {
      par_opt <- list(alpha = opt$optimo$alpha)
    } else if (metodo == "holt") {
      par_opt <- list(alpha = opt$optimo$alpha,
                      beta  = opt$optimo$beta)
    } else {
      par_opt <- list()
    }
  } else {
    opt     <- NULL
    par_opt <- if (!is.null(par_tend)) par_tend else list()
  }
  
  ajuste <- switch(metodo,
                   "media"     = ajustar_media(y_est),
                   "mm"        = ajustar_mm(y_est, k = par_opt$k),
                   "ses"       = ajustar_ses(y_est, alpha = par_opt$alpha),
                   "dmm"       = ajustar_dmm(y_est, k = par_opt$k),
                   "tendencia" = ajustar_tendencia(y_est,
                                                   tipo           = par_opt$tipo,
                                                   corregir_sesgo = par_opt$corregir_sesgo),
                   "holt"      = ajustar_holt(y_est,
                                              alpha = par_opt$alpha,
                                              beta  = par_opt$beta)
  )
  
  # ---------- 5. Medidas de error -------------------------------------
  yhat_est <- ajuste$yhat
  med_est  <- medidas(y_est, yhat_est, m = m_ing)
  
  pron_val <- ajuste$pronosticar(h)
  med_val  <- medidas(y_val, pron_val, m = m_ing)
  
  # Referente ingenuo
  s_ing    <- if (s > 1) s else 1L
  pron_ing <- pronostico_ingenuo(y_est, h = h, s = s_ing)
  med_ing_val <- medidas(y_val, pron_ing, m = m_ing)
  
  # MASE del metodo y del ingenuo
  e_val       <- y_val - pron_val
  e_ing_val   <- y_val - pron_ing
  mase_metodo <- mase(e_val,     y_est, s = s_ing)
  mase_ing    <- mase(e_ing_val, y_est, s = s_ing)
  
  # Tabla comparativa
  comparacion <- data.frame(
    referente = c("Metodo", "Ingenuo"),
    MSE  = c(med_val$MSE,  med_ing_val$MSE),
    MAD  = c(med_val$MAD,  med_ing_val$MAD),
    MAPE = c(med_val$MAPE, med_ing_val$MAPE),
    MASE = c(mase_metodo,  mase_ing),
    N    = c(med_val$N,    med_ing_val$N)
  )
  
  # ---------- 6. Validacion de errores --------------------------------
  e_est    <- y_est - yhat_est
  p_metodo <- switch(metodo,
                     "media"     = 1L,
                     "mm"        = 1L,
                     "ses"       = 1L,
                     "dmm"       = 2L,
                     "tendencia" = if (!is.null(par_opt$tipo) &&
                                       par_opt$tipo == "cuadratica") 3L else 2L,
                     "holt"      = 2L
  )
  
  val_errores <- validar_errores(e_est, p = p_metodo)
  
  # ---------- 6b. Pruebas con los seis elementos ----------------------
  k_prima <- switch(metodo,
                    "media"     = 0L,
                    "mm"        = 1L,
                    "ses"       = 1L,
                    "dmm"       = 2L,
                    "holt"      = 2L,
                    "tendencia" = if (!is.null(par_opt$tipo) &&
                                      par_opt$tipo == "cuadratica") 2L else 1L
  )
  e_clean <- val_errores$errores
  cotas <- obtener_cotas_dw(T = length(e_clean), k_prima = k_prima)
  
  r_est <- acf_muestral(e_clean, lag.max = length(val_errores$acf))
  
  pruebas <- list(
    lb_serie      = prueba_ljung_box(r_serie, T = T_full, m = m_acf, p = 0),
    t_media       = prueba_t_media_cero(e_clean),
    lb_errores    = prueba_ljung_box(
      r_est,
      T = length(e_clean),
      m = length(r_est),
      p = p_metodo
    ),
    jarque_bera   = prueba_jarque_bera(e_clean),
    durbin_watson = prueba_durbin_watson(
      e_clean,
      dL = as.numeric(cotas["dL"]),
      dU = as.numeric(cotas["dU"])
    )
  )
  
  
  # ---------- 7. Grafico final ----------------------------------------
  df_graf <- data.frame(
    t    = seq_len(T_full),
    y    = y_full,
    yhat = NA_real_,
    tipo = "estimacion"
  )
  df_graf$yhat[seq_len(T_full - h)]     <- yhat_est
  df_graf$yhat[(T_full - h + 1):T_full] <- pron_val
  df_graf$tipo[(T_full - h + 1):T_full] <- "validacion"
  
  g_final <- ggplot(df_graf, aes(x = t)) +
    geom_line(aes(y = y), color = "grey30", linewidth = 0.4) +
    geom_line(aes(y = yhat), color = "steelblue", linewidth = 0.7) +
    geom_vline(xintercept = T_full - h + 0.5,
               linetype = "dashed", color = "red") +
    labs(title = paste("Serie", nombre_serie, "-", metodo),
         x = "t", y = unidad) +
    theme_minimal()
  
  # ---------- Salida --------------------------------------------------
  list(
    info        = info,
    datos       = datos,
    g_serie     = g_serie,
    g_acf       = g_acf,
    r_serie     = r_serie,
    lb_serie    = lb_serie,
    h           = h,
    y_est       = y_est,
    y_val       = y_val,
    rejilla     = rejilla,
    opt         = opt,
    par_opt     = par_opt,
    ajuste      = ajuste,
    med_est     = med_est,
    med_val     = med_val,
    med_ing_val = med_ing_val,
    mase_metodo = mase_metodo,
    mase_ing    = mase_ing,
    comparacion = comparacion,
    val_errores = val_errores,
    g_final     = g_final,
    pron_val    = pron_val,
    pron_ing    = pron_ing,
    pruebas     = pruebas
  )
}

# =========================================================================
# Rejillas de optimizacion
# =========================================================================
rej_mm   <- data.frame(k = 2:12)
rej_dmm  <- data.frame(k = 2:12)
rej_ses  <- data.frame(alpha = seq(0.02, 0.98, by = 0.02))
rej_holt <- expand.grid(alpha = seq(0.05, 0.95, by = 0.05),
                        beta  = seq(0.05, 0.95, by = 0.05))

# =========================================================================
# Los ocho ejemplos del literal 4
# =========================================================================

# Ejemplo 1: media simple
ej1 <- ejecutar_ejemplo(
  serie_ts     = sunspot.year,
  nombre_serie = "sunspot.year",
  unidad       = "numero de manchas solares",
  metodo       = "media"
)
# Ejemplo 2: media movil
ej2 <- ejecutar_ejemplo(
  serie_ts     = discoveries,
  nombre_serie = "discoveries",
  unidad       = "numero de descubrimientos",
  metodo       = "mm",
  rejilla      = rej_mm
)

# Ejemplo 3: suavizamiento exponencial simple
ej3 <- ejecutar_ejemplo(
  serie_ts     = Nile,
  nombre_serie = "Nile",
  unidad       = "10^8 m^3",
  metodo       = "ses",
  rejilla      = rej_ses
)

# Ejemplo 4: doble media movil
ej4 <- ejecutar_ejemplo(
  serie_ts     = discoveries,
  nombre_serie = "discoveries",
  unidad       = "numero de descubrimientos",
  metodo       = "dmm",
  rejilla      = rej_dmm
)

# Ejemplo 5: tendencia lineal
ej5 <- ejecutar_ejemplo(
  serie_ts     = airmiles,
  nombre_serie = "airmiles",
  unidad       = "millas pasajero-ano",
  metodo       = "tendencia",
  par_tend     = list(tipo = "lineal")
)

# Ejemplo 6: tendencia cuadratica
ej6 <- ejecutar_ejemplo(
  serie_ts     = co2,
  nombre_serie = "co2",
  unidad       = "ppm",
  metodo       = "tendencia",
  par_tend     = list(tipo = "cuadratica")
)

# Ejemplo 7: tendencia exponencial
ej7 <- ejecutar_ejemplo(
  serie_ts     = JohnsonJohnson,
  nombre_serie = "JohnsonJohnson",
  unidad       = "dolares por accion",
  metodo       = "tendencia",
  par_tend     = list(tipo = "exponencial", corregir_sesgo = TRUE),
  m_ing        = 4
)

# Ejemplo 8: Holt lineal
ej8 <- ejecutar_ejemplo(
  serie_ts     = airmiles,
  nombre_serie = "airmiles",
  unidad       = "millas pasajero-ano",
  metodo       = "holt",
  rejilla      = rej_holt
)

# =========================================================================
# Contraejemplo del literal 4(e)
# =========================================================================
ej_ce <- ejecutar_ejemplo(
  serie_ts     = AirPassengers,
  nombre_serie = "AirPassengers (contraejemplo)",
  unidad       = "miles de pasajeros",
  metodo       = "holt",
  rejilla      = rej_holt,
  m_ing        = 12
)
# =========================================================================
# Guardar figuras en figs/
# =========================================================================
guardar_fig <- function(ej, prefijo) {
  ggsave(file.path(raiz, "figs", paste0(prefijo, "-serie.png")),
         ej$g_serie, width = 8, height = 4, dpi = 100)
  ggsave(file.path(raiz, "figs", paste0(prefijo, "-acf.png")),
         ej$g_acf, width = 8, height = 6, dpi = 100)
  ggsave(file.path(raiz, "figs", paste0(prefijo, "-final.png")),
         ej$g_final, width = 8, height = 4, dpi = 100)
}

guardar_fig(ej1,   "ej1-media")
guardar_fig(ej2,   "ej2-mm")
guardar_fig(ej3,   "ej3-ses")
guardar_fig(ej4,   "ej4-dmm")
guardar_fig(ej5,   "ej5-tend-lineal")
guardar_fig(ej6,   "ej6-tend-cuadratica")
guardar_fig(ej7,   "ej7-tend-exp")
guardar_fig(ej8,   "ej8-holt")
guardar_fig(ej_ce, "ej-contraejemplo")

# =========================================================================
# Verificacion contra stats
# =========================================================================
cat("\n=== Verificacion ACF contra stats::acf() ===\n")
print(verificar_acf(as.numeric(AirPassengers), lag.max = 24))

cat("\n=== Verificacion Qm contra stats::Box.test() ===\n")
print(verificar_ljung_box(as.numeric(AirPassengers), lag.max = 24, p = 0))

# =========================================================================
# Tabla resumen para el README
# =========================================================================
tabla_readme <- do.call(rbind, lapply(
  list(ej1, ej2, ej3, ej4, ej5, ej6, ej7, ej8),
  function(ej) {
    data.frame(
      serie       = ej$info$nombre,
      parametros  = paste(names(ej$par_opt), unlist(ej$par_opt),
                          sep = "=", collapse = ", "),
      MASE_metodo = round(ej$mase_metodo, 4),
      MASE_ing    = round(ej$mase_ing, 4)
    )
  }))

cat("\n=== Tabla resumen ===\n")
print(tabla_readme)

# =========================================================================
# sesion-info.txt
# =========================================================================
writeLines(capture.output(sessionInfo()), file.path(raiz, "sesion-info.txt"))

cat("\nListo. Ocho ejemplos + contraejemplo ejecutados.\n")
