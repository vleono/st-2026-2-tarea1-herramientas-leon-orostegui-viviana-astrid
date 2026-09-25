# =========================================================================
# 03-evaluacion.R
# Medidas de error, validacion de residuos y pruebas de hipotesis.
#
# Contiene:
#   1. Medidas de error (literal 4a)
#   2. Pruebas de calculo puro (literal 2d)
#   3. Auxiliares para las pruebas (t, validacion)
#   4. Pruebas con los seis elementos (literal 4b)
#   5. Referente ingenuo y MASE (literal 4c)
#   6. Tabla de cotas Durbin-Watson y helper
# =========================================================================


# =========================================================================
# 1. MEDIDAS DE ERROR (literal 4a)
# =========================================================================

# -------------------------------------------------------------------------
# medidas(y, yhat, m)
#
# y:    valores observados.
# yhat: valores ajustados o pronosticados (mismo length que y).
# m:    periodo estacional para el MASE (1 si no estacional, NULL si
#       no se quiere MASE).
#
# Devuelve MSE, MAD, MAPE, MASE y el numero de errores no NA.
# -------------------------------------------------------------------------
medidas <- function(y, yhat, m = NULL) {
  stopifnot(length(y) == length(yhat))
  ok <- !is.na(yhat)
  if (sum(ok) == 0) {
    return(list(MSE = NA, MAD = NA, MAPE = NA, MASE = NA, N = 0L))
  }
  e   <- y[ok] - yhat[ok]
  N   <- length(e)
  MSE <- mean(e^2)
  MAD <- mean(abs(e))
  MAPE <- if (all(y[ok] != 0)) mean(abs(e / y[ok])) * 100 else NA_real_
  
  if (!is.null(m) && m >= 1) {
    difs  <- diff(y, lag = m)
    denom <- if (length(difs) > 0) mean(abs(difs), na.rm = TRUE) else NA_real_
    MASE  <- if (is.finite(denom) && denom > 0) MAD / denom else NA_real_
  } else {
    MASE <- NA_real_
  }
  list(MSE = MSE, MAD = MAD, MAPE = MAPE, MASE = MASE, N = N)
}


# =========================================================================
# 2. PRUEBAS DE CALCULO PURO (literal 2d)
# =========================================================================

# -------------------------------------------------------------------------
# ljung_box(r, T, m, p)
#
# Q_m = T(T+2) sum_{h=1}^m r_h^2 / (T - h), referido a chi^2_{m-p}.
# r: vector de autocorrelaciones muestrales r_1, ..., r_m.
# T: numero de observaciones.
# m: numero de rezagos.
# p: numero de parametros estimados (0 si sobre la serie).
#
# Devuelve: estadistico, gl, critico al 5 %, valor p.
# -------------------------------------------------------------------------
ljung_box <- function(r, T, m, p = 0) {
  stopifnot(is.numeric(r), length(r) >= m, m >= 1, T >= 1, p >= 0)
  Q  <- T * (T + 2) * sum(r[1:m]^2 / (T - (1:m)))
  gl <- m - p
  if (gl <= 0) stop("Los grados de libertad m - p deben ser positivos.")
  crit <- qchisq(0.95, df = gl)
  pval <- 1 - pchisq(Q, df = gl)
  list(estadistico = Q, gl = gl, critico = crit, p_valor = pval)
}
# -------------------------------------------------------------------------
# verificar_ljung_box(y, lag.max, p)
#
# Bloque de verificacion: contrasta Qm calculado a mano contra el de
# Box.test() de stats.
#
# y:       vector numerico (la serie, o los errores de un metodo).
# lag.max: numero de rezagos m.
# p:       numero de parametros estimados (0 si sobre la serie).
#
# Devuelve la maxima diferencia absoluta en Qm y si cumple < 1e-10.
# -------------------------------------------------------------------------
verificar_ljung_box <- function(y, lag.max = NULL, p = 0) {
  stopifnot(is.numeric(y), length(y) >= 5, !anyNA(y))
  
  T <- length(y)
  if (is.null(lag.max)) lag.max <- min(floor(T / 4), 24)
  
  # Qm a mano
  r      <- acf_muestral(y, lag.max = lag.max)
  Q_prop <- ljung_box(r, T = T, m = lag.max, p = p)$estadistico
  
  # Qm de stats::Box.test()
  bt <- stats::Box.test(y, lag = lag.max, type = "Ljung-Box",
                        fitdf = p)
  Q_stats <- as.numeric(bt$statistic)
  
  dif_max <- abs(Q_prop - Q_stats)
  
  list(
    Q_propio = Q_prop,
    Q_stats  = Q_stats,
    dif_max  = dif_max,
    cumple   = dif_max < 1e-10
  )
}
# -------------------------------------------------------------------------
# jarque_bera(e)
#
# JB = N/6 [A^2 + (K-3)^2/4], referido a chi^2_2.
# e: vector de errores.
#
# Devuelve: estadistico, gl, critico al 5 %, valor p, N, asimetria,
# curtosis.
# -------------------------------------------------------------------------
jarque_bera <- function(e) {
  e <- e[!is.na(e)]
  N <- length(e)
  if (N < 4) {
    return(list(estadistico = NA_real_, gl = 2L,
                critico = qchisq(0.95, 2), p_valor = NA_real_,
                N = N, asimetria = NA_real_, curtosis = NA_real_))
  }
  m2 <- mean((e - mean(e))^2)
  m3 <- mean((e - mean(e))^3)
  m4 <- mean((e - mean(e))^4)
  A  <- m3 / m2^(3/2)
  K  <- m4 / m2^2
  JB <- N / 6 * (A^2 + (K - 3)^2 / 4)
  
  list(estadistico = JB, gl = 2L,
       critico = qchisq(0.95, 2),
       p_valor = 1 - pchisq(JB, df = 2),
       N = N, asimetria = A, curtosis = K)
}

# -------------------------------------------------------------------------
# durbin_watson(e)
#
# d = sum_{t=2}^T (e_t - e_{t-1})^2 / sum_{t=1}^T e_t^2.
# Se compara contra dL y dU de la tabla, no contra un chi^2.
#
# Devuelve: estadistico, gl (NA), critico (NA), valor p (NA).
# -------------------------------------------------------------------------
durbin_watson <- function(e) {
  e <- e[!is.na(e)]
  d <- sum(diff(e)^2) / sum(e^2)
  list(estadistico = d, gl = NA_integer_,
       critico = NA_real_, p_valor = NA_real_)
}


# =========================================================================
# 3. AUXILIARES
# =========================================================================

# -------------------------------------------------------------------------
# t_media_cero(e)
#
# t = mean(e) / (sd(e) / sqrt(N)), referido a t_{N-1}.
# -------------------------------------------------------------------------
t_media_cero <- function(e) {
  e  <- e[!is.na(e)]
  N  <- length(e)
  tt <- mean(e) / (sd(e) / sqrt(N))
  gl <- N - 1
  list(t = tt, gl = gl,
       critico = qt(0.975, df = gl),
       p_valor = 2 * pt(-abs(tt), df = gl),
       N = N)
}

# -------------------------------------------------------------------------
# validar_errores(e, m, p, m_lb)
#
# Reune las pruebas de validacion de los errores de un paso:
#   t de media cero, Ljung-Box (m - p gl), Jarque-Bera, Durbin-Watson,
#   y la ACF a mano.
# -------------------------------------------------------------------------
validar_errores <- function(e, m = NULL, p = 0, m_lb = NULL) {
  e <- e[!is.na(e)]
  N <- length(e)
  if (is.null(m_lb)) m_lb <- min(floor(N / 4), 12)
  
  r <- acf_muestral(e, lag.max = m_lb)
  
  list(
    errores       = e,
    N             = N,
    t_media       = t_media_cero(e),
    ljung_box     = ljung_box(r, T = N, m = m_lb, p = p),
    jarque_bera   = jarque_bera(e),
    durbin_watson = durbin_watson(e),
    acf           = r,
    m_lb          = m_lb,
    p             = p
  )
}


# =========================================================================
# 4. PRUEBAS CON LOS SEIS ELEMENTOS (literal 4b)
# =========================================================================

# -------------------------------------------------------------------------
# prueba_rh_individual(r, T, h)
#
# H0: rho_h = 0 vs Ha: rho_h != 0.
# Z = r_h / (1/sqrt(T)), asintotica N(0,1) bajo H0.
# -------------------------------------------------------------------------
prueba_rh_individual <- function(r, T, h) {
  stopifnot(h >= 1, h <= length(r))
  rh   <- r[h]
  sd_r <- 1 / sqrt(T)
  z    <- rh / sd_r
  crit <- qnorm(0.975)
  pval <- 2 * pnorm(-abs(z))
  decision <- if (abs(z) > crit) "Se rechaza H0" else "No se rechaza H0"
  lectura <- if (abs(z) > crit) {
    paste0("El rezago h=", h, " tiene autocorrelacion muestral ",
           "significativa al 5 %.")
  } else {
    paste0("El rezago h=", h, " no es distinguible de cero al 5 %.")
  }
  list(
    nombre        = paste0("Contraste individual de r_", h),
    H0            = paste0("rho_", h, " = 0"),
    Ha            = paste0("rho_", h, " != 0"),
    estadistico   = paste0("Z = r_", h, " / (1/sqrt(T))"),
    distribucion  = "N(0,1) bajo H0, asintotica",
    gl            = NA_integer_,
    valor_obs     = z,
    valor_p       = pval,
    region        = paste0("|Z| > ", round(crit, 4)),
    decision      = decision,
    lectura       = lectura
  )
}

# -------------------------------------------------------------------------
# prueba_ljung_box(r, T, m, p)
#
# H0: rho_1 = ... = rho_m = 0 vs Ha: algun rho_h != 0.
# Q_m ~ chi^2_{m - p} bajo H0.
# -------------------------------------------------------------------------
prueba_ljung_box <- function(r, T, m, p = 0) {
  res  <- ljung_box(r, T, m, p = p)
  gl   <- m - p
  crit <- qchisq(0.95, df = gl)
  decision <- if (res$estadistico > crit) "Se rechaza H0" else "No se rechaza H0"
  lectura <- if (res$estadistico > crit) {
    paste0("Los primeros ", m, " rezagos contienen autocorrelacion ",
           "conjunta significativa. Queda estructura lineal sin capturar.")
  } else {
    paste0("No hay evidencia de autocorrelacion conjunta en los primeros ",
           m, " rezagos.")
  }
  list(
    nombre        = if (p == 0) "Ljung-Box sobre la serie" else "Ljung-Box sobre los errores",
    H0            = paste0("rho_1 = rho_2 = ... = rho_", m, " = 0"),
    Ha            = paste0("Algun rho_h != 0, h <= ", m),
    estadistico   = "Q_m = T(T+2) sum_{h=1}^m r_h^2 / (T-h)",
    distribucion  = paste0("chi^2_", gl, " bajo H0"),
    gl            = gl,
    valor_obs     = res$estadistico,
    valor_p       = res$p_valor,
    region        = paste0("Q_m > ", round(crit, 4)),
    decision      = decision,
    lectura       = lectura
  )
}

# -------------------------------------------------------------------------
# prueba_t_media_cero(e)
#
# H0: E(e_t) = 0 vs Ha: E(e_t) != 0.
# t = mean(e) / (sd(e)/sqrt(N)), t_{N-1} bajo H0.
# -------------------------------------------------------------------------
prueba_t_media_cero <- function(e) {
  res  <- t_media_cero(e)
  crit <- qt(0.975, df = res$gl)
  decision <- if (abs(res$t) > crit) "Se rechaza H0" else "No se rechaza H0"
  lectura <- if (abs(res$t) > crit) {
    "El error medio es distinto de cero; hay sesgo sistematico en el pronostico."
  } else {
    "No hay evidencia de sesgo sistematico en el pronostico."
  }
  list(
    nombre        = "t de media cero sobre los errores",
    H0            = "E(e_t) = 0",
    Ha            = "E(e_t) != 0",
    estadistico   = "t = mean(e) / (sd(e)/sqrt(N))",
    distribucion  = paste0("t_", res$gl, " bajo H0"),
    gl            = res$gl,
    valor_obs     = res$t,
    valor_p       = res$p_valor,
    region        = paste0("|t| > ", round(crit, 4)),
    decision      = decision,
    lectura       = lectura
  )
}

# -------------------------------------------------------------------------
# prueba_jarque_bera(e)
#
# H0: A = 0 y K = 3 (normalidad) vs Ha: A != 0 o K != 3.
# JB ~ chi^2_2 bajo H0. Se advierte si N < 20.
# -------------------------------------------------------------------------
prueba_jarque_bera <- function(e) {
  res  <- jarque_bera(e)
  crit <- qchisq(0.95, df = 2)
  decision <- if (is.na(res$estadistico)) "No calculable (N < 4)" else
    if (res$estadistico > crit) "Se rechaza H0" else "No se rechaza H0"
  lectura <- paste0(
    "N = ", res$N, ". ",
    if (res$N < 20) {
      "Con menos de 20 errores, JB es asintotica y su decision es solo indicativa. "
    } else "",
    if (is.na(res$estadistico)) "No calculable." else
      if (res$estadistico > crit) "Los errores no son normales." else
        "No hay evidencia contra la normalidad."
  )
  list(
    nombre        = "Jarque-Bera",
    H0            = "A = 0 y K = 3 (normalidad)",
    Ha            = "A != 0 o K != 3",
    estadistico   = "JB = N/6 [A^2 + (K-3)^2/4]",
    distribucion  = "chi^2_2 bajo H0",
    gl            = 2L,
    valor_obs     = res$estadistico,
    valor_p       = res$p_valor,
    region        = paste0("JB > ", round(crit, 4)),
    decision      = decision,
    lectura       = lectura
  )
}

# -------------------------------------------------------------------------
# prueba_durbin_watson(e, dL, dU)
#
# H0: rho = 0 vs Ha: rho > 0.
# Se compara d contra dL y dU de la tabla.
# -------------------------------------------------------------------------
prueba_durbin_watson <- function(e, dL, dU) {
  d <- durbin_watson(e)$estadistico
  if (d < dL) {
    decision <- "Se rechaza H0 (autocorrelacion positiva)"
    lectura  <- paste0("d = ", round(d, 4), " < dL = ", dL,
                       ". Los errores tienen autocorrelacion positiva; ",
                       "los valores p ordinarios no son interpretables.")
  } else if (d > dU) {
    decision <- "No se rechaza H0"
    lectura  <- paste0("d = ", round(d, 4), " > dU = ", dU,
                       ". No hay evidencia de autocorrelacion positiva.")
  } else {
    decision <- "Zona de indecision"
    lectura  <- paste0("d = ", round(d, 4), " esta entre dL = ", dL,
                       " y dU = ", dU, ". La prueba no concluye.")
  }
  list(
    nombre        = "Durbin-Watson",
    H0            = "rho = 0 (sin autocorrelacion de primer orden)",
    Ha            = "rho > 0",
    estadistico   = "d = sum_{t=2}^T (e_t - e_{t-1})^2 / sum_{t=1}^T e_t^2",
    distribucion  = paste0("Cotas: dL = ", dL, ", dU = ", dU, " al 5 %"),
    gl            = NA_integer_,
    valor_obs     = d,
    valor_p       = NA_real_,
    region        = paste0("d < ", dL),
    decision      = decision,
    lectura       = lectura
  )
}

# -------------------------------------------------------------------------
# prueba_t_coeficiente(nombre_coef, estimacion, se_robusto, gl, H0_texto)
#
# Contraste t sobre un coeficiente de tendencia con EE robusto HAC.
# -------------------------------------------------------------------------
prueba_t_coeficiente <- function(nombre_coef, estimacion, se_robusto, gl,
                                 H0_texto = NULL) {
  if (is.null(H0_texto)) H0_texto <- paste0(nombre_coef, " = 0")
  t_obs <- estimacion / se_robusto
  crit  <- qt(0.975, df = gl)
  pval  <- 2 * pt(-abs(t_obs), df = gl)
  decision <- if (abs(t_obs) > crit) "Se rechaza H0" else "No se rechaza H0"
  lectura <- if (abs(t_obs) > crit) {
    paste0("El coeficiente ", nombre_coef,
           " es significativo al 5 % con EE robusto.")
  } else {
    paste0("El coeficiente ", nombre_coef,
           " no es distinguible de cero al 5 % con EE robusto.")
  }
  list(
    nombre        = paste0("t sobre ", nombre_coef, " (EE robusto)"),
    H0            = H0_texto,
    Ha            = sub("= 0", "!= 0", H0_texto),
    estadistico   = paste0("t = ", nombre_coef, "_hat / se_robusto"),
    distribucion  = paste0("t_", gl, " bajo H0"),
    gl            = gl,
    valor_obs     = t_obs,
    valor_p       = pval,
    region        = paste0("|t| > ", round(crit, 4)),
    decision      = decision,
    lectura       = lectura
  )
}

# -------------------------------------------------------------------------
# imprimir_prueba(p)
#
# Imprime una prueba con sus seis elementos en orden.
# -------------------------------------------------------------------------
imprimir_prueba <- function(p) {
  cat("===", p$nombre, "===\n")
  cat("1. H0:", p$H0, "\n")
  cat("   Ha:", p$Ha, "\n")
  cat("2. Estadistico:", p$estadistico, "\n")
  cat("   Distribucion:", p$distribucion, "\n")
  if (!is.na(p$gl)) cat("   gl =", p$gl, "\n")
  cat("3. Region de rechazo al 5 %:", p$region, "\n")
  cat("4. Valor observado:", round(p$valor_obs, 4),
      "| valor p:", round(p$valor_p, 4), "\n")
  cat("5. Decision:", p$decision, "\n")
  cat("6. Lectura:", p$lectura, "\n\n")
  invisible(p)
}


# =========================================================================
# 5. REFERENTE INGENUO Y MASE (literal 4c)
# =========================================================================

# -------------------------------------------------------------------------
# pronostico_ingenuo(y_est, h, s)
#
# y_est: tramo de estimacion.
# h:     horizonte.
# s:     periodo estacional (1 si no estacional).
# -------------------------------------------------------------------------
pronostico_ingenuo <- function(y_est, h, s = 1L) {
  stopifnot(is.numeric(y_est), length(y_est) >= s, h >= 1, s >= 1,
            s == as.integer(s))
  if (s == 1L) {
    rep(y_est[length(y_est)], h)
  } else {
    base <- y_est[(length(y_est) - s + 1):length(y_est)]
    rep(base, length.out = h)
  }
}

# -------------------------------------------------------------------------
# mase(e_val, y_est, s)
#
# MASE = MAD(e_val) / MAD_ingenuo(y_est).
# El denominador se calcula DENTRO del tramo de estimacion.
# -------------------------------------------------------------------------
mase <- function(e_val, y_est, s = 1L) {
  stopifnot(is.numeric(e_val), is.numeric(y_est), s >= 1)
  if (s == 1L) {
    dif_est <- diff(y_est)
  } else {
    dif_est <- y_est[(s + 1):length(y_est)] - y_est[1:(length(y_est) - s)]
  }
  mad_ing_est <- mean(abs(dif_est), na.rm = TRUE)
  mad_val     <- mean(abs(e_val),   na.rm = TRUE)
  mad_val / mad_ing_est
}


# =========================================================================
# 6. TABLA DE COTAS DURBIN-WATSON Y HELPER
# =========================================================================

# Cotas al 5 % para k' = 1 (un regresor sin intercepto) y k' = 2.
dw_cotas <- data.frame(
  T     = c(15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200),
  dL_k1 = c(1.08, 1.20, 1.29, 1.35, 1.44, 1.50, 1.55, 1.61, 1.65, 1.72, 1.76),
  dU_k1 = c(1.36, 1.41, 1.45, 1.49, 1.54, 1.59, 1.62, 1.66, 1.69, 1.75, 1.78),
  dL_k2 = c(0.95, 1.10, 1.21, 1.28, 1.39, 1.46, 1.51, 1.59, 1.63, 1.71, 1.75),
  dU_k2 = c(1.54, 1.54, 1.55, 1.57, 1.60, 1.63, 1.65, 1.69, 1.72, 1.77, 1.80)
)

# -------------------------------------------------------------------------
# obtener_cotas_dw(T, k_prima)
#
# Interpola dL y dU segun T y el numero de regresores sin intercepto.
# k_prima = 0 (media simple), 1 (lineal, exponencial, mm, ses) o 2
# (cuadratica, dmm, holt).
# -------------------------------------------------------------------------
obtener_cotas_dw <- function(T, k_prima) {
  stopifnot(is.numeric(T), T >= 1, k_prima %in% c(0L, 1L, 2L))
  
  if (k_prima == 0L) {
    # Sin regresores: usar cotas de k' = 1 como aproximacion conservadora.
    dL_col <- dw_cotas$dL_k1
    dU_col <- dw_cotas$dU_k1
  } else if (k_prima == 1L) {
    dL_col <- dw_cotas$dL_k1
    dU_col <- dw_cotas$dU_k1
  } else {
    dL_col <- dw_cotas$dL_k2
    dU_col <- dw_cotas$dU_k2
  }
  
  if (T <= min(dw_cotas$T)) {
    return(c(dL = dL_col[1], dU = dU_col[1]))
  }
  if (T >= max(dw_cotas$T)) {
    return(c(dL = tail(dL_col, 1), dU = tail(dU_col, 1)))
  }
  c(dL = approx(dw_cotas$T, dL_col, xout = T)$y,
    dU = approx(dw_cotas$T, dU_col, xout = T)$y)
}
