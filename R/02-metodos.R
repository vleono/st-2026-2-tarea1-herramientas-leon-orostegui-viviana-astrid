
# =========================================================================
# 02-metodos.R
# Los ocho metodos del literal 3 y la funcion de optimizacion.
#
# No usa funciones de forecast, fable, smooth, TTR, zoo ni de otros
# paquetes de pronostico. Tampoco usa HoltWinters(), filter(),
# decompose() ni stl(). Solo R base.
#
# Cada funcion:
#   - Recibe un vector numerico y, sin NA, ordenado en el tiempo.
#   - Devuelve una lista con yhat, pronosticar y parametros.
#   - Valida sus entradas con stopifnot().
#   - No lee objetos globales.
# =========================================================================


# =========================================================================
# 3(a) METODOS DE NIVEL
# =========================================================================

# -------------------------------------------------------------------------
# ajustar_media(y)
#
# Media simple con actualizacion recursiva.
#
# Convencion (tabla del literal 3a y Definicion 10 de la Clase 3):
#   - Inicializacion: \hat{Y}_2 = Y_1.
#   - Dentro de muestra: \hat{Y}_{t+1} = t^{-1} \sum_{i<=t} Y_i,
#     version recursiva, no la media de toda la muestra.
#   - Extramuestral: \hat{Y}_{T+h} = \bar{Y}_T para todo h.
#   - Calentamiento: yhat[1] = NA.
# -------------------------------------------------------------------------
ajustar_media <- function(y) {
  stopifnot(
    is.numeric(y),
    length(y) >= 1,
    !anyNA(y)
  )
  
  T <- length(y)
  yhat <- rep(NA_real_, T)
  
  if (T >= 2) {
    suma <- 0
    for (t in 2:T) {
      suma <- suma + y[t - 1]
      yhat[t] <- suma / (t - 1)
    }
  }
  
  media_final <- sum(y) / T
  
  pronosticar <- function(h) {
    stopifnot(
      is.numeric(h), length(h) == 1,
      h >= 1, h == as.integer(h)
    )
    rep(media_final, h)
  }
  
  parametros <- list(
    media_final = media_final,
    T           = T,
    suma_total  = sum(y)
  )
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = parametros
  )
}


# -------------------------------------------------------------------------
# ajustar_mm(y, k)
#
# Media movil de orden k.
#
# Convencion (tabla del literal 3a y Definicion 11 de la Clase 3):
#   - MM_t(k) = k^{-1} \sum_{i=0}^{k-1} Y_{t-i}, para t >= k.
#   - Dentro de muestra: \hat{Y}_{t+1} = MM_t(k).
#   - Extramuestral: \hat{Y}_{T+h} = MM_T(k) para todo h.
#   - Calentamiento: primeros k valores NA (primer pronostico en yhat[k+1]).
# -------------------------------------------------------------------------
ajustar_mm <- function(y, k) {
  stopifnot(
    is.numeric(y),
    length(y) >= 1,
    !anyNA(y),
    is.numeric(k),
    length(k) == 1,
    k >= 2,
    k == as.integer(k),
    k <= length(y)
  )
  
  T  <- length(y)
  cs <- cumsum(y)
  yhat <- rep(NA_real_, T)
  
  if (T > k) {
    for (t in k:(T - 1)) {
      suma <- cs[t] - if (t - k >= 1) cs[t - k] else 0
      yhat[t + 1] <- suma / k
    }
  }
  
  suma_final <- cs[T] - if (T - k >= 1) cs[T - k] else 0
  mm_final   <- suma_final / k
  
  pronosticar <- function(h) {
    stopifnot(
      is.numeric(h), length(h) == 1,
      h >= 1, h == as.integer(h)
    )
    rep(mm_final, h)
  }
  
  parametros <- list(
    k        = k,
    mm_final = mm_final,
    T        = T
  )
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = parametros
  )
}


# -------------------------------------------------------------------------
# ajustar_ses(y, alpha)
#
# Suavizamiento exponencial simple.
#
# Convencion (tabla del literal 3a y Definicion 12 de la Clase 3):
#   - Inicializacion: \hat{Y}_2 = Y_1.
#   - Correccion de error: \hat{Y}_{t+1} = \hat{Y}_t + alpha * e_t,
#     con e_t = Y_t - \hat{Y}_t. Equivale a
#     \hat{Y}_{t+1} = alpha Y_t + (1 - alpha) \hat{Y}_t.
#   - Extramuestral: pronostico plano
#     \hat{Y}_{T+h} = alpha Y_T + (1 - alpha) \hat{Y}_T.
#   - Calentamiento: yhat[1] = NA, yhat[2] = Y_1.
#   - Verifica numericamente la equivalencia de las dos formas.
# -------------------------------------------------------------------------
ajustar_ses <- function(y, alpha) {
  stopifnot(
    is.numeric(y),
    length(y) >= 1,
    !anyNA(y),
    is.numeric(alpha),
    length(alpha) == 1,
    alpha > 0,
    alpha < 1
  )
  
  T <- length(y)
  yhat <- rep(NA_real_, T)
  
  if (T >= 2) {
    yhat[2] <- y[1]
    if (T >= 3) {
      for (t in 2:(T - 1)) {
        yhat[t + 1] <- yhat[t] + alpha * (y[t] - yhat[t])
      }
    }
  }
  
  f_final <- if (T == 1) {
    y[1]
  } else {
    yhat[T] + alpha * (y[T] - yhat[T])
  }
  
  # Verificacion: correccion de error vs. promedio ponderado
  if (T >= 2) {
    yhat_pond <- rep(NA_real_, T)
    yhat_pond[2] <- y[1]
    if (T >= 3) {
      for (t in 2:(T - 1)) {
        yhat_pond[t + 1] <- alpha * y[t] + (1 - alpha) * yhat_pond[t]
      }
    }
    dif <- max(abs(yhat - yhat_pond), na.rm = TRUE)
    if (dif >= 1e-10) {
      stop("La forma de correccion de error no coincide con la forma ",
           "de promedio ponderado (diferencia = ", dif, ").")
    }
  }
  
  pronosticar <- function(h) {
    stopifnot(
      is.numeric(h), length(h) == 1,
      h >= 1, h == as.integer(h)
    )
    rep(f_final, h)
  }
  
  parametros <- list(
    alpha       = alpha,
    nivel_final = f_final,
    T           = T
  )
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = parametros
  )
}


# -------------------------------------------------------------------------
# ajustar_dmm(y, k)
#
# Doble media movil de orden k.
#
# Convencion (tabla del literal 3a):
#   - MM_t(k) = media movil simple de orden k.
#   - MM2_t(k) = media movil de orden k aplicada a MM_t(k).
#   - \hat{E}_t = 2 MM_t(k) - MM2_t(k).
#   - \hat{\beta}_1(t) = 2/(k-1) [MM_t(k) - MM2_t(k)].
#   - \hat{Y}_{t+1} = \hat{E}_t + \hat{\beta}_1(t), para t >= 2k-1.
#   - Extramuestral: \hat{Y}_{T+h} = \hat{E}_T + \hat{\beta}_1(T) h.
#   - Calentamiento: primeros 2k-1 valores NA (primer pronostico en
#     yhat[2k]).
# -------------------------------------------------------------------------
ajustar_dmm <- function(y, k) {
  stopifnot(
    is.numeric(y),
    length(y) >= 1,
    !anyNA(y),
    is.numeric(k),
    length(k) == 1,
    k >= 2,
    k == as.integer(k),
    length(y) >= 2 * k - 1
  )
  
  T  <- length(y)
  cs <- cumsum(y)
  
  # MM_t(k)
  mm <- rep(NA_real_, T)
  for (t in k:T) {
    suma   <- cs[t] - if (t - k >= 1) cs[t - k] else 0
    mm[t]  <- suma / k
  }
  
  # Cumsum auxiliar para MM2_t(k), tratando NA como 0
  mm_cs <- cumsum(ifelse(is.na(mm), 0, mm))
  
  mm2   <- rep(NA_real_, T)
  E     <- rep(NA_real_, T)
  beta1 <- rep(NA_real_, T)
  
  for (t in (2 * k - 1):T) {
    suma_mm  <- mm_cs[t] - if (t - k >= 1) mm_cs[t - k] else 0
    mm2[t]   <- suma_mm / k
    E[t]     <- 2 * mm[t] - mm2[t]
    beta1[t] <- 2 / (k - 1) * (mm[t] - mm2[t])
  }
  
  yhat <- rep(NA_real_, T)
  if (T > 2 * k - 1) {
    for (t in (2 * k - 1):(T - 1)) {
      yhat[t + 1] <- E[t] + beta1[t]
    }
  }
  
  E_T     <- E[T]
  beta1_T <- beta1[T]
  
  pronosticar <- function(h) {
    stopifnot(
      is.numeric(h), length(h) == 1,
      h >= 1, h == as.integer(h)
    )
    E_T + beta1_T * seq_len(h)
  }
  
  parametros <- list(
    k       = k,
    E       = E,
    beta1   = beta1,
    E_T     = E_T,
    beta1_T = beta1_T,
    T       = T
  )
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = parametros
  )
}


# =========================================================================
# 3(b) TENDENCIAS POR MINIMOS CUADRADOS
# =========================================================================

# -------------------------------------------------------------------------
# ajustar_tendencia(y, tipo, corregir_sesgo)
#
# Tendencias lineal, cuadratica y exponencial por ecuaciones normales:
#   beta_hat = solve(crossprod(X), crossprod(X, y)).
#
# Convencion (tabla del literal 3a):
#   - Calentamiento: ninguno (yhat son los valores ajustados del modelo).
#   - Lineal:      Y_t = b0 + b1 t.
#   - Cuadratica:  Y_t = b0 + b1 t + b2 t^2.
#   - Exponencial: Y_t = b0 b1^t, estimada sobre log(Y_t):
#                  \hat{Y}_{T+h} = exp(a + theta (T+h)).
#     Documenta que el pronostico exp(a + theta t) estima la MEDIANA
#     condicional y no la media. Con corregir_sesgo = TRUE multiplica
#     por exp(sigma2 / 2), que es la media condicional bajo normalidad
#     de los errores en log.
#
# Devuelve en parametros la tabla de coeficientes con:
#   estimacion, se_ordinario, se_robusto (HAC Bartlett),
#   t_ordinario, p_ordinario, t_robusto, p_robusto,
#   junto con R2, sigma2, durbin_watson, L, p, n.
# -------------------------------------------------------------------------
ajustar_tendencia <- function(y,
                              tipo = c("lineal", "cuadratica", "exponencial"),
                              corregir_sesgo = FALSE) {
  tipo <- match.arg(tipo)
  if (is.null(corregir_sesgo)) corregir_sesgo <- FALSE
  
  stopifnot(
    is.numeric(y),
    length(y) >= 1,
    !anyNA(y),
    is.logical(corregir_sesgo),
    length(corregir_sesgo) == 1
  )
  
  n     <- length(y)
  t_idx <- seq_len(n)
  
  if (tipo == "lineal") {
    if (n < 2) stop("La tendencia lineal requiere al menos 2 observaciones.")
    X <- cbind(1, t_idx)
    y_ajuste <- y
    p <- 2L
    nombres <- c("beta0", "beta1")
  } else if (tipo == "cuadratica") {
    if (n < 3) stop("La tendencia cuadratica requiere al menos 3 observaciones.")
    X <- cbind(1, t_idx, t_idx^2)
    y_ajuste <- y
    p <- 3L
    nombres <- c("beta0", "beta1", "beta2")
  } else {
    if (any(y <= 0)) {
      stop("La tendencia exponencial requiere valores estrictamente positivos; ",
           "la serie tiene valores no positivos.")
    }
    if (n < 2) stop("La tendencia exponencial requiere al menos 2 observaciones.")
    X <- cbind(1, t_idx)
    y_ajuste <- log(y)
    p <- 2L
    nombres <- c("a", "theta")
  }
  
  XtX      <- crossprod(X)
  Xty      <- crossprod(X, y_ajuste)
  beta_hat <- solve(XtX, Xty)
  
  y_hat_ajuste <- as.numeric(X %*% beta_hat)
  e            <- y_ajuste - y_hat_ajuste
  gl           <- n - p
  sigma2       <- sum(e^2) / gl
  
  XtX_inv <- solve(XtX)
  se_ord  <- sqrt(diag(sigma2 * XtX_inv))
  
  # HAC Newey-West con nucleo de Bartlett, L = floor(4 (n/100)^(2/9))
  L <- floor(4 * (n / 100)^(2/9))
  meat <- crossprod(X * e)
  if (L >= 1) {
    for (h in 1:L) {
      w_h    <- 1 - h / (L + 1)
      idx_t  <- (h + 1):n
      idx_s  <- 1:(n - h)
      x_t    <- X[idx_t, , drop = FALSE]
      x_s    <- X[idx_s, , drop = FALSE]
      Gamma_h <- crossprod(x_t * e[idx_t], x_s * e[idx_s])
      meat   <- meat + w_h * (Gamma_h + t(Gamma_h))
    }
  }
  se_rob <- sqrt(diag(XtX_inv %*% meat %*% XtX_inv))
  
  t_ord <- as.numeric(beta_hat) / se_ord
  p_ord <- 2 * pt(-abs(t_ord), df = gl)
  t_rob <- as.numeric(beta_hat) / se_rob
  p_rob <- 2 * pt(-abs(t_rob), df = gl)
  
  tabla <- data.frame(
    parametro    = nombres,
    estimacion   = as.numeric(beta_hat),
    se_ordinario = as.numeric(se_ord),
    se_robusto   = as.numeric(se_rob),
    t_ordinario  = t_ord,
    p_ordinario  = p_ord,
    t_robusto    = t_rob,
    p_robusto    = p_rob,
    row.names    = NULL
  )
  
  R2 <- 1 - sum(e^2) / sum((y_ajuste - mean(y_ajuste))^2)
  dw <- sum(diff(e)^2) / sum(e^2)
  
  if (tipo == "exponencial") {
    yhat <- exp(y_hat_ajuste)
    if (corregir_sesgo) yhat <- yhat * exp(sigma2 / 2)
  } else {
    yhat <- y_hat_ajuste
  }
  
  if (tipo == "lineal") {
    b0 <- beta_hat[1]; b1 <- beta_hat[2]
    pronosticar <- function(h) {
      stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == as.integer(h))
      b0 + b1 * (n + seq_len(h))
    }
  } else if (tipo == "cuadratica") {
    b0 <- beta_hat[1]; b1 <- beta_hat[2]; b2 <- beta_hat[3]
    pronosticar <- function(h) {
      stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == as.integer(h))
      tt <- n + seq_len(h)
      b0 + b1 * tt + b2 * tt^2
    }
  } else {
    a     <- beta_hat[1]
    theta <- beta_hat[2]
    factor <- if (corregir_sesgo) exp(sigma2 / 2) else 1
    pronosticar <- function(h) {
      stopifnot(is.numeric(h), length(h) == 1, h >= 1, h == as.integer(h))
      factor * exp(a + theta * (n + seq_len(h)))
    }
  }
  
  parametros <- list(
    tipo           = tipo,
    coeficientes   = tabla,
    R2             = R2,
    sigma2         = sigma2,
    durbin_watson  = dw,
    L              = L,
    p              = p,
    n              = n,
    corregir_sesgo = corregir_sesgo
  )
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = parametros
  )
}


# =========================================================================
# 3(c) HOLT LINEAL
# =========================================================================

# -------------------------------------------------------------------------
# ajustar_holt(y, alpha, beta)
#
# Holt lineal.
#
# Convencion (tabla del literal 3a y notas de la Clase 3):
#   - Inicializacion: L_1 = Y_1, \hat{T}_1 = 0.
#   - Nivel:    L_t = alpha Y_t + (1 - alpha) (L_{t-1} + \hat{T}_{t-1}).
#   - Pendiente:\hat{T}_t = beta (L_t - L_{t-1}) + (1 - beta) \hat{T}_{t-1}.
#   - Dentro de muestra: \hat{Y}_{t+1} = L_t + \hat{T}_t.
#   - Extramuestral: \hat{Y}_{T+h} = L_T + \hat{T}_T h.
#   - Calentamiento: yhat[1] = NA; yhat[2] = L_1 + \hat{T}_1 = Y_1.
#   - Verifica numericamente la forma de correccion de error:
#       L_t = L_{t-1} + \hat{T}_{t-1} + alpha e_t,
#       \hat{T}_t = \hat{T}_{t-1} + alpha beta e_t,
#     con e_t = Y_t - \hat{Y}_t.
#   - Verifica numericamente la forma de correccion de error con una
#     tolerancia de 1e-10, adecuada para acumulacion de punto flotante
#     en recursiones largas.
# -------------------------------------------------------------------------
ajustar_holt <- function(y, alpha, beta) {
  stopifnot(
    is.numeric(y),
    length(y) >= 2,
    !anyNA(y),
    is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1,
    is.numeric(beta),  length(beta)  == 1, beta  > 0, beta  < 1
  )
  
  T <- length(y)
  L     <- rep(NA_real_, T)
  Trend <- rep(NA_real_, T)   # no se llama T para no chocar con length(y)
  yhat  <- rep(NA_real_, T)
  
  L[1]     <- y[1]
  Trend[1] <- 0
  
  if (T >= 2) yhat[2] <- L[1] + Trend[1]
  
  if (T >= 3) {
    for (t in 2:(T - 1)) {
      L[t]     <- alpha * y[t] + (1 - alpha) * (L[t - 1] + Trend[t - 1])
      Trend[t] <- beta  * (L[t] - L[t - 1]) + (1 - beta) * Trend[t - 1]
      yhat[t + 1] <- L[t] + Trend[t]
    }
  }
  
  if (T >= 2) {
    L[T]     <- alpha * y[T] + (1 - alpha) * (L[T - 1] + Trend[T - 1])
    Trend[T] <- beta  * (L[T] - L[T - 1]) + (1 - beta) * Trend[T - 1]
  }
  
  # --- Verificacion: forma de correccion de error ----------------------
  L_ce     <- rep(NA_real_, T)
  T_ce     <- rep(NA_real_, T)
  L_ce[1]  <- y[1]
  T_ce[1]  <- 0
  
  if (T >= 3) {
    for (t in 2:(T - 1)) {
      e_t      <- y[t] - (L_ce[t - 1] + T_ce[t - 1])
      L_ce[t]  <- L_ce[t - 1] + T_ce[t - 1] + alpha * e_t
      T_ce[t]  <- T_ce[t - 1] + alpha * beta * e_t
    }
  }
  if (T >= 2) {
    e_T      <- y[T] - (L_ce[T - 1] + T_ce[T - 1])
    L_ce[T]  <- L_ce[T - 1] + T_ce[T - 1] + alpha * e_T
    T_ce[T]  <- T_ce[T - 1] + alpha * beta * e_T
  }
  
  dif_L <- max(abs(L     - L_ce), na.rm = TRUE)
  dif_T <- max(abs(Trend - T_ce), na.rm = TRUE)
  if (dif_L >= 1e-8 || dif_T >= 1e-8) {
    stop("La forma de correccion de error no coincide con la forma estandar ",
         "(dif L = ", dif_L, ", dif T = ", dif_T, ").")
  }
  
  L_T     <- L[T]
  Trend_T <- Trend[T]
  
  pronosticar <- function(h) {
    stopifnot(
      is.numeric(h), length(h) == 1,
      h >= 1, h == as.integer(h)
    )
    L_T + Trend_T * seq_len(h)
  }
  
  parametros <- list(
    alpha        = alpha,
    beta         = beta,
    L            = L,
    T_tray       = Trend,
    L_T          = L_T,
    T_T          = Trend_T,
    verificacion = list(dif_L = dif_L, dif_T = dif_T),
    T_obs        = T
  )
  
  list(
    yhat        = yhat,
    pronosticar = pronosticar,
    parametros  = parametros
  )
}


# =========================================================================
# 3(d) OPTIMIZACION DE CONSTANTES Y VENTANAS
# =========================================================================

# -------------------------------------------------------------------------
# optimizar(y, metodo, rejilla)
#
# Evaluates MSE de un paso dentro del tramo de estimacion sobre una rejilla
# declarada. Devuelve la rejilla completa con su MSE y el optimo.
#
# y:       vector numerico (tramo de ESTIMACION).
# metodo:  "mm", "dmm", "ses" o "holt".
# rejilla: data.frame con columnas:
#            - k              para mm y dmm
#            - alpha          para ses
#            - alpha, beta    para holt
# -------------------------------------------------------------------------
optimizar <- function(y, metodo, rejilla) {
  stopifnot(
    is.numeric(y),
    length(y) >= 2,
    !anyNA(y),
    is.character(metodo),
    length(metodo) == 1
  )
  
  metodo <- match.arg(metodo, c("mm", "dmm", "ses", "holt"))
  
  if (!is.data.frame(rejilla)) rejilla <- as.data.frame(rejilla)
  if (nrow(rejilla) < 1) stop("La rejilla esta vacia.")
  
  if (metodo %in% c("mm", "dmm")) {
    if (!"k" %in% names(rejilla))
      stop("La rejilla para '", metodo, "' debe tener columna 'k'.")
    if (any(rejilla$k < 2) || any(rejilla$k != as.integer(rejilla$k)))
      stop("Los valores de k deben ser enteros >= 2.")
    if (any(rejilla$k > length(y)))
      stop("Hay valores de k mayores que la longitud de la serie.")
  } else if (metodo == "ses") {
    if (!"alpha" %in% names(rejilla))
      stop("La rejilla para 'ses' debe tener columna 'alpha'.")
    if (any(rejilla$alpha <= 0) || any(rejilla$alpha >= 1))
      stop("Los valores de alpha deben estar estrictamente en (0, 1).")
  } else {
    if (!all(c("alpha", "beta") %in% names(rejilla)))
      stop("La rejilla para 'holt' debe tener columnas 'alpha' y 'beta'.")
    if (any(rejilla$alpha <= 0) || any(rejilla$alpha >= 1) ||
        any(rejilla$beta  <= 0) || any(rejilla$beta  >= 1))
      stop("Los valores de alpha y beta deben estar estrictamente en (0, 1).")
  }
  
  mse_uno_paso <- function(ajuste, y) {
    yhat    <- ajuste$yhat
    validos <- !is.na(yhat)
    if (sum(validos) == 0) return(NA_real_)
    mean((y[validos] - yhat[validos])^2)
  }
  
  mse <- numeric(nrow(rejilla))
  for (i in seq_len(nrow(rejilla))) {
    ajuste <- switch(metodo,
                     "mm"   = ajustar_mm(y,  k = rejilla$k[i]),
                     "dmm"  = ajustar_dmm(y, k = rejilla$k[i]),
                     "ses"  = ajustar_ses(y, alpha = rejilla$alpha[i]),
                     "holt" = ajustar_holt(y, alpha = rejilla$alpha[i],
                                           beta  = rejilla$beta[i])
    )
    mse[i] <- mse_uno_paso(ajuste, y)
  }
  
  resultado         <- rejilla
  resultado$mse     <- mse
  resultado$optimo  <- FALSE
  
  idx_opt <- which.min(mse)
  resultado$optimo[idx_opt] <- TRUE
  optimo <- resultado[idx_opt, , drop = FALSE]
  
  tol <- 1e-10
  en_borde <- FALSE
  if (metodo %in% c("mm", "dmm")) {
    k_min <- min(rejilla$k); k_max <- max(rejilla$k)
    en_borde <- (optimo$k == k_min) || (optimo$k == k_max)
  } else if (metodo == "ses") {
    a_min <- min(rejilla$alpha); a_max <- max(rejilla$alpha)
    en_borde <- (abs(optimo$alpha - a_min) < tol) ||
      (abs(optimo$alpha - a_max) < tol)
  } else {
    a_min <- min(rejilla$alpha); a_max <- max(rejilla$alpha)
    b_min <- min(rejilla$beta);  b_max <- max(rejilla$beta)
    en_borde <- (abs(optimo$alpha - a_min) < tol) ||
      (abs(optimo$alpha - a_max) < tol) ||
      (abs(optimo$beta  - b_min) < tol) ||
      (abs(optimo$beta  - b_max) < tol)
  }
  
  list(
    rejilla        = resultado,
    optimo         = optimo,
    metodo         = metodo,
    n_evaluaciones = nrow(rejilla),
    en_borde       = en_borde
  )
}