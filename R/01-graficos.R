# -------------------------------------------------------------------------
# correlograma(datos, m)
#
# datos: tibble con columna y (salida de leer_serie()) O vector numerico
#        (por ejemplo, los errores de un metodo).
# m:     numero de rezagos. Por defecto min(floor(T/4), 24).
#
# ACF a mano (divisor T) arriba, PACF con pacf() abajo.
# Banda de ruido blanco: +-qnorm(0.975)/sqrt(n) con n el numero de
# observaciones de la sucesion graficada.
#
# Devuelve un patchwork con dos paneles.
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# acf_muestral(y, lag.max)
# -------------------------------------------------------------------------
acf_muestral <- function(y, lag.max = NULL) {
  stopifnot(is.numeric(y), length(y) >= 5, !anyNA(y))
  T <- length(y)
  if (is.null(lag.max)) lag.max <- min(floor(T / 4), 24)
  if (lag.max >= T) stop("lag.max debe ser menor que T.")
  if (lag.max < 1)  stop("lag.max debe ser al menos 1.")
  
  ybar  <- mean(y)
  denom <- sum((y - ybar)^2)
  r     <- numeric(lag.max)
  for (h in seq_len(lag.max)) {
    num   <- sum((y[(h + 1):T] - ybar) * (y[1:(T - h)] - ybar))
    r[h]  <- num / denom
  }
  r
}

# -------------------------------------------------------------------------
# graficar_serie(datos, titulo)
# -------------------------------------------------------------------------
graficar_serie <- function(datos, titulo = "Serie de tiempo") {
  stopifnot(
    inherits(datos, "data.frame"),
    all(c("t", "fecha", "y") %in% names(datos))
  )
  frec   <- attr(datos, "frecuencia")
  fuente <- attr(datos, "fuente")
  unidad <- attr(datos, "unidad")
  n      <- nrow(datos)
  
  ggplot2::ggplot(datos, ggplot2::aes(x = fecha, y = y)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.5) +
    ggplot2::labs(
      title   = titulo,
      x       = "Fecha",
      y       = paste0("Y (", unidad, ")"),
      caption = paste0("Fuente: ", fuente, " | n = ", n,
                       " obs | frecuencia = ", frec)
    ) +
    ggplot2::theme_minimal()
}
# -------------------------------------------------------------------------
correlograma <- function(datos, m = NULL) {
  
  # ---------- Aceptar tibble o vector numerico ------------------------
  if (inherits(datos, "data.frame")) {
    if (!"y" %in% names(datos)) {
      stop("Si 'datos' es un data.frame, debe tener una columna 'y'.")
    }
    y <- datos$y
  } else if (is.numeric(datos)) {
    y <- datos
  } else {
    stop("'datos' debe ser un tibble con columna 'y' o un vector numerico.")
  }
  
  stopifnot(
    length(y) >= 5,
    !anyNA(y)
  )
  
  n.used <- length(y) 
  if (is.null(m)) m <- min(floor(T / 4), 24)
  if (m >= n.used) stop("m debe ser menor que T.")
  if (m < 1)  stop("m debe ser al menos 1.")
  
  r <- acf_muestral(y, lag.max = m)
  
  ci    <- 0.95
  clim0 <- qnorm((1 + ci) / 2) / sqrt(n.used)
  
  df_acf <- data.frame(h = seq_len(m), r = r)
  
  pp <- stats::pacf(y, lag.max = m, plot = FALSE)
  df_pacf <- data.frame(h = seq_len(m), r = as.numeric(pp$acf))
  
  g_acf <- ggplot2::ggplot(df_acf, ggplot2::aes(x = h, y = r)) +
    ggplot2::geom_segment(ggplot2::aes(xend = h, yend = 0),
                          linewidth = 0.6, color = "steelblue") +
    ggplot2::geom_hline(yintercept = c(-clim0, clim0),
                        linetype = "dashed", color = "red") +
    ggplot2::geom_hline(yintercept = 0, color = "grey40") +
    ggplot2::labs(title = "ACF",
                  x = "Rezago h",
                  y = expression(r[h])) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_minimal()
  
  g_pacf <- ggplot2::ggplot(df_pacf, ggplot2::aes(x = h, y = r)) +
    ggplot2::geom_segment(ggplot2::aes(xend = h, yend = 0),
                          linewidth = 0.6, color = "darkorange") +
    ggplot2::geom_hline(yintercept = c(-clim0, clim0),
                        linetype = "dashed", color = "red") +
    ggplot2::geom_hline(yintercept = 0, color = "grey40") +
    ggplot2::labs(title = "PACF",
                  x = "Rezago h",
                  y = expression(phi[hh])) +
    ggplot2::ylim(-1, 1) +
    ggplot2::theme_minimal()
  
  patchwork::wrap_plots(g_acf, g_pacf, ncol = 1)
}
# -------------------------------------------------------------------------
# verificar_acf(y, lag.max)
#
# Bloque de verificacion: contraste contra acf() de stats.
# Devuelve la maxima diferencia absoluta y si cumple < 1e-12.
# -------------------------------------------------------------------------
verificar_acf <- function(y, lag.max = NULL) {
  if (is.null(lag.max)) lag.max <- min(floor(length(y) / 4), 24)
  
  r_propio <- acf_muestral(y, lag.max = lag.max)
  r_stats  <- as.numeric(stats::acf(y, lag.max = lag.max,
                                    plot = FALSE)$acf)[-1]
  dif_max  <- max(abs(r_propio - r_stats))
  
  list(
    dif_max = dif_max,
    cumple  = dif_max < 1e-12
  )
}