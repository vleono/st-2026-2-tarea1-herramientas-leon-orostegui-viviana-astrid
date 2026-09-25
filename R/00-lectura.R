# =========================================================================
# 00-lectura.R
# Lectura y estandarizacion de una serie de tiempo univariada.
# =========================================================================

# -------------------------------------------------------------------------
# leer_serie(x, fuente, unidad)
#
# x: un objeto ts UNIVARIADO, o la ruta a un csv con columnas
#    'fecha' y 'valor'.
# fuente: cadena con la fuente original tal como la documenta help().
# unidad: cadena con la unidad de la serie.
#
# Devuelve un tibble con columnas t (entero desde 1), fecha (Date), y,
# mas los atributos frecuencia, fuente y unidad.
#
# Verifica que las fechas son crecientes y equiespaciadas segun la
# frecuencia. Se detiene si la serie es multivariada o si el csv no tiene
# las columnas requeridas.
# -------------------------------------------------------------------------

leer_serie <- function(x, fuente, unidad) {
  
  stopifnot(
    is.character(fuente), length(fuente) == 1,
    is.character(unidad), length(unidad) == 1
  )
  
  # =====================================================================
  # Caso 1: objeto ts
  # =====================================================================
  if (inherits(x, "ts")) {
    
    # --- Verificacion: la serie debe ser univariada ------------------
    if (!is.null(dim(x)) && ncol(x) > 1) {
      stop("leer_serie() espera un ts univariado. ",
           "Extrae una columna antes de llamarla, por ejemplo: ",
           "serie <- Seatbelts[, \"DriversKilled\"]")
    }
    
    y    <- as.numeric(x)
    frec <- frequency(x)
    ini  <- start(x)
    n    <- length(y)
    
    if (n < 2) stop("La serie debe tener al menos 2 observaciones.")
    
    # --- Construccion de fechas segun la frecuencia ------------------
    if (frec == 12) {
      year  <- ini[1]
      month <- ini[2]
      fecha <- seq(as.Date(sprintf("%04d-%02d-01", year, month)),
                   by = "month", length.out = n)
      
    } else if (frec == 4) {
      year  <- ini[1]
      q     <- ini[2]
      mes   <- (q - 1) * 3 + 1
      fecha <- seq(as.Date(sprintf("%04d-%02d-01", year, mes)),
                   by = "3 months", length.out = n)
      
    } else if (frec == 1) {
      year  <- ini[1]
      fecha <- seq(as.Date(sprintf("%04d-01-01", year)),
                   by = "year", length.out = n)
      
    } else if (frec == 52) {
      # Frecuencia semanal: partir del inicio del ano del primer dato
      start_year <- floor(ini[1])
      start_week <- round((ini[1] - start_year) * 52) + 1
      fecha <- seq(as.Date(sprintf("%04d-01-01", start_year)) +
                     (start_week - 1) * 7,
                   by = "week", length.out = n)
      
    } else {
      # Frecuencia no estandar: usar posicion entera desde 1970-01-01
      fecha <- as.Date("1970-01-01") + (seq_len(n) - 1)
    }
    
    tt <- seq_len(n)
  }
  
  # =====================================================================
  # Caso 2: ruta a un archivo csv
  # =====================================================================
  else if (is.character(x) && length(x) == 1) {
    
    if (!file.exists(x)) stop("No existe el archivo: ", x)
    
    df <- utils::read.csv(x, stringsAsFactors = FALSE)
    
    if (!all(c("fecha", "valor") %in% names(df))) {
      stop("El csv debe tener columnas 'fecha' y 'valor'.")
    }
    
    fecha <- as.Date(df$fecha)
    y     <- as.numeric(df$valor)
    n     <- length(y)
    tt    <- seq_len(n)
    
    if (n < 2) stop("La serie debe tener al menos 2 observaciones.")
    
    # --- Inferencia de la frecuencia a partir de la separacion ------
    difs <- as.numeric(diff(fecha), units = "days")
    
    if (length(difs) == 0) {
      frec <- 1
    } else if (all(abs(difs - 1)   <= 1)) {
      frec <- 365
    } else if (all(abs(difs - 7)   <= 1)) {
      frec <- 52
    } else if (all(abs(difs - 30)  <= 3)) {
      frec <- 12
    } else if (all(abs(difs - 91)  <= 3)) {
      frec <- 4
    } else if (all(abs(difs - 365) <= 3)) {
      frec <- 1
    } else {
      frec <- 1
    }
  }
  
  # =====================================================================
  # Caso 3: tipo no soportado
  # =====================================================================
  else {
    stop("leer_serie() espera un objeto ts univariado o una ruta a un csv.")
  }
  
  # =====================================================================
  # Verificacion de fechas y valores
  # =====================================================================
  if (any(is.na(fecha))) stop("Hay fechas no validas.")
  if (any(is.na(y)))     stop("Hay valores NA en la serie.")
  
  if (is.unsorted(fecha, strictly = TRUE)) {
    stop("Las fechas no son estrictamente crecientes.")
  }
  
  # --- Verificacion de equiespaciado segun la frecuencia --------------
  if (length(fecha) >= 3) {
    if (frec == 12) {
      # Mensual: dia 1 del mes, mes avanza en 1 (o retrocede 11 en el
      # salto de diciembre a enero)
      if (any(as.integer(format(fecha, "%d")) != 1)) {
        stop("Las fechas mensuales deben caer en el dia 1.")
      }
      meses <- as.integer(format(fecha, "%m"))
      dif_m <- diff(meses)
      if (!all(dif_m == 1 | dif_m == -11)) {
        stop("Las fechas no estan equiespaciadas segun la frecuencia mensual.")
      }
      
    } else if (frec == 4) {
      # Trimestral: dia 1 del mes, mes avanza en 3 (o retrocede 9)
      if (any(as.integer(format(fecha, "%d")) != 1)) {
        stop("Las fechas trimestrales deben caer en el dia 1.")
      }
      meses <- as.integer(format(fecha, "%m"))
      dif_m <- diff(meses)
      if (!all(dif_m == 3 | dif_m == -9)) {
        stop("Las fechas no estan equiespaciadas segun la frecuencia trimestral.")
      }
      
    } else {
      # Otras frecuencias: chequeo en dias
      difs <- as.numeric(diff(fecha), units = "days")
      tol <- switch(as.character(frec),
                    "365" = 1,    # diaria
                    "52"  = 1,    # semanal
                    "1"   = 2,    # anual
                    5             # otro: tolerancia amplia
      )
      ref <- switch(as.character(frec),
                    "365" = 1,
                    "52"  = 7,
                    "1"   = 365,
                    mean(difs)
      )
      if (any(abs(difs - ref) > tol)) {
        stop("Las fechas no estan equiespaciadas segun la frecuencia.")
      }
    }
  }
  
  # =====================================================================
  # Construccion del tibble y atributos
  # =====================================================================
  out <- tibble::tibble(
    t     = tt,
    fecha = fecha,
    y     = y
  )
  
  attr(out, "frecuencia") <- frec
  attr(out, "fuente")     <- fuente
  attr(out, "unidad")     <- unidad
  
  out
}

