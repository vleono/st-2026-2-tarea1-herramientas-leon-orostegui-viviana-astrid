https://github.com/vleono/st-2026-2-tarea1-herramientas-leon-orostegui-viviana-astrid

# Tarea 1 - Herramientas de pronostico

Series de Tiempo Univariadas, codigo 3009297. Semestre 2026-II.
Universidad Nacional de Colombia, sede Medellin.
Estudiante: Viviana Astrid León Orostegui.

## Que contiene el repositorio

| Archivo | Contenido | Dependencias |
|---------|-----------|--------------|
| `R/00-lectura.R` | `leer_serie()` | R base, tibble |
| `R/01-graficos.R` | `acf_muestral()`, `graficar_serie()`, `correlograma()`, `verificar_acf()` | R base, ggplot2, patchwork |
| `R/02-metodos.R` | Los ocho metodos y `optimizar()` | R base |
| `R/03-evaluacion.R` | Medidas, pruebas con seis elementos, referente ingenuo, MASE | R base |
| `ejemplos/ejemplos.R` | Corre los ocho ejemplos y el contraejemplo | ggplot2, patchwork, tibble |
| `informe/informe.qmd` | Analisis detallado | Quarto |
| `informe/informe.html` | Informe renderizado | - |
| `figs/` | Figuras generadas | - |
| `sesion-info.txt` | Salida de sessionInfo() | - |

## Como se corre

Desde la raiz del repositorio:

```r
source("ejemplos/ejemplos.R")
```

Tiempo de ejecucion observado: 45.36263 segundos en un equipo con R 4.5.1 sobre Windows.

Para renderizar el informe:

```bash
quarto render informe/informe.qmd
```

## Como se usan las funciones

Ejemplo minimo con Holt lineal:

```r
source("R/02-metodos.R")


source("R/02-metodos.R")

# Serie de prueba: los primeros 30 meses de AirPassengers
y   <- as.numeric(AirPassengers)[1:30]

# Ajuste con Holt lineal
res <- ajustar_holt(y, alpha = 0.3, beta = 0.1)

# Valores ajustados (con NA en el calentamiento)
head(res$yhat)
# [1] NA 112.0000 113.9800 ...

# Pronóstico extramuestral a 5 pasos
res$pronosticar(5)

# Parámetros y estados finales
res$parametros$alpha
res$parametros$beta
res$parametros$L_T
res$parametros$T_T

# Verificación numérica de la forma de corrección de error
res$parametros$verificacion
# $dif_L  [1] 1.70e-13
# $dif_T  [1] 1.64e-14
```

## Convenciones que fijan los numeros

- Media simple: inicializacion Y_2 = Y_1; recursiva; calentamiento 1 periodo.
- Media movil: MM_t(k) para t >= k; calentamiento k periodos.
- Doble media movil: calentamiento 2k-1 periodos.
- Suavizamiento exponencial simple: Y_2 = Y_1; calentamiento 1 periodo.
- Tendencias por mínimos cuadrados: ecuaciones normales; sin calentamiento.
- Holt lineal: L_1 = Y_1, T_1 = 0; calentamiento 1 periodo.
- ACF con divisor unico T segun Definicion 2.5. de las notas de Clase 3.
- Banda del correlograma: la de plot.acf(), clim0 = qnorm((1+ci)/2)/sqrt(n.used) con ci = 0.95 y n.used el número de observaciones de la sucesión graficada.
- MASE: MAD de los errores en validación dividido por el MAD del ingenuo dentro del tramo de estimación.
## Resumen de resultados

| Ejemplo | Serie | Parametros | MASE metodo | MASE ingenuo |
|---------|-------|------------|-------------|---------------|
| 1 | sunspot.year | - | 2.87 | 3.88 |
| 2 | discoveries | k = 9 | 0.59 | 1.14 |
| 3 | Nile | alpha = 0.24 | 0.81 | 0.83 |
| 4 | discoveries | k = 9 | 0.58 | 1.14 |
| 5 | airmiles | tendencia lineal | 6.39 | 4.50 |
| 6 | co2 | tendencia cuadratica | 1.61 | 0.89 |
| 7 | JohnsonJohnson | tendencia exponencial | 3.59 | 6.53 |
| 8 | airmiles | alpha = 0.95, beta = 0.05 | 1.24 | 4.50 |
| 4(e) | AirPassengers (contraejemplo) | alpha = 0.95, beta = 0.05 | 2.43 | 1.57 |

## Declaracion de uso de IA

Se utilizo un asistente de inteligencia artificial DeepSeek durante el desarrollo de esta tarea.

**Lo que se pidio.** Apoyo para redactar encabezados de funciones, plantillas del README y estructura del informe.

**Lo que se recibio.** Borradores de codigo R, fragmentos de texto y sugerencias de estructura.

**Lo que se verifico por cuenta propia.** Se verifico la coincidencia de la ACF a mano con acf() hasta 1e-12 (diferencia maxima: 5.55e-16). Se verifico Qm de Ljung-Box contra Box.test() (diferencia: 6.82e-13). Se verifico la forma de correccion de error de Holt y SES. Todo el codigo se entiende y se puede explicar en la sustentacion oral.
