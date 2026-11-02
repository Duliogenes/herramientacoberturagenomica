
# ============================================================
# HERRAMIENTA DE COBERTURA GENOMICA
# MODULO 03: COMPARACION DE COBERTURA POR REGION
# KOLMOGOROV-SMIRNOV + MEDIA + FDR
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. FUNCION DE COMPARACION KS
# ------------------------------------------------------------

comparar_grupos_KS <- function(x, y, minimo = 10L) {

  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  n_x <- length(x)
  n_y <- length(y)

  media_x <- if (n_x > 0) mean(x) else NA_real_
  media_y <- if (n_y > 0) mean(y) else NA_real_

  diferencia <- media_y - media_x

  if (n_x < minimo || n_y < minimo) {
    return(list(
      n_menores_60 = n_x,
      n_60_o_mas = n_y,
      media_menores_60 = media_x,
      media_60_o_mas = media_y,
      diferencia_media = diferencia,
      P = NA_real_
    ))
  }

  p_valor <- tryCatch(
    suppressWarnings(
      ks.test(
        x,
        y,
        alternative = "two.sided",
        exact = FALSE
      )$p.value
    ),
    error = function(e) NA_real_
  )

  list(
    n_menores_60 = n_x,
    n_60_o_mas = n_y,
    media_menores_60 = media_x,
    media_60_o_mas = media_y,
    diferencia_media = diferencia,
    P = p_valor
  )
}


# ------------------------------------------------------------
# 2. COMPARACION POR REGION SUBTELOMERICA
# ------------------------------------------------------------

analizar_cobertura_por_region <- function(
    archivo_subset,
    claves,
    tabla_analisis,
    tamano_bloque = 1000L,
    minimo = 10L
) {

  encabezado <- names(
    fread(
      archivo_subset,
      nrows = 0,
      check.names = FALSE
    )
  )

  ids_muestras <- as.character(
    tabla_analisis$id_muestra
  )

  indices_muestras <- match(
    ids_muestras,
    encabezado
  )

  if (anyNA(indices_muestras)) {
    stop("Hay muestras que no aparecen en la matriz.")
  }

  indices_columnas <- c(1:6, indices_muestras)

  nombres_seleccionados <- encabezado[
    indices_columnas
  ]

  meta_regiones <- copy(
    as.data.table(claves)
  )

  setorder(meta_regiones, fila_matriz)

  info_muestras <- as.data.table(
    copy(tabla_analisis)
  )

  info_muestras <- info_muestras[
    match(ids_muestras, id_muestra)
  ]

  grupo_joven <- which(
    info_muestras$grupo_edad == "<60 años"
  )

  grupo_mayor <- which(
    info_muestras$grupo_edad == "60 años o más"
  )

  if (length(grupo_joven) == 0 ||
      length(grupo_mayor) == 0) {
    stop("Uno de los grupos de edad está vacío.")
  }

  n_regiones <- nrow(meta_regiones)

  n_bloques <- ceiling(
    n_regiones / tamano_bloque
  )

  resultados_lista <- vector(
    "list",
    n_bloques
  )

  cat("\n===== COMPARACION POR REGION =====\n")
  cat("Prueba: Kolmogorov-Smirnov\n")
  cat("Descriptivo: media\n")
  cat("Regiones:", n_regiones, "\n")

  for (b in seq_len(n_bloques)) {

    inicio <- (b - 1L) * tamano_bloque

    cantidad <- min(
      tamano_bloque,
      n_regiones - inicio
    )

    bloque <- fread(
      archivo_subset,
      skip = inicio + 1L,
      nrows = cantidad,
      header = FALSE,
      col.names = nombres_seleccionados,
      select = indices_columnas,
      check.names = FALSE,
      showProgress = FALSE
    )

    if (nrow(bloque) != cantidad) {
      stop("Número inesperado de filas en bloque ", b)
    }

    matriz_bloque <- as.matrix(
      bloque[, -(1:6), with = FALSE]
    )

    storage.mode(matriz_bloque) <- "numeric"

    resultados_bloque <- lapply(
      seq_len(nrow(matriz_bloque)),
      function(i) {

        comparar_grupos_KS(
          matriz_bloque[i, grupo_joven],
          matriz_bloque[i, grupo_mayor],
          minimo
        )
      }
    )

    resumen_bloque <- rbindlist(
      lapply(
        resultados_bloque,
        as.data.table
      )
    )

    meta_bloque <- meta_regiones[
      seq.int(
        inicio + 1L,
        inicio + cantidad
      ),
      .(
        chromosome,
        start,
        end,
        gene_name,
        exon_rank,
        transcript_id,
        extremo
      )
    ]

    resultados_lista[[b]] <- cbind(
      meta_bloque,
      resumen_bloque
    )

    cat(
      "Bloque", b, "de", n_bloques,
      "| Regiones:", inicio + cantidad,
      "de", n_regiones, "\n"
    )

    rm(
      bloque,
      matriz_bloque,
      resultados_bloque,
      resumen_bloque
    )

    gc(verbose = FALSE)
  }

  resultados_regiones <- rbindlist(
    resultados_lista,
    fill = TRUE
  )

  # Corrección por comparaciones múltiples
  resultados_regiones[
    ,
    FDR := p.adjust(P, method = "BH")
  ]

  setorder(
    resultados_regiones,
    FDR,
    P,
    na.last = TRUE
  )

  cat("\nAnálisis regional terminado.\n")
  cat("Regiones analizadas:", nrow(resultados_regiones), "\n")
  cat(
    "Regiones con FDR < 0.05:",
    sum(resultados_regiones$FDR < 0.05, na.rm = TRUE),
    "\n"
  )

  return(resultados_regiones)
}
