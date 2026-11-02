
# ============================================================
# HERRAMIENTA DE COBERTURA GENOMICA
# MODULO 05: ANALISIS DE COBERTURA POR GEN
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
# 2. FUNCION PRINCIPAL: ANALISIS POR GEN
# ------------------------------------------------------------

analizar_cobertura_por_gen <- function(
    archivo_subset,
    tabla_analisis,
    minimo = 10L
) {

  # ----------------------------------------------------------
  # 2.1 PREPARAR INFORMACION DE PARTICIPANTES
  # ----------------------------------------------------------

  info_muestras <- as.data.table(
    copy(tabla_analisis)
  )

  info_muestras <- info_muestras[
    !is.na(id_muestra) &
      grupo_edad %in% c(
        "<60 años",
        "60 años o más"
      )
  ]

  # Evitar muestras duplicadas
  if (anyDuplicated(info_muestras$id_muestra)) {
    stop("Hay identificadores de muestra duplicados.")
  }

  # ----------------------------------------------------------
  # 2.2 IDENTIFICAR COLUMNAS DE LA MATRIZ
  # ----------------------------------------------------------

  columnas_anotacion <- c(
    "chromosome",
    "start",
    "end",
    "gene_name",
    "exon_rank",
    "transcript_id"
  )

  encabezado <- names(
    fread(
      archivo_subset,
      nrows = 0,
      check.names = FALSE
    )
  )

  if (!all(columnas_anotacion %in% encabezado)) {
    stop("Faltan columnas de anotación en la matriz.")
  }

  ids_muestras <- intersect(
    info_muestras$id_muestra,
    encabezado
  )

  if (length(ids_muestras) == 0) {
    stop("No se encontraron muestras coincidentes.")
  }

  info_muestras <- info_muestras[
    match(ids_muestras, id_muestra)
  ]

  stopifnot(
    identical(
      info_muestras$id_muestra,
      ids_muestras
    )
  )

  grupo_joven <- which(
    info_muestras$grupo_edad == "<60 años"
  )

  grupo_mayor <- which(
    info_muestras$grupo_edad == "60 años o más"
  )

  cat("\n===== MUESTRAS INCLUIDAS =====\n")
  cat("Total:", length(ids_muestras), "\n")
  cat("Menores de 60:", length(grupo_joven), "\n")
  cat("60 años o más:", length(grupo_mayor), "\n")

  if (length(grupo_joven) < minimo ||
      length(grupo_mayor) < minimo) {
    stop("Uno de los grupos tiene menos muestras que el mínimo requerido.")
  }

  # ----------------------------------------------------------
  # 2.3 CARGAR ANOTACIONES Y COBERTURA
  # ----------------------------------------------------------

  datos <- fread(
    archivo_subset,
    select = c(columnas_anotacion, ids_muestras),
    check.names = FALSE,
    showProgress = TRUE
  )

  cat("\nRegiones cargadas:", nrow(datos), "\n")

  # Eliminar genes sin anotación
  datos <- datos[
    !is.na(gene_name) &
      trimws(gene_name) != "" &
      gene_name != "NOANNOTATION"
  ]

  # Eliminar anotaciones duplicadas
  datos <- unique(
    datos,
    by = c(
      "chromosome",
      "start",
      "end",
      "gene_name"
    )
  )

  genes <- unique(datos$gene_name)

  cat("Regiones anotadas:", nrow(datos), "\n")
  cat("Genes identificados:", length(genes), "\n")

  # ----------------------------------------------------------
  # 2.4 CALCULAR MEDIA DE COBERTURA POR GEN
  # ----------------------------------------------------------

  resultados_lista <- vector(
    "list",
    length(genes)
  )

  for (i in seq_along(genes)) {

    gen_actual <- genes[i]

    filas_gen <- which(
      datos$gene_name == gen_actual
    )

    # Matriz: regiones del gen x participantes
    mat <- as.matrix(
      datos[filas_gen, ..ids_muestras]
    )

    storage.mode(mat) <- "numeric"

    # Media de cobertura por participante
    cobertura_media <- colMeans(
      mat,
      na.rm = TRUE
    )

    cobertura_media[
      !is.finite(cobertura_media)
    ] <- NA_real_

    # Separar grupos
    x <- cobertura_media[grupo_joven]
    y <- cobertura_media[grupo_mayor]

    # Comparación estadística KS
    comparacion <- comparar_grupos_KS(
      x,
      y,
      minimo = minimo
    )

    # Cromosomas asociados al gen
    cromosomas <- unique(
      datos[filas_gen, chromosome]
    )

    resultados_lista[[i]] <- data.table(
      gene_name = gen_actual,
      chromosome = paste(
        cromosomas,
        collapse = ";"
      ),
      n_regiones = length(filas_gen),
      n_menores_60 = comparacion$n_menores_60,
      n_60_o_mas = comparacion$n_60_o_mas,
      media_menores_60 = comparacion$media_menores_60,
      media_60_o_mas = comparacion$media_60_o_mas,
      diferencia_media = comparacion$diferencia_media,
      P = comparacion$P
    )

    if (i %% 100 == 0 || i == length(genes)) {
      cat(
        "Genes analizados:",
        i,
        "de",
        length(genes),
        "\n"
      )
    }
  }

  # ----------------------------------------------------------
  # 2.5 CORRECCION POR COMPARACIONES MULTIPLES
  # ----------------------------------------------------------

  resultados_genes <- rbindlist(
    resultados_lista,
    fill = TRUE
  )

  resultados_genes[
    ,
    FDR := p.adjust(
      P,
      method = "BH"
    )
  ]

  setorder(
    resultados_genes,
    FDR,
    P,
    na.last = TRUE
  )

  # ----------------------------------------------------------
  # 2.6 IDENTIFICAR GENES SIGNIFICATIVOS
  # ----------------------------------------------------------

  genes_significativos <- resultados_genes[
    !is.na(FDR) & FDR < 0.05
  ]

  # ----------------------------------------------------------
  # 2.7 RESUMEN
  # ----------------------------------------------------------

  cat("\n========================================\n")
  cat("ANALISIS POR GEN COMPLETADO\n")
  cat("========================================\n")

  cat("Método: Kolmogorov-Smirnov\n")
  cat("Descriptivo: media\n")
  cat("Genes analizados:", nrow(resultados_genes), "\n")

  cat(
    "Genes con FDR < 0.05:",
    nrow(genes_significativos),
    "\n"
  )

  cat("\nPrimeros 20 genes:\n")
  print(head(resultados_genes, 20))

  cat("\nGenes significativos:\n")
  print(genes_significativos)

  # ----------------------------------------------------------
  # 2.8 DEVOLVER RESULTADOS
  # ----------------------------------------------------------

  return(
    list(
      resultados_genes = resultados_genes,
      genes_significativos = genes_significativos
    )
  )
}


# ------------------------------------------------------------
# 3. EJEMPLO DE EJECUCION
# ------------------------------------------------------------

# resultado_genes <- analizar_cobertura_por_gen(
#   archivo_subset = "cobertura_subtelomerica.tsv",
#   tabla_analisis = tabla_analisis
# )
#
# resultados_genes <- resultado_genes$resultados_genes
# genes_significativos <- resultado_genes$genes_significativos
