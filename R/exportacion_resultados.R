
# ============================================================
# HERRAMIENTA DE COBERTURA GENOMICA
# MODULO 05: EXPORTACION Y AGRUPACION DE RESULTADOS
# GENES SIGNIFICATIVOS POR CROMOSOMA
# COMPARACION: MEDIA + KOLMOGOROV-SMIRNOV
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. FUNCION PRINCIPAL
# ------------------------------------------------------------

exportar_genes_por_cromosoma <- function(
    resultados_genes,
    datos,
    directorio_salida = "resultados",
    umbral_fdr = 0.05
) {

  # ----------------------------------------------------------
  # 2. VALIDAR INFORMACION DE ENTRADA
  # ----------------------------------------------------------

  resultados_genes <- as.data.table(
    copy(resultados_genes)
  )

  datos <- as.data.table(copy(datos))

  columnas_resultados <- c(
    "gene_name",
    "n_regiones",
    "media_menores_60",
    "media_60_o_mas",
    "diferencia_media",
    "P",
    "FDR"
  )

  columnas_anotacion <- c(
    "gene_name",
    "chromosome"
  )

  if (!all(columnas_resultados %in% names(resultados_genes))) {
    stop("Faltan columnas en resultados_genes.")
  }

  if (!all(columnas_anotacion %in% names(datos))) {
    stop("Faltan columnas de anotacion en datos.")
  }

  # ----------------------------------------------------------
  # 3. CREAR DIRECTORIO DE RESULTADOS
  # ----------------------------------------------------------

  if (!dir.exists(directorio_salida)) {
    dir.create(
      directorio_salida,
      recursive = TRUE
    )
  }

  # ----------------------------------------------------------
  # 4. ASOCIAR CADA GEN CON SU CROMOSOMA
  # ----------------------------------------------------------

  mapa_genes <- unique(
    datos[
      !is.na(gene_name) &
        !is.na(chromosome),
      .(
        gene_name,
        chromosome
      )
    ]
  )

  resultados_completos <- merge(
    resultados_genes,
    mapa_genes,
    by = "gene_name",
    all.x = TRUE
  )

  # ----------------------------------------------------------
  # 5. FILTRAR GENES SIGNIFICATIVOS
  # ----------------------------------------------------------

  genes_significativos <- resultados_completos[
    !is.na(P) &
      !is.na(FDR) &
      FDR < umbral_fdr
  ]

  setorder(
    genes_significativos,
    chromosome,
    FDR,
    P
  )

  # ----------------------------------------------------------
  # 6. EXPORTAR TODOS LOS GENES SIGNIFICATIVOS
  # ----------------------------------------------------------

  fwrite(
    genes_significativos,
    file.path(
      directorio_salida,
      "genes_significativos_por_cromosoma.csv"
    )
  )

  # ----------------------------------------------------------
  # 7. RESUMEN DE GENES POR CROMOSOMA
  # ----------------------------------------------------------

  resumen_cromosomas <- genes_significativos[
    ,
    .(
      numero_genes_significativos = uniqueN(gene_name),
      P_minimo = min(P, na.rm = TRUE),
      FDR_minimo = min(FDR, na.rm = TRUE)
    ),
    by = chromosome
  ]

  setorder(
    resumen_cromosomas,
    chromosome
  )

  fwrite(
    resumen_cromosomas,
    file.path(
      directorio_salida,
      "resumen_genes_significativos_por_cromosoma.csv"
    )
  )

  # ----------------------------------------------------------
  # 8. CREAR LISTA NUMERADA POR CROMOSOMA
  # ----------------------------------------------------------

  lista_genes <- genes_significativos[
    ,
    .(
      numero = seq_len(.N),
      gen = gene_name,
      regiones = n_regiones,
      media_menores_60,
      media_60_o_mas,
      diferencia_media,
      P,
      FDR
    ),
    by = chromosome
  ]

  fwrite(
    lista_genes,
    file.path(
      directorio_salida,
      "lista_genes_significativos_por_cromosoma.csv"
    )
  )

  # ----------------------------------------------------------
  # 9. AGRUPAR UNICAMENTE LOS NOMBRES DE LOS GENES
  # ----------------------------------------------------------

  lista_nombres <- genes_significativos[
    order(FDR, P),
    .(
      genes = paste(
        unique(gene_name),
        collapse = ", "
      )
    ),
    by = chromosome
  ]

  setorder(
    lista_nombres,
    chromosome
  )

  fwrite(
    lista_nombres,
    file.path(
      directorio_salida,
      "lista_nombres_genes_por_cromosoma.csv"
    )
  )

  # ----------------------------------------------------------
  # 10. MOSTRAR RESULTADOS
  # ----------------------------------------------------------

  cat("\n========================================\n")
  cat("EXPORTACION DE GENES SIGNIFICATIVOS\n")
  cat("========================================\n")

  cat("Comparacion: media + KS\n")

  cat(
    "Umbral FDR:",
    umbral_fdr,
    "\n"
  )

  cat(
    "Total de genes significativos:",
    uniqueN(genes_significativos$gene_name),
    "\n"
  )

  cat(
    "Cromosomas con genes significativos:",
    uniqueN(genes_significativos$chromosome),
    "\n"
  )

  cat("\nResumen por cromosoma:\n")

  print(resumen_cromosomas)

  cat("\nLista de genes significativos:\n")

  print(lista_genes, nrows = Inf)

  cat("\nNombres agrupados por cromosoma:\n")

  print(lista_nombres)

  cat("\nArchivos exportados en:", directorio_salida, "\n")

  # ----------------------------------------------------------
  # 11. DEVOLVER RESULTADOS
  # ----------------------------------------------------------

  return(
    list(
      resultados_completos = resultados_completos,
      genes_significativos = genes_significativos,
      resumen_cromosomas = resumen_cromosomas,
      lista_genes = lista_genes,
      lista_nombres = lista_nombres
    )
  )
}

# ------------------------------------------------------------
# 12. EJEMPLO DE EJECUCION
# ------------------------------------------------------------

# resultado_exportacion <- exportar_genes_por_cromosoma(
#   resultados_genes = resultados_genes,
#   datos = datos,
#   directorio_salida = "resultados",
#   umbral_fdr = 0.05
# )
