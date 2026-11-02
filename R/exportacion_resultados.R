
# ============================================================
# HERRAMIENTA DE COBERTURA GENOMICA
# MODULO 05: EXPORTACION DE RESULTADOS
# AGRUPACION DE GENES SIGNIFICATIVOS POR CROMOSOMA
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. FUNCION PRINCIPAL
# ------------------------------------------------------------

exportar_genes_por_cromosoma <- function(
    resultados_genes,
    directorio_salida = "resultados",
    umbral_fdr = 0.05
) {

  resultados_genes <- as.data.table(
    copy(resultados_genes)
  )

  # ----------------------------------------------------------
  # 2. VALIDAR COLUMNAS
  # ----------------------------------------------------------

  columnas_requeridas <- c(
    "gene_name",
    "chromosome",
    "n_regiones",
    "media_menores_60",
    "media_60_o_mas",
    "diferencia_media",
    "P",
    "FDR"
  )

  faltantes <- setdiff(
    columnas_requeridas,
    names(resultados_genes)
  )

  if (length(faltantes) > 0) {
    stop(
      "Faltan columnas: ",
      paste(faltantes, collapse = ", ")
    )
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
  # 4. FILTRAR GENES SIGNIFICATIVOS
  # ----------------------------------------------------------

  genes_significativos <- resultados_genes[
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
  # 5. EXPORTAR GENES SIGNIFICATIVOS
  # ----------------------------------------------------------

  fwrite(
    genes_significativos,
    file.path(
      directorio_salida,
      "genes_significativos_por_cromosoma.csv"
    )
  )

  # ----------------------------------------------------------
  # 6. RESUMEN DE GENES POR CROMOSOMA
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
  # 7. CREAR LISTA NUMERADA POR CROMOSOMA
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
  # 8. AGRUPAR SOLAMENTE LOS NOMBRES DE LOS GENES
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
  # 9. MOSTRAR RESULTADOS
  # ----------------------------------------------------------

  cat("\n========================================\n")
  cat("EXPORTACION DE RESULTADOS COMPLETADA\n")
  cat("========================================\n")

  cat("Método estadístico: Kolmogorov-Smirnov\n")
  cat("Estadístico descriptivo: media\n")
  cat("Umbral FDR:", umbral_fdr, "\n")

  cat(
    "Total de genes significativos:",
    uniqueN(genes_significativos$gene_name),
    "\n"
  )

  cat(
    "Cromosomas identificados:",
    uniqueN(genes_significativos$chromosome),
    "\n"
  )

  cat("\nRESUMEN POR CROMOSOMA:\n")
  print(resumen_cromosomas)

  cat("\nLISTA DE GENES SIGNIFICATIVOS:\n")
  print(lista_genes, nrows = Inf)

  cat("\nNOMBRES AGRUPADOS POR CROMOSOMA:\n")
  print(lista_nombres)

  cat("\nArchivos guardados en:", directorio_salida, "\n")

  # ----------------------------------------------------------
  # 10. DEVOLVER RESULTADOS
  # ----------------------------------------------------------

  return(
    list(
      genes_significativos = genes_significativos,
      resumen_cromosomas = resumen_cromosomas,
      lista_genes = lista_genes,
      lista_nombres = lista_nombres
    )
  )
}


# ------------------------------------------------------------
# 11. EJEMPLO DE EJECUCION
# ------------------------------------------------------------

# resultado_exportacion <- exportar_genes_por_cromosoma(
#   resultados_genes = resultados_genes,
#   directorio_salida = "resultados",
#   umbral_fdr = 0.05
# )
