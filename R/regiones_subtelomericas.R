
# ============================================================
# HERRAMIENTA DE COBERTURA GENOMICA
# MODULO 02: IDENTIFICACION DE REGIONES SUBTELOMERICAS
# REFERENCIA: GRCh38
# VENTANA: 500 kb
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. CONFIGURACION
# ------------------------------------------------------------

ventana_telomero <- 500000

# Archivo de referencia GRCh38 de UCSC
url_chrom <- paste0(
  "https://hgdownload.soe.ucsc.edu/goldenPath/",
  "hg38/bigZips/hg38.chrom.sizes"
)

# ------------------------------------------------------------
# 2. OBTENER LONGITUDES CROMOSOMICAS
# ------------------------------------------------------------

obtener_longitudes_grch38 <- function(
    archivo_chrom = "hg38.chrom.sizes"
) {

  if (!file.exists(archivo_chrom)) {
    download.file(
      url_chrom,
      destfile = archivo_chrom,
      mode = "wb"
    )
  }

  tamanos_chr <- fread(
    archivo_chrom,
    sep = "\t",
    header = FALSE,
    col.names = c("chromosome", "longitud")
  )

  tamanos_chr <- tamanos_chr[
    grepl(
      "^chr([1-9]|1[0-9]|2[0-2]|X|Y)$",
      chromosome
    )
  ]

  tamanos_chr[
    ,
    chromosome := sub("^chr", "", chromosome)
  ]

  return(tamanos_chr)
}

# ------------------------------------------------------------
# 3. CONSTRUIR VENTANAS SUBTELOMERICAS
# ------------------------------------------------------------

crear_ventanas_subtelomericas <- function(
    tamanos_chr,
    ventana_telomero = 500000
) {

  ventanas <- rbindlist(
    lapply(seq_len(nrow(tamanos_chr)), function(i) {

      chr <- tamanos_chr$chromosome[i]
      longitud <- tamanos_chr$longitud[i]

      data.table(
        chromosome = chr,
        extremo = c("p", "q"),
        inicio = c(
          1,
          max(1, longitud - ventana_telomero + 1)
        ),
        fin = c(
          min(ventana_telomero, longitud),
          longitud
        )
      )
    })
  )

  return(ventanas)
}

# ------------------------------------------------------------
# 4. IDENTIFICAR REGIONES SUBTELOMERICAS
# ------------------------------------------------------------

identificar_regiones_subtelomericas <- function(
    anotaciones,
    ventanas_subtelomericas
) {

  anotaciones <- as.data.table(copy(anotaciones))

  anotaciones[
    ,
    chromosome := sub("^chr", "", as.character(chromosome))
  ]

  anotaciones[, start := as.numeric(start)]
  anotaciones[, end := as.numeric(end)]

  # Conservar número de fila original de la matriz.
  # Si ya existe fila_matriz, se conserva.
  if (!"fila_matriz" %in% names(anotaciones)) {
    anotaciones[
      ,
      fila_matriz := .I
    ]
  }

  regiones <- rbindlist(
    lapply(seq_len(nrow(ventanas_subtelomericas)), function(i) {

      w <- ventanas_subtelomericas[i]

      anotaciones[
        chromosome == w$chromosome &
          start <= w$fin &
          end >= w$inicio,
        .(
          fila_matriz,
          chromosome,
          start,
          end,
          gene_name,
          exon_rank,
          transcript_id
        )
      ][
        ,
        `:=`(
          extremo = w$extremo,
          ventana_inicio = w$inicio,
          ventana_fin = w$fin
        )
      ]
    }),
    use.names = TRUE,
    fill = TRUE
  )

  setorder(
    regiones,
    chromosome,
    extremo,
    start,
    end
  )

  # Resumen por cromosoma y extremo
  resumen <- regiones[
    ,
    .(
      regiones = .N,
      genes = uniqueN(gene_name)
    ),
    by = .(chromosome, extremo)
  ][order(chromosome, extremo)]

  cat("\n===== REGIONES SUBTELOMERICAS =====\n")

  cat(
    "Regiones encontradas:",
    nrow(regiones),
    "\n"
  )

  cat(
    "Genes distintos:",
    uniqueN(regiones$gene_name),
    "\n\n"
  )

  print(resumen)

  return(
    list(
      regiones = regiones,
      resumen = resumen
    )
  )
}

# ------------------------------------------------------------
# 5. EJECUCION DEL MODULO
# ------------------------------------------------------------

# Ejemplo de uso:
#
# tamanos_chr <- obtener_longitudes_grch38()
#
# ventanas_subtelomericas <- crear_ventanas_subtelomericas(
#   tamanos_chr,
#   ventana_telomero = ventana_telomero
# )
#
# resultado_subtel <- identificar_regiones_subtelomericas(
#   anotaciones_check,
#   ventanas_subtelomericas
# )
#
# regiones_subtelomericas <- resultado_subtel$regiones
# resumen_subtel <- resultado_subtel$resumen
