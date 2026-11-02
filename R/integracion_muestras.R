
# ============================================================
# HERRAMIENTA DE COBERTURA GENOMICA
# MODULO 01: INTEGRACION DE MUESTRAS
# ============================================================

library(data.table)

# ------------------------------------------------------------
# 1. INTEGRAR MUESTRAS GENOMICAS CON PARTICIPANTES
# ------------------------------------------------------------

integrar_muestras <- function(
    archivo_cuestionario,
    archivo_cobertura,
    columnas_anotacion = 6
) {

  # Verificar archivos
  if (!file.exists(archivo_cuestionario)) {
    stop("No se encontro el cuestionario.")
  }

  if (!file.exists(archivo_cobertura)) {
    stop("No se encontro la matriz de cobertura.")
  }

  # ----------------------------------------------------------
  # 2. CARGAR IDENTIFICADORES Y EDADES DEL CUESTIONARIO
  # ----------------------------------------------------------

  datos <- fread(
    archivo_cuestionario,
    sep = "\t",
    select = c(
      "codParticipante",
      "PARTICIPANTE.EDAD_REGISTRO"
    )
  )

  tabla_edad <- data.table(
    codParticipante = trimws(
      as.character(datos$codParticipante)
    ),
    edad = suppressWarnings(
      as.numeric(datos$PARTICIPANTE.EDAD_REGISTRO)
    )
  )

  # Eliminar identificadores vacios
  tabla_edad <- tabla_edad[
    !is.na(codParticipante) &
    codParticipante != ""
  ]

  # Verificar identificadores duplicados
  duplicados <- tabla_edad[
    ,
    .N,
    by = codParticipante
  ][N > 1]

  if (nrow(duplicados) > 0) {
    stop("Hay identificadores duplicados en el cuestionario.")
  }

  # Eliminar edades faltantes o invalidas
  tabla_edad <- tabla_edad[
    !is.na(edad) & edad >= 0
  ]

  # ----------------------------------------------------------
  # 3. IDENTIFICAR MUESTRAS EN LA MATRIZ DE COBERTURA
  # ----------------------------------------------------------

  encabezado <- names(
    fread(
      archivo_cobertura,
      nrows = 0,
      check.names = FALSE
    )
  )

  if (length(encabezado) <= columnas_anotacion) {
    stop("No se encontraron columnas de muestras.")
  }

  ids_muestras <- encabezado[
    -(seq_len(columnas_anotacion))
  ]

  tabla_muestras <- data.table(
    id_muestra = ids_muestras,
    codParticipante = sub(
      "^.*_",
      "",
      ids_muestras
    )
  )

  # ----------------------------------------------------------
  # 4. VINCULAR MUESTRAS Y PARTICIPANTES
  # ----------------------------------------------------------

  tabla_analisis <- merge(
    tabla_muestras,
    tabla_edad,
    by = "codParticipante",
    all = FALSE,
    sort = FALSE
  )

  # ----------------------------------------------------------
  # 5. CLASIFICAR GRUPOS DE EDAD
  # ----------------------------------------------------------

  tabla_analisis[
    ,
    grupo_edad := fifelse(
      edad < 60,
      "<60 años",
      "60 años o más"
    )
  ]

  # ----------------------------------------------------------
  # 6. RESUMEN DE INTEGRACION
  # ----------------------------------------------------------

  cat("\n===== INTEGRACION DE MUESTRAS =====\n")

  cat("Participantes:", nrow(tabla_edad), "\n")

  cat("Muestras en matriz:", length(ids_muestras), "\n")

  cat("Muestras vinculadas:", nrow(tabla_analisis), "\n")

  cat(
    "Muestras sin coincidencia:",
    length(ids_muestras) - nrow(tabla_analisis),
    "\n"
  )

  cat("\nParticipantes por grupo de edad:\n")

  print(
    tabla_analisis[
      ,
      .N,
      by = grupo_edad
    ]
  )

  return(tabla_analisis)
}
