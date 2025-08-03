# ---- CONFIGURAR DIRECTORIO DE TRABAJO ----
# Carpeta base donde se encuentra la carpeta DWCA y se guardarán los resultados
# Cambiar "D:/IAvH/M1" por la ruta correspondiente en su equipo
DIR_BASE <- "D:/IAvH/M1"
DIR_DWCA <- file.path(DIR_BASE, "dwca-ictiofauna_colombiana_dulceacuicola-v2.17")
if (!dir.exists(DIR_BASE)) stop("El directorio base no existe: ", DIR_BASE)
if (!dir.exists(DIR_DWCA)) stop("No se encontró la carpeta DWCA en: ", DIR_DWCA)
setwd(DIR_BASE)
message("Directorio de trabajo: ", getwd())

# ---- CARGAR PAQUETES Y CONFIGURAR MULTIPROCESAMIENTO ----
pkgs <- c("rgbif", "data.table", "future.apply", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

# Configuración de threads y multiprocesamiento
cores <- parallel::detectCores()
setDTthreads(cores)
future::plan(multisession, workers = max(1, cores - 1))
message("Plan de future establecido con ", max(1, cores - 1), " workers")

# ---- DEFINIR RUTAS Y CREAR DIRECTORIOS ----
PATH_TAXON <- file.path(DIR_DWCA, "taxon.txt")
PATH_DISTRIBUTION <- file.path(DIR_DWCA, "distribution.txt")
PATH_OUTPUT_DIR <- file.path(DIR_BASE, "output")
PATH_RAW <- file.path(PATH_OUTPUT_DIR, "gbif_presencia_raw.csv")
PATH_FINAL <- file.path(PATH_OUTPUT_DIR, "gbif_presencia_filtrada.csv")

for (d in unique(c(PATH_OUTPUT_DIR))) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---- FUNCIÓN ASSERT ----
auto_assert <- function(cond, msg) if (!cond) stop(msg)

# ---- LEER LISTA DE ESPECIES NO ENDEMICAS ----
auto_assert(file.exists(PATH_TAXON), paste("Archivo taxon no encontrado:", PATH_TAXON))
auto_assert(file.exists(PATH_DISTRIBUTION), paste("Archivo distribution no encontrado:", PATH_DISTRIBUTION))

# Leer archivos
taxon_df <- fread(PATH_TAXON, sep = "\t", select = c("id", "scientificName"))
dist_df  <- fread(PATH_DISTRIBUTION, sep = "\t", select = c("id", "establishmentMeans"))

# Filtrar ids no endémicos (incluye NA)
non_endemic_ids <- dist_df[is.na(establishmentMeans) | establishmentMeans != "Endémico", id]

# Mantener solo especies no endémicas
spec_df <- taxon_df[id %in% non_endemic_ids]
auto_assert(nrow(spec_df) > 0, "No se encontraron especies no endémicas")

taxon_names <- unique(na.omit(spec_df$scientificName))
auto_assert(length(taxon_names) > 0, "Lista de nombres científicos vacía")
message("Especies no endémicas a procesar: ", length(taxon_names))

# ---- RESOLVER taxonKey EN PARALELO ----
df_list <- future_lapply(taxon_names, function(nm) {
  res <- rgbif::name_backbone(name = nm)
  if (!is.null(res$usageKey))
    data.frame(scientificName = nm, taxonKey = res$usageKey, stringsAsFactors = FALSE)
  else
    NULL
}, future.seed = TRUE)

df_valid <- do.call(rbind, Filter(Negate(is.null), df_list))
auto_assert(nrow(df_valid) > 0, "No se resolvieron taxonKeys válidos en GBIF")
message("TaxonKeys válidos obtenidos: ", nrow(df_valid))

# ---- PARÁMETROS OPCIONALES ----
# Valores posibles para basisOfRecord (GBIF):
# "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN", "LIVING_SPECIMEN",
# "HUMAN_OBSERVATION", "MACHINE_OBSERVATION", "MATERIAL_SAMPLE",
# "LITERATURE", "UNKNOWN"
BASIS_OF_RECORD <- NULL        # e.g., c("HUMAN_OBSERVATION")
YEAR_MIN <- NULL               # e.g., 2000
YEAR_MAX <- NULL               # e.g., 2024
BATCH_SIZE <- 100              # Número de especies por lote

# ---- DESCARGAR REGISTROS POR LOTES ----
num_batches <- ceiling(nrow(df_valid) / BATCH_SIZE)
if (file.exists(PATH_RAW)) file.remove(PATH_RAW)  # comenzar desde cero

for (b in seq_len(num_batches)) {
  idx_start <- (b - 1) * BATCH_SIZE + 1
  idx_end   <- min(b * BATCH_SIZE, nrow(df_valid))
  df_batch  <- df_valid[idx_start:idx_end, ]
  message("\n=== Procesando lote ", b, " de ", num_batches, " (", nrow(df_batch), " especies) ===")
  registros <- list()

  for (i in seq_len(nrow(df_batch))) {
    key <- df_batch$taxonKey[i]
    name <- df_batch$scientificName[i]
    cat("🔍 Descargando registros para:", name, "(taxonKey:", key, ")\n")

    res <- tryCatch({
      occ_search(
        taxonKey = key,
        hasCoordinate = TRUE,
        limit = 100000,
        basisOfRecord = BASIS_OF_RECORD,
        year = if (!is.null(YEAR_MIN) && !is.null(YEAR_MAX)) paste(YEAR_MIN, YEAR_MAX, sep = ",") else NULL
      )$data
    }, error = function(e) {
      cat("  ⚠️  Error al descargar para:", name, "\n")
      NULL
    })

    if (!is.null(res) && nrow(res) > 0) {
      if (!is.null(BASIS_OF_RECORD)) res <- dplyr::filter(res, basisOfRecord %in% BASIS_OF_RECORD)
      if (!is.null(YEAR_MIN))        res <- dplyr::filter(res, !is.na(year) & year >= YEAR_MIN)
      if (!is.null(YEAR_MAX))        res <- dplyr::filter(res, !is.na(year) & year <= YEAR_MAX)
      if (nrow(res) > 0) {
        res$scientificNameQueried <- name
        registros[[length(registros) + 1]] <- res
      }
    }
  }

  if (length(registros) > 0) {
    registros_batch <- dplyr::bind_rows(registros)
    data.table::fwrite(
      registros_batch,
      PATH_RAW,
      append = file.exists(PATH_RAW),
      col.names = !file.exists(PATH_RAW)
    )
    cat("📦 Lote", b, "guardado con", nrow(registros_batch), "registros\n")
  }
}

# ---- POST-PROCESO: UNIR Y RESUMEN ----
if (file.exists(PATH_RAW)) {
  registros_todos <- data.table::fread(PATH_RAW)
  message("📦 Total de registros descargados: ", nrow(registros_todos))
}

# ---- GENERAR RESUMEN PARA DATOS CRUDOS ----
crear_resumen_txt <- function(df, nombre_archivo, taxon_names, df_valid, tipo = "RAW") {
  resumen <- c()
  resumen <- c(resumen, paste0("=== RESUMEN DE DESCARGA GBIF - ", tipo, " ==="), "")
  resumen <- c(resumen, ">>> Parámetros de búsqueda:")
  resumen <- c(resumen, " - hasCoordinate = TRUE")
  if (!is.null(BASIS_OF_RECORD)) resumen <- c(resumen, paste0(" - basisOfRecord = ", paste(BASIS_OF_RECORD, collapse = ",")))
  if (!is.null(YEAR_MIN) || !is.null(YEAR_MAX)) resumen <- c(resumen, paste0(" - year = ", YEAR_MIN, "-", YEAR_MAX))
  resumen <- c(resumen, "")

  resumen <- c(resumen, ">>> Especies:")
  resumen <- c(resumen, paste0(" - Total especies en lista: ", length(taxon_names)))
  resumen <- c(resumen, paste0(" - Con taxonKey válido: ", nrow(df_valid)))
  resumen <- c(resumen, paste0(" - Sin taxonKey válido: ", length(taxon_names) - nrow(df_valid)))
  resumen <- c(resumen, "")

  resumen <- c(resumen, ">>> Registros:")
  resumen <- c(resumen, paste0(" - Total de registros: ", nrow(df)))
  total_especies <- if ("scientificNameQueried" %in% names(df)) length(unique(df$scientificNameQueried)) else length(unique(df$scientificName))
  resumen <- c(resumen, paste0(" - Total de especies únicas: ", total_especies))
  resumen <- c(resumen, "")

  if (tipo == "CRUDO") {
    sin_key <- setdiff(taxon_names, df_valid$scientificName)
    resumen <- c(resumen, ">>> Especies sin taxonKey:")
    if (length(sin_key) == 0) resumen <- c(resumen, " - Ninguna") else resumen <- c(resumen, paste0(" - ", sin_key))
    resumen <- c(resumen, "")
  }

  writeLines(resumen, con = nombre_archivo)
  message("📝 Resumen escrito en: ", nombre_archivo)
}

if (exists("registros_todos")) {
  crear_resumen_txt(
    df = registros_todos,
    nombre_archivo = file.path(PATH_OUTPUT_DIR, "resumen_raw_gbif.txt"),
    taxon_names = taxon_names,
    df_valid = df_valid,
    tipo = "CRUDO"
  )
}

