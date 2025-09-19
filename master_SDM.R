# =========================================================
# sdm_master.R  —  Flujo completo SDM (conciliación → GBIF → MaxEnt → riqueza)
# =========================================================
# Requisitos: R 4.x; Java instalado; dismo con maxent.jar disponible.
# ---------------------------------------------------------

# -------------------- PARAMS (EDITA AQUÍ) ----------------
DIR_BASE   <- "C:/Users/henry.garcia/Desktop/Primates"

# Entradas
PRIMATES_LISTA_CSV <- file.path(DIR_BASE, "primates_lista.csv")           # lista inicial de especies
PRESENT_TIF        <- "C:/Users/henry.garcia/Desktop/Primates/predictors/climate/predictors_vif/present_VIFsel.tif"
FUTURES            <- list(
  list(tag="ssp585_2040",
       tif="C:/Users/henry.garcia/Desktop/Primates/predictors/climate/predictors_vif/future_ssp585_2040_VIFsel.tif")
)

# Salidas y rutas intermedias
PATH_CONC          <- file.path(DIR_BASE, "primates_backbone_conciliacion.csv")
GBIF_DIR           <- file.path(DIR_BASE, "output_gbif")
GBIF_RAW_CSV       <- file.path(GBIF_DIR, "gbif_presencia_raw.csv")
GBIF_FINAL_CSV     <- file.path(GBIF_DIR, "gbif_presencia_filtrada.csv")
MODELS_DIR         <- file.path(DIR_BASE, "models_maxent_ALL_with_future") # usado también por 'riqueza'
RICH_OUT_DIR       <- file.path(MODELS_DIR, "richness_maps")

# Parámetros GBIF
YEAR_MIN <- 1970
YEAR_MAX <- as.integer(format(Sys.Date(), "%Y"))
BASIS_OF_RECORD <- NULL   # p.ej., c("HUMAN_OBSERVATION","PRESERVED_SPECIMEN")

# Ejecutar pasos
RUN_STEP_1_CONCILIACION <- TRUE
RUN_STEP_2_GBIF         <- TRUE
RUN_STEP_3_MAXENT       <- TRUE
RUN_STEP_4_RIQUEZA      <- TRUE
# ---------------------------------------------------------

set.seed(123)
options(java.parameters = "-Xmx4g")  # ajusta según tu RAM

# Paquetes base
pkgs <- c(
  "readr","dplyr","data.table","stringr","purrr",
  "rgbif","terra","sf","raster","dismo"
)
new <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# Utilidades
`%||%` <- function(a,b) if (!is.null(a)) a else b
auto_assert <- function(cond, msg) if (!cond) stop(msg)

dir.create(DIR_BASE, showWarnings = FALSE, recursive = TRUE)
dir.create(GBIF_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(MODELS_DIR, showWarnings = FALSE, recursive = TRUE)

# =========================================================
# STEP 1 — Conciliación taxonómica (GBIF backbone)
# =========================================================
.clean_sciname <- function(x) {
  x %>% as.character() %>% stringr::str_squish() %>%
    stringr::str_replace("\\s*\\([^)]*\\)\\s*$", "") %>%
    stringr::str_replace_all("(?i)\\b(cf\\.|aff\\.|sp\\.?|spp\\.?|s\\.l\\.|s\\.s\\.|\\?)\\b", "") %>%
    stringr::str_squish()
}
backbone_conciliar <- function(df, name_col = "scientificName", enforce_primates = TRUE) {
  stopifnot(name_col %in% names(df))
  nombres <- .clean_sciname(df[[name_col]]) %>% unique()
  pb <- utils::txtProgressBar(min = 0, max = length(nombres), style = 3); i <- 0L
  resolver_uno <- function(sci) {
    on.exit({ i <<- i + 1L; utils::setTxtProgressBar(pb, i) }, add = TRUE)
    out <- tryCatch({
      parts <- strsplit(sci, "\\s+")[[1]]
      bb <- rgbif::name_backbone(name = sci, rank = "SPECIES")
      if (is.null(bb$usageKey) && length(parts)>=2)
        bb <- rgbif::name_backbone(name = paste(parts[1],parts[2]), rank="SPECIES")
      if (is.null(bb$usageKey) && length(parts)>=1)
        bb <- rgbif::name_backbone(name = parts[1], rank="GENUS")
      if (is.null(bb$usageKey)) {
        return(tibble::tibble(scientificName_input=sci, gbif_scientificName=NA_character_, acceptedUsageKey=NA_integer_))
      }
      if (isTRUE(enforce_primates) && !identical(bb$order,"Primates"))
        return(tibble::tibble(scientificName_input=sci, gbif_scientificName=NA_character_, acceptedUsageKey=NA_integer_))
      key <- bb$acceptedUsageKey %||% bb$usageKey %||% NA_integer_
      nm  <- bb$scientificName %||% NA_character_
      tibble::tibble(scientificName_input=sci, gbif_scientificName=nm, acceptedUsageKey=as.integer(key))
    }, error = function(e) tibble::tibble(scientificName_input=sci, gbif_scientificName=NA_character_, acceptedUsageKey=NA_integer_))
    if (length(nombres) > 100) Sys.sleep(0.05)
    out
  }
  res <- purrr::map_dfr(nombres, resolver_uno); close(pb); res
}

if (RUN_STEP_1_CONCILIACION) {
  auto_assert(file.exists(PRIMATES_LISTA_CSV),
              paste("No existe el CSV de entrada con la lista de especies:", PRIMATES_LISTA_CSV))
  sp <- readr::read_csv(PRIMATES_LISTA_CSV, show_col_types = FALSE)
  auto_assert("scientificName" %in% names(sp),
              "Se espera una columna 'scientificName' en primates_lista.csv")
  
  conc <- backbone_conciliar(sp, name_col = "scientificName", enforce_primates = TRUE)
  conc_export <- conc %>% dplyr::select(scientificName_input, gbif_scientificName, acceptedUsageKey)
  readr::write_csv(conc_export, PATH_CONC)
  
  n_ok   <- sum(!is.na(conc_export$acceptedUsageKey))
  n_fail <- sum(is.na(conc_export$acceptedUsageKey))
  cat("STEP1 ✅ Conciliación ->", PATH_CONC, "| OK:", n_ok, "| Sin match:", n_fail, "\n")
} else {
  cat("STEP1 ⏭ Saltado (se asume PATH_CONC existente)\n")
}

# =========================================================
# STEP 2 — Descarga/curaduría de ocurrencias (GBIF API)
# =========================================================
fetch_occ_all <- function(taxonKey, hasCoordinate=TRUE, basisOfRecord=NULL, year=NULL,
                          page_size=300, max_records=1e6) {
  offset <- 0L; out <- list()
  repeat {
    dat <- tryCatch({
      rgbif::occ_search(
        taxonKey = taxonKey, hasCoordinate = hasCoordinate,
        basisOfRecord = basisOfRecord, year = year,
        limit = page_size, start = offset
      )$data
    }, error=function(e) NULL)
    if (is.null(dat) || nrow(dat)==0) break
    out[[length(out)+1L]] <- dat
    got <- nrow(dat); offset <- offset + got
    if (got < page_size || offset >= max_records) break
  }
  if (length(out)) data.table::rbindlist(out, fill=TRUE) else NULL
}
fetch_with_retries <- function(taxonKey, hasCoordinate=TRUE, basisOfRecord=NULL, year=NULL, tries=4) {
  res <- NULL
  for (att in 0:(tries-1)) {
    res <- tryCatch(fetch_occ_all(taxonKey, hasCoordinate, basisOfRecord, year), error=function(e) NULL)
    if (!is.null(res) && nrow(res)>0) return(res)
    Sys.sleep(min(60, (2^att) + stats::runif(1,0,0.7)))
  }
  res
}

if (RUN_STEP_2_GBIF) {
  auto_assert(file.exists(PATH_CONC), paste("No existe el CSV de conciliación:", PATH_CONC))
  conc <- readr::read_csv(PATH_CONC, show_col_types = FALSE)
  req_cols <- c("scientificName_input","gbif_scientificName","acceptedUsageKey")
  auto_assert(all(req_cols %in% names(conc)), "Faltan columnas requeridas en conciliación.")
  
  conc <- conc %>%
    mutate(
      taxonKey = as.integer(acceptedUsageKey),
      scientificName = ifelse(!is.na(gbif_scientificName) & nzchar(gbif_scientificName),
                              gbif_scientificName, scientificName_input)
    )
  df_valid <- conc %>% filter(!is.na(taxonKey)) %>% distinct(taxonKey, scientificName, .keep_all = FALSE)
  auto_assert(nrow(df_valid)>0, "df_valid vacío: no hay acceptedUsageKey válidos.")
  
  if (file.exists(GBIF_RAW_CSV)) file.remove(GBIF_RAW_CSV)
  failed_fetch <- character(0); total_records <- 0L
  BATCH_SIZE <- 100
  num_batches <- ceiling(nrow(df_valid)/BATCH_SIZE)
  
  for (b in seq_len(num_batches)) {
    idx  <- ((b-1)*BATCH_SIZE+1):min(b*BATCH_SIZE, nrow(df_valid))
    lote <- df_valid[idx, ]
    message("\nSTEP2 ▶ Lote ", b, "/", num_batches, " (", nrow(lote), " spp)")
    
    registros <- vector("list", nrow(lote))
    for (i in seq_len(nrow(lote))) {
      key  <- lote$taxonKey[i]; name <- lote$scientificName[i]
      res <- fetch_with_retries(
        taxonKey = key, hasCoordinate = TRUE,
        basisOfRecord = BASIS_OF_RECORD,
        year = paste(YEAR_MIN, YEAR_MAX, sep = ","),
        tries = 4
      )
      if (!is.null(res) && nrow(res)>0) {
        if (!is.null(BASIS_OF_RECORD)) res <- dplyr::filter(res, basisOfRecord %in% BASIS_OF_RECORD)
        if (!is.null(YEAR_MIN))        res <- dplyr::filter(res, !is.na(year) & year >= YEAR_MIN)
        if (!is.null(YEAR_MAX))        res <- dplyr::filter(res, !is.na(year) & year <= YEAR_MAX)
        if (nrow(res)>0) { res$scientificNameQueried <- name; registros[[i]] <- res }
      } else failed_fetch <- c(failed_fetch, paste0(name,"|",key))
    }
    
    regs_non_null <- Filter(Negate(is.null), registros)
    if (length(regs_non_null)) {
      registros_batch <- dplyr::bind_rows(regs_non_null) %>% as.data.frame()
      for (nn in c("gbifID","occurrenceID","eventDate"))
        if (nn %in% names(registros_batch)) registros_batch[[nn]] <- as.character(registros_batch[[nn]])
      
      columnas_importantes <- c(
        "scientificNameQueried","scientificName","verbatimScientificName",
        "kingdom","phylum","class","order","family","genus","species","taxonRank","taxonKey",
        "gbifID","occurrenceID","occurrenceStatus","basisOfRecord",
        "eventDate","year","month","day",
        "decimalLatitude","decimalLongitude","coordinatePrecision",
        "coordinateUncertaintyInMeters","geodeticDatum",
        "elevation","elevationAccuracy","depth","depthAccuracy",
        "country","countryCode","stateProvince","county","municipality","locality",
        "datasetKey","publishingOrgKey","license",
        "institutionCode","collectionCode","catalogNumber",
        "recordedBy","identifiedBy","samplingProtocol","samplingEffort",
        "issues"
      )
      for (col in columnas_importantes) if (!col %in% names(registros_batch)) registros_batch[[col]] <- NA
      keep <- intersect(columnas_importantes, names(registros_batch))
      registros_batch <- registros_batch[, keep, drop = FALSE]
      
      # aplanar list-cols si aparecen
      list_cols <- names(registros_batch)[vapply(registros_batch, is.list, logical(1))]
      if (length(list_cols)) {
        registros_batch[list_cols] <- lapply(registros_batch[list_cols], function(v) {
          vapply(v, function(x) if (is.null(x) || length(x) == 0) NA_character_ else paste0(as.character(unlist(x)), collapse = ";"),
                 character(1))
        })
      }
      
      data.table::setDT(registros_batch)
      if ("gbifID" %in% names(registros_batch)) registros_batch <- unique(registros_batch, by = "gbifID")
      total_records <- total_records + nrow(registros_batch)
      data.table::fwrite(registros_batch, GBIF_RAW_CSV,
                         append = file.exists(GBIF_RAW_CSV),
                         col.names = !file.exists(GBIF_RAW_CSV))
      cat("   📦 Lote", b, "->", nrow(registros_batch), "registros guardados\n")
    } else {
      cat("   ⚠️ Lote", b, "sin registros válidos\n")
    }
    rm(registros); gc()
  }
  
  # CSV final compacto
  if (file.exists(GBIF_RAW_CSV)) {
    registros_todos <- readr::read_csv(GBIF_RAW_CSV, show_col_types = FALSE, guess_max = 100000)
    if ("gbifID" %in% names(registros_todos)) {
      data.table::setDT(registros_todos)
      registros_todos <- unique(registros_todos, by = "gbifID")
      registros_todos <- as.data.frame(registros_todos)
    }
    cols_final <- c(
      "scientificNameQueried","scientificName","taxonKey","gbifID","occurrenceID","basisOfRecord",
      "eventDate","year","month","day","decimalLatitude","decimalLongitude",
      "coordinateUncertaintyInMeters","country","countryCode","stateProvince","municipality","locality",
      "datasetKey","license","institutionCode","collectionCode","catalogNumber","recordedBy","identifiedBy","issues"
    )
    keep_final <- intersect(cols_final, names(registros_todos))
    readr::write_csv(registros_todos[, keep_final, drop = FALSE], GBIF_FINAL_CSV)
    message("STEP2 ✅ GBIF final -> ", GBIF_FINAL_CSV, " | Registros: ", nrow(registros_todos))
  } else {
    stop("STEP2: No se generó GBIF_RAW_CSV; revisa conexión/queries.")
  }
} else {
  cat("STEP2 ⏭ Saltado (se asume GBIF_FINAL_CSV existente)\n")
}

# =========================================================
# STEP 3 — Modelación MaxEnt (presente + futuros)
# =========================================================
safe_dir <- function(x) gsub("[^A-Za-z0-9_\\-]", "_", x)
haversine_km <- function(lon1, lat1, lon2, lat2){
  rad <- pi/180; dlat <- (lat2-lat1)*rad; dlon <- (lon2-lon1)*rad
  a <- sin(dlat/2)^2 + cos(lat1*rad)*cos(lat2*rad)*sin(dlon/2)^2
  6371*2*atan2(sqrt(a), sqrt(1-a))
}
make_MCP_buffer_robusto <- function(pres_df, buffer_km=25, trim_q=0.995){
  c_lon <- median(pres_df$lon); c_lat <- median(pres_df$lat)
  d_km  <- haversine_km(pres_df$lon, pres_df$lat, c_lon, c_lat)
  thr   <- as.numeric(quantile(d_km, trim_q, na.rm=TRUE))
  keep  <- d_km <= thr
  pres_trim <- pres_df[keep, , drop=FALSE]
  pts_sf <- sf::st_as_sf(pres_trim, coords=c("lon","lat"), crs=4326)
  mcp    <- pts_sf |> sf::st_union() |> sf::st_convex_hull()
  mcp_buf <- mcp |> sf::st_transform(3857) |> sf::st_buffer(buffer_km*1000) |> sf::st_transform(4326)
  list(M=mcp_buf, pres_trim=pres_trim, removed=sum(!keep), thr_km=thr, mcp=mcp)
}
normalize_names <- function(nm){
  nm <- gsub("[\u00A0\u2007\u202F]", " ", nm); nm <- gsub("\\.+", " ", nm); gsub("\\s+", " ", trimws(nm))
}
numify <- function(x){
  v <- suppressWarnings(as.numeric(x))
  if (any(is.na(v))) v <- suppressWarnings(readr::parse_number(x, locale=readr::locale(decimal_mark=",", grouping_mark=".")))
  if (any(is.na(v))) v <- suppressWarnings(readr::parse_number(x, locale=readr::locale(decimal_mark=".", grouping_mark=",")))
  v
}
getcol_exact <- function(df, label){
  hits <- which(tolower(names(df)) == tolower(label)); if (length(hits)) names(df)[hits[1]] else character(0)
}
fc_to_args <- function(fc){
  base <- c("autofeature=false","linear=false","quadratic=false","product=false","threshold=false","hinge=false")
  c(base, if (grepl("L",fc)) "linear=true",
    if (grepl("Q",fc)) "quadratic=true",
    if (grepl("H",fc)) "hinge=true",
    if (grepl("P",fc)) "product=true",
    if (grepl("T",fc)) "threshold=true")
}
guess_col <- function(df, candidates) {
  hits <- names(df)[tolower(names(df)) %in% tolower(candidates)]
  if (length(hits)) hits[1] else stop("No se encontraron columnas: ", paste(candidates, collapse=", "))
}

if (RUN_STEP_3_MAXENT) {
  # Checks Java/MaxEnt
  maxent_jar <- file.path(system.file("java", package = "dismo"), "maxent.jar")
  auto_assert(file.exists(maxent_jar), paste("No se encontró maxent.jar en dismo:", maxent_jar, "\nInstala Java/dismo correctamente."))
  auto_assert(file.exists(PRESENT_TIF), paste("No existe stack de predictores presente:", PRESENT_TIF))
  dir.create(MODELS_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Occs (usa el CSV final del STEP2)
  auto_assert(file.exists(GBIF_FINAL_CSV), paste("No existe GBIF_FINAL_CSV:", GBIF_FINAL_CSV))
  occ_raw <- readr::read_csv(GBIF_FINAL_CSV, show_col_types = FALSE)
  col_species <- guess_col(occ_raw, c("species_bm","species","scientificName","scientificname","scientificNameQueried"))
  col_lon     <- guess_col(occ_raw, c("longitude","decimalLongitude","lon"))
  col_lat     <- guess_col(occ_raw, c("latitude","decimalLatitude","lat"))
  
  occ <- occ_raw %>%
    transmute(
      species = .data[[col_species]],
      lon     = suppressWarnings(as.numeric(.data[[col_lon]])),
      lat     = suppressWarnings(as.numeric(.data[[col_lat]]))
    ) %>%
    filter(is.finite(lon), is.finite(lat), !is.na(species)) %>%
    distinct(species, lon, lat, .keep_all = TRUE)
  
  sp_counts   <- occ %>% count(species, name="n_pres") %>% filter(n_pres > 10) %>% arrange(species)
  species_all <- sp_counts$species
  auto_assert(length(species_all)>0, "No hay especies con >10 presencias.")
  
  env_full <- terra::rast(PRESENT_TIF); names(env_full) <- paste0("v", seq_len(nlyr(env_full)))
  FC_SET <- c("LQ","LQH","LQHP"); RM_SET <- c(0.5,1,2); REPS <- 5
  
  run_one_species <- function(sp) {
    out_dir <- file.path(MODELS_DIR, safe_dir(sp))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    stats_csv <- file.path(out_dir, "model_stats.csv")
    if (file.exists(stats_csv)) return(readr::read_csv(stats_csv, show_col_types = FALSE, n_max = 1))
    
    pres_all <- occ %>% filter(species == sp) %>% select(lon, lat)
    if (nrow(pres_all) <= 10) { writeLines("<=10 presencias", con=file.path(out_dir,"ERROR.txt")); return(NULL) }
    
    # Área M
    Mres   <- make_MCP_buffer_robusto(pres_all, buffer_km=25, trim_q=0.995)
    M_poly <- Mres$M
    area_km2 <- as.numeric(sf::st_area(sf::st_transform(M_poly, 3857))) / 1e6
    # Plot M
    png(file.path(out_dir, "area_M_plot.png"), width=1200, height=900, res=150)
    plot(sf::st_geometry(M_poly), col=adjustcolor("orange",0.25), border="orange", lwd=2,
         main=sprintf("Área M ( %s )", sp), reset=FALSE)
    plot(sf::st_geometry(Mres$mcp), add=TRUE, col=NA, border="red", lwd=2, lty=2)
    in_idx <- paste(pres_all$lon, pres_all$lat) %in% paste(Mres$pres_trim$lon, Mres$pres_trim$lat)
    if (any(in_idx))  points(pres_all$lon[in_idx],  pres_all$lat[in_idx],  pch=20, col="blue")
    if (any(!in_idx)) points(pres_all$lon[!in_idx], pres_all$lat[!in_idx], pch=4,  col="darkorange", lwd=1.5)
    legend("topleft", legend=c("MCP + 25 km","MCP","Pres (usadas)","Outliers MCP"),
           lwd=c(2,2,NA,NA), lty=c(1,2,NA,NA), pch=c(NA,NA,20,4),
           col=c("orange","red","blue","darkorange"), bty="n")
    invisible(dev.off())
    sf::st_write(M_poly, file.path(out_dir, "area_M.gpkg"), delete_dsn=TRUE, quiet=TRUE)
    
    env_M_t   <- terra::mask(terra::crop(env_full, terra::vect(M_poly)), terra::vect(M_poly))
    vals_pres <- terra::extract(env_M_t, as.matrix(Mres$pres_trim))
    keep      <- complete.cases(vals_pres)
    pres_use  <- Mres$pres_trim[keep, , drop=FALSE]
    if (nrow(pres_use) <= 10) { writeLines("<=10 presencias válidas (NA en predictores).", con=file.path(out_dir,"ERROR.txt")); return(NULL) }
    
    tmp_tif  <- file.path(out_dir, "env_M.tif")
    terra::writeRaster(env_M_t, tmp_tif, overwrite=TRUE)
    env_M_rs <- raster::stack(tmp_tif)
    
    # Tuning
    tuning_rows <- list()
    for (FC in FC_SET) for (RM in RM_SET) {
      cv_dir <- file.path(out_dir, sprintf("maxent_cv_FC-%s_RM-%s", FC, RM))
      dir.create(cv_dir, showWarnings = FALSE)
      args_cv <- c("outputformat=logistic",
                   paste0("replicates=", REPS), "replicatetype=crossvalidate", "randomtestpoints=0",
                   fc_to_args(FC), paste0("betamultiplier=", RM))
      invisible(dismo::maxent(x=env_M_rs, p=as.matrix(pres_use[,c("lon","lat")]), path=cv_dir, args=args_cv))
      mx_csv <- file.path(cv_dir, "maxentResults.csv"); if (!file.exists(mx_csv)) next
      mx <- suppressWarnings(readr::read_csv(mx_csv, show_col_types = FALSE,
                                             locale=readr::locale(decimal_mark=".", grouping_mark=",")))
      names(mx) <- normalize_names(names(mx))
      spec_col <- names(mx)[1]; avg_row <- which(grepl("\\(average\\)\\s*$", tolower(mx[[spec_col]])))
      if (!length(avg_row)) next
      col_tr <- getcol_exact(mx,"Training AUC"); col_ts <- getcol_exact(mx,"Test AUC")
      thr_labels <- c("Maximum test sensitivity plus specificity Logistic threshold",
                      "Equal test sensitivity and specificity Logistic threshold",
                      "Maximum training sensitivity plus specificity Logistic threshold",
                      "Equal training sensitivity and specificity Logistic threshold")
      col_thr <- character(0); for (lbl in thr_labels){ cnd<-getcol_exact(mx,lbl); if (length(cnd)){col_thr<-cnd; break}}
      omit_labels <- c("Maximum test sensitivity plus specificity test omission",
                       "Equal test sensitivity and specificity test omission")
      col_omit <- character(0); for (lbl in omit_labels){ cnd<-getcol_exact(mx,lbl); if (length(cnd)){col_omit<-cnd; break}}
      tuning_rows[[length(tuning_rows)+1]] <- data.frame(
        species=sp, FC=FC, RM=RM,
        AUC_train_mean = if (length(col_tr)) round(numify(mx[[col_tr]][avg_row]),5) else NA,
        AUC_test_mean  = if (length(col_ts)) round(numify(mx[[col_ts]][avg_row]),5) else NA,
        thr_logistic   = if (length(col_thr)) round(numify(mx[[col_thr]][avg_row]),5) else NA,
        test_omission  = if (length(col_omit)) round(numify(mx[[col_omit]][avg_row]),5) else NA,
        stringsAsFactors=FALSE
      )
    }
    tuning_df <- if (length(tuning_rows)) dplyr::bind_rows(tuning_rows) else NULL
    if (is.null(tuning_df) || nrow(tuning_df)==0) { writeLines("Tuning vacío.", con=file.path(out_dir,"ERROR.txt")); return(NULL) }
    readr::write_csv(tuning_df, file.path(out_dir, "tuning_summary.csv"))
    
    ord <- with(tuning_df, order(-AUC_test_mean, ifelse(is.na(test_omission), Inf, test_omission), -RM))
    best <- tuning_df[ord[1], , drop=FALSE]
    best_FC <- best$FC[1]; best_RM <- best$RM[1]; best_thr <- best$thr_logistic[1]; best_auc_ts <- best$AUC_test_mean[1]
    writeLines(sprintf("Best combo: FC=%s RM=%s | Test AUC=%.3f | thr=%.4f",
                       best_FC, best_RM, best_auc_ts, best_thr),
               con = file.path(out_dir, "best_FC_RM.txt"))
    if (is.na(best_thr)) { writeLines("Sin umbral MSS (TEST) para la combinación ganadora.", con=file.path(out_dir,"ERROR.txt")); return(NULL) }
    
    args_final <- c("outputformat=logistic", fc_to_args(best_FC), paste0("betamultiplier=", best_RM))
    final_dir  <- file.path(out_dir, "maxent_final"); dir.create(final_dir, showWarnings = FALSE)
    me_final   <- dismo::maxent(x=env_M_rs, p=as.matrix(pres_use[,c("lon","lat")]), path=final_dir, args=args_final)
    
    r_cont  <- raster::predict(env_M_rs, me_final, progress="text")
    r_bin01 <- raster::calc(r_cont, function(x) as.integer(x >= best_thr))
    stem <- file.path(out_dir, gsub(" ", "_", sp))
    raster::writeRaster(r_cont,  paste0(stem,"_present_continuous.tif"),
                        overwrite=TRUE, datatype="FLT4S", options=c("COMPRESS=LZW"))
    raster::writeRaster(r_bin01, paste0(stem,"_present_binary_spec_sens.tif"),
                        overwrite=TRUE, datatype="INT1U", options=c("COMPRESS=LZW"))
    
    km2_pres <- sum(raster::values(raster::area(r_bin01)) * (raster::values(r_bin01)==1), na.rm=TRUE)
    
    fut_cols <- list()
    if (length(FUTURES)){
      for (sc in FUTURES){
        tag <- sc$tag; fut_path <- sc$tif
        km2_fut <- NA_real_; delta_pct <- NA_real_
        try({
          env_F0_t <- terra::rast(fut_path)
          if (terra::nlyr(env_F0_t) != terra::nlyr(env_full))
            stop(sprintf("Capas futuras (%d) != presentes (%d)", terra::nlyr(env_F0_t), terra::nlyr(env_full)))
          names(env_F0_t) <- names(env_full)
          env_F_t  <- terra::mask(terra::crop(env_F0_t, terra::vect(M_poly)), terra::vect(M_poly))
          tmpF     <- file.path(out_dir, paste0("env_F_", tag, ".tif"))
          terra::writeRaster(env_F_t, tmpF, overwrite=TRUE)
          env_F_rs <- raster::stack(tmpF)
          
          r_cont_fut <- raster::predict(env_F_rs, me_final, progress="text",
                                        args=c("outputformat=logistic","extrapolate=false","doclamp=false"))
          r_bin_fut  <- raster::calc(r_cont_fut, function(x) as.integer(x >= best_thr))
          
          raster::writeRaster(r_cont_fut, paste0(stem,"_future_continuous_",tag,".tif"),
                              overwrite=TRUE, datatype="FLT4S", options=c("COMPRESS=LZW"))
          raster::writeRaster(r_bin_fut,  paste0(stem,"_future_binary_spec_sens_",tag,".tif"),
                              overwrite=TRUE, datatype="INT1U", options=c("COMPRESS=LZW"))
          
          km2_fut   <- sum(raster::values(raster::area(r_bin_fut)) * (raster::values(r_bin_fut)==1), na.rm=TRUE)
          delta_pct <- if (is.finite(km2_pres) && km2_pres>0) 100*(km2_fut - km2_pres)/km2_pres else NA_real_
        })
        fut_cols[[paste0("area_km2_future_", tag)]] <- round(km2_fut, 2)
        fut_cols[[paste0("delta_pct_", tag)]]       <- round(delta_pct, 2)
      }
    }
    
    stats_df <- data.frame(
      species             = sp,
      n_pres_total        = nrow(pres_all),
      n_pres_used_in_M    = nrow(pres_use),
      n_predictors        = raster::nlayers(env_M_rs),
      predictors          = paste(names(env_M_rs), collapse=";"),
      area_M_km2          = round(area_km2, 2),
      area_km2_present    = round(km2_pres, 2),
      replicates          = 5,
      replicatetype       = "crossvalidate",
      best_FC             = best_FC,
      best_RM             = best_RM,
      AUC_test_mean       = round(best_auc_ts, 4),
      threshold_TEST_MSS  = round(best_thr, 4),
      stringsAsFactors    = FALSE
    )
    if (length(fut_cols)) stats_df <- cbind(stats_df, as.data.frame(fut_cols, check.names = FALSE))
    readr::write_csv(stats_df, file.path(out_dir, "model_stats.csv"))
    stats_df
  }
  
  all_rows <- list()
  for (sp in species_all) {
    cat("\nSTEP3 ▶", sp, "\n")
    res <- try(run_one_species(sp), silent = TRUE)
    if (inherits(res, "try-error") || is.null(res)) {
      out_dir <- file.path(MODELS_DIR, safe_dir(sp))
      writeLines(as.character(attr(res,"condition")), con = file.path(out_dir, "ERROR.txt"))
      message("ERROR en especie: ", sp)
    } else {
      all_rows[[length(all_rows)+1]] <- res
    }
    gc()
  }
  if (length(all_rows)) {
    all_df <- dplyr::bind_rows(all_rows)
    readr::write_csv(all_df, file.path(MODELS_DIR, "all_species_summary.csv"))
    message("STEP3 ✅ Resumen global -> all_species_summary.csv")
  } else {
    message("STEP3 ⚠ Sin modelos exitosos; revisa ERROR.txt por especie")
  }
} else {
  cat("STEP3 ⏭ Saltado (se asume MODELS_DIR con modelos listos)\n")
}

# =========================================================
# STEP 4 — Mapas de riqueza (presente + futuros)
# =========================================================
safe_name <- function(x) {
  x %>%
    gsub("[^A-Za-z0-9_\\-]+", "_", .) %>%
    gsub("_+", "_", .) %>%
    gsub("^_|_$", "", .)
}
if (RUN_STEP_4_RIQUEZA) {
  base_dir <- MODELS_DIR
  PAT_PRESENT <- "_present_binary_spec_sens\\.tif$"
  FUTS <- lapply(FUTURES, function(z) list(tag=z$tag, pat=paste0("_future_binary_spec_sens_", z$tag, "\\.tif$")))
  out_dir <- RICH_OUT_DIR; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  species_dirs_all <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  species_dirs     <- species_dirs_all[species_dirs_all != out_dir & dir.exists(species_dirs_all)]
  find_any_present <- function(sp_dir) {
    fs <- list.files(sp_dir, pattern = PAT_PRESENT, full.names = TRUE, ignore.case = TRUE)
    if (length(fs)) fs[1] else NA_character_
  }
  present_files <- lapply(species_dirs, find_any_present)
  names(present_files) <- basename(species_dirs)
  present_files <- present_files[!is.na(unlist(present_files))]
  auto_assert(length(present_files)>0, paste("No se encontraron TIF de presente con el patrón:", PAT_PRESENT))
  
  first_present <- try(terra::rast(present_files[[1]]), silent = TRUE)
  auto_assert(!inherits(first_present,"try-error"), "No pude abrir el primer raster presente.")
  if (is.na(terra::crs(first_present))) terra::crs(first_present) <- "epsg:4326"
  
  ext_union <- terra::ext(first_present)
  for (p in present_files) {
    r <- try(terra::rast(p), silent = TRUE); if (inherits(r, "try-error")) next
    if (!identical(terra::crs(r), terra::crs(first_present))) {
      r <- try(terra::project(r, first_present, method = "near"), silent = TRUE); if (inherits(r, "try-error")) next
    }
    ext_union <- terra::union(ext_union, terra::ext(r))
  }
  canvas <- terra::rast(extent = ext_union, resolution = terra::res(first_present), crs = terra::crs(first_present))
  
  present_layers <- list(); band_index_present <- list()
  for (sp in names(present_files)) {
    f <- present_files[[sp]]
    r <- try(terra::rast(f), silent = TRUE); if (inherits(r,"try-error")) next
    if (!identical(terra::crs(r), terra::crs(canvas))) {
      r <- try(terra::project(r, canvas, method = "near"), silent = TRUE); if (inherits(r,"try-error")) next
    }
    rA <- terra::resample(r, canvas, method = "near"); rA <- rA >= 1; rA <- terra::ifel(is.na(rA), 0, rA)
    nm <- safe_name(sp); names(rA) <- nm
    present_layers[[length(present_layers)+1]] <- rA
    band_index_present[[length(band_index_present)+1]] <- data.frame(band=length(present_layers), species=sp, filename=f)
  }
  auto_assert(length(present_layers)>0, "No quedaron capas presentes apiladas.")
  present_stack <- if (length(present_layers)==1) present_layers[[1]] else do.call(c, present_layers)
  present_rich  <- sum(present_stack); names(present_rich) <- "richness"
  present_final <- c(present_stack, present_rich)
  present_out_tif <- file.path(out_dir, "richness_present_stack.tif")
  terra::writeRaster(present_final, present_out_tif, overwrite=TRUE,
                     wopt=list(datatype="INT2U", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER")))
  present_bands_df <- dplyr::bind_rows(band_index_present)
  present_bands_df <- rbind(present_bands_df,
                            data.frame(band=nlyr(present_final), species="richness", filename=present_out_tif))
  readr::write_csv(present_bands_df, file.path(out_dir, "band_index_present.csv"))
  cat("STEP4 ✅ Stack presente -> ", present_out_tif, "\n", sep="")
  
  present_species_stack <- present_final[[1:(nlyr(present_final)-1)]]
  
  for (sc in FUTS) {
    tag <- sc$tag; pat <- sc$pat
    cat("\nSTEP4 ▶ Escenario futuro:", tag, "\n")
    future_files <- list()
    for (spdir in species_dirs) {
      files <- list.files(spdir, pattern = pat, full.names = TRUE, ignore.case = TRUE)
      if (length(files)) future_files[[basename(spdir)]] <- files[1]
    }
    if (!length(future_files)) { message("⚠ Sin archivos futuros para patrón: ", pat); next }
    
    fut_layers <- list(); band_index_future <- list()
    for (sp in names(future_files)) {
      f <- future_files[[sp]]
      r <- try(terra::rast(f), silent = TRUE); if (inherits(r,"try-error")) next
      if (!identical(terra::crs(r), terra::crs(canvas))) {
        r <- try(terra::project(r, canvas, method = "near"), silent = TRUE); if (inherits(r,"try-error")) next
      }
      rA <- terra::resample(r, canvas, method = "near"); rA <- rA >= 1; rA <- terra::ifel(is.na(rA), 0, rA)
      nm <- safe_name(sp); names(rA) <- nm
      fut_layers[[length(fut_layers)+1]] <- rA
      band_index_future[[length(band_index_future)+1]] <- data.frame(band=length(fut_layers), species=sp, filename=f)
    }
    if (!length(fut_layers)) { message("⚠ No quedaron capas futuras apiladas para ", tag); next }
    future_stack <- if (length(fut_layers)==1) fut_layers[[1]] else do.call(c, fut_layers)
    future_rich  <- sum(future_stack); names(future_rich) <- "richness"
    future_final <- c(future_stack, future_rich)
    fut_out_tif  <- file.path(out_dir, paste0("richness_future_stack_", tag, ".tif"))
    terra::writeRaster(future_final, fut_out_tif, overwrite=TRUE,
                       wopt=list(datatype="INT2U", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER")))
    future_bands_df <- dplyr::bind_rows(band_index_future)
    future_bands_df <- rbind(future_bands_df,
                             data.frame(band=nlyr(future_final), species="richness", filename=fut_out_tif))
    readr::write_csv(future_bands_df, file.path(out_dir, paste0("band_index_future_", tag, ".csv")))
    cat("STEP4 ✅ Futuro (", tag, ") -> ", fut_out_tif, "\n", sep="")
    
    # Cambios de riqueza (especies comunes)
    future_species_stack <- future_final[[1:(nlyr(future_final)-1)]]
    common_sp <- intersect(names(present_species_stack), names(future_species_stack))
    if (length(common_sp) >= 1) {
      P_common <- present_species_stack[[common_sp]]
      F_common <- future_species_stack[[common_sp]]
      rich_P_common <- sum(P_common); rich_F_common <- sum(F_common)
      rich_delta_abs <- rich_F_common - rich_P_common
      names(rich_delta_abs) <- paste0("richness_delta_abs_", tag)
      rich_delta_pct <- (rich_F_common - rich_P_common) * 100
      rich_delta_pct <- rich_delta_pct / terra::ifel(rich_P_common > 0, rich_P_common, NA)
      names(rich_delta_pct) <- paste0("richness_delta_pct_", tag)
      gain_stack <- (F_common > 0) & (P_common == 0)
      loss_stack <- (P_common > 0) & (F_common == 0)
      gain_count <- sum(gain_stack); names(gain_count) <- paste0("richness_gain_only_", tag)
      loss_count <- sum(loss_stack); names(loss_count) <- paste0("richness_loss_only_", tag)
      turnover <- gain_count + loss_count; names(turnover) <- paste0("richness_turnover_", tag)
      
      terra::writeRaster(rich_delta_abs, file.path(out_dir, paste0("richness_change_abs_common_", tag, ".tif")),
                         overwrite=TRUE, wopt=list(datatype="INT2S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER")))
      terra::writeRaster(rich_delta_pct, file.path(out_dir, paste0("richness_change_pct_common_", tag, ".tif")),
                         overwrite=TRUE, wopt=list(datatype="FLT4S", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER")))
      terra::writeRaster(gain_count, file.path(out_dir, paste0("richness_gain_only_common_", tag, ".tif")),
                         overwrite=TRUE, wopt=list(datatype="INT2U", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER")))
      terra::writeRaster(loss_count, file.path(out_dir, paste0("richness_loss_only_common_", tag, ".tif")),
                         overwrite=TRUE, wopt=list(datatype="INT2U", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER")))
      terra::writeRaster(turnover, file.path(out_dir, paste0("richness_turnover_common_", tag, ".tif")),
                         overwrite=TRUE, wopt=list(datatype="INT2U", gdal=c("COMPRESS=LZW","BIGTIFF=IF_SAFER")))
      readr::write_csv(data.frame(species = common_sp),
                       file.path(out_dir, paste0("richness_change_common_species_", tag, ".csv")))
      cat("STEP4 ✓ Cambio de riqueza (", tag, ") guardado.\n", sep="")
    } else {
      message("STEP4 ⚠ Sin especies comunes entre presente y futuro para ", tag)
    }
  }
} else {
  cat("STEP4 ⏭ Saltado\n")
}

cat("\n✔ Flujo SDM completo (master) finalizado.\n")
