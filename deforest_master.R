# ============================================================
# deforest_master.R  —  Pipeline completo (Steps 1–6)
# Autor: tú :)
# Descripción:
#   (1) Reclasificación MapBiomas → bosque binario (1/0)
#   (2) Distancia a borde bosque–no bosque (m)
#   (3) Máscaras de transición (deforestación 1→0; bosque estable 1→1)
#   (4) Muestreo supervisado desde transiciones + estáticos
#   (5) XGBoost: entrenamiento, validación, prueba y exportes
#   (6) Inferencia a rejilla (probabilidad y binario; acumulada a 2040)
# Convenciones:
#   1 = bosque, 0 = no-bosque; NA = sin dato
# ============================================================

suppressPackageStartupMessages({
  need <- c("terra","rnaturalearth","data.table","dplyr","xgboost","pROC","PRROC")
  miss <- setdiff(need, rownames(installed.packages()))
  if (length(miss)) install.packages(miss, dependencies = TRUE)
  lapply(need, require, character.only = TRUE)
})

# ---------------- PARAMS (EDITA AQUÍ) ----------------
root        <- "C:/Users/henry.garcia/Desktop/Primates/predictors/DefRisk"
in_dir_mb   <- "C:/Users/henry.garcia/Desktop/Primates/mapbiomas"   # classification_YYYY.tif
crs_metric  <- "EPSG:6933"    # CRS métrico para distancias si el input está en lon/lat

# Ventanas de transición (train/valid/test)
win_train   <- list(c(2000,2005), c(2005,2010), c(2010,2015))
win_valid   <- c(2015,2020)
win_test    <- c(2020,2023)

# Tamaños de muestra por clase (por ventana)
npos <- 15000; nneg <- 15000

# Rásteres estáticos para muestreo (STEP 4)
dir_static  <- file.path(root, "rasters", "static")
static_paths <- file.path(dir_static, c(
  "elev.tif","slope.tif","dist_runap_crsfix.tif",
  "dist_road_crsfix.tif","dist_river_crsfix.tif","dist_town_crsfix.tif",
  "soil_C.tif"
))
stopifnot(all(file.exists(static_paths)))

# Inferencia (STEP 6)
base_year       <- 2023
res_tag         <- "90m"
window_years    <- 3      # p.ej., 3 años si valid/test fueron 2020–2023
target_year_cum <- 2040   # acumulada hasta 2040 (opcional)
# Capas de inferencia — NOMBRES deben coincidir con las features del modelo
paths_infer <- c(
  elev              = file.path(root, "rasters/static_90m/elev_90m.tif"),
  slope             = file.path(root, "rasters/static_90m/slope_90m.tif"),
  dist_runap_crsfix = file.path(root, "rasters/static_90m/dist_runap_crsfix_90m.tif"),
  dist_road_crsfix  = file.path(root, "rasters/static_90m/dist_road_crsfix_90m.tif"),
  dist_river_crsfix = file.path(root, "rasters/static_90m/dist_river_crsfix_90m.tif"),
  dist_town_crsfix  = file.path(root, "rasters/static_90m/dist_town_crsfix_90m.tif"),
  soil_C            = file.path(root, "rasters/static_90m/soil_C_90m.tif"),
  dist_edge         = file.path(root, sprintf("rasters/dist_edge_90m/dist_edge_%d_90m.tif", base_year))
)
stopifnot(all(file.exists(paths_infer)))
forest_mask_file <- file.path(root, "rasters", "forest01", sprintf("forest01_%d_%s.tif", base_year, res_tag))
use_forest_mask  <- file.exists(forest_mask_file)

# Control de ejecución
RUN_STEP1_RECLASS  <- TRUE
RUN_STEP2_DISTEDGE <- TRUE
RUN_STEP3_TRANS    <- TRUE
RUN_STEP4_SAMPLE   <- TRUE
RUN_STEP5_XGB      <- TRUE
RUN_STEP6_INFER    <- TRUE
# -----------------------------------------------------

# ===== Opciones IO / rendimiento
Sys.setenv(GDAL_CACHEMAX = "1024")
terraOptions(todisk = TRUE, memfrac = 0.7, progress = 2)
# terraOptions(tmpdir = "D:/terra_tmp")  # opcional (SSD)

# ===== Carpetas de trabajo
dir_forest    <- file.path(root, "rasters", "forest01")
dir_distedge  <- file.path(root, "rasters", "dist_edge")
dir_trans     <- file.path(root, "rasters", "transitions")
dir_samples   <- file.path(root, "samples")
dir_models    <- file.path(root, "models")
dir_riskmaps  <- file.path(root, "risk_maps")
for (d in c(dir_forest, dir_distedge, dir_trans, dir_samples, dir_models, dir_riskmaps)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

log <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse=" ")))
write_tif <- function(x, path, wopt=list(), overwrite=TRUE){
  writeRaster(x, path, overwrite=overwrite,
              wopt = modifyList(list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER")), wopt))
}

# ============================================================
# STEP 1 — Reclasificación MapBiomas → bosque binario (1/0)
# ============================================================
if (RUN_STEP1_RECLASS) {
  log("STEP 1: Reclasificación MapBiomas → bosque 1/0")
  # AOI Colombia
  co <- rnaturalearth::ne_countries(country = "Colombia", scale = "medium", returnclass = "sf") |> terra::vect()
  # Buscar archivos
  fs  <- list.files(in_dir_mb, pattern = "classification_\\d{4}\\.tif$", full.names = TRUE)
  stopifnot("No se encontraron classification_YYYY.tif" = length(fs) > 0)
  yrs <- as.integer(gsub(".*classification_(\\d{4})\\.tif$", "\\1", fs))
  ord <- order(yrs); fs <- fs[ord]; yrs <- yrs[ord]
  
  reclass_to_forest01 <- function(r) {
    r1 <- ifel((r == 0) | (r == 27), NA, r)                           # 0/27 → NA
    r2 <- ifel((r1 == 3) | (r1 == 5) | (r1 == 6) | (r1 == 49), 1, r1) # bosque → 1
    ifel(is.na(r2), NA, ifel(r2 == 1, 1, 0))                          # resto → 0
  }
  
  summary_list <- vector("list", length(fs))
  for (i in seq_along(fs)) {
    log(sprintf("[%d/%d] %s", i, length(fs), basename(fs[i])))
    r <- rast(fs[i])
    co_proj <- project(co, crs(r))
    r_crop  <- crop(r, co_proj, snap = "out")
    r_mask  <- mask(r_crop, co_proj)
    rb <- reclass_to_forest01(r_mask)
    out_i <- file.path(dir_forest, sprintf("forest01_%d.tif", yrs[i]))
    write_tif(rb, out_i, wopt=list(datatype="INT1U", NAflag=255))
    # resumen
    n_total <- ncell(rb)
    n_na    <- global(is.na(rb), "sum", na.rm=TRUE)[1,1]
    n_ones  <- global(rb == 1, "sum", na.rm=TRUE)[1,1]
    n_zeros <- global(rb == 0, "sum", na.rm=TRUE)[1,1]
    summary_list[[i]] <- data.frame(
      year=yrs[i], zeros=n_zeros, ones=n_ones, nas=n_na,
      pct_0=round(100*n_zeros/n_total,4),
      pct_1=round(100*n_ones /n_total,4),
      pct_na=round(100*n_na  /n_total,4)
    )
    rm(r, r_crop, r_mask, rb); gc()
  }
  summary_df <- do.call(rbind, summary_list)
  data.table::fwrite(summary_df, file.path(dir_forest, "forest_bin_summary_by_year.csv"))
  log("STEP 1 ✅")
} else log("STEP 1 ⏭")

# ============================================================
# STEP 2 — Distancia a borde (metros) por año
# ============================================================
if (RUN_STEP2_DISTEDGE) {
  log("STEP 2: Distancia a borde bosque–no bosque (m)")
  fs_forest <- list.files(dir_forest, pattern="^forest01_\\d{4}\\.tif$", full.names=TRUE)
  stopifnot("No hay forest01_YYYY.tif" = length(fs_forest) > 0)
  yrs <- as.integer(gsub(".*forest01_(\\d{4})\\.tif$", "\\1", fs_forest))
  ord <- order(yrs); fs_forest <- fs_forest[ord]; yrs <- yrs[ord]
  
  for (i in seq_along(fs_forest)) {
    in_path <- fs_forest[i]
    yr      <- yrs[i]
    out_path <- file.path(dir_distedge, sprintf("dist_edge_%d.tif", yr))
    if (file.exists(out_path)) { log("• Ya existe:", basename(out_path)); next }
    r <- rast(in_path)
    # reproyección temporal si lon/lat
    need_reproj <- is.lonlat(r)
    r_m <- if (need_reproj) project(r, crs_metric, method="near") else r
    # fuente de distancia = NO-bosque (0→1); bosque (1→NA)
    src_nonforest <- classify(r_m, rcl = rbind(c(0,1), c(1,NA)), include.lowest = TRUE)
    d_all <- distance(src_nonforest)
    forest_mask <- classify(r_m, rcl = rbind(c(1,1), c(0,NA)), include.lowest = TRUE)
    dist_edge_m <- d_all * forest_mask; names(dist_edge_m) <- "dist_edge_m"
    # devolver al grid original si reproyectaste
    if (need_reproj) dist_edge_m <- project(dist_edge_m, r, method = "bilinear")
    write_tif(dist_edge_m, out_path, wopt=list(datatype="FLT4S"))
    log("✔ Guardado:", out_path)
    rm(r, r_m, src_nonforest, d_all, forest_mask, dist_edge_m); gc()
  }
  log("STEP 2 ✅")
} else log("STEP 2 ⏭")

# ============================================================
# STEP 3 — Máscaras de transición (defor 1→0; estable 1→1)
# ============================================================
if (RUN_STEP3_TRANS) {
  log("STEP 3: Transiciones 1→0 y 1→1")
  write_transition_masks <- function(t_prev, t_curr) {
    f_prev <- rast(file.path(dir_forest,   sprintf("forest01_%d.tif", t_prev)))
    f_curr <- rast(file.path(dir_forest,   sprintf("forest01_%d.tif", t_curr)))
    de_prev<- rast(file.path(dir_distedge, sprintf("dist_edge_%d.tif", t_prev)))
    if (!same.crs(f_curr, f_prev)) f_curr <- project(f_curr, y=f_prev, method="near")
    f_curr  <- resample(f_curr, f_prev, method="near")
    if (!same.crs(de_prev, f_prev)) de_prev <- project(de_prev, y=f_prev, method="bilinear")
    de_prev <- resample(de_prev, f_prev, method="bilinear")
    y_pos_r <- (f_prev == 1) & (f_curr == 0)
    y_neg_r <- (f_prev == 1) & (f_curr == 1)
    valid_de <- !is.na(de_prev)
    pos_mask <- ifel(y_pos_r & valid_de, 1, NA)
    neg_mask <- ifel(y_neg_r & valid_de, 1, NA)
    wp <- list(gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"))
    writeRaster(pos_mask, file.path(dir_trans, sprintf("trans_pos_%d_%d.tif", t_prev, t_curr)), overwrite=TRUE, wopt=wp)
    writeRaster(neg_mask, file.path(dir_trans, sprintf("trans_neg_%d_%d.tif", t_prev, t_curr)), overwrite=TRUE, wopt=wp)
    log("✔ Trans:", t_prev, "->", t_curr)
  }
  # construir desde forest disponibles
  years <- sort(as.integer(gsub(".*forest01_(\\d{4})\\.tif$", "\\1", list.files(dir_forest, full.names = TRUE))))
  windows_all <- c(lapply(win_train, identity), list(win_valid), list(win_test))
  for (w in windows_all) write_transition_masks(w[1], w[2])
  log("STEP 3 ✅")
} else log("STEP 3 ⏭")

# ============================================================
# STEP 4 — Muestreo desde transiciones + estáticos
# ============================================================
if (RUN_STEP4_SAMPLE) {
  log("STEP 4: Muestreo supervisado")
  
  save_samples <- function(dt, out_csv) {
    if ("dist_edge" %in% names(dt)) dt <- dt[!is.na(dist_edge)]
    data.table::fwrite(dt, out_csv); log("Guardado:", out_csv, "| n =", nrow(dt))
  }
  sample_cells_mask <- function(mask, n, seed = 42, batch = NULL, max_tries = 50L) {
    set.seed(seed)
    n_total <- as.integer(global(!is.na(mask), "sum", na.rm = TRUE)[1,1])
    if (!is.finite(n_total) || n_total == 0) return(integer(0))
    n <- min(n, n_total); if (n == 0) return(integer(0))
    if (is.null(batch)) batch <- max(2L * n, 10000L)
    picked <- integer(0); tries <- 0L; nc <- ncell(mask)
    while (length(picked) < n && tries < max_tries) {
      tries <- tries + 1L
      cand <- sample.int(nc, size = batch, replace = FALSE)
      xy   <- xyFromCell(mask, cand)
      val  <- tryCatch(terra::extract(mask, vect(xy, crs = crs(mask)), ID = FALSE)[,1],
                       error = function(e) rep(NA_real_, nrow(xy)))
      ok   <- !is.na(val)
      if (any(ok)) {
        good_cells <- cand[ok]
        picked <- unique(c(picked, good_cells))
        if (length(picked) > n) picked <- picked[seq_len(n)]
      }
    }
    picked
  }
  extract_static_by_pts <- function(paths, pts_in_ref_crs, method_map = NULL) {
    out <- vector("list", length(paths))
    nm  <- gsub("\\.tif$", "", basename(paths))
    for (i in seq_along(paths)) {
      r <- rast(paths[i])
      pts_curr <- if (!same.crs(r, pts_in_ref_crs)) project(pts_in_ref_crs, crs(r)) else pts_in_ref_crs
      m <- if (!is.null(method_map) && nm[i] %in% names(method_map)) method_map[[nm[i]]] else "bilinear"
      out[[i]] <- tryCatch(extract(r, pts_curr, ID = FALSE, method = m)[,1],
                           error = function(e) rep(NA_real_, nrow(pts_curr)))
    }
    names(out) <- nm; as.data.frame(out)
  }
  sample_from_transition <- function(t_prev, t_curr, out_csv, npos=10000, nneg=10000, seed=42) {
    pos_path <- file.path(dir_trans,    sprintf("trans_pos_%d_%d.tif", t_prev, t_curr))
    neg_path <- file.path(dir_trans,    sprintf("trans_neg_%d_%d.tif", t_prev, t_curr))
    de_path  <- file.path(dir_distedge, sprintf("dist_edge_%d.tif",     t_prev))
    stopifnot(file.exists(pos_path), file.exists(neg_path), file.exists(de_path))
    pos_mask <- rast(pos_path); neg_mask <- rast(neg_path); de_prev  <- rast(de_path)
    if (!same.crs(de_prev, pos_mask)) de_prev <- project(de_prev, y = pos_mask, method = "bilinear")
    de_prev <- resample(de_prev, pos_mask, method = "bilinear")
    npos_avail <- as.integer(global(!is.na(pos_mask), "sum", na.rm=TRUE)[1,1])
    nneg_avail <- as.integer(global(!is.na(neg_mask), "sum", na.rm=TRUE)[1,1])
    stopifnot(npos_avail > 0, nneg_avail > 0)
    npos_eff <- min(npos, npos_avail); nneg_eff <- min(nneg, nneg_avail)
    pos_cells <- sample_cells_mask(pos_mask, npos_eff, seed = seed,      batch = max(4L*npos_eff, 20000L))
    neg_cells <- sample_cells_mask(neg_mask, nneg_eff, seed = seed + 17, batch = max(4L*nneg_eff, 20000L))
    stopifnot(length(pos_cells) > 0, length(neg_cells) > 0)
    xy_pos <- xyFromCell(pos_mask, pos_cells); pts_pos <- vect(xy_pos, crs = crs(pos_mask))
    xy_neg <- xyFromCell(neg_mask, neg_cells); pts_neg <- vect(xy_neg, crs = crs(neg_mask))
    de_pos <- extract(de_prev, pts_pos, ID = FALSE)[,1]
    de_neg <- extract(de_prev, pts_neg, ID = FALSE)[,1]
    method_map <- NULL
    Xs_pos <- extract_static_by_pts(static_paths, pts_pos, method_map)
    Xs_neg <- extract_static_by_pts(static_paths, pts_neg, method_map)
    dt_pos <- data.table::as.data.table(Xs_pos); dt_pos[, dist_edge := de_pos][, y := 1L][, `:=`(t_prev=t_prev, t_curr=t_curr)]
    dt_neg <- data.table::as.data.table(Xs_neg); dt_neg[, dist_edge := de_neg][, y := 0L][, `:=`(t_prev=t_prev, t_curr=t_curr)]
    dt <- data.table::rbindlist(list(dt_pos, dt_neg), use.names = TRUE, fill = TRUE)
    save_samples(dt, out_csv); invisible(TRUE)
  }
  
  # Ejecutar
  summary_rows <- list()
  for (w in win_train) {
    out_csv <- file.path(dir_samples, sprintf("train_%d_%d.csv", w[1], w[2]))
    if (!file.exists(out_csv)) try(sample_from_transition(w[1], w[2], out_csv, npos, nneg), silent = TRUE)
    if (file.exists(out_csv)) {
      df <- data.table::fread(out_csv, nThread = 1L)
      summary_rows[[length(summary_rows)+1]] <- data.frame(split="train", t_prev=w[1], t_curr=w[2],
                                                           n=nrow(df), pos=sum(df$y==1L), neg=sum(df$y==0L))
    }
  }
  valid_csv <- file.path(dir_samples, sprintf("valid_%d_%d.csv", win_valid[1], win_valid[2]))
  if (!file.exists(valid_csv)) try(sample_from_transition(win_valid[1], win_valid[2], valid_csv, npos, nneg), silent = TRUE)
  if (file.exists(valid_csv)) {
    df <- data.table::fread(valid_csv, nThread = 1L)
    summary_rows[[length(summary_rows)+1]] <- data.frame(split="valid", t_prev=win_valid[1], t_curr=win_valid[2],
                                                         n=nrow(df), pos=sum(df$y==1L), neg=sum(df$y==0L))
  }
  test_csv <- file.path(dir_samples, sprintf("test_%d_%d.csv", win_test[1], win_test[2]))
  if (!file.exists(test_csv)) try(sample_from_transition(win_test[1], win_test[2], test_csv, npos, nneg), silent = TRUE)
  if (file.exists(test_csv)) {
    df <- data.table::fread(test_csv, nThread = 1L)
    summary_rows[[length(summary_rows)+1]] <- data.frame(split="test", t_prev=win_test[1], t_curr=win_test[2],
                                                         n=nrow(df), pos=sum(df$y==1L), neg=sum(df$y==0L))
  }
  if (length(summary_rows)) {
    sampling_summary <- dplyr::bind_rows(summary_rows)
    data.table::fwrite(sampling_summary, file.path(dir_samples, "sampling_summary.csv"))
    log("Resumen -> samples/sampling_summary.csv")
  }
  log("STEP 4 ✅")
} else log("STEP 4 ⏭")

# ============================================================
# STEP 5 — XGBoost: entrenamiento/evaluación/exportes
# ============================================================
if (RUN_STEP5_XGB) {
  log("STEP 5: XGBoost")
  fpath <- function(tag) file.path(dir_samples, tag)
  read_split <- function(path) {
    dt <- data.table::fread(path); stopifnot("y" %in% names(dt))
    dt[, y := as.integer(y)]
    drop <- intersect(c("t_prev","t_curr"), names(dt)); if (length(drop)) dt[, (drop) := NULL]
    feats <- setdiff(names(dt), "y"); for (v in feats) dt[[v]] <- as.numeric(dt[[v]])
    dt
  }
  # Cargar splits
  train_tags <- c("train_2000_2005.csv","train_2005_2010.csv","train_2010_2015.csv")
  valid_tag  <- "valid_2015_2020.csv"
  test_tag   <- "test_2020_2023.csv"
  dtrain <- data.table::rbindlist(lapply(train_tags, \(t) read_split(fpath(t))), use.names=TRUE, fill=TRUE)
  dvalid <- read_split(fpath(valid_tag))
  dtest  <- read_split(fpath(test_tag))
  log("Train:", nrow(dtrain), "| Valid:", nrow(dvalid), "| Test:", nrow(dtest))
  features <- Reduce(intersect, list(setdiff(names(dtrain),"y"), setdiff(names(dvalid),"y"), setdiff(names(dtest),"y")))
  stopifnot(length(features) > 0); features <- sort(features)
  # Matrices
  Xtr <- as.matrix(dtrain[, ..features]); ytr <- dtrain$y
  Xva <- as.matrix(dvalid[, ..features]); yva <- dvalid$y
  Xte <- as.matrix(dtest[,  ..features]); yte <- dtest$y
  # Pesos por desbalance
  pos_ratio <- sum(ytr == 0) / max(1, sum(ytr == 1))
  spw <- if (abs(pos_ratio - 1) < 0.25) 1 else pos_ratio
  params <- list(
    objective = "binary:logistic", eval_metric = "auc",
    eta = 0.05, max_depth = 6, min_child_weight = 5,
    subsample = 0.8, colsample_bytree = 0.8,
    lambda = 1.0, alpha = 0.0, scale_pos_weight = spw
  )
  dtr <- xgb.DMatrix(Xtr, label=ytr, missing=NA)
  dva <- xgb.DMatrix(Xva, label=yva, missing=NA)
  dte <- xgb.DMatrix(Xte, label=yte, missing=NA)
  # Entrenamiento con early stopping
  set.seed(42)
  bst <- xgb.train(params=params, data=dtr, nrounds=5000,
                   watchlist=list(train=dtr, eval=dva),
                   early_stopping_rounds=100, verbose=1)
  best_iter <- bst$best_iteration; best_score <- bst$best_score
  log("best_iteration =", best_iter, "| best_valid_AUC =", round(best_score, 4))
  # Predicciones / métricas
  pred_va <- predict(bst, dva); pred_te <- predict(bst, dte)
  auc_va <- as.numeric(pROC::auc(response=yva, predictor=pred_va))
  auc_te <- as.numeric(pROC::auc(response=yte, predictor=pred_te))
  pra_va <- PRROC::pr.curve(scores.class0=pred_va[yva==1], scores.class1=pred_va[yva==0])$auc.integral
  pra_te <- PRROC::pr.curve(scores.class0=pred_te[yte==1], scores.class1=pred_te[yte==0])$auc.integral
  brier  <- function(p, y) mean((p - y)^2)
  brier_va <- brier(pred_va, yva); brier_te <- brier(pred_te, yte)
  # Umbral óptimo (F1) en valid
  f1_at_thr <- function(p, y, thr) {
    yhat <- as.integer(p >= thr)
    tp <- sum(yhat==1 & y==1); fp <- sum(yhat==1 & y==0); fn <- sum(yhat==0 & y==1)
    prec <- if ((tp+fp)==0) 0 else tp/(tp+fp); rec <- if ((tp+fn)==0) 0 else tp/(tp+fn)
    if (prec+rec==0) return(0); 2*prec*rec/(prec+rec)
  }
  thr_grid <- seq(0.05, 0.95, 0.01)
  thr_opt  <- thr_grid[ which.max(sapply(thr_grid, \(t) f1_at_thr(pred_va, yva, t))) ]
  # Métricas de test con thr_opt
  yhat_te <- as.integer(pred_te >= thr_opt)
  tp <- sum(yhat_te==1 & yte==1); fp <- sum(yhat_te==1 & yte==0)
  tn <- sum(yhat_te==0 & yte==0); fn <- sum(yhat_te==0 & yte==1)
  precision <- if ((tp+fp)==0) 0 else tp/(tp+fp)
  recall    <- if ((tp+fn)==0) 0 else tp/(tp+fn)
  f1_te     <- if ((precision+recall)==0) 0 else 2*precision*recall/(precision+recall)
  acc_te    <- (tp+tn)/max(1,(tp+tn+fp+fn))
  log(sprintf("VALID | AUC=%.4f  PR-AUC=%.4f  Brier=%.4f", auc_va, pra_va, brier_va))
  log(sprintf("TEST  | AUC=%.4f  PR-AUC=%.4f  Brier=%.4f  Acc=%.3f  F1=%.3f (thr=%.2f)",
              auc_te, pra_te, brier_te, acc_te, f1_te, thr_opt))
  # Exportes
  saveRDS(bst, file.path(dir_models, "xgb_defrisk_model.rds"))
  data.table::fwrite(data.frame(threshold = thr_opt), file.path(dir_models, "xgb_threshold.csv"))
  imp <- xgb.importance(model=bst, feature_names=features)
  data.table::fwrite(imp, file.path(dir_models, "feature_importance.csv"))
  metrics_df <- data.frame(
    split=c("valid","test"), auc=c(auc_va,auc_te), pr_auc=c(pra_va,pra_te),
    brier=c(brier_va,brier_te), acc=c(NA,acc_te), f1=c(NA,f1_te), thr=c(thr_opt,thr_opt),
    tp=c(NA,tp), fp=c(NA,fp), tn=c(NA,tn), fn=c(NA,fn), precision=c(NA,precision),
    recall=c(NA,recall)
  )
  data.table::fwrite(metrics_df, file.path(dir_models, "metrics_valid_test.csv"))
  # Reentrenar final (train+valid) con nrounds = best_iter
  dtrain_full <- data.table::rbindlist(list(dtrain, dvalid), use.names=TRUE, fill=TRUE)
  Xtr_full <- as.matrix(dtrain_full[, ..features]); ytr_full <- dtrain_full$y
  dtr_full <- xgb.DMatrix(Xtr_full, label=ytr_full, missing=NA)
  set.seed(42)
  bst_final <- xgb.train(params=params, data=dtr_full, nrounds=best_iter, verbose=1)
  saveRDS(bst_final, file.path(dir_models, "xgb_defrisk_model_final.rds"))
  data.table::fwrite(data.frame(threshold = thr_opt), file.path(dir_models, "xgb_threshold_final.csv"))
  imp_final <- xgb.importance(model=bst_final, feature_names=features)
  data.table::fwrite(imp_final, file.path(dir_models, "feature_importance_final.csv"))
  log("STEP 5 ✅  (modelo y umbral en 'models/')")
  
} else log("STEP 5 ⏭")

# ============================================================
# STEP 6 — Inferencia a rejilla (probabilidad/binario/acumulada)
# ============================================================
if (RUN_STEP6_INFER) {
  log("STEP 6: Inferencia a rejilla")
  model_path <- file.path(dir_models, "xgb_defrisk_model_final.rds")
  stopifnot(file.exists(model_path))
  bst_final <- readRDS(model_path)
  thr_file_final <- file.path(dir_models, "xgb_threshold_final.csv")
  thr_file_fbk   <- file.path(dir_models, "xgb_threshold.csv")
  thr_opt <- if (file.exists(thr_file_final)) {
    as.numeric(data.table::fread(thr_file_final)$threshold[1])
  } else if (file.exists(thr_file_fbk)) {
    as.numeric(data.table::fread(thr_file_fbk)$threshold[1])
  } else 0.5
  if (!is.finite(thr_opt)) thr_opt <- 0.5
  # Validar features esperadas
  fi_path <- file.path(dir_models, "feature_importance_final.csv")
  if (!file.exists(fi_path)) fi_path <- file.path(dir_models, "feature_importance.csv")
  features_expected <- if (file.exists(fi_path)) unique(data.table::fread(fi_path)$Feature) else names(paths_infer)
  if (!setequal(names(paths_infer), features_expected)) {
    missing_in_paths <- setdiff(features_expected, names(paths_infer))
    extra_in_paths   <- setdiff(names(paths_infer), features_expected)
    stop("Desajuste de features para inferencia.\nFaltan en 'paths_infer': ", paste(missing_in_paths, collapse=", "),
         "\nSobran en 'paths_infer': ", paste(extra_in_paths, collapse=", "))
  }
  paths_infer <- paths_infer[features_expected]  # orden exacto
  # Construir stack predictor
  pred_stack <- rast(unname(paths_infer)); names(pred_stack) <- names(paths_infer)
  if (use_forest_mask) {
    tpl <- rast(forest_mask_file)
    if (!identical(crs(tpl), crs(pred_stack))) tpl <- project(tpl, pred_stack, method="near")
    tpl <- resample(tpl, pred_stack, method="near")
    pred_stack <- mask(pred_stack, ifel(tpl==1, 1, NA))
  }
  # Predict por bloques
  xgb_fun <- function(model, data, ...) { data <- as.matrix(data); as.numeric(predict(model, data)) }
  wp_prob <- list(datatype = "FLT4S",
                  gdal = c("COMPRESS=ZSTD","NUM_THREADS=ALL_CPUS","TILED=YES","BIGTIFF=IF_SAFER","BLOCKXSIZE=512","BLOCKYSIZE=512"))
  wp_bin  <- list(datatype = "INT1U", gdal = c("COMPRESS=LZW","TILED=YES","BIGTIFF=IF_SAFER"))
  prob_tif <- file.path(dir_riskmaps, sprintf("defrisk_prob_%d_%s.tif", base_year, res_tag))
  bin_tif  <- file.path(dir_riskmaps, sprintf("defrisk_bin_%d_%s_thr.tif",  base_year, res_tag))
  log("→ Prediciendo probabilidad…")
  prob <- terra::predict(pred_stack, model=bst_final, fun=xgb_fun,
                         filename=prob_tif, overwrite=TRUE, wopt=wp_prob, progress=1)
  names(prob) <- "risk_prob"; log("✔", prob_tif)
  # Binario
  risk_bin <- ifel(prob >= thr_opt, 1, ifel(is.na(prob), NA, 0))
  writeRaster(risk_bin, bin_tif, overwrite=TRUE, wopt=wp_bin); log("✔", bin_tif)
  # QA
  n_total <- ncell(prob)
  n_na    <- global(is.na(prob), "sum", na.rm = TRUE)[1,1]
  n_ones  <- global(risk_bin == 1, "sum", na.rm = TRUE)[1,1]
  n_zeros <- global(risk_bin == 0, "sum", na.rm = TRUE)[1,1]
  p_mean  <- global(prob, "mean", na.rm = TRUE)[1,1]
  p_med   <- global(prob, "median", na.rm = TRUE)[1,1]
  p_q90   <- global(prob, "quantile", p=0.90, na.rm = TRUE)[1,1]
  qa <- data.frame(base_year=base_year, res=res_tag, n_total=n_total, n_na=n_na, n_ones=n_ones, n_zeros=n_zeros,
                   prob_mean=round(p_mean,6), prob_median=round(p_med,6), prob_q90=round(p_q90,6), thr=thr_opt)
  data.table::fwrite(qa, file.path(dir_riskmaps, sprintf("inference_summary_%d.csv", base_year)))
  # Acumulada hasta target_year_cum
  if (!is.null(target_year_cum) && target_year_cum > base_year) {
    n_steps <- ceiling((target_year_cum - base_year) / window_years)
    log(sprintf("→ Probabilidad acumulada hasta %d (n=%d pasos de %d años)…", target_year_cum, n_steps, window_years))
    prob_cum <- 1 - (1 - prob)^n_steps; names(prob_cum) <- "risk_prob_cum"
    prob_cum_tif <- file.path(dir_riskmaps, sprintf("defrisk_prob_cum_to_%d_%s.tif", target_year_cum, res_tag))
    writeRaster(prob_cum, prob_cum_tif, overwrite=TRUE, wopt=wp_prob)
    bin_cum_tif <- file.path(dir_riskmaps, sprintf("defrisk_bin_cum_to_%d_%s_thr.tif", target_year_cum, res_tag))
    risk_bin_cum <- ifel(prob_cum >= thr_opt, 1, ifel(is.na(prob_cum), NA, 0))
    writeRaster(risk_bin_cum, bin_cum_tif, overwrite=TRUE, wopt=wp_bin)
    log("✔", prob_cum_tif); log("✔", bin_cum_tif)
  }
  log("STEP 6 ✅")
} else log("STEP 6 ⏭")

log("✓ Pipeline completado.")
