# ==============================================================================
# deforisk_master Ajuste y evaluación temporal del modelo XGBoost
# ==============================================================================
#
# PROPÓSITO
# ---------
# Ajustar un modelo de clasificación XGBoost para estimar el riesgo relativo de
# deforestación, evaluar su transferencia temporal y guardar el modelo final que
# posteriormente podrá aplicarse sobre la cobertura de bosque del año base.
#
# Este script corresponde a una versión independiente, revisada y documentada
# del bloque "STEP 5 — XGBoost" del script maestro `master_deforisk`.
#
# QUÉ NECESITA
# ------------
# 1. R y los paquetes:
#      data.table, xgboost, pROC y PRROC.
#
# 2. Una carpeta de muestras con un archivo CSV por ventana temporal:
#
#      samples/train_2000_2005.csv
#      samples/train_2005_2010.csv
#      samples/train_2010_2015.csv
#      samples/valid_2015_2020.csv
#      samples/test_2020_2024.csv
#
#    Los nombres se construyen automáticamente a partir de los años declarados
#    en la sección PARÁMETROS EDITABLES.
#
# 3. Cada CSV debe contener:
#      - `y`: variable respuesta binaria;
#             1 = transición de bosque a no bosque (deforestación)
#             0 = permanencia de bosque (bosque estable)
#      - una columna numérica por predictor;
#      - opcionalmente `t_prev` y `t_curr`, que se eliminan antes del ajuste.
#
#    Predictores esperados por defecto:
#      elevation, soil_carbon, dist_river, dist_road,
#      dist_edge, dist_town y dist_runap.
#
#    Estos nombres pueden modificarse en `expected_features`, siempre que sean
#    exactamente iguales a los nombres de las columnas presentes en los CSV.
#
# SUPUESTOS METODOLÓGICOS
# -----------------------
# - Las muestras siguen un diseño caso-control aproximadamente balanceado.
# - Por esta razón, la salida `binary:logistic` se interpreta como un ÍNDICE DE
#   RIESGO RELATIVO y no como una probabilidad absoluta calibrada respecto a la
#   prevalencia nacional de la deforestación.
# - La selección del número de iteraciones se realiza únicamente con el conjunto
#   de validación mediante detención temprana sobre ROC-AUC.
# - El umbral que maximiza F1 se selecciona únicamente con validación.
# - La prueba temporal independiente no interviene en el ajuste, la selección de
#   iteraciones ni la selección del umbral.
# - Después de la evaluación, el modelo final se reajusta con entrenamiento más
#   validación, usando el número de iteraciones previamente seleccionado.
#
# QUÉ HACE
# --------
# 1. Lee y verifica las muestras de entrenamiento, validación y prueba.
# 2. Comprueba la presencia de la respuesta y los predictores requeridos.
# 3. Une las ventanas históricas destinadas al entrenamiento.
# 4. Ajusta XGBoost con validación temporal y detención temprana.
# 5. Calcula ROC-AUC, PR-AUC y Brier score en validación y prueba.
# 6. Selecciona en validación el umbral que maximiza F1.
# 7. Aplica ese umbral a la prueba temporal independiente.
# 8. Exporta métricas, umbral, importancia de variables y modelos.
# 9. Reajusta el modelo final con entrenamiento + validación.
#
# QUÉ PRODUCE
# -----------
# En la carpeta `models/`:
#
#   xgb_defrisk_model.rds
#       Modelo utilizado para evaluación temporal. Conserva el ajuste realizado
#       con entrenamiento y validación empleada como watchlist.
#
#   xgb_threshold.csv
#       Umbral seleccionado exclusivamente en validación mediante F1.
#
#   metrics_valid_test.csv
#       ROC-AUC, PR-AUC, Brier score y métricas dependientes del umbral.
#
#   feature_importance.csv
#       Importancia de variables del modelo de evaluación según Gain, Cover y
#       Frequency, calculadas por XGBoost.
#
#   xgb_defrisk_model_final.rds
#       Modelo reajustado con entrenamiento + validación. Este es el modelo que
#       debe utilizarse posteriormente para generar el ráster continuo de riesgo.
#
#   xgb_threshold_final.csv
#       Copia del umbral de validación asociada al modelo final. Se conserva para
#       reproducibilidad y evaluación; no controla la demanda de pérdida forestal.
#
#   feature_importance_final.csv
#       Importancia de variables del modelo final.
#
#   model_features.csv
#       Orden exacto de los predictores utilizados. Debe conservarse durante la
#       inferencia espacial.
#
#   model_settings.csv
#       Años, número óptimo de iteraciones, umbral y parámetros del modelo.
#
# ESTE SCRIPT NO REALIZA
# ----------------------
# - preparación de rásteres o variables predictoras;
# - muestreo de píxeles;
# - inferencia espacial sobre toda Colombia;
# - estimación o asignación ecorregional de la pérdida forestal;
# - calibración de probabilidades respecto a la prevalencia nacional.
#
# EJECUCIÓN
# ---------
# Ubique este archivo en la raíz del proyecto o defina la variable de entorno
# `DEFORISK_PROJECT_ROOT`. Luego ejecute:
#
#   source("02_modelo.R")
#
# ==============================================================================

suppressPackageStartupMessages({
  required_packages <- c("data.table", "xgboost", "pROC", "PRROC")
  missing_packages <- setdiff(required_packages, rownames(installed.packages()))

  if (length(missing_packages) > 0L) {
    stop(
      "Faltan los siguientes paquetes: ",
      paste(missing_packages, collapse = ", "),
      ". Instálelos antes de ejecutar el modelo."
    )
  }

  invisible(lapply(required_packages, library, character.only = TRUE))
})

# ==============================================================================
# 1. PARÁMETROS EDITABLES
# ==============================================================================

# Raíz del proyecto. Por defecto se usa el directorio de trabajo actual.
# También puede definirse externamente mediante DEFORISK_PROJECT_ROOT.
project_root <- normalizePath(
  Sys.getenv("DEFORISK_PROJECT_ROOT", unset = getwd()),
  winslash = "/",
  mustWork = TRUE
)

# Carpetas relativas de entrada y salida.
dir_samples <- file.path(project_root, "samples")
dir_models  <- file.path(project_root, "models")
dir.create(dir_models, recursive = TRUE, showWarnings = FALSE)

# Partición temporal.
win_train <- list(
  c(2000L, 2005L),
  c(2005L, 2010L),
  c(2010L, 2015L)
)
win_valid <- c(2015L, 2020L)
win_test  <- c(2020L, 2024L)

# Nombres exactos de las variables predictoras en los archivos CSV.
expected_features <- c(
  "elevation",
  "soil_carbon",
  "dist_river",
  "dist_road",
  "dist_edge",
  "dist_town",
  "dist_runap"
)

# Reproducibilidad y control del entrenamiento.
seed_model            <- 42L
nrounds_max            <- 5000L
early_stopping_rounds <- 100L
threshold_grid         <- seq(0.05, 0.95, by = 0.01)
verbose_training       <- 1L

# Hiperparámetros establecidos en la metodología.
xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  eta = 0.05,
  max_depth = 6L,
  min_child_weight = 5,
  subsample = 0.8,
  colsample_bytree = 0.8,
  lambda = 1.0,
  alpha = 0.0,
  scale_pos_weight = 1.0
)

# ==============================================================================
# 2. FUNCIONES AUXILIARES SENCILLAS
# ==============================================================================

log_message <- function(...) {
  cat(sprintf(
    "[%s] %s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    paste0(..., collapse = " ")
  ))
}

window_file <- function(prefix, window) {
  file.path(
    dir_samples,
    sprintf("%s_%d_%d.csv", prefix, window[1], window[2])
  )
}

read_sample_file <- function(path, expected_features) {
  if (!file.exists(path)) {
    stop("No se encontró el archivo de muestras: ", path)
  }

  dt <- data.table::fread(path, showProgress = FALSE)

  if (!"y" %in% names(dt)) {
    stop("El archivo no contiene la columna respuesta `y`: ", path)
  }

  optional_columns <- intersect(c("t_prev", "t_curr"), names(dt))
  if (length(optional_columns) > 0L) {
    dt[, (optional_columns) := NULL]
  }

  missing_features <- setdiff(expected_features, names(dt))
  if (length(missing_features) > 0L) {
    stop(
      "Faltan predictores en ", basename(path), ": ",
      paste(missing_features, collapse = ", ")
    )
  }

  # Se conservan únicamente la respuesta y los predictores definidos.
  dt <- dt[, c("y", expected_features), with = FALSE]
  dt[, y := as.integer(y)]

  invalid_response <- setdiff(unique(stats::na.omit(dt$y)), c(0L, 1L))
  if (length(invalid_response) > 0L) {
    stop("La columna `y` contiene valores distintos de 0 y 1 en: ", path)
  }

  if (anyNA(dt$y)) {
    stop("La columna `y` contiene valores NA en: ", path)
  }

  for (feature in expected_features) {
    dt[[feature]] <- as.numeric(dt[[feature]])

    if (all(is.na(dt[[feature]]))) {
      stop(
        "El predictor `", feature,
        "` está completamente vacío en: ", path
      )
    }
  }

  if (!all(c(0L, 1L) %in% unique(dt$y))) {
    stop("El archivo no contiene ambas clases de respuesta: ", path)
  }

  dt
}

brier_score <- function(probability, observed) {
  mean((probability - observed)^2)
}

classification_metrics <- function(probability, observed, threshold) {
  predicted <- as.integer(probability >= threshold)

  tp <- sum(predicted == 1L & observed == 1L)
  fp <- sum(predicted == 1L & observed == 0L)
  tn <- sum(predicted == 0L & observed == 0L)
  fn <- sum(predicted == 0L & observed == 1L)

  precision <- if ((tp + fp) == 0L) NA_real_ else tp / (tp + fp)
  recall    <- if ((tp + fn) == 0L) NA_real_ else tp / (tp + fn)
  accuracy  <- (tp + tn) / length(observed)

  f1 <- if (
    is.na(precision) || is.na(recall) || (precision + recall) == 0
  ) {
    0
  } else {
    2 * precision * recall / (precision + recall)
  }

  data.frame(
    threshold = threshold,
    tp = tp,
    fp = fp,
    tn = tn,
    fn = fn,
    precision = precision,
    recall = recall,
    accuracy = accuracy,
    f1 = f1
  )
}

# ==============================================================================
# 3. IDENTIFICACIÓN Y COMPROBACIÓN DE LOS ARCHIVOS DE MUESTRA
# ==============================================================================

train_files <- vapply(
  win_train,
  function(window) window_file("train", window),
  FUN.VALUE = character(1)
)
valid_file <- window_file("valid", win_valid)
test_file  <- window_file("test", win_test)

required_files <- c(train_files, valid_file, test_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Faltan archivos de muestras:\n",
    paste(missing_files, collapse = "\n")
  )
}

log_message("Archivos de entrada verificados:", length(required_files))

# ==============================================================================
# 4. LECTURA Y ORGANIZACIÓN DE ENTRENAMIENTO, VALIDACIÓN Y PRUEBA
# ==============================================================================

train_list <- lapply(
  train_files,
  read_sample_file,
  expected_features = expected_features
)

dtrain <- data.table::rbindlist(
  train_list,
  use.names = TRUE,
  fill = FALSE
)
dvalid <- read_sample_file(valid_file, expected_features)
dtest  <- read_sample_file(test_file, expected_features)

log_message(
  "Muestras | entrenamiento:", nrow(dtrain),
  "| validación:", nrow(dvalid),
  "| prueba:", nrow(dtest)
)

class_summary <- data.frame(
  split = c("train", "validation", "test"),
  n = c(nrow(dtrain), nrow(dvalid), nrow(dtest)),
  positives = c(sum(dtrain$y == 1L), sum(dvalid$y == 1L), sum(dtest$y == 1L)),
  negatives = c(sum(dtrain$y == 0L), sum(dvalid$y == 0L), sum(dtest$y == 0L))
)
class_summary$positive_fraction <- class_summary$positives / class_summary$n

data.table::fwrite(
  class_summary,
  file.path(dir_models, "sample_class_summary.csv")
)

# El muestreo metodológico es balanceado. Se reporta una advertencia si la
# diferencia entre clases supera 10 %, pero se conserva scale_pos_weight = 1.
if (abs(class_summary$positive_fraction[1] - 0.5) > 0.10) {
  warning(
    "El conjunto de entrenamiento no está aproximadamente balanceado. ",
    "Revise el muestreo; `scale_pos_weight` se mantiene en 1 según el diseño ",
    "metodológico establecido."
  )
}

# ==============================================================================
# 5. CONSTRUCCIÓN DE MATRICES XGBOOST
# ==============================================================================

X_train <- as.matrix(dtrain[, ..expected_features])
y_train <- dtrain$y

X_valid <- as.matrix(dvalid[, ..expected_features])
y_valid <- dvalid$y

X_test <- as.matrix(dtest[, ..expected_features])
y_test <- dtest$y

matrix_train <- xgboost::xgb.DMatrix(X_train, label = y_train, missing = NA)
matrix_valid <- xgboost::xgb.DMatrix(X_valid, label = y_valid, missing = NA)
matrix_test  <- xgboost::xgb.DMatrix(X_test, label = y_test, missing = NA)

# ==============================================================================
# 6. AJUSTE DEL MODELO DE EVALUACIÓN
# ==============================================================================

set.seed(seed_model)

model_evaluation <- xgboost::xgb.train(
  params = xgb_params,
  data = matrix_train,
  nrounds = nrounds_max,
  watchlist = list(
    train = matrix_train,
    validation = matrix_valid
  ),
  early_stopping_rounds = early_stopping_rounds,
  verbose = verbose_training
)

best_iteration <- model_evaluation$best_iteration
best_validation_auc <- model_evaluation$best_score

if (is.null(best_iteration) || !is.finite(best_iteration)) {
  stop("XGBoost no devolvió un número óptimo de iteraciones válido.")
}

log_message(
  "Iteración óptima:", best_iteration,
  "| AUC de validación:", round(best_validation_auc, 4)
)

saveRDS(
  model_evaluation,
  file.path(dir_models, "xgb_defrisk_model.rds")
)

# ==============================================================================
# 7. PREDICCIONES Y MÉTRICAS INDEPENDIENTES DEL UMBRAL
# ==============================================================================

pred_valid <- predict(model_evaluation, matrix_valid)
pred_test  <- predict(model_evaluation, matrix_test)

if (
  any(!is.finite(pred_valid)) || any(pred_valid < 0 | pred_valid > 1) ||
  any(!is.finite(pred_test))  || any(pred_test < 0 | pred_test > 1)
) {
  stop("Se obtuvieron predicciones no válidas fuera del intervalo [0, 1].")
}

auc_valid <- as.numeric(
  pROC::auc(response = y_valid, predictor = pred_valid, quiet = TRUE)
)
auc_test <- as.numeric(
  pROC::auc(response = y_test, predictor = pred_test, quiet = TRUE)
)

pr_auc_valid <- PRROC::pr.curve(
  scores.class0 = pred_valid[y_valid == 1L],
  scores.class1 = pred_valid[y_valid == 0L]
)$auc.integral

pr_auc_test <- PRROC::pr.curve(
  scores.class0 = pred_test[y_test == 1L],
  scores.class1 = pred_test[y_test == 0L]
)$auc.integral

brier_valid <- brier_score(pred_valid, y_valid)
brier_test  <- brier_score(pred_test, y_test)

# ==============================================================================
# 8. SELECCIÓN DEL UMBRAL EN VALIDACIÓN
# ==============================================================================

threshold_results <- data.table::rbindlist(
  lapply(
    threshold_grid,
    function(threshold) {
      classification_metrics(pred_valid, y_valid, threshold)
    }
  )
)

best_threshold_row <- threshold_results[which.max(f1)]
optimal_threshold <- best_threshold_row$threshold

if (!is.finite(optimal_threshold)) {
  stop("No fue posible seleccionar un umbral válido mediante F1.")
}

data.table::fwrite(
  threshold_results,
  file.path(dir_models, "threshold_selection_validation.csv")
)

data.table::fwrite(
  data.frame(threshold = optimal_threshold),
  file.path(dir_models, "xgb_threshold.csv")
)

# ==============================================================================
# 9. MÉTRICAS DEPENDIENTES DEL UMBRAL
# ==============================================================================

metrics_valid_threshold <- classification_metrics(
  pred_valid,
  y_valid,
  optimal_threshold
)

metrics_test_threshold <- classification_metrics(
  pred_test,
  y_test,
  optimal_threshold
)

metrics_df <- data.table::rbindlist(list(
  data.frame(
    split = "validation",
    auc = auc_valid,
    pr_auc = pr_auc_valid,
    brier = brier_valid,
    metrics_valid_threshold
  ),
  data.frame(
    split = "test",
    auc = auc_test,
    pr_auc = pr_auc_test,
    brier = brier_test,
    metrics_test_threshold
  )
), fill = TRUE)

data.table::fwrite(
  metrics_df,
  file.path(dir_models, "metrics_valid_test.csv")
)

log_message(sprintf(
  "VALIDACIÓN | AUC=%.4f | PR-AUC=%.4f | Brier=%.4f | F1=%.4f",
  auc_valid,
  pr_auc_valid,
  brier_valid,
  metrics_valid_threshold$f1
))

log_message(sprintf(
  "PRUEBA | AUC=%.4f | PR-AUC=%.4f | Brier=%.4f | F1=%.4f | umbral=%.2f",
  auc_test,
  pr_auc_test,
  brier_test,
  metrics_test_threshold$f1,
  optimal_threshold
))

# ==============================================================================
# 10. IMPORTANCIA DE VARIABLES DEL MODELO DE EVALUACIÓN
# ==============================================================================

importance_evaluation <- xgboost::xgb.importance(
  feature_names = expected_features,
  model = model_evaluation
)

data.table::fwrite(
  importance_evaluation,
  file.path(dir_models, "feature_importance.csv")
)

# ==============================================================================
# 11. REAJUSTE DEL MODELO FINAL CON ENTRENAMIENTO + VALIDACIÓN
# ==============================================================================

# La prueba temporal se mantiene excluida del modelo final para conservar su
# independencia como evaluación externa del procedimiento de modelamiento.
train_full <- data.table::rbindlist(
  list(dtrain, dvalid),
  use.names = TRUE,
  fill = FALSE
)

X_train_full <- as.matrix(train_full[, ..expected_features])
y_train_full <- train_full$y

matrix_train_full <- xgboost::xgb.DMatrix(
  X_train_full,
  label = y_train_full,
  missing = NA
)

set.seed(seed_model)

model_final <- xgboost::xgb.train(
  params = xgb_params,
  data = matrix_train_full,
  nrounds = best_iteration,
  verbose = verbose_training
)

saveRDS(
  model_final,
  file.path(dir_models, "xgb_defrisk_model_final.rds")
)

data.table::fwrite(
  data.frame(threshold = optimal_threshold),
  file.path(dir_models, "xgb_threshold_final.csv")
)

importance_final <- xgboost::xgb.importance(
  feature_names = expected_features,
  model = model_final
)

data.table::fwrite(
  importance_final,
  file.path(dir_models, "feature_importance_final.csv")
)

# ==============================================================================
# 12. METADATOS PARA REPRODUCIBILIDAD E INFERENCIA
# ==============================================================================

data.table::fwrite(
  data.frame(
    position = seq_along(expected_features),
    feature = expected_features
  ),
  file.path(dir_models, "model_features.csv")
)

settings_table <- data.frame(
  parameter = c(
    "train_windows",
    "validation_window",
    "test_window",
    "best_iteration",
    "best_validation_auc",
    "optimal_threshold_f1",
    "seed",
    names(xgb_params)
  ),
  value = c(
    paste(vapply(win_train, paste, collapse = "-", FUN.VALUE = character(1)), collapse = "; "),
    paste(win_valid, collapse = "-"),
    paste(win_test, collapse = "-"),
    best_iteration,
    best_validation_auc,
    optimal_threshold,
    seed_model,
    unlist(xgb_params, use.names = FALSE)
  ),
  stringsAsFactors = FALSE
)

data.table::fwrite(
  settings_table,
  file.path(dir_models, "model_settings.csv")
)

log_message("Modelo XGBoost completado correctamente.")
log_message("Productos guardados en:", dir_models)
