# ============================================================
# POSTPROCESO Y ANÁLISIS INTEGRADO
# - Mapas comparativos por especie (2×2: binario/continuo, 2025/2040)
# - Mapas de riqueza (1×3: 2025 | 2040 | Δ)
# - Estadística zonal por ecorregión (resúmenes, boxplots, pruebas pareadas)
# - Riesgo de deforestación: mapa, proyección acumulada y áreas por ecorregión
# Autor: [Tu nombre]
# Fecha: [YYYY-MM-DD]
# Notas:
#   * Se mantienen rutas absolutas originales.
#   * No se remuestrean rásteres salvo donde es estrictamente necesario
#     (alinear 2040 a la rejilla 2025 o similar).
#   * Se eliminaron bloques duplicados y se homogeneizó la documentación.
# ============================================================


# ============================================================
# A) ESPECIE EJEMPLO — Plecturocebus caquetensis (Panel 2×2)
#    Arriba: mapas binarios (gris=no apto, negro=apto)
#    Abajo: mapas continuos (viridis)
#    Izq: 2025 | Der: 2040
#    Ecorregiones WWF (Colombia) en contorno
# ============================================================

# Paquetes
pkgs <- c("terra","sf","ggplot2","tidyterra","patchwork")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, dependencies = TRUE)
library(terra); library(sf); library(ggplot2); library(tidyterra); library(patchwork)

# Rutas
f_bin_2025  <- "C:/Users/henry.garcia/Desktop/PC_bin_2025.tif"
f_bin_2040  <- "C:/Users/henry.garcia/Desktop/PC_bin_2040.tif"
f_cont_2025 <- "C:/Users/henry.garcia/Desktop/PC_cont_2025.tif"
f_cont_2040 <- "C:/Users/henry.garcia/Desktop/PC_cont_2040.tif"
eco_shp     <- "C:/Users/henry.garcia/Desktop/Primates/ecoregions/wwf_terr_ecos_col.shp"
outfile     <- "C:/Users/henry.garcia/Desktop/PCaquetensis_2x2_eco_2025_2040.png"

for (p in c(f_bin_2025,f_bin_2040,f_cont_2025,f_cont_2040)) if (!file.exists(p)) stop("No existe: ", p)
if (!file.exists(eco_shp)) stop("No existe shapefile: ", eco_shp)

# Leer rásteres (SIN alinear ni remuestrear)
r_bin_2025  <- rast(f_bin_2025)
r_bin_2040  <- rast(f_bin_2040)
r_cont_2025 <- rast(f_cont_2025)
r_cont_2040 <- rast(f_cont_2040)

# Ecorregiones base
eco0 <- st_read(eco_shp, quiet = TRUE) |> st_make_valid()

# Adaptar ecorregiones al CRS/extent del ráster (sin modificar el ráster)
prep_eco_for <- function(r){
  eco <- suppressWarnings(st_transform(eco0, crs = st_crs(terra::crs(r, proj = TRUE))))
  sv  <- terra::vect(eco)
  svc <- tryCatch(terra::crop(sv, r), error = function(e) sv)
  st_as_sf(svc)
}

# Binarios a factor 0/1 con etiquetas
to_factor01 <- function(r){
  r01 <- terra::ifel(r > 0, 1, 0)  # >0 -> 1 (Apto); <=0 -> 0 (No apto); NA preserva
  rf  <- as.factor(r01)
  lv  <- data.frame(ID = c(0, 1), label = c("No apto", "Apto"))
  if (is.list(levels(rf))) levels(rf) <- list(lv) else levels(rf) <- lv
  rf
}
rb25f <- to_factor01(r_bin_2025)
rb40f <- to_factor01(r_bin_2040)

# Tema
theme_map <- theme_void(base_size = 11) +
  theme(plot.title = element_text(face="bold", margin=margin(b=3)),
        legend.title = element_text(),
        legend.position = "right")

# BINARIO 2025 (arriba-izq)
eco_b25 <- prep_eco_for(rb25f)
p_bin_2025 <- ggplot() +
  geom_spatraster(data = rb25f, na.rm = TRUE) +
  scale_fill_manual(
    values = c("No apto" = "grey85", "Apto" = "black"),
    name   = "Binario",
    breaks = c("No apto","Apto"),
    limits = c("No apto","Apto"),
    drop   = FALSE
  ) +
  geom_sf(data = eco_b25, fill = NA, color = "black", linewidth = 0.25) +
  coord_sf(expand = FALSE, datum = NA) +
  labs(title = "Binario 2025") + theme_map

# BINARIO 2040 (arriba-der)
eco_b40 <- prep_eco_for(rb40f)
p_bin_2040 <- ggplot() +
  geom_spatraster(data = rb40f, na.rm = TRUE) +
  scale_fill_manual(
    values = c("No apto" = "grey85", "Apto" = "black"),
    name   = "Binario",
    breaks = c("No apto","Apto"),
    limits = c("No apto","Apto"),
    drop   = FALSE
  ) +
  geom_sf(data = eco_b40, fill = NA, color = "black", linewidth = 0.25) +
  coord_sf(expand = FALSE, datum = NA) +
  labs(title = "Binario 2040 (SSP585)") + theme_map

# CONTINUO 2025 (abajo-izq)
eco_c25 <- prep_eco_for(r_cont_2025)
p_cont_2025 <- ggplot() +
  geom_spatraster(data = r_cont_2025, na.rm = TRUE) +
  scale_fill_viridis_c(name = "Idoneidad", na.value = NA) +
  geom_sf(data = eco_c25, fill = NA, color = "black", linewidth = 0.25) +
  coord_sf(expand = FALSE, datum = NA) +
  labs(title = "Continuo 2025") + theme_map

# CONTINUO 2040 (abajo-der)
eco_c40 <- prep_eco_for(r_cont_2040)
p_cont_2040 <- ggplot() +
  geom_spatraster(data = r_cont_2040, na.rm = TRUE) +
  scale_fill_viridis_c(name = "Idoneidad", na.value = NA) +
  geom_sf(data = eco_c40, fill = NA, color = "black", linewidth = 0.25) +
  coord_sf(expand = FALSE, datum = NA) +
  labs(title = "Continuo 2040 (SSP585)") + theme_map

# Panel 2×2 y guardado
panel <- (p_bin_2025 + p_bin_2040) /
  (p_cont_2025 + p_cont_2040) +
  plot_annotation(
    title = "Plecturocebus caquetensis — 2025 vs 2040",
    subtitle = "Arriba: binarios; Abajo: continuos — contorno: ecorregiones WWF (Colombia)",
    theme = theme(plot.title = element_text(face="bold", size=14),
                  plot.subtitle = element_text(size=11))
  )
ggsave(outfile, panel, width = 12, height = 10, dpi = 300)
message("✅ Figura guardada en: ", outfile)



# ============================================================
# B) RIQUEZA DE PRIMATES — Panel 1×3 (2025 | 2040 | Δ=2040−2025)
#    - Alinea 2040 a la rejilla de 2025 (baseline)
#    - Ecorregiones en contorno
# ============================================================

pkgs <- c("terra","sf","ggplot2","tidyterra","patchwork")
new  <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(new)) install.packages(new, dependencies = TRUE)
library(terra); library(sf); library(ggplot2); library(tidyterra); library(patchwork)

# Rutas
f_rich_2025 <- "C:/Users/henry.garcia/Desktop/Primates/models_maxent_ALL_ecoregM/0.richness_maps/riqueza_2025_col.tif"
f_rich_2040 <- "C:/Users/henry.garcia/Desktop/Primates/models_maxent_ALL_ecoregM/0.richness_maps/riqueza_2040_col.tif"
f_rich_delta<- "C:/Users/henry.garcia/Desktop/Primates/models_maxent_ALL_ecoregM/0.richness_maps/riqueza_dellta_col.tif"  # si existe, se usa
eco_shp     <- "C:/Users/henry.garcia/Desktop/Primates/ecoregions/wwf_terr_ecos_col.shp"

for (p in c(f_rich_2025,f_rich_2040)) if (!file.exists(p)) stop("No existe: ", p)
if (!file.exists(eco_shp)) stop("No existe shapefile: ", eco_shp)

# Leer
r25    <- rast(f_rich_2025)
r40raw <- rast(f_rich_2040)

# Alinear 2040 a 2025
r40 <- if (!compareGeom(r25, r40raw, stopOnError = FALSE)) {
  project(r40raw, r25, method = "bilinear")
} else r40raw

# Delta (usar archivo si existe; si no, calcular)
if (file.exists(f_rich_delta)) {
  rD_raw <- rast(f_rich_delta)
  rD <- if (!compareGeom(r25, rD_raw, stopOnError = FALSE)) {
    project(rD_raw, r25, method = "bilinear")
  } else rD_raw
} else {
  rD <- r40 - r25
}

# Ecorregiones a cada ráster
eco0 <- st_read(eco_shp, quiet = TRUE) |> st_make_valid()
prep_eco_for <- function(r){
  rcrs <- crs(r, proj = TRUE)
  eco  <- if (nzchar(rcrs)) suppressWarnings(st_transform(eco0, st_crs(rcrs))) else eco0
  sv   <- vect(eco); svc <- tryCatch(crop(sv, r), error = function(e) sv)
  st_as_sf(svc)
}

# Escalas
vmax <- max(as.numeric(global(r25, "max", na.rm = TRUE)),
            as.numeric(global(r40, "max", na.rm = TRUE)), na.rm = TRUE)
dmin <- as.numeric(global(rD, "min", na.rm = TRUE))
dmax <- as.numeric(global(rD, "max", na.rm = TRUE))

theme_map <- theme_void(base_size = 11) +
  theme(plot.title = element_text(face="bold", margin=margin(b=3)),
        legend.title = element_text(),
        legend.position = "right")
coord_map <- coord_sf(expand = FALSE, datum = NA)

# Mapas y panel
eco25 <- prep_eco_for(r25)
p25 <- ggplot() +
  geom_spatraster(data = r25, na.rm = TRUE) +
  scale_fill_viridis_c(name = "Riqueza", limits = c(0, vmax), oob = scales::squish) +
  geom_sf(data = eco25, fill = NA, color = "black", linewidth = 0.25) +
  coord_map + labs(title = "Riqueza 2025") + theme_map

eco40 <- prep_eco_for(r40)
p40 <- ggplot() +
  geom_spatraster(data = r40, na.rm = TRUE) +
  scale_fill_viridis_c(name = "Riqueza", limits = c(0, vmax), oob = scales::squish) +
  geom_sf(data = eco40, fill = NA, color = "black", linewidth = 0.25) +
  coord_map + labs(title = "Riqueza 2040 (SSP585)") + theme_map

ecoD <- prep_eco_for(rD)
pD <- ggplot() +
  geom_spatraster(data = rD, na.rm = TRUE) +
  scale_fill_viridis_c(name = "Δ riqueza (2040 − 2025)",
                       limits = c(dmin, dmax), oob = scales::squish) +
  geom_sf(data = ecoD, fill = NA, color = "black", linewidth = 0.25) +
  coord_map + labs(title = "Δ Riqueza (2040 − 2025)") + theme_map

panel <- p25 + p40 + pD +
  plot_annotation(
    title    = "Riqueza de primates — 2025 | 2040 | Δ 2040−2025 (viridis)",
    subtitle = "Ecorregiones WWF (Colombia) en contorno",
    theme = theme(plot.title = element_text(face="bold", size=14),
                  plot.subtitle = element_text(size = 11))
  )

outfile <- "C:/Users/henry.garcia/Desktop/richness_1x3_viridis_2025_2040_delta.png"
ggsave(outfile, panel, width = 16, height = 6, dpi = 300)
message("✅ Figura guardada en: ", outfile)



# ============================================================
# C) RIQUEZA PROMEDIO POR ECORREGIÓN Y BOXPLOTS 2025 vs 2040
#    - Estadística zonal (media) + n celdas válidas
#    - Boxplots por periodo y (opcional) Δ por ecorregión
# ============================================================

pkgs <- c("terra","sf","dplyr","ggplot2","tidyr","readr")
new  <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)
library(terra); library(sf); library(dplyr); library(ggplot2); library(tidyr); library(readr)

# Rutas
f_rich_2025 <- "C:/Users/henry.garcia/Desktop/Primates/models_maxent_ALL_ecoregM/0.richness_maps/riqueza_2025_col.tif"
f_rich_2040 <- "C:/Users/henry.garcia/Desktop/Primates/models_maxent_ALL_ecoregM/0.richness_maps/riqueza_2040_col.tif"
eco_shp     <- "C:/Users/henry.garcia/Desktop/Primates/ecoregions/wwf_terr_ecos_col.shp"
stopifnot(file.exists(f_rich_2025), file.exists(f_rich_2040), file.exists(eco_shp))

# Leer y alinear 2040 a 2025
r25    <- rast(f_rich_2025)
r40raw <- rast(f_rich_2040)
r40    <- if (!compareGeom(r25, r40raw, stopOnError = FALSE)) project(r40raw, r25, method = "bilinear") else r40raw

# Ecorregiones al CRS de r25 + ID y rasterización
eco_sf <- st_read(eco_shp, quiet = TRUE) |> st_make_valid() |> st_transform(st_crs(crs(r25, proj = TRUE)))
stopifnot("ECO_NAME" %in% names(eco_sf))
eco_sf <- eco_sf |> mutate(ECO_ID = as.integer(factor(ECO_NAME)))
eco_v  <- vect(eco_sf)
eco_r  <- rasterize(eco_v, r25, field = "ECO_ID", touches = TRUE)

# Zonal (medias y nº celdas válidas)
z_mean_25 <- zonal(r25, eco_r, fun = "mean", na.rm = TRUE) |> as.data.frame()
z_mean_40 <- zonal(r40, eco_r, fun = "mean", na.rm = TRUE) |> as.data.frame()
z_n_25    <- zonal(!is.na(r25), eco_r, fun = "sum",  na.rm = TRUE) |> as.data.frame()
z_n_40    <- zonal(!is.na(r40), eco_r, fun = "sum",  na.rm = TRUE) |> as.data.frame()

names(z_mean_25) <- c("ECO_ID","mean_2025")
names(z_mean_40) <- c("ECO_ID","mean_2040")
names(z_n_25)    <- c("ECO_ID","n25")
names(z_n_40)    <- c("ECO_ID","n40")

eco_names <- eco_sf |> st_drop_geometry() |> select(ECO_ID, ECO_NAME) |> distinct()

df <- eco_names %>%
  left_join(z_mean_25, by = "ECO_ID") %>%
  left_join(z_mean_40, by = "ECO_ID") %>%
  left_join(z_n_25,    by = "ECO_ID") %>%
  left_join(z_n_40,    by = "ECO_ID") %>%
  filter((n25 > 0) | (n40 > 0)) %>%
  mutate(delta = mean_2040 - mean_2025)

# Boxplot comparando promedios por ecorregión
df_long <- df |>
  pivot_longer(cols = c(mean_2025, mean_2040),
               names_to = "periodo", values_to = "riqueza_media") |>
  mutate(periodo = recode(periodo, mean_2025 = "2025", mean_2040 = "2040"))

p_box <- ggplot(df_long, aes(x = periodo, y = riqueza_media, fill = periodo)) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.15, alpha = 0.65, size = 1.6, aes(color = periodo)) +
  scale_fill_viridis_d(end = 0.9, guide = "none") +
  scale_color_viridis_d(end = 0.9, guide = "none") +
  labs(
    title = "Riqueza promedio por ecorregión (ECO_NAME)",
    subtitle = sprintf("n ecorregiones: %d | Promedio global: 2025=%.2f, 2040=%.2f",
                       dplyr::n_distinct(df$ECO_ID),
                       mean(df$mean_2025, na.rm = TRUE),
                       mean(df$mean_2040, na.rm = TRUE)),
    x = NULL, y = "Riqueza media (promedio por ecorregión)"
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

out_png <- "C:/Users/henry.garcia/Desktop/box_riqueza_promedio_ecoregiones_2025_2040.png"
out_csv <- "C:/Users/henry.garcia/Desktop/riqueza_promedio_por_ecoregion_2025_2040.csv"
ggsave(out_png, p_box, width = 8.8, height = 6.0, dpi = 300)
write_csv(df |> arrange(ECO_NAME), out_csv)
message("✅ Boxplot: ", out_png)
message("✅ Tabla (por ecorregión): ", out_csv)

# Boxplot del Δ por ecorregión (opcional)
p_delta <- ggplot(df, aes(y = delta, x = "Δ 2040−2025", fill = "Δ 2040−2025")) +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width = 0.12, alpha = 0.65, size = 1.6, color = "grey20") +
  scale_fill_viridis_d(end = 0.9, guide = "none") +
  labs(title = "Cambio en riqueza promedio por ecorregión (2040−2025)",
       x = NULL, y = "Δ riqueza (promedio por ecorregión)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))
out_png_delta <- "C:/Users/henry.garcia/Desktop/box_delta_riqueza_promedio_ecoregiones.png"
ggsave(out_png_delta, p_delta, width = 6.2, height = 6.0, dpi = 300)
message("✅ Box Δ: ", out_png_delta)



# ============================================================
# D) PRUEBAS PAREADAS POR ECORREGIÓN: 2025 vs 2040 (RIQUEZA)
#    - t pareada sobre Δ si normalidad no se rechaza;
#      si no, Wilcoxon; si muchos empates, Sign test
#    - Exporta resumen global y por ecorregión
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(data.table)
})

# Verificar insumo 'df' (proveniente del bloque C)
if (!exists("df")) stop("No existe el objeto 'df'. Ejecuta primero el bloque C (zonal y tabla 'df').")

dd <- df %>%
  dplyr::filter(!is.na(mean_2025), !is.na(mean_2040)) %>%
  dplyr::mutate(delta = mean_2040 - mean_2025)
if (nrow(dd) < 3) stop("Muy pocas ecorregiones con datos en ambos años para prueba pareada.")

# Normalidad orientativa de Δ (global, solo reportada)
sw <- tryCatch(shapiro.test(dd$delta), error = function(e) NULL)
tt <- t.test(dd$mean_2040, dd$mean_2025, paired = TRUE, alternative = "two.sided")
dz <- mean(dd$delta, na.rm = TRUE) / stats::sd(dd$delta, na.rm = TRUE)
wx <- suppressWarnings(
  wilcox.test(dd$mean_2040, dd$mean_2025,
              paired = TRUE, exact = FALSE,
              conf.int = TRUE, conf.level = 0.95)
)

summary_df <- data.frame(
  n_ecoregiones = nrow(dd),
  media_2025 = mean(dd$mean_2025, na.rm = TRUE),
  media_2040 = mean(dd$mean_2040, na.rm = TRUE),
  delta_media = mean(dd$delta, na.rm = TRUE),
  delta_sd = sd(dd$delta, na.rm = TRUE),
  t_stat = unname(tt$statistic),
  gl = unname(tt$parameter),
  p_t = tt$p.value,
  ci_delta_low = tt$conf.int[1],
  ci_delta_high = tt$conf.int[2],
  cohen_dz = dz,
  W_wilcoxon = unname(wx$statistic),
  p_wilcoxon = wx$p.value,
  wilcox_ci_low = if (!is.null(wx$conf.int)) wx$conf.int[1] else NA_real_,
  wilcox_ci_high = if (!is.null(wx$conf.int)) wx$conf.int[2] else NA_real_
)

out_sum <- "C:/Users/henry.garcia/Desktop/paired_test_riqueza_summary.csv"
out_per <- "C:/Users/henry.garcia/Desktop/paired_test_riqueza_por_ecoregion.csv"
utils::write.csv(summary_df, out_sum, row.names = FALSE, fileEncoding = "UTF-8")
utils::write.csv(dd[, c("ECO_NAME","mean_2025","mean_2040","delta")],
                 out_per, row.names = FALSE, fileEncoding = "UTF-8")

cat("============================================================\n")
cat(sprintf("Ecorregiones incluidas: %d\n", nrow(dd)))
if (!is.null(sw)) {
  cat(sprintf("Shapiro–Wilk (Δ): W=%.3f, p=%.4g -> %s\n",
              sw$statistic, sw$p.value,
              if (sw$p.value >= 0.05) "normalidad no rechazada (t pareada OK)"
              else "normalidad rechazada (reportar también Wilcoxon)"))
}
cat(sprintf("t pareada: t=%.3f (gl=%d), p=%.5g\n", tt$statistic, as.integer(tt$parameter), tt$p.value))
cat(sprintf("Δ media = %.3f (IC95%%: %.3f a %.3f), Cohen's dz = %.3f\n",
            mean(dd$delta), tt$conf.int[1], tt$conf.int[2], dz))
cat(sprintf("Wilcoxon pareada: W=%g, p=%.5g", wx$statistic, wx$p.value))
if (!is.null(wx$conf.int)) cat(sprintf(", IC mediana Δ [%.3f, %.3f]\n", wx$conf.int[1], wx$conf.int[2])) else cat("\n")
cat("Archivos:\n - Resumen: ", out_sum, "\n - Por ecorregión: ", out_per, "\n")
cat("============================================================\n")


# ============================================================
# E) POST HOC POR ECORREGIÓN (selección automática de prueba)
#    - t(Δ) si normalidad OK; si no, Wilcoxon; con ≥30% empates → Sign test
#    - Exporta CSV y gráfico de Δ coloreado por significancia (FDR BH)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# Reconstruir insumos mínimos si faltan en memoria
if (!exists("r25") || !exists("r40")) {
  f_rich_2025 <- "C:/Users/henry.garcia/Desktop/Primates/models_maxent_ALL_ecoregM/0.richness_maps/riqueza_2025_col.tif"
  f_rich_2040 <- "C:/Users/henry.garcia/Desktop/Primates/models_maxent_ALL_ecoregM/0.richness_maps/riqueza_2040_col.tif"
  stopifnot(file.exists(f_rich_2025), file.exists(f_rich_2040))
  r25    <- rast(f_rich_2025)
  r40raw <- rast(f_rich_2040)
  r40    <- if (!compareGeom(r25, r40raw, stopOnError = FALSE)) project(r40raw, r25, method = "bilinear") else r40raw
}
if (!exists("eco_r") || !exists("eco_sf")) {
  eco_shp <- "C:/Users/henry.garcia/Desktop/Primates/ecoregions/wwf_terr_ecos_col.shp"
  stopifnot(file.exists(eco_shp))
  eco_sf <- st_read(eco_shp, quiet = TRUE) |> st_make_valid() |>
    st_transform(st_crs(crs(r25, proj = TRUE)))
  stopifnot("ECO_NAME" %in% names(eco_sf))
  eco_sf <- eco_sf |> mutate(ECO_ID = as.integer(factor(ECO_NAME)))
  eco_r  <- rasterize(vect(eco_sf), r25, field = "ECO_ID", touches = TRUE)
}

# Extraer pares por píxel y ecorregión
v25 <- values(r25, mat = FALSE)
v40 <- values(r40, mat = FALSE)
zid <- values(eco_r, mat = FALSE)
ok  <- !is.na(zid) & !is.na(v25) & !is.na(v40)
DT  <- data.table(ECO_ID = zid[ok], v25 = v25[ok], v40 = v40[ok], delta = v40[ok] - v25[ok])

# Utilidades de prueba
shapiro_smart <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3) return(list(p.value = NA_real_, W = NA_real_))
  if (n > 5000) x <- sample(x, 5000)
  tryCatch(shapiro.test(x), error = function(e) list(p.value = NA_real_, W = NA_real_))
}
sign_test <- function(d) {
  d <- d[is.finite(d)]
  nz <- d[d != 0]
  n  <- length(nz)
  if (n == 0) return(list(p.value = 1, stat = 0, n = 0))
  kpos <- sum(nz > 0)
  bt   <- binom.test(kpos, n, p = 0.5, alternative = "two.sided")
  list(p.value = bt$p.value, stat = kpos, n = n)
}

# Analítica por ecorregión (selección de prueba)
analyze_zone <- function(df_zone) {
  n <- nrow(df_zone)
  mean25 <- as.numeric(mean(df_zone$v25, na.rm = TRUE))
  mean40 <- as.numeric(mean(df_zone$v40, na.rm = TRUE))
  d <- df_zone$delta
  delta_mean <- as.numeric(mean(d, na.rm = TRUE))
  if (n < 10) {
    return(data.frame(test = "NA", n = as.integer(n), mean25 = mean25, mean40 = mean40,
                      delta_mean = delta_mean, stat = NA_real_, p = NA_real_,
                      method_note = "n<10", stringsAsFactors = FALSE))
  }
  prop_ties <- mean(d == 0, na.rm = TRUE)
  if (prop_ties >= 0.30) {
    nz <- d[d != 0]
    if (length(nz) == 0) {
      return(data.frame(test = "Sign (binomial)", n = as.integer(n), mean25 = mean25, mean40 = mean40,
                        delta_mean = delta_mean, stat = 0, p = 1,
                        method_note = sprintf("ties=%.1f%%; sin datos no nulos", 100*prop_ties),
                        stringsAsFactors = FALSE))
    }
    kpos <- sum(nz > 0)
    bt <- binom.test(kpos, length(nz), p = 0.5, alternative = "two.sided")
    return(data.frame(test = "Sign (binomial)", n = as.integer(n), mean25 = mean25, mean40 = mean40,
                      delta_mean = delta_mean, stat = as.numeric(kpos), p = as.numeric(bt$p.value),
                      method_note = sprintf("ties=%.1f%%", 100*prop_ties), stringsAsFactors = FALSE))
  }
  sw <- shapiro_smart(d)
  if (!is.na(sw$p.value) && sw$p.value >= 0.05) {
    tt <- t.test(d, mu = 0)
    sd_d <- sd(d, na.rm = TRUE)
    dz <- if (is.finite(sd_d) && sd_d > 0) mean(d, na.rm = TRUE)/sd_d else NA_real_
    data.frame(test = "t pareada (Δ)", n = as.integer(n), mean25 = mean25, mean40 = mean40,
               delta_mean = delta_mean, stat = as.numeric(tt$statistic), p = as.numeric(tt$p.value),
               method_note = sprintf("Shapiro p=%.3g; dz=%.3f", sw$p.value, dz), stringsAsFactors = FALSE)
  } else {
    wx <- suppressWarnings(wilcox.test(df_zone$v40, df_zone$v25, paired = TRUE, exact = FALSE))
    data.frame(test = "Wilcoxon pareada", n = as.integer(n), mean25 = mean25, mean40 = mean40,
               delta_mean = delta_mean, stat = as.numeric(wx$statistic), p = as.numeric(wx$p.value),
               method_note = if (is.na(sw$p.value)) "Shapiro NA" else sprintf("Shapiro p=%.3g", sw$p.value),
               stringsAsFactors = FALSE)
  }
}

# Ejecutar por ecorregión
res_list <- lapply(split(DT, DT$ECO_ID), analyze_zone)
RES <- data.table::rbindlist(res_list, idcol = "ECO_ID", use.names = TRUE, fill = TRUE)
RES[, ECO_ID := as.integer(ECO_ID)]

# Añadir nombres de ecorregión
if (exists("eco_sf")) {
  name_map <- sf::st_drop_geometry(eco_sf)[, c("ECO_ID","ECO_NAME")] |> dplyr::distinct()
  RES <- merge(RES, data.table::as.data.table(name_map), by = "ECO_ID", all.x = TRUE)
} else {
  RES[, ECO_NAME := paste0("ECO_", ECO_ID)]
}

# Ajuste FDR (BH), ordenar y exportar
RES <- RES[order(p)]
RES[, p_adj_BH := p.adjust(p, method = "BH")]
RES[, decision := ifelse(!is.na(p_adj_BH) & p_adj_BH < 0.05, "Significativo (FDR<0.05)", "ns")]

out_csv <- "C:/Users/henry.garcia/Desktop/posthoc_pareado_por_ecoregion.csv"
data.table::fwrite(RES[, .(ECO_ID, ECO_NAME, n, mean25, mean40, delta_mean,
                           test, stat, p, p_adj_BH, decision, method_note)],
                   out_csv)
message("✅ Resultados por ecorregión guardados en: ", out_csv)

# Gráfico Δ por ecorregión coloreado por significancia
p_sig <- ggplot(RES, aes(x = reorder(ECO_NAME, delta_mean), y = delta_mean, fill = decision)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Significativo (FDR<0.05)" = "#21918c", "ns" = "grey80")) +
  labs(title = "Δ riqueza (2040−2025) por ecorregión",
       x = NULL, y = "Δ riqueza promedio", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave("C:/Users/henry.garcia/Desktop/posthoc_delta_bar_sig.png", p_sig, width = 9, height = 10, dpi = 300)
message("✅ Gráfico Δ por ecorregión: C:/Users/henry.garcia/Desktop/posthoc_delta_bar_sig.png")



# ============================================================
# F) RIESGO DE DEFORESTACIÓN (90 m)
#    - Mapa de probabilidad (colores V-A-R) en bosque 2023
#    - Probabilidad acumulada P = 1 - (1 - p)^n y forest01_2040 (umbral 0.45)
#    - Áreas de bosque (ha) por ecorregión: 2023 vs 2040 + barras
# ============================================================

pkgs <- c("terra","sf","ggplot2","tidyterra","data.table","dplyr","tidyr")
new  <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)
library(terra); library(sf); library(ggplot2); library(tidyterra)
library(data.table); library(dplyr); library(tidyr)

# Rutas
risk_path  <- "C:/Users/henry.garcia/Desktop/Primates/outputs/defrisk_prob_2025_2030_forestOnly_90m.tif"
f2023_path <- "C:/Users/henry.garcia/Desktop/Primates/predictors/DefRisk/rasters/forest01/forest01_2023_90m.tif"
eco_shp    <- "C:/Users/henry.garcia/Desktop/Primates/ecoregions/wwf_terr_ecos_col.shp"
stopifnot(file.exists(risk_path), file.exists(f2023_path), file.exists(eco_shp))

# Leer y alinear riesgo a malla 2023 si hace falta
risk  <- rast(risk_path)   # p en [0,1] sobre bosque; fuera = NA
f2023 <- rast(f2023_path)  # 1=bosque, 0=no bosque
if (!compareGeom(f2023, risk, stopOnError = FALSE)) {
  risk <- project(risk, f2023, method = "bilinear")
}

# Ecorregiones (contorno y raster para zonal)
eco_sf0 <- st_read(eco_shp, quiet = TRUE) |> st_make_valid()
eco_sf  <- suppressWarnings(st_transform(eco_sf0, st_crs(crs(f2023, proj = TRUE))))
stopifnot("ECO_NAME" %in% names(eco_sf))
eco_sf <- eco_sf |> mutate(ECO_ID = as.integer(factor(ECO_NAME)))
eco_r  <- rasterize(vect(eco_sf), f2023, field = "ECO_ID", touches = TRUE)

# (1) Mapa de riesgo (verde→amarillo→rojo) solo en bosque
risk_plot <- terra::mask(risk, f2023, maskvalue = 0)
p_risk <- ggplot() +
  geom_spatraster(data = risk_plot, na.rm = TRUE) +
  scale_fill_gradientn(
    name   = "Prob. deforestación\n(2025–2030)",
    colours = c("#2CA02C", "#FFFF00", "#D62728"),
    values  = scales::rescale(c(0, 0.5, 1)),
    limits  = c(0, 1),
    oob     = scales::squish,
    na.value = NA
  ) +
  geom_sf(data = eco_sf, fill = NA, color = "black", linewidth = 0.25) +
  coord_sf(expand = FALSE, datum = NA) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA)) +
  labs(title = "Riesgo de deforestación (90 m) — V-A-R")
ggsave("C:/Users/henry.garcia/Desktop/map_defrisk_VYR.png", p_risk, width = 9, height = 7, dpi = 300)

# (2) Probabilidad acumulada y forest01_2040 con umbral 0.45
n_steps <- 3
p <- risk
p <- terra::ifel(p < 0, 0, p)
p <- terra::ifel(p > 1, 1, p)
p_cum <- 1 - (1 - p)^n_steps
thr   <- 0.45
lost        <- p_cum >= thr
forest2040  <- ifel(f2023 == 1 & lost, 0, f2023)  # 1=bosque, 0=no bosque
if (all(is.na(values(forest2040)))) warning("forest2040 resultó todo NA; revisa insumos/umbral.")
writeRaster(forest2040, "C:/Users/henry.garcia/Desktop/forest01_2040_90m.tif", overwrite = TRUE)

# (3) Máscaras de bosque (0 -> NA) para cálculo de área
f2023_mask <- ifel(f2023 == 1, 1, NA)
f2040_mask <- ifel(forest2040 == 1, 1, NA)

# (4) Área de BOSQUE por ecorregión (ha) 2023 vs 2040
px_ha <- (90 * 90) / 10000  # 0.81 ha por píxel
z_2023 <- zonal(f2023_mask, eco_r, fun = "sum", na.rm = TRUE) |> as.data.frame()
z_2040 <- zonal(f2040_mask, eco_r, fun = "sum", na.rm = TRUE) |> as.data.frame()
names(z_2023) <- c("ECO_ID","n_pix_bosque_2023")
names(z_2040) <- c("ECO_ID","n_pix_bosque_2040")

areas <- dplyr::full_join(z_2023, z_2040, by = "ECO_ID") |>
  dplyr::mutate(across(starts_with("n_pix_bosque_"), ~ tidyr::replace_na(., 0))) |>
  dplyr::mutate(area_ha_bosque_2023 = n_pix_bosque_2023 * px_ha,
                area_ha_bosque_2040 = n_pix_bosque_2040 * px_ha,
                delta_ha = area_ha_bosque_2040 - area_ha_bosque_2023) |>
  dplyr::left_join(st_drop_geometry(eco_sf)[, c("ECO_ID","ECO_NAME")] |> dplyr::distinct(), by = "ECO_ID") |>
  dplyr::relocate(ECO_ID, ECO_NAME)

# Guardar tabla y barras
out_csv <- "C:/Users/henry.garcia/Desktop/area_bosque_por_ecoregion_2023_2040.csv"
data.table::fwrite(areas, out_csv)

areas_long <- areas |>
  dplyr::select(ECO_NAME, area_ha_bosque_2023, area_ha_bosque_2040) |>
  tidyr::pivot_longer(cols = starts_with("area_ha_bosque_"),
                      names_to = "periodo", values_to = "area_ha") |>
  dplyr::mutate(periodo = dplyr::recode(periodo,
                                        area_ha_bosque_2023 = "2023",
                                        area_ha_bosque_2040 = "2040"))
ord <- areas |> dplyr::arrange(dplyr::desc(area_ha_bosque_2040)) |> dplyr::pull(ECO_NAME)
areas_long$ECO_NAME <- factor(areas_long$ECO_NAME, levels = ord)

p_bosque <- ggplot(areas_long, aes(x = ECO_NAME, y = area_ha, fill = periodo)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  scale_fill_viridis_d(end = 0.9, name = NULL) +
  labs(title = "Área de BOSQUE por ecorregión (90 m)",
       subtitle = sprintf("Prob. acumulada con n=%d; umbral P≥%.2f para pérdida", n_steps, thr),
       x = NULL, y = "Área de bosque (ha)") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave("C:/Users/henry.garcia/Desktop/barras_area_bosque_ecoregion_2023_2040.png",
       p_bosque, width = 12, height = max(6, 0.35*length(ord)+2), dpi = 300)

message("✅ Mapa riesgo: C:/Users/henry.garcia/Desktop/map_defrisk_VYR.png")
message("✅ forest01_2040 (1=bosque): C:/Users/henry.garcia/Desktop/forest01_2040_90m.tif")
message("✅ Tabla áreas (ha): ", out_csv)
message("✅ Barras área bosque: C:/Users/henry.garcia/Desktop/barras_area_bosque_ecoregion_2023_2040.png")
