# 1. Install and load packages -------------------------------------------------
packages <- c(
  "sf", "terra", "exactextractr", "viridis", "maptiles",
  "tidyterra", "ggspatial", "readr", "ggplot2", "here")

# install any missing packages
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install) > 0) {
  install.packages(to_install)
}

invisible(lapply(packages, library, character.only = TRUE))

# 2. Load and prepare data -----------------------------------------------------
# use relative paths for better portability
# (edit these paths if your directory structure is different)

sinap <- st_read(here::here("data", "RUNAP_clip.shp"))
sinap$zoneID <- seq_len(nrow(sinap))

landuse <- st_read(here::here("data", "Coberturas_clip_reclass.shp"))

# create binary field: 1 for forest (N) and 0 otherwise
landuse$bin <- ifelse(landuse$area_type == "N", 1, 0)
landuse_vect <- vect(landuse)

# create base raster at 1 km resolution
res_km <- 1000
raster_template <- rast(ext(landuse_vect), resolution = res_km, crs = crs(landuse_vect))

# rasterize using binary field
landuse_raster <- rasterize(landuse_vect, raster_template, field = "bin", touches = TRUE)

# 3. Compute mean forest area density (FAD) ------------------------------------
# a 7 x 7 window is used but you can modify `win_size` if needed
win_size <- 7
w <- matrix(1, nrow = win_size, ncol = win_size)

forest_sum <- focal(landuse_raster, w = w, fun = sum, na.rm = TRUE, fillvalue = 0)
forest_density <- forest_sum / (win_size * win_size)

# 4. Compute mean FAD per protected area using only forest pixels -------------

sinap <- st_transform(sinap, crs = crs(landuse_raster))
sinap_vect <- vect(sinap)

# mask density raster so only forest pixels contribute
forest_density_only <- mask(forest_density, landuse_raster, maskvalues = 0, updatevalue = NA)

ext_df <- extract(forest_density_only, sinap_vect, fun = mean, na.rm = TRUE, exact = TRUE,
                  touches = TRUE, ID = TRUE)

names(ext_df)[names(ext_df) == "focal"] <- "mean_FAD"

sinap$mean_FAD <- ext_df$mean_FAD[match(sinap$zoneID, ext_df$ID)]

# 5. Export results ------------------------------------------------------------
output_df <- sinap[, c("zoneID", "nombre", "mean_FAD")]
write_csv(output_df, here::here("results", "FAD_results.csv"))

# 6. Visualise results ---------------------------------------------------------

basemap <- get_tiles(sinap, provider = "Esri.WorldImagery", crop = TRUE)

p <- ggplot() +
  geom_spatraster_rgb(data = basemap) +
  geom_sf(data = sinap, aes(fill = mean_FAD), color = NA) +
  scale_fill_gradientn(
    colours = c("red", "#FFD301", "#006400"),
    values = scales::rescale(c(0, 0.5, 1)),
    limits = c(0, 1),
    na.value = "grey90",
    name = "Mean FAD"
  ) +
  guides(fill = guide_colorbar()) +
  coord_sf(datum = NA) +
  theme_void() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
  labs(title = "Densidad media de bosque (FAD) por área protegida")

# 7. Save plot ----------------------------------------------------------------

ggsave(here::here("results", "FAD_map.png"), plot = p, width = 10, height = 8, dpi = 300)

