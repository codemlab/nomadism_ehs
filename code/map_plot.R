## install.packages(c("elevatr", "ggspatial", "rnaturalearth", "rnaturalearthdata"))
## tidyterra is NOT required — hillshade is plotted via geom_raster on a plain data frame

library(elevatr)
library(terra)
library(sf)
library(ggplot2)
library(ggrepel)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# ---- Waypoints ----
waypoints <- data.frame(
  name = c("Tagendoust", "Amellagou", "Aït Hani", "Tinghir",
           "Jbel Saghro", "Boulmane\nDades", "Msemrir", "Agoudal"),
  lat  = c(32.336398, 32.011425, 31.778305, 31.512597,
           31.150331, 31.378035, 31.705714, 32.023067),
  lon  = c(-4.820442, -4.996739, -5.455544, -5.520666,
           -5.658711, -6.007226, -5.810226, -5.479298),
  season = c("summer", "summer", "village", "town",
             "winter", "village", "village", "summer")
)

# ---- Routes ----
# Spring: south to north through Todgha corridor
# Autumn: north to south via western route
spring_names <- c("Jbel Saghro", "Tinghir", "Aït Hani", "Amellagou", "Tagendoust")
autumn_names <- c("Tagendoust", "Agoudal", "Msemrir", "Boulmane\nDades", "Jbel Saghro")

make_route <- function(names, wp) {
  pts <- wp[match(names, wp$name), ]
  st_sf(geometry = st_sfc(st_linestring(as.matrix(pts[, c("lon","lat")])), crs = 4326))
}

spring_sf <- make_route(spring_names, waypoints)
autumn_sf <- make_route(autumn_names, waypoints)
pts_sf    <- st_as_sf(waypoints, coords = c("lon", "lat"), crs = 4326)

# ---- DEM and hillshade ----
loc_sf <- st_as_sf(
  data.frame(x = c(-6.5, -4.4), y = c(30.8, 32.6)),
  coords = c("x", "y"), crs = 4326
)
dem   <- get_elev_raster(locations = loc_sf, z = 9, clip = "bbox", neg_to_na = FALSE)
dem_t <- rast(dem)
hill  <- shade(
  terrain(dem_t, "slope",  unit = "radians"),
  terrain(dem_t, "aspect", unit = "radians"),
  angle = 40, direction = 315  # northwest illumination
)

# ---- Blend elevation colors + hillshade into a single color per cell ----
# Stack DEM and hillshade so rows align perfectly
combined    <- c(dem_t, hill)
combined_df <- as.data.frame(combined, xy = TRUE)
colnames(combined_df)[3:4] <- c("elev", "hill")

# Elevation color palette: sandy lowlands → ochre slopes → brown ridges → grey/white peaks
elev_pal  <- colorRampPalette(c("#e8d8b0", "#d4b87a", "#b09050",
                                 "#8a6840", "#b8b0a0", "#e0e0e0", "#ffffff"))
n_cols    <- 512
elev_rng  <- range(combined_df$elev, na.rm = TRUE)
na_cells  <- is.na(combined_df$elev) | is.na(combined_df$hill)

elev_idx  <- round((combined_df$elev - elev_rng[1]) / diff(elev_rng) * (n_cols - 1)) + 1
elev_idx  <- pmax(1L, pmin(n_cols, replace(elev_idx, is.na(elev_idx), 1L)))
pal_cols  <- elev_pal(n_cols)
elev_rgb  <- col2rgb(pal_cols[elev_idx])

# Hillshade factor — keep 35% ambient light so colour survives in shadows
hill_norm    <- (combined_df$hill - min(combined_df$hill, na.rm = TRUE)) /
                diff(range(combined_df$hill, na.rm = TRUE))
hill_norm    <- replace(hill_norm, is.na(hill_norm), 0)
shade_factor <- 0.35 + 0.65 * hill_norm

# Multiply RGB channels by shade, clamp to 0-255
r_b <- pmin(255L, round(elev_rgb[1L, ] * shade_factor))
g_b <- pmin(255L, round(elev_rgb[2L, ] * shade_factor))
b_b <- pmin(255L, round(elev_rgb[3L, ] * shade_factor))

combined_df$color <- rgb(r_b, g_b, b_b, maxColorValue = 255)
combined_df$color[na_cells] <- "#ffffff"

# ---- Main map ----
main_map <- ggplot() +
  # Blended terrain + hillshade — I() passes hex colors directly, no scale needed
  geom_raster(
    data = combined_df,
    aes(x = x, y = y, fill = I(color)),
    interpolate = TRUE, inherit.aes = FALSE
  ) +
  # Spring route (solid, warm)
  geom_sf(
    data = spring_sf, colour = "#E69F00", linewidth = 1.3,
    arrow = arrow(length = unit(0.35, "cm"), type = "closed", ends = "last")
  ) +
  # Autumn route (dashed, cool)
  geom_sf(
    data = autumn_sf, colour = "#56B4E9", linewidth = 1.3,
    linetype = "dashed",
    arrow = arrow(length = unit(0.35, "cm"), type = "closed", ends = "last")
  ) +
  # Waypoints
  geom_sf(
    data = pts_sf, size = 3, shape = 21,
    fill = "white", colour = "grey20", stroke = 0.8
  ) +
  # Labels
  geom_label_repel(
    data = waypoints,
    aes(x = lon, y = lat, label = name),
    size = 3.2, label.padding = unit(0.15, "lines"),
    box.padding = unit(0.5, "lines"),
    fill = alpha("white", 0.88), label.size = 0.2,
    segment.colour = "grey30", segment.size = 0.4,
    seed = 42
  ) +
  # Scale bar and north arrow
  annotation_scale(location = "br", width_hint = 0.22, style = "ticks",
                   text_cex = 0.75) +
  annotation_north_arrow(
    location = "tr", which_north = "true",
    style = north_arrow_minimal(text_size = 9),
    height = unit(1, "cm"), width = unit(0.7, "cm")
  ) +
  coord_sf(xlim = c(-6.35, -4.5), ylim = c(30.9, 32.5), expand = FALSE) +
  labs(
    x = NULL, y = NULL,
    caption = paste0(
      "——  Spring migration (south → north)",
      "          - - -  Autumn migration (north → south)"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid      = element_line(colour = alpha("white", 0.25), linewidth = 0.3),
    axis.text       = element_text(size = 8, colour = "grey40"),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.caption    = element_text(size = 8.5, hjust = 0.5, colour = "grey30",
                                   margin = margin(t = 6))
  )

# ---- Inset: Morocco + Western Sahara locator ----
mor_ws <- ne_countries(
  country = c("Morocco", "Western Sahara"),
  returnclass = "sf", scale = "medium"
)

inset <- ggplot() +
  geom_sf(data = mor_ws, fill = "grey80", colour = "white", linewidth = 0.3) +
  annotate("rect",
    xmin = -6.5, xmax = -4.4, ymin = 30.8, ymax = 32.6,
    fill = NA, colour = "#D55E00", linewidth = 0.9
  ) +
  coord_sf(xlim = c(-13.5, -1.0), ylim = c(20.5, 36.0), expand = FALSE) +
  theme_void() +
  theme(
    panel.border    = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
    plot.background = element_rect(fill = "white", colour = NA)
  )

# ---- Combine ----
main_map + inset_element(inset, left = 0.0, bottom = 0.0, right = 0.20, top = 0.22)

ggsave("manuscripts/nomadism_ehs/figures/migration_map.png",
       width = 8, height = 8, dpi = 300, bg = "white")
