PATH = "C:/Users/utente/Nextcloud/01_work/01_research/01_projects/02_SPATIAL_COX/package/"
setwd(PATH)
{
  # Package names
  packages <- c(
    "ggplot2", "survival", "fdaPDE", "tidytable", "ggspatial", "INLA",
    "Matrix","fields","dplyr","sf","scales","rmapshaper", "lwgeom", "ggspatial"
  )
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages],repos='http://cran.us.r-project.org')
  }
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
} #LOAD PACKAGES

# Custom formatters for longitude and latitude
lon_formatter <- function(x) {
  paste0(abs(x), "°", ifelse(x < 0, "W", "E"))
}
lat_formatter <- function(x) {
  paste0(abs(x), "°", ifelse(x < 0, "S", "N"))
}

## --- 1) Breslow baseline cumulative hazard (handles ties) --------------------
breslow_cumhaz_at_time <- function(time, status, risk) {
  stopifnot(length(time)==length(status), length(status)==length(risk))
  ord <- order(time, decreasing = FALSE)
  t <- time[ord]; d <- status[ord]; r <- risk[ord]
  
  # unique event/censor times in order
  ut <- unique(t)
  
  # index of first occurrence of each unique time in the ordered vectors
  first_idx <- match(ut, t)
  
  # risk set sum at each unique time u: sum of risk for subjects with time >= u
  # Using reverse cumulative sum and grabbing the value at first_idx
  risk_cumrev <- rev(cumsum(rev(r)))
  denom_ut <- risk_cumrev[first_idx]
  
  # number of events at each unique time
  d_ut <- as.numeric(tapply(d, t, sum))[match(ut, sort(unique(t)))]
  
  # Breslow hazard increment at each unique time
  inc_ut <- ifelse(denom_ut > 0, d_ut / denom_ut, 0)
  
  # cumulative baseline hazard evaluated at each individual's time
  cumhaz_ut <- cumsum(inc_ut)
  cumhaz_i  <- cumhaz_ut[match(t, ut)]
  
  # return on original order
  cumhaz <- numeric(length(time))
  cumhaz[ord] <- cumhaz_i
  cumhaz
}

eval_FEM_basis_all <- function(FEMbasis, locations) {
  stopifnot(is.matrix(locations) || is.data.frame(locations))
  locations <- as.matrix(locations)
  
  nbasis <- FEMbasis$nbasis
  if (is.null(nbasis)) stop("FEMbasis$nbasis not found. Is this a valid FEMbasis?")
  
  res <- matrix(NA_real_, nrow = nrow(locations), ncol = nbasis)
  
  # Pre-allocate a coefficient vector and reuse it
  coeffs <- numeric(nbasis)
  
  for (j in seq_len(nbasis)) {
    coeffs[] <- 0
    coeffs[j] <- 1
    FEMfun <- fdaPDE::FEM(coeffs, FEMbasis)
    # eval.FEM returns a vector of length nrow(locations)
    res[, j] <- fdaPDE::eval.FEM(FEMfun, locations)
  }
  
  colnames(res) <- paste0("phi", seq_len(nbasis))
  rownames(res) <- paste0("x", seq_len(nrow(locations)))
  res
}

active = read.csv(paste(PATH, "active.csv", sep = ""))
trigger = read.csv(paste(PATH, "trigger_mod.csv", sep = ""))
head(active)
head(trigger)
options(digits.secs = 3)
trigger$timestamp = as.POSIXct(trigger$timestamp, format = "%Y-%m-%d %H:%M:%OS", tz = "UTC")
library(tidytable)
df <- active %>%
  left_join(trigger %>% dplyr::select(user_id, timestamp), by = "user_id")
min(df$timestamp,na.rm = TRUE)
max(df$timestamp,na.rm = TRUE)
df$status = as.numeric(!is.na(df$timestamp))
df[df$status==0, c("timestamp")] = max(df$timestamp,na.rm = TRUE)
head(df)
heartquake_timestamp= as.POSIXct("18-02-2025 02:22:19", format = "%d-%m-%Y %H:%M:%OS", tz = "UTC")
# constants
eq_lat <- 40.8313
eq_lon <- 14.1423
eq_depth <- 2 # km
R <- 6371     # Earth radius in meters

# haversine formula for surface distance (arc length)
haversine <- function(lat1, lon1, lat2, lon2, R) {
  dlat <- (lat2 - lat1) * pi/180
  dlon <- (lon2 - lon1) * pi/180
  lat1 <- lat1 * pi/180
  lat2 <- lat2 * pi/180
  
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  return(c)  # central angle (radians)
}
df$time = as.numeric(df$timestamp - heartquake_timestamp)
# compute distances
df$dist <- sapply(1:nrow(df), function(i) {
  # central angle
  delta_sigma <- haversine(df$lat[i], df$lon[i], eq_lat, eq_lon, R)
  
  # apply formula
  sqrt(eq_depth^2 + (2 * R * sin(delta_sigma / 2))^2)
})
max_dist_trigger = max(df%>% filter(status == 1) %>% pull(dist))
df = df %>% filter(dist <= 5)#max_dist_trigger
nrow(df)
1-sum(df$status)/nrow(df)
hist(df$dist)
hist(df %>% filter(status == 1) %>% pull(time))

naples <- st_read(paste(PATH, "SITRC_COMUNI_CAMPANIA.shp", sep = ""))
st_crs(naples) <- 32633
# Do NOT force a CRS unless you're 100% sure. Just transform to 4326:
naples_ll <- st_transform(naples, 4326)

# --- 2) Points (df) in lon/lat -> sf (x=lon, y=lat) with CRS 4326 ---
df_sf <- st_as_sf(df, coords = c("lon","lat"), crs = 4326)
df_ll <- df_sf

# --- 3) Bbox of your points (+ small buffer in degrees) ---
points_bbox <- st_as_sfc(st_bbox(df_ll), crs = 4326)
points_bbox <- st_buffer(points_bbox, dist = 0.05)  # ~5–6 km at these lats

# --- 4) Clip shapefile to speed up plotting ---
naples_clipped <- st_intersection(naples_ll, points_bbox)
naples_boundary = st_boundary(st_union(naples_clipped))
naples_bd_simp = rmapshaper::ms_simplify(naples_boundary, keep = 0.025) # 0.05
bd_coords = st_coordinates(naples_bd_simp)[,-3]
# avoid repetead nodes
bd_coords = bd_coords[-1,]
# check --- 
plot(bd_coords)
text(bd_coords[,1], bd_coords[,2], labels = 1:nrow(bd_coords), pos = 3, cex = 0.7)
points(df$lon, df$lat, col = "red", pch = 16, cex = 0.4)

bd_coords = rbind(bd_coords[1:21,],
                  c(14.08779, 40.85),
                  bd_coords[22:52,],
                  c(14.19672, 40.85),
                  bd_coords[53:59,])
plot(bd_coords)
text(bd_coords[,1], bd_coords[,2], labels = 1:nrow(bd_coords), pos = 3, cex = 0.7)
points(df$lon, df$lat, col = "red", pch = 16, cex = 0.4)

bd_coords = bd_coords[-c(55:59,61,1:14, 16:21),]
plot(bd_coords)
text(bd_coords[,1], bd_coords[,2], labels = 1:nrow(bd_coords), pos = 3, cex = 0.7)
points(bd_coords[1,1], bd_coords[1,2], col="red", pch=16)
points(bd_coords[2,1], bd_coords[2,2], col="blue", pch=16)
points(bd_coords[nrow(bd_coords),1], bd_coords[nrow(bd_coords),2], pch=16)

segments = cbind(1:(nrow(bd_coords)-1), 2:nrow(bd_coords))
segments =  rbind(segments, c(nrow(bd_coords), 1))
mesh_2 = fdaPDE::create.mesh.2D(nodes = bd_coords, segments = segments)
plot(mesh_2)


# install.packages(c("sf", "ggplot2", "rnaturalearth", "rnaturalearthdata", "lwgeom"))
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(lwgeom)   # for st_polygonize (more robust if segments cross/touch)

# --- 0) Assumptions ---
# - bd_coords is an ordered matrix/data.frame of boundary coords in lon/lat (CRS 4326)
# - df has columns lon/lat (optional, only if you want to show your points)

# --- 1) Build a polygon from your boundary coords ---
bdm <- as.matrix(bd_coords)

# Ensure the ring is closed
if (!all(bdm[1, ] == bdm[nrow(bdm), ])) {
  bdm <- rbind(bdm, bdm[1, ])
}

# Simple polygon from the ring
final_area <- st_sfc(st_polygon(list(bdm)), crs = 4326) |>
  st_make_valid()

# If your boundary was created from segments and might have topological quirks,
# uncomment this more robust route using polygonize:
# segments <- cbind(1:(nrow(bd_coords)-1), 2:nrow(bd_coords))
# segments <- rbind(segments, c(nrow(bd_coords), 1))
# lines_list <- lapply(seq_len(nrow(segments)), function(i) as.matrix(bd_coords[segments[i, ], ]))
# ml <- st_sfc(st_multilinestring(lines_list), crs = 4326)
# final_area <- lwgeom::st_polygonize(ml) |>
#   st_collection_extract("POLYGON") |>
#   st_make_valid()

# --- 2) Get Italy outline ---
italy <-ne_countries(scale = "medium", country = "Italy", returnclass = "sf")

# Align CRS just in case
final_area <- st_transform(final_area, st_crs(italy))

# --- 3) (Optional) Points as sf to plot too ---
# df_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)

# --- 4) Plot: Italy base + highlighted area ---
p1 <- ggplot() +
  geom_sf(data = italy, fill = "grey95", color = "grey70", linewidth = 0.3) +
  geom_sf(data = final_area, fill = "#E41A1C", color = "#E41A1C", linewidth = 0.25, alpha = 0.95) +
  # geom_sf(data = df_sf, size = 0.6, alpha = 0.7) +  # <- uncomment if you want to show your points
  annotate("rect",
           xmin = 13, xmax = 15, ymin = 40, ymax = 41.5,
           fill = NA, color = "black", linewidth = 0.6) +
  coord_sf(xlim = st_bbox(italy)[c("xmin", "xmax")],
           ylim = st_bbox(italy)[c("ymin", "ymax")],
           expand = FALSE) +
  theme_void()

p1


p2 <- ggplot() +
  geom_sf(data = italy, fill = "grey95", color = "grey70", linewidth = 0.3) +
  geom_sf(data = final_area,
          fill = "#E41A1C",      # vivid red
          color = "black",        # dark outline for contrast
          linewidth = 0.5,
          alpha = 0.9) +
  coord_sf(xlim = c(13, 15), ylim = c(40, 41.5), expand = FALSE) + # zoom on southern Italy
  theme_void() + 
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))


p2
library(gridExtra)
p = grid.arrange(p1, p2, ncol=2)


# # --- 1) Build an sf rectangle matching p2 limits (CRS 4326) ---
# zoom_bbox <- st_as_sfc(st_bbox(c(xmin = 13, xmax = 15,
#                                  ymin = 40, ymax = 41.5), 
#                                crs = 4326))
# zoom_bbox <- st_transform(zoom_bbox, st_crs(italy))  # match map CRS
# 
# # --- 2) Add the box to p1 (double-stroke for visibility) ---
# p1_box <- ggplot() +
#   geom_sf(data = italy, fill = "grey95", color = "grey70", linewidth = 0.3) +
#   # your highlighted area
#   geom_sf(data = final_area, fill = "#E41A1C", color = "#E41A1C", linewidth = 0.25, alpha = 0.95) +
#   # thick white halo under the box to pop over the basemap
#   geom_sf(data = zoom_bbox, fill = NA, color = "white", linewidth = 1.6) +
#   # black outline box
#   geom_sf(data = zoom_bbox, fill = NA, color = "black", linewidth = 0.7) +
#   coord_sf(xlim = st_bbox(italy)[c("xmin", "xmax")],
#            ylim = st_bbox(italy)[c("ymin", "ymax")],
#            expand = FALSE) +
#   theme_void()
# 
# grid.arrange(p1_box, p2, ncol=2)
# --- 5) Save image ---
ggsave("italy_final_area_highlight.png", p, width = 7, height = 8, dpi = 300)
