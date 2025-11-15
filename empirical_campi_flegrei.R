PATH = ""
setwd(PATH)
{
  # Package names
  packages <- c(
    "ggplot2", "survival", "fdaPDE", "tidytable", "ggspatial","doRNG",
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
alpha = .45
gamma = .55
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


#### load boundary
library(sf)
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
mesh_2_boundary = fdaPDE::create.mesh.2D(nodes = bd_coords, segments = segments)
plot(mesh_2_boundary)

# Refine mesh based on sample size
mesh_2_pts <- fdaPDE::refine.mesh.2D(
  mesh_2_boundary,
  minimum_angle = 30,
  maximum_area = 0.000001
)

plot(mesh_2_pts)

# Compute total area
int_triangles <- mesh_2_pts$triangles
int_coords_mesh <- mesh_2_pts$nodes
int_triangle_area <- function(tri) {
  p1 <- int_coords_mesh[tri[1], ]
  p2 <- int_coords_mesh[tri[2], ]
  p3 <- int_coords_mesh[tri[3], ]
  0.5 * abs(det(matrix(c(p2 - p1, p3 - p1), ncol = 2)))
}
int_triangle_barycenter <- function(tri) {
  pts <- int_coords_mesh[tri, , drop = FALSE]
  colMeans(pts)  # (p1 + p2 + p3) / 3
}
int_areas <- apply(int_triangles, 1, int_triangle_area)
int_barycenters <- t(apply(int_triangles, 1, int_triangle_barycenter))

integrate_over_omega_from_barycenters=function(f_vals){
  sum(int_areas*f_vals)
}


(AREA =  integrate_over_omega_from_barycenters(rep(1,nrow(int_barycenters))))


# controllare "area totale" della regione prima di scegliere 'maximum_area'
mesh_2 = fdaPDE::refine.mesh.2D(mesh_2_boundary, minimum_angle = 30, maximum_area = 0.1*AREA/(nrow(df)^alpha))
plot(mesh_2, pch=".")
points(df$lon, df$lat, col = "red", pch = 16, cex = 0.4)


########################### creare boundary from segments
pts_xy <- as.matrix(bd_coords)

ls_sfc <- st_sfc(
  lapply(seq_len(nrow(segments)), function(k) {
    st_linestring(pts_xy[segments[k, ], , drop = FALSE])
  }),
  crs = 4326  # set your CRS if needed
)

# dissolve tiny gaps and connect edges
ls_union <- st_union(ls_sfc)

# make polygons from the line network
polys <- st_polygonize(ls_union)

polys_only <- st_collection_extract(polys, "POLYGON")  # or use 6 instead of "POLYGON"

# now take the polygon boundary (as MULTILINESTRING)
bnd <- st_boundary(polys_only)

plot(polys, col = NA, border = 'grey')
plot(bnd, add = TRUE, lwd = 2)
plot(bnd)
points(mesh_2$nodes, col = "red", pch = 16, cex = 0.4)

# mesh_2_inla = inla.mesh.create(
#   loc = mesh_2$nodes, 
#   tv = mesh_2$triangles,
#   crs = 4326
# )
# plot(mesh_2_inla)
# pts <- sf::st_sample(polys, 10000, type = "regular")
# points(st_coordinates(pts), col = "red", pch = 16, cex = 0.4)




# --- nodes ---
nodes_df <- as.data.frame(mesh_2$nodes)
names(nodes_df)[1:2] <- c("x","y")

# --- ALL mesh edges from mesh_2$edges ---
edges_df <- as.data.frame(mesh_2$edges)
# first two columns are the node indices (works whether there are 2 or 3+ cols)
names(edges_df)[1:2] <- c("v1","v2")

edges_sfc <- lapply(seq_len(nrow(edges_df)), function(i){
  idx <- as.integer(edges_df[i, c("v1","v2")])
  st_linestring(as.matrix(nodes_df[idx, c("x","y")]))
}) |> st_sfc()

mesh_edges_sf <- st_sf(geometry = edges_sfc)

# --- (optional) boundary in bold from mesh_2$segments ---
segs_df <- as.data.frame(mesh_2$segments)
names(segs_df)[1:2] <- c("v1","v2")
segs_sfc <- lapply(seq_len(nrow(segs_df)), function(i){
  idx <- as.integer(segs_df[i, c("v1","v2")])
  st_linestring(as.matrix(nodes_df[idx, c("x","y")]))
}) |> st_sfc()
boundary_sf <- st_sf(geometry = segs_sfc)


# --- plot ---
p_mesh <- ggplot() +
  # all mesh edges (internal + boundary)
  geom_sf(data = mesh_edges_sf, color = "grey60", linewidth = 0.5) +
  # your points, blue if status==1 else red
  geom_point(
    data = df,
    aes(x = lon, y = lat, color = status == 0),
    size = 2
  ) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "blue"),
    name   = "Censored"
  ) +
  geom_sf(data = boundary_sf, color = "black", linewidth = 1) +
  geom_point(aes(x = eq_lon, y = eq_lat),
             shape = 42, fill = "black", color = "black",
             size = 10, stroke = 1.2, inherit.aes = FALSE) +
  coord_sf(expand=FALSE)+
  # correct way to label axes
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )
p_mesh

######################## NEW PLOT MESH

st_crs(mesh_edges_sf) <- 4326
st_crs(boundary_sf)   <- 4326

# also convert your points to sf in WGS84
pts_sf <- st_as_sf(df, coords = c("lon","lat"), crs = 4326)
eq_sf  <- st_sf(geometry = st_sfc(st_point(c(eq_lon, eq_lat)), crs = 4326))

# 2) Choose a **metric** CRS near your study area (UTM zone by centroid)
ctr <- st_coordinates(st_centroid(st_union(boundary_sf)))
utm_zone <- floor((ctr[1] + 180) / 6) + 1
epsg <- if (ctr[2] >= 0) 32600 + utm_zone else 32700 + utm_zone  # 326xx for N, 327xx for S
target_crs <- st_crs(epsg)

# 3) Reproject everything
mesh_edges_m <- st_transform(mesh_edges_sf, target_crs)
boundary_m   <- st_transform(boundary_sf,   target_crs)
pts_m        <- st_transform(pts_sf,        target_crs)
eq_m         <- st_transform(eq_sf,         target_crs)

# 4) Plot in meters (axis units are now meters; aspect is true)
p_mesh <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = status == 0), size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "Censored") +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  geom_sf(data = eq_m, shape = 42, fill = "black", color = "black", size = 20, stroke = 1.2) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  ) +
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr")

p_mesh


########################  add radiis
radii_km <- c(1,2)              # tweak as you like
radii_m  <- radii_km * 1000

# 1) make circular buffers around the epicenter (polygons, in meters)
rings_poly <- st_buffer(eq_m, dist = radii_m)

# 2) turn those polygons into linework so they plot as outlines only
# Base R (senza purrr)
rings_line <- do.call(rbind, lapply(seq_along(radii_m), function(i){
  g <- st_boundary(st_buffer(eq_m, dist = radii_m[i]))   # una geometria per volta
  st_sf(radius_km = radii_km[i], geometry = st_geometry(g))
}))

# (optional) quick label positions due east of the epicenter
eq_xy <- st_coordinates(eq_m)[1, ]
labels_df <- data.frame(
  x = eq_xy[1] + radii_m,   # easting + radius
  y = eq_xy[2],             # same northing
  lab = paste0(radii_km, " km")
)
labels_sf <- st_as_sf(labels_df, coords = c("x","y"), crs = st_crs(eq_m))

# # 3) add to your plot
# p_mesh +
#   geom_sf(data = rings_line, linetype = "solid", linewidth = 1, alpha = 0.8) +
#   geom_sf_text(data = labels_sf, aes(label = lab), nudge_y = 0, size = 3)

ggsave(
  filename = "francesco_mesh2.png",
  plot     = p_mesh,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

########################## time 

p_time <- ggplot() +
  # all mesh edges (internal + boundary)
  geom_sf(data = mesh_edges_sf, color = "grey60", linewidth = 0.5) +
  # boundary emphasized (optional)
  # your points, blue if status==1 else red
  geom_point(data = df,
             aes(x = lon, y = lat, color = time),
             size = 2) +
  geom_sf(data = boundary_sf, color = "black", linewidth = 1) +
  geom_point(aes(x = eq_lon, y = eq_lat),
             shape = 42, fill = "black", color = "black",
             size = 20, stroke = 1.2, inherit.aes = FALSE) +
  coord_sf(expand=FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_viridis_c(option = "plasma", name = "Seconds")
# annotate("rect",
#          xmin = min(df$lon)-.001, xmax = max(df$lon)+.001,
#          ymin = min(df$lat)-.001, ymax = max(df$lat)+.001,
#          color = "black", fill = NA, linewidth = 1)+
#theme(legend.position = "none")

#ggsave("mesh_full_with_points.png", p, width = 7, height = 6, dpi = 300)
p_time
ggsave(
  filename = "francesco_time.png",
  plot     = p_time,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

######################## NEW PLOT TIME


p_time <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = time), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  geom_sf(data = eq_m, shape = 42, fill = "black", color = "black", size = 20, stroke = 1.2) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_viridis_c(option = "inferno", name = "Seconds") +
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr")+
  geom_sf(data = rings_line, linetype = "solid", linewidth = 1, alpha = 0.8) +
  geom_sf_text(data = labels_sf, aes(label = lab), nudge_y = 0, size = 3)


p_time

ggsave(
  filename = "francesco_time2.png",
  plot     = p_time,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

# p_time_status <- ggplot() +
#   # all mesh edges (internal + boundary)
#   geom_sf(data = mesh_edges_sf, color = "grey60", linewidth = 0.5) +
#   
#   # points: color by time, shape by status
#   geom_point(
#     data = df,
#     aes(x = lon, y = lat, color = time, shape = factor(status)),
#     size = 3
#   ) +
#   
#   # boundary emphasized
#   geom_sf(data = boundary_sf, color = "black", linewidth = 1) +
#   geom_point(aes(x = eq_lon, y = eq_lat),
#              shape = 42, fill = "black", color = "black",
#              size = 10, stroke = 1.2, inherit.aes = FALSE) +
#   
#   coord_sf(expand = FALSE) +
#   labs(x = "Longitude", y = "Latitude") +
#   scale_x_continuous(labels = lon_formatter) +
#   scale_y_continuous(labels = lat_formatter) +
#   
#   # scales
#   scale_color_viridis_c(option = "plasma", name = "Seconds") +
#   scale_shape_manual(
#     name = "Status",
#     values = c("0" = 1, "1" = 16),    # open circle for censored, filled circle for uncensored
#     labels = c("0" = "Censored", "1" = "Uncensored")
#   ) +
#   
#   # theme
#   theme_minimal() +
#   theme(
#     axis.title = element_text(size = 12, face = "bold"),
#     axis.text  = element_text(size = 10)
#   )
# p_time_status
# ggsave(
#   filename = "francesco_time_status.png",
#   plot     = p_time_status,
#   width    = 10,
#   height   = 8,
#   dpi      = 300,
#   units    = "in"
# )
################ PLOT F TRUE
# find the global min/max across all three variables



FEMbasis = fdaPDE::create.FEM.basis(mesh = mesh_2)
# mass
M0 = fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)
# stiff
M1 = fdaPDE:::CPP_get.FEM.Stiff.Matrix(FEMbasis)
# penalty ?
K4 = M1%*% solve(M0) %*% M1

# fem    <- INLA::inla.mesh.fem(mesh_2_inla, order = 2)
# M0     <- fem$c0
# K4     <- fem$g2

# Sum-to-zero constraint via mass matrix
m_vec <- rowSums(M0)
Z     <- pracma::null(t(m_vec))
Kred  <- t(Z) %*% K4 %*% Z

coords = df[,c("lon", "lat")]
## 4. Design matrices (link data to basis) -----------------------
A    <- eval_FEM_basis_all(FEMbasis, coords)
#A    <- INLA::inla.spde.make.A(mesh_2_inla, loc = as.matrix(coords))
Phi  <- A %*% Z                                  # n × (K-1)  constrained basis

Xcov <- as.matrix(df[, c("dist")])         # add your covariates here


cumexp <- function(eta) rev(cumsum(rev(exp(eta))))  # ∑_{j ∈ R(t_i)} e^{η_j}

negPL_pen <- function(par, X, Phi, Kred, time, status, λ) {
  p      <- ncol(X)
  beta   <- par[seq_len(p)]
  theta  <- par[-seq_len(p)]
  eta    <- as.vector(X %*% beta + Phi %*% theta)
  
  ord    <- order(time)
  eta    <- eta[ord]; status <- status[ord]
  
  logrisk <- log(cumexp(eta))
  ll      <- sum(status * (eta - logrisk))
  
  pen     <- -0.5 * λ * as.numeric(crossprod(theta, Kred %*% theta))  # minus = add in objective
  as.numeric(-ll - pen)   # NEGATIVE objective for minimiser
}

negPL_pen_grad <- function(par, X, Phi, Kred, time, status, λ) {
  p      <- ncol(X)
  beta   <- par[seq_len(p)]
  theta  <- par[-seq_len(p)]
  eta    <- as.vector(X %*% beta + Phi %*% theta)
  
  ord    <- order(time)
  eta    <- eta[ord]; status <- status[ord]
  Xo     <- X[ord, , drop = FALSE]
  Phio   <- Phi[ord, , drop = FALSE]
  
  risk   <- cumexp(eta)
  
  # Efficient risk-set means using reverse cumsums
  rev_cumsum <- function(mat) apply(mat[nrow(mat):1, , drop = FALSE], 2, cumsum)[nrow(mat):1, ]
  Xbar   <- rev_cumsum(Xo * exp(eta)) / risk
  Pbar   <- rev_cumsum(Phio * exp(eta)) / risk
  
  sco_beta  <- colSums(status * (Xo - Xbar))
  sco_theta <- colSums(status * (Phio - Pbar)) - λ * drop(Kred %*% theta)
  
  -c(sco_beta, sco_theta)   # NEGATIVE gradient for minimiser
}

# ------- CV helpers (PL + CV-PLD selection) -------

# (normalized) log partial likelihood, matching the definition in the text
logPL_norm <- function(beta, theta, X, Phi, time, status) {
  n   <- nrow(X)
  eta <- as.vector(X %*% beta + Phi %*% theta)
  
  ord <- order(time)
  eta <- eta[ord]; status <- status[ord]
  
  logrisk <- log(cumexp(eta))               # log ∑_{j∈R(t_i)} e^{η_j}
  # note: log((1/n)∑) = logrisk - log(n)
  ll_vec  <- status * (eta - (logrisk - log(n)))
  (1/n) * sum(ll_vec)
}

# Fit penalized Cox (wrapper around your objective/gradient)
fit_pen_cox <- function(X, Phi, Kred, time, status, lambda, start = NULL, control = list(), hessian = FALSE) {
  pX <- ncol(X); pPhi <- ncol(Phi)
  if (is.null(start)) start <- rep(0, pX + pPhi)
  
  opt <- optim(
    par    = start,
    fn     = negPL_pen,
    gr     = negPL_pen_grad,
    method = "BFGS",
    X      = X,
    Phi    = Phi,
    Kred   = Kred,
    time   = time,
    status = status,
    λ      = lambda,
    control = modifyList(list(maxit = 500, reltol = 1e-9), control),
    hessian = hessian
  )
  list(
    beta  = opt$par[seq_len(pX)],
    theta = opt$par[-seq_len(pX)],
    value = opt$value,
    convergence = opt$convergence,
    hessian = opt$hessian
  )
}

# Choose lambda by K-fold CV-PLD (difference of PL on full vs retained data)
select_lambda_cvpld <- function(
    X, Phi, Kred, time, status,
    K = 5,
    lambda_grid = NULL,
    seed = 123
) {
  set.seed(seed)
  
  n <- nrow(X)
  # default: log-spaced grid around n^{-1/2}
  if (is.null(lambda_grid)) {
    lam0 <- n^(-0.5)
    #lambda_grid = seq(0.01, 0.99, length.out = 10)*lam0
    lambda_grid <- seq(0.1*lam0, 10*lam0, length.out = 10)
  }
  
  # random folds
  folds <- sample(rep(1:K, length.out = n))
  
  cvpld <- numeric(length(lambda_grid))
  
  for (g in seq_along(lambda_grid)) {
    print(g)
    lam <- lambda_grid[g]
    score_g <- 0
    warm <- NULL
    
    for (k in 1:K) {
      idx_tr <- which(folds != k)
      
      fit_k <- fit_pen_cox(
        X   = X[idx_tr, , drop = FALSE],
        Phi = Phi[idx_tr, , drop = FALSE],
        Kred = Kred,
        time = time[idx_tr],
        status = status[idx_tr],
        lambda = lam,
        start = warm
      )
      warm <- c(fit_k$beta, fit_k$theta)  # warm start across folds
      
      # Evaluate on full data (D) and retained data (D_{-k}) with the SAME estimates
      ll_full <- logPL_norm(fit_k$beta, fit_k$theta, X, Phi, time, status)
      ll_tr   <- logPL_norm(
        fit_k$beta, fit_k$theta,
        X[idx_tr, , drop = FALSE],
        Phi[idx_tr, , drop = FALSE],
        time[idx_tr], status[idx_tr]
      )
      
      score_g <- score_g + (-2) * (ll_full - ll_tr)
    }
    
    cvpld[g] <- score_g
  }
  
  list(
    lambda_opt = lambda_grid[which.min(cvpld)],
    grid       = lambda_grid,
    cvpld      = cvpld
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ===========================================================
# Fast, parallel CV for penalized Cox with selectable criterion
#   criterion = "val"  -> sum_k  -2 * logPL(D_k | θ̂_{-k})
#   criterion = "diff" -> sum_k  -2 * [ logPL(D | θ̂_{-k}) - logPL(D_{-k} | θ̂_{-k}) ]
# ===========================================================

select_lambda_cvpld_fast <- function(
    X, Phi, Kred, time, status,
    K = 5,
    lambda_grid = NULL,
    seed = 123,
    workers = max(1, parallel::detectCores() - 1),
    pkgs = NULL,             # e.g. c("survival","Matrix")
    source_file = NULL,      # e.g. "path/to/defs.R"
    extra_exports = NULL,    # character vector of extra symbols to export
    criterion = c("val","diff")  # "val" = validation-only, "diff" = match original
) {
  criterion <- match.arg(criterion)
  stopifnot(is.matrix(X), is.matrix(Phi), is.matrix(Kred))
  stopifnot(nrow(X) == length(time), length(status) == length(time))
  
  # required functions present?
  req_fns <- c("fit_pen_cox","logPL_norm","negPL_pen","negPL_pen_grad","cumexp")
  miss <- req_fns[!vapply(req_fns, exists, logical(1), mode = "function")]
  if (length(miss)) stop("Missing functions: ", paste(miss, collapse = ", "))
  
  set.seed(seed)
  n <- nrow(X)
  if (is.null(lambda_grid)) {
    lam0 <- n^(-0.5)
    lambda_grid <- seq(0.1 * lam0, 10 * lam0, length.out = 10)
  } else {
    lambda_grid <- as.numeric(lambda_grid)
  }
  
  # fixed folds
  folds <- sample(rep(seq_len(K), length.out = n))
  idx_tr_list  <- lapply(seq_len(K), function(k) which(folds != k))
  idx_val_list <- lapply(seq_len(K), function(k) which(folds == k))
  
  # parallel setup
  use_parallel <- workers > 1
  if (use_parallel) {
    if (!requireNamespace("foreach", quietly = TRUE) ||
        !requireNamespace("doParallel", quietly = TRUE) ||
        !requireNamespace("doRNG", quietly = TRUE))
      stop("Need packages: foreach, doParallel, doRNG.")
    
    cl <- parallel::makeCluster(workers)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    doParallel::registerDoParallel(cl)
    parallel::clusterSetRNGStream(cl, iseed = seed)
    
    to_export <- unique(c(
      "X","Phi","Kred","time","status",
      "fit_pen_cox","logPL_norm","negPL_pen","negPL_pen_grad","cumexp",
      "criterion", extra_exports
    ))
    parallel::clusterExport(cl, varlist = to_export, envir = environment())
    if (!is.null(source_file)) parallel::clusterEvalQ(cl, source(source_file, local = TRUE))
    if (!is.null(pkgs) && length(pkgs))
      parallel::clusterEvalQ(cl, invisible(lapply(pkgs, require, character.only = TRUE)))
  }
  
  warm_list <- vector("list", K)
  cvpld <- numeric(length(lambda_grid))
  
  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    
    one_fold <- function(k) {
      idx_tr  <- idx_tr_list[[k]]
      idx_val <- idx_val_list[[k]]
      
      fit_k <- fit_pen_cox(
        X      = X[idx_tr, , drop = FALSE],
        Phi    = Phi[idx_tr, , drop = FALSE],
        Kred   = Kred,
        time   = time[idx_tr],
        status = status[idx_tr],
        lambda = lam,
        start  = warm_list[[k]]
      )
      
      if (criterion == "val") {
        # Validation-only PL
        ll_val <- logPL_norm(
          fit_k$beta, fit_k$theta,
          X[idx_val, , drop = FALSE],
          Phi[idx_val, , drop = FALSE],
          time[idx_val], status[idx_val]
        )
        crit <- (-2) * ll_val
      } else { # "diff" -> reproduce original scoring exactly
        ll_full <- logPL_norm(fit_k$beta, fit_k$theta, X, Phi, time, status)
        ll_tr   <- logPL_norm(
          fit_k$beta, fit_k$theta,
          X[idx_tr, , drop = FALSE],
          Phi[idx_tr, , drop = FALSE],
          time[idx_tr], status[idx_tr]
        )
        crit <- (-2) * (ll_full - ll_tr)
      }
      
      list(crit = crit, warm = c(fit_k$beta, fit_k$theta))
    }
    
    fold_res <- if (use_parallel) {
      doRNG::`%dorng%`(
        foreach::foreach(
          k = seq_len(K),
          .export   = c("fit_pen_cox","logPL_norm","negPL_pen","negPL_pen_grad","cumexp","criterion"),
          .packages = pkgs %||% character()
        ),
        one_fold(k)
      )
    } else {
      lapply(seq_len(K), one_fold)
    }
    
    cvpld[i]   <- sum(vapply(fold_res, `[[`, numeric(1), "crit"))
    warm_list  <- lapply(fold_res, `[[`, "warm")
  }
  
  list(
    lambda_opt = lambda_grid[which.min(cvpld)],
    grid       = lambda_grid,
    cvpld      = cvpld
  )
}



## 7. Fit with optim() (BFGS) ------------------------------------
# ---------- Standard PH (no spatial) ----------
fit_std  <- survival::coxph(survival::Surv(time, status) ~ dist, data = df)
bhat_std <- stats::coef(fit_std)


#(lambda_grid =  exp(seq(log(0.01), log(5), length.out = 10))*(nrow(df)^-1))
#(lambda_grid <- exp(seq(log(0.05), log(50), length.out = 10)) * (nrow(df)^-gamma) * AREA)  #exp(seq(log(10^-5), log(10^5), length.out = 10))
(lambda_grid <- exp(seq(log(0.05), log(50), length.out = 10)) * (nrow(df)^-gamma) * AREA) 

#run again only if you want to check that lambda = lambda_grid[10] is selected
# cvsel <- select_lambda_cvpld(
#   X = Xcov,
#   Phi = Phi,
#   Kred = Kred,
#   time = df$time,
#   status = df$status,
#   K = 5,                        # adjust if desired
#   lambda_grid = lambda_grid,           # or provide numeric vector
#   seed = 1234   # CV fold seed
# )
# (lambda <- cvsel$lambda_opt)

# RUN TO CHECK THAT THE PICKED ONE IS lambda_grid[10]
# start_time = Sys.time()
# cvsel <- select_lambda_cvpld_fast(
#   X = Xcov, Phi = Phi, Kred = as.matrix(Kred),
#   time = df$time, status = df$status,
#   K = 5,
#   lambda_grid = lambda_grid,
#   seed = 1000,
#   workers = min(10,parallel::detectCores() - 1),   # parallel
#   criterion = "diff"
# )
# Sys.time()-start_time
# (lambda <- cvsel$lambda_opt)

lambda = lambda_grid[10]
# ----- Final fit at selected lambda -----
start.time <- Sys.time()
start <- rep(0, ncol(Xcov) + ncol(Phi))
fit   <- fit_pen_cox(
  X = Xcov, Phi = Phi, Kred = Kred,
  time = df$time, status = df$status,
  lambda = lambda, start = start, 
  hessian = TRUE
)
print(Sys.time() - start.time)
beta_hat  <- fit$beta
theta_hat    <- fit$theta

p <- ncol(Xcov)

H <- fit$hessian  
# Standard errors for beta
std_beta1 <- sqrt(1/H[1:p, 1:p])
std_beta2 <- sqrt(solve(H)[1:p, 1:p])

# Name them
names(beta_hat) <- colnames(Xcov)
names(std_beta1) <- colnames(Xcov)
names(std_beta2) <- colnames(Xcov)

# View result
print(beta_hat)
print(std_beta1)
print(std_beta2)


# Also write output to a text file
sink("data_francesco9_results.txt")  # Start writing output to a file
cat("beta_hat:\n")
print(beta_hat)
cat("\nstd_beta1 (extract and invert):\n")
print(std_beta1)
cat("\nstd_beta2 (invert and extract):\n")
print(std_beta2)
cat("\nbhat_std (Cox model coefficients):\n")
print(bhat_std)
cat("\nSummary of Cox model:\n")
print(summary(fit_std))
sink()  # stop writing to file

poly <- boundary_sf |>
  st_union() |>
  st_line_merge() |>
  st_polygonize() |>
  st_make_valid()



# mesh_2 = fdaPDE::create.mesh.2D(nodes = bd_coords, segments = segments)
# # Refine mesh based on sample size
# mesh_2 <- fdaPDE::refine.mesh.2D(
#   mesh_2,
#   minimum_angle = 30,
#   maximum_area = 0.001
# )
# 
# plot(mesh_2)



pts <- sf::st_sample(poly, 100000, type = "regular") 
plot(mesh_2, pch=".")
points(sf::st_coordinates(pts), col="black", pch='.')
grid <- as.data.frame(sf::st_coordinates(pts))
rm(pts)
names(grid) = c("x", "y")
Agrid <- eval_FEM_basis_all(FEMbasis, grid[,c("x", "y")])

#Agrid = INLA::inla.spde.make.A(mesh_2_inla, loc = as.matrix(grid[, c("x", "y")]))
grid$f_hat = as.vector(Agrid %*% (Z %*% theta_hat))

p_hat <- ggplot(grid, aes(x = x, y = y, color = f_hat)) +
  geom_point(size = 2) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    name = "Value"
  ) +
  geom_sf(data = boundary_sf, inherit.aes = FALSE,
          color = "black", linewidth = 1) +
  geom_point(aes(x = eq_lon, y = eq_lat),
             shape = 42, fill = "black", color = "black",
             size = 20, stroke = 1.2, inherit.aes = FALSE) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )
p_hat

ggsave(
  filename = "francesco_hat.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)


######################## NEW PLOT HAT

pts_sf <- st_as_sf(grid, coords = c("x","y"), crs = 4326)
eq_sf  <- st_sf(geometry = st_sfc(st_point(c(eq_lon, eq_lat)), crs = 4326))
pts_m        <- st_transform(pts_sf,        target_crs)
range_f <- round(range(grid$f_hat, na.rm = TRUE),1) + c(-.1,+.1)

p_hat <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, aes(color = f_hat), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
#  geom_sf(data = eq_m, shape = 42, fill = "black", color = "black", size = 20, stroke = 1.2) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = range_f,               # ensure full range is shown
    breaks = c(range_f[1], 0, range_f[2]),  # show min, mid, max on legend
    name = "Value"
  )  + 
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr")
  # geom_sf(data = rings_line, linetype = "solid", linewidth = 1, alpha = 0.8) +
  # geom_sf_text(data = labels_sf, aes(label = lab), nudge_y = 0, size = 3)


p_hat

ggsave(
  filename = "francesco_hat2.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)





brks <- c(-1.5, -1.0, -0.5, 0, 0.5, 1.0, 1.5)  # example breaks
brks <- round(seq(-1.5, 1.5, length.out = 20),2)

p_hat <- ggplot() + 
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, aes(color = f_hat), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  #geom_sf(data = eq_m, shape = 42, fill = "black", color = "black", size = 20, stroke = 1.2) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_steps2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    breaks = brks,
    limits = range(brks),
    name = "Value"
  )  + 
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr")
  
 

p_hat =  p_hat + guides(
  color = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(10, "cm"),  # increase to separate labels vertically
    barwidth = unit(0.6, "cm"),
    label.theme = element_text(size = 9)
  )
)

ggsave(
  filename = "francesco_hat3.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)



################ NEW PLOT RISK

grid$dist =  sapply(1:nrow(grid), function(i) {
  # central angle
  delta_sigma <- haversine(grid$y[i], grid$x[i], eq_lat, eq_lon, R)
  
  # apply formula
  sqrt(eq_depth^2 + (2 * R * sin(delta_sigma / 2))^2)
})

grid$risk = exp(grid$dist*beta_hat + grid$f_hat)

pts_sf <- st_as_sf(grid, coords = c("x","y"), crs = 4326)
eq_sf  <- st_sf(geometry = st_sfc(st_point(c(eq_lon, eq_lat)), crs = 4326))
pts_m        <- st_transform(pts_sf,        target_crs)

p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, aes(color = risk), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  geom_sf(data = eq_m, shape = 42, fill = "black", color = "black", size = 20, stroke = 1.2) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_viridis_c(option = "inferno", name = "Risk") +
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr")+
  geom_sf(data = rings_line, linetype = "solid", linewidth = 1, alpha = 0.8) +
  geom_sf_text(data = labels_sf, aes(label = lab), nudge_y = 0, size = 3)


p_risk

ggsave(
  filename = "francesco_risk2.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)



# define breaks as before (adjust range and step size as appropriate for 'risk')
brks_risk <- round(seq(0,
                       max(pts_m$risk, na.rm = TRUE),
                       length.out = 20), 2)

p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, aes(color = risk), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  geom_sf(data = eq_m, shape = 42, fill = "black", color = "black", size = 20, stroke = 1.2) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  ) +
  scale_color_stepsn(
    colors = viridis::viridis(20, option = "inferno"),
    breaks = brks_risk,
    limits = range(brks_risk),
    name = "Risk"
  ) +
  annotation_scale(location = "bl", unit_category = "metric") +
  annotation_north_arrow(location = "tr") +
  geom_sf(data = rings_line, linetype = "solid", linewidth = 1, alpha = 0.8) +
  geom_sf_text(data = labels_sf, aes(label = lab), nudge_y = 0, size = 3)

# optional: adjust colorbar styling like before
p_risk <- p_risk + guides(
  color = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(10, "cm"),
    barwidth = unit(0.6, "cm"),
    label.theme = element_text(size = 9)
  )
)


ggsave(
  filename = "francesco_risk3.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

