#https://data.sfgov.org/Public-Safety/Fire-Department-and-Emergency-Medical-Services-Dis/nuek-vuh3/about_data
PATH <- "C:/Users/utente/Downloads/"
#PATH <- "C:/Users/utente/Nextcloud/01_work/01_research/01_projects/02_SPATIAL_COX/10_scandinavian/"

setwd(PATH)
# --- 1. SETUP & LIBRARIES ---
{
  packages <- c(
    "tidyr", "dplyr", "lubridate", "stringr", "ggplot2", "sf", "tigris", 
    "osmdata", "rmapshaper", "fdaPDE", "smoothr", "survival", "jsonlite", 
    "timeDate", "Matrix", "ggspatial", "scales", "parallel", "Rcpp", 
    "RcppArmadillo", "foreach", "doParallel", "doRNG"
  )
  new_vars <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_vars)) install.packages(new_vars, repos = "http://cran.us.r-project.org")
  invisible(lapply(packages, library, character.only = TRUE))
  
  # Load custom optimization functions
  if (file.exists("script_functions_download.R")) {
    source("script_functions_download.R")
  } else {
    stop("Please ensure 'script_functions_down.R' is in the working directory.")
  }
}

#PATH = "/home/ltedesco/spcoxnew/"
# --- 1. Get San Francisco County Boundary ---
# cb = TRUE roughly clips to the shoreline, removing most of the "Sea"/Bay automatically.
sf_boundary <- counties(state = "CA", cb = TRUE, class = "sf") %>%
  filter(NAME == "San Francisco") %>%
  st_transform(26910) # NAD83 / UTM zone 10N (Meters)



# --- 2. Get WATER Polygons from OpenStreetMap ---
# Instead of parks, we query for water (Lakes, Reservoirs, Ponds)
water_query <- opq(bbox = st_bbox(st_transform(sf_boundary, 4326))) %>%
  add_osm_feature(key = "natural", value = "water") %>%
  osmdata_sf()

# Extract polygons and transform to same CRS
# We create a union of all water bodies to punch holes later
water_polys <- water_query$osm_polygons %>%
  st_transform(26910) %>%
  st_make_valid() %>% 
  st_union()

# --- 3. Subtract Water from the City Boundary ---
# This cuts holes for lakes (e.g., Lake Merced, Stow Lake)
# The 'sea' is mostly handled by the county boundary, but this ensures internal water is gone.
urban_area <- st_difference(sf_boundary, water_polys)




# 4. Plot to verify
plot(st_geometry(urban_area), col = "grey", main = "SF Urban Area (No Parks)")


# --- 4. Smoothing (As before) ---
urban_smooth <- smoothr::smooth(urban_area, method = "ksmooth", smoothness = 2)


# --- 5. Clean Islands and Holes ---

polys <- st_cast(urban_smooth, "POLYGON") %>%
  st_transform(26910)


# A. FILTER ISLANDS: Keep only the Top 2 (Mainland + Treasure Island)
# We calculate area, sort largest-to-smallest, and keep the top 2.
# This guarantees the Farallon Islands (which are much smaller and 3rd/4th) are dropped.
polys <- polys %>%
  mutate(area = as.numeric(st_area(.))) %>%
  arrange(desc(area)) %>% 
  slice(1:2) 

# B. FILTER HOLES: Fill small ponds/pools
# Keep this as is to handle internal lakes like Merced
MIN_HOLE_AREA <- 50000 
polys_clean <- fill_holes(polys, threshold = MIN_HOLE_AREA)

# --- Final Plot ---
plot(st_geometry(polys_clean), col = "lightblue", border = "blue", main = "SF: Land Only (Mainland + Treasure Island)")








# --- 1. PREPARE DATA ---
# Remove empty geometries immediately after simplification
simple_polys <- ms_simplify(polys_clean, keep = 0.01, keep_shapes = TRUE) %>%
  st_cast("POLYGON") 
simple_polys <- simple_polys[!st_is_empty(simple_polys), ]
#plot(simple_polys)
# --- 2. ROBUST EXTRACTION FUNCTION ---
sf_to_mesh_final <- function(poly_sf) {
  
  all_pts <- list()
  all_segs <- list()
  all_holes <- list()
  
  node_cursor <- 0
  
  # Get the list of geometries (each item is a list of rings)
  geom_list <- st_geometry(poly_sf)
  
  for (i in seq_along(geom_list)) {
    rings <- geom_list[[i]]
    
    # Iterate over RINGS (1=Outer, 2+=Holes)
    for (j in seq_along(rings)) {
      
      ring_mat <- rings[[j]]
      
      # --- FIX FOR THE ERROR ---
      # Check if the ring is valid (must be a matrix and have >3 points to be a loop)
      if (!is.matrix(ring_mat) || nrow(ring_mat) < 4) {
        next # Skip this degenerate ring
      }
      
      # A. Extract Boundary Coordinates
      # Remove the repeated last point
      pts <- ring_mat[1:(nrow(ring_mat)-1), , drop=FALSE]
      n <- nrow(pts)
      
      # B. Store Nodes
      all_pts[[length(all_pts) + 1]] <- pts
      
      # C. Build Segments
      idx <- 1:n
      seg_from <- idx
      seg_to <- c(idx[-1], 1)
      all_segs[[length(all_segs) + 1]] <- cbind(seg_from, seg_to) + node_cursor
      
      node_cursor <- node_cursor + n
      
      # D. DETECT HOLES
      if (j > 1) {
        # Create a polygon of just this hole ring
        # We wrap this in tryCatch because collapsed holes might fail polygon creation
        try({
          hole_poly <- st_polygon(list(ring_mat))
          # st_point_on_surface finds a point *inside* the hole
          hole_pt <- st_coordinates(st_point_on_surface(hole_poly))
          
          # Only add if we actually got a point back
          if(nrow(hole_pt) > 0) {
            all_holes[[length(all_holes) + 1]] <- hole_pt[1, 1:2]
          }
        }, silent = TRUE)
      }
    }
  }
  
  # Check if we found any nodes at all
  if (length(all_pts) == 0) stop("No valid geometry found after cleaning.")
  
  # Combine lists
  raw_nodes <- do.call(rbind, all_pts)
  raw_segments <- do.call(rbind, all_segs)
  if(length(all_holes) > 0) {
    raw_holes <- do.call(rbind, all_holes)
  } else {
    raw_holes <- matrix(0, nrow=0, ncol=2)
  }
  
  # --- CLEAN DUPLICATES ---
  # 1. Round coordinates to merge floating-point noise
  nodes_round <- round(raw_nodes, 1) 
  
  # 2. Identify unique nodes
  node_ids <- paste(nodes_round[,1], nodes_round[,2], sep="_")
  unique_ids <- unique(node_ids)
  
  # 3. Create Clean Node List (Original coordinates, but filtered)
  # match(unique_ids, node_ids) gives the index of the *first* occurrence of each unique node
  clean_nodes <- raw_nodes[match(unique_ids, node_ids), ]
  
  # 4. Remap Segments to new unique IDs
  mapping <- match(node_ids, unique_ids)
  clean_segments <- matrix(mapping[raw_segments], ncol=2)
  
  # 5. Remove degenerate segments (Start == End)
  # These happen if a tiny edge got rounded into a single point
  keep_segs <- clean_segments[,1] != clean_segments[,2]
  clean_segments <- clean_segments[keep_segs, ]
  
  return(list(nodes = clean_nodes, segments = clean_segments, holes = raw_holes))
}

# --- 3. GENERATE MESH ---
mesh_input <- sf_to_mesh_final(simple_polys)

# Create the mesh with holes
mesh <- create.mesh.2D(
  nodes = mesh_input$nodes, 
  segments = mesh_input$segments,
  holes = mesh_input$holes, 
  order = 1
)

# --- 4. VERIFY ---
plot(mesh, main = "Final SF Mesh")
# Optional: Plot the hole seeds to see where the parks were cut
if (nrow(mesh_input$holes) > 0) {
  points(mesh_input$holes, col="red", pch=19, cex=0.5)
}

final_mesh2 = 
  fdaPDE::refine.mesh.2D(mesh, 
                         minimum_angle = 20, maximum_area = 300000)
plot(final_mesh2)
nrow(final_mesh2$nodes)


df <- read.cnodesdf <- read.csv(
  paste0(PATH, "Fire_Department_and_Emergency_Medical_Services_Dispatched_Calls_for_Service_20260120.csv")
)
df$Call.Number
table(df$Call.Type)
table(df$Original.Priority)
table(df$Unit.sequence.in.call.dispatch)
# coordinates
df <- df %>%
  mutate(
    lon = as.numeric(str_extract(case_location, "-?\\d+\\.\\d+(?=\\s)")),
    lat = as.numeric(str_extract(case_location, "(?<=\\s)-?\\d+\\.\\d+"))
  )

#plot(df$lat,df$lon)
# 1. Clean and Prepare the Points
# Ensure we drop rows where extraction failed (NAs)
df <- df %>%
  filter(!is.na(lon) & !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) # keep original columns
df_projected <- st_transform(df, st_crs(simple_polys))
df_inside_mesh <- st_filter(df_projected, simple_polys)


# B. Count Check
cat("Total calls:", nrow(df), "\n")
cat("Calls inside mesh:", nrow(df_inside_mesh), "\n")
cat("Dropped (Outside domain):", nrow(df) - nrow(df_inside_mesh), "\n")

df = df_inside_mesh
rm(df_projected, df_inside_mesh)
# parse datetimes
df <- df %>%
  mutate(
    received  = parse_date_time(Received.DtTm,  "Y b d I:M:S p"),
    dispatch  = parse_date_time(Dispatch.DtTm,  "Y b d I:M:S p"),
    response  = parse_date_time(Response.DtTm,  "Y b d I:M:S p"),
    on_scene  = parse_date_time(On.Scene.DtTm,  "Y b d I:M:S p")
  )





# Round your timestamps to the nearest hour (to match API data structure)
df$received_hour <- round_date(df$received, unit = "hour")

# 3. FETCH WEATHER DATA
# ----------------------------------------
# We define the range of dates we need based on your data
lat_query <- mean(df$lat)
lon_query <- mean(df$lon)
start_date <- as.Date(min(df$received))
end_date   <- as.Date(max(df$received))

# Construct the API URL manually
# We request: Temperature, Rain, Wind Speed, Weather Code
base_url <- "https://archive-api.open-meteo.com/v1/archive"
query_params <- paste0(
  "?latitude=", lat_query,
  "&longitude=", lon_query,
  "&start_date=", start_date,
  "&end_date=", end_date,
  "&hourly=temperature_2m,precipitation,windspeed_10m,weathercode",
  "&timezone=UTC"
)

# Fetch and convert JSON to Dataframe
weather_raw <- fromJSON(paste0(base_url, query_params))
weather_df  <- as.data.frame(weather_raw$hourly)
weather_df$time <- as.POSIXct(weather_df$time, format="%Y-%m-%dT%H:%M", tz="UTC")
df <- df %>%
  left_join(weather_df, by = c("received_hour" = "time"))

names(df)
head(df)
# Definition of time-to-arrival and event indicator:
#
# The analysis is framed as a time-to-event (survival) problem, where the
# event of interest is a unit arriving on scene.
#
# Time origin (t = 0):
#   - Received.DtTm: time when the call was received by dispatch.
#
# Event:
#   - event = 1 if On.Scene.DtTm is observed (unit arrived on scene).
#   - event = 0 if On.Scene.DtTm is missing (unit did not arrive during
#     the observed time window, e.g. cancelled or diverted).
#
# Time variable (time_min):
#   - If event = 1: time_min is the exact arrival time
#       (On.Scene.DtTm - Received.DtTm).
#   - If event = 0: time_min is the right-censoring time, defined as the
#     last known timestamp for the unit (Response.DtTm, Dispatch.DtTm,
#     or Received.DtTm if no later information is available).
#
# Thus, time_min always represents the total time the unit was under
# observation, and the event indicator specifies whether arrival occurred
# within that time. This construction encodes right censoring explicitly
# and allows standard survival analysis methods to be applied.

df <- df %>%
  mutate(
    event = if_else(!is.na(on_scene), 1L, 0L),
    
    end_time = case_when(
      !is.na(on_scene) ~ on_scene,
      !is.na(response) ~ response,
      !is.na(dispatch) ~ dispatch,
      TRUE             ~ received
    ),
    
    time_min = as.numeric(difftime(end_time, received, units = "mins"))
  ) %>%
  filter(time_min >= 0)

MAX_MIN <- 60  # administrative censoring threshold (minutes)
df <- df %>%
  mutate(
    # administratively censored time
    time_min_admin = pmin(time_min, MAX_MIN),
    
    # event occurs only if arrival happened before the cutoff
    event_admin = if_else(event == 1 & time_min <= MAX_MIN, 1L, 0L)
  )




table(df$Call.Type)


# --- A. STRICT FILTERING ---
# 1. Structure Fires ONLY
# 2. First Arriving Unit ONLY (Sequence == 1)
df_first_response <- df %>%
  filter(
    Original.Priority %in% c(2, 3),
    Call.Type == "Medical Incident",#"Structure Fire / Smoke in Building",#
    Unit.sequence.in.call.dispatch == 1,
    time_min_admin > 0
    #ALS.Unit == "true"
  )

df_first_response$emergency = as.numeric(df_first_response$Original.Priority == 3)
hist(df_first_response$time_min_admin)
nrow(df_first_response)
sum(df_first_response$time_min_admin == 60)/nrow(df_first_response)
cat("Analysis Set: Structure Fires (First Response Only)\n")
cat("Total Rows:", nrow(df_first_response), "\n")


# convert (remove the trailing " UTC" then parse as UTC)
df_first_response$received_dt <- as.POSIXct(
  sub(" UTC$", "", df_first_response$received),
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)

# keep only 2025 rows
df_first_response <- df_first_response[format(df_first_response$received_dt, "%Y") == "2025", ]
df_first_response <- df_first_response[format(df_first_response$received_dt, "%m") == "02", ]
nrow(df_first_response)
table(df_first_response$event_admin)/nrow(df_first_response)
#df_first_response = df_first_response %>% filter(Unit.Type == "MEDIC")
# --- A. FEATURE ENGINEERING ---
df_analysis <- df_first_response %>%
  mutate(
    # 1. TIME: Decimal hour for Harmonics
    time_numeric = hour(received) + minute(received)/60 + second(received)/3600,
    
    # Ensure Priority is numeric if it isn't already (assuming 1,2,3 etc)
    # If Priority is "A", "B", "C", change this to as.factor()
    priority_numeric = as.numeric(as.character(Original.Priority)), 
    
    # 4. HOLIDAYS: US Holidays (SF is in US)
    # Convert call date to Date object
    call_date_obj = as.Date(received),
    is_holiday = ifelse(isHoliday(as.timeDate(call_date_obj), holidayNYSE()), 1, 0),
    
    # 5. WEATHER: Ensure numeric types
    temp = as.numeric(temperature_2m),
    precip = as.numeric(precipitation),
    wind = as.numeric(windspeed_10m)
    # Note: weathercode is usually categorical. If you want to use it, 
    # uncomment the line below, but it may create too many dummy columns.
    # weather_code = as.factor(weathercode)
  ) %>%
  # Filter missing critical data
  filter(!is.na(time_numeric), !is.na(time_min_admin))

# Check the distribution
table(df_analysis$is_holiday)
write.csv(df_analysis, "df_san_francisco.csv", row.names = FALSE)
# --- B. GENERATE HARMONICS (K=2) ---
# K1 (24h) + K2 (12h)
# We calculate this separately to keep the matrix clean
time_mat <- matrix(NA, nrow = nrow(df_analysis), ncol = 4)
colnames(time_mat) <- c("sin_24h", "cos_24h", "sin_12h", "cos_12h")

for(i in 1:2) {
  omega <- i * 2 * pi * df_analysis$time_numeric / 24
  time_mat[, (i*2)-1] <- sin(omega)
  time_mat[, (i*2)]   <- cos(omega)
}

# --- C. BUILD DESIGN MATRIX (X) ---
# We use model.matrix to handle factors (Unit Type) and dummies automatically.
# We include: Priority, Severity, Unit Type, Weather, Holiday
X_main <- model.matrix(
  ~ emergency + 
    temp + 
    precip + 
    wind + 
    is_holiday - 1,  # -1 removes the Intercept (Standard for Cox/glmnet)
  data = df_analysis
)
# --- D. MERGE AND CLEANUP ---
# Combine main covariates with Time Harmonics
Xcov <- cbind(X_main, time_mat)

# Print final column names to verify
# You should see: "unit_classOther", "unit_classTruck" (Engine is hidden as reference)
print(colnames(Xcov))

# Preview first few rows
head(Xcov)
dim(Xcov)
# Prepare Vectors
time_obs  <- df_analysis$time_min_admin
event_obs <- df_analysis$event_admin
coords = as.matrix(df_analysis[,c("lon", "lat")])[,c(1,2)]
head(coords)

table(event_obs)/nrow(Xcov)
hist(time_obs)


#n = nrow(Xcov)
#time_obs = time_obs[1:n]
#event_obs = event_obs[1:n]
#Xcov = Xcov[1:n,]
#coords = coords[1:n,]
#n



FEMbasis = fdaPDE::create.FEM.basis(mesh = final_mesh2)
# mass
M0 = fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)
# stiff
M1 = fdaPDE:::CPP_get.FEM.Stiff.Matrix(FEMbasis)
# penalty ?
K4 = M1%*% solve(M0) %*% M1


# Sum-to-zero constraint via mass matrix
m_vec <- rowSums(M0)
Z     <- pracma::null(t(m_vec))
Kred  <- t(Z) %*% K4 %*% Z


# 1. Coordinate System Alignment (Critical!)
# Re-create an sf object from your analysis data (Lat/Lon)
points_sf <- st_as_sf(df_analysis, coords = c("lon", "lat"), crs = 4326)

# Transform points to the EXACT same CRS as your mesh/polygons
# (Assuming 'simple_polys' was used to create 'final_mesh2')
points_projected <- st_transform(points_sf, st_crs(simple_polys))

# 2. Extract Coordinates as a Pure Matrix
# st_coordinates returns a matrix, which fixes the "REAL() can only be applied to a 'numeric'" error
coords_mat <- st_coordinates(points_projected)

# Double check dimensions and type
print(class(coords_mat)) # Should be "matrix" "array"
print(head(coords_mat))  # Should now look like mesh nodes (large numbers)
print(head(final_mesh2$nodes))


## 4. Design matrices (link data to basis) -----------------------
tmp = fdaPDE:::CPP_search_points(final_mesh2, coords_mat)
sum(tmp == -1)
A <- fdaPDE:::CPP_get.psi.Matrix(FEMbasis, coords_mat)
#A    <- eval_FEM_basis_all(FEMbasis, coords_mat)
dim(A)

sum(is.na(A))
Phi  <- A %*% Z                                  # n × (K-1)  constrained basis

df_analysis = as.data.frame(cbind(Xcov,time_obs,event_obs,coords))
names(df_analysis) = c(colnames(Xcov), "time", "status", "lon", "lat")
head(df_analysis)
df_analysis[] <- lapply(df_analysis, function(x) as.numeric(as.character(x)))
fit_cox = coxph(Surv(time,status)~
                  emergency + 
                  temp + 
                  precip + 
                  wind + 
                  is_holiday + 
                  sin_24h + 
                  cos_24h + 
                  sin_12h + 
                  cos_12h, data = df_analysis)
(bhat_cox <- stats::coef(fit_cox))

# mesh_2 is your mesh created via create.mesh.2D()
triangles <- final_mesh2$triangles
nodes <- final_mesh2$nodes

# Function to compute area of each triangle
triangle_area <- function(tri) {
  p1 <- nodes[tri[1], ]
  p2 <- nodes[tri[2], ]
  p3 <- nodes[tri[3], ]
  0.5 * abs(det(matrix(c(p2 - p1, p3 - p1), ncol = 2)))
}

# Apply over all triangles
areas <- apply(triangles, 1, triangle_area)
(AREA <- sum(areas))

gamma = .55
(lambda_grid <- exp(seq(log(0.05), log(50), length.out = 10)) * (nrow(df_analysis)^-gamma) * AREA ) 

dim(Xcov) 
# ## 7. Fit with optim() (BFGS) ------------------------------------


start_time <- Sys.time()
cv_res <- select_lambda_cvpld_opt(
  X = Xcov, 
  Phi = as.matrix(Phi), 
  Kred =as.matrix(Kred), 
  time = time_obs, 
  status = event_obs,
  lambda_grid = lambda_grid,
  K = 5,
  workers = parallel::detectCores(logical = FALSE),  # Adjust based on your cores
  criterion = "val"
)
end_time <- Sys.time()
print(end_time-start_time)
(choosen_index = which(lambda_grid == cv_res$lambda_opt))
choosen_index = 5
lambda = lambda_grid[choosen_index]
# ----- Final fit at selected lambda -----
start.time <- Sys.time()
print(start.time)
start <- rep(0, ncol(Xcov) + ncol(Phi))
fit <- fit_pen_cox_opt(
  X = Xcov, Phi = as.matrix(Phi), Kred = as.matrix(Kred),
  time = time_obs, status = event_obs,
  lambda = lambda, start = start, 
  hessian = TRUE
)
print(Sys.time() - start.time)
beta_hat  <- fit$beta
theta_hat    <- fit$theta

p <- ncol(Xcov)

H <- fit$hessian  


# Standard errors for beta
std_beta1 <- sqrt(diag(solve(H[1:p, 1:p])))
std_beta2 <- sqrt(diag(solve(H)[1:p, 1:p]))

# Name them
names(beta_hat) <- colnames(Xcov)
names(std_beta1) <- colnames(Xcov)
names(std_beta2) <- colnames(Xcov)

# View result
print(bhat_cox)
print(beta_hat)
print(std_beta1)
print(std_beta2)


############# PLOT



# --- 1. GENERATE PREDICTION GRID ---
# We sample points directly inside your 'simple_polys' (SF Boundary)
# simple_polys is already in UTM (EPSG:26910), so the grid will be too.
grid_sf <- st_sample(simple_polys, size = 50000, type = "regular") %>%
  st_as_sf()

# Extract Matrix coordinates for fdaPDE evaluation
grid_coords <- st_coordinates(grid_sf)

# --- 2. FAST FIELD EVALUATION ---
# Reconstruct the full coefficient vector on the mesh nodes
# theta_hat are the estimated reduced coefficients.
# Z is the constraint matrix used to map them back to the full mesh.
coeffs_full <- as.numeric(Z %*% theta_hat)

# Create a FEM functional object with these coefficients
# This allows us to evaluate the field without building the massive A matrix
fem_obj <- fdaPDE::FEM(coeffs_full, FEMbasis)

# Evaluate the spatial field at the grid points
grid_values <- fdaPDE::eval.FEM(fem_obj, grid_coords)

# Assign values back to the sf object
grid_sf$f_hat <- as.vector(grid_values)

# --- 3. PREPARE PLOT LIMITS ---
# Symmetrical limits are best for diverging (blue-white-red) scales
max_val <- max(abs(range(grid_sf$f_hat, na.rm = TRUE)))
limits_f <- c(-max_val, max_val)

# --- 4. PLOT ---
p_hat <- ggplot() +
  # A. The Spatial Field (Grid Points)
  # Using size=0.8 or smaller creates a smooth "raster-like" look
  geom_sf(data = grid_sf, aes(color = f_hat), size = 0.5, stroke = 0) +
  
  # B. The Mesh Edges (Optional: Low alpha to just show structure)
  # Convert mesh segments to SF for plotting if desired, or skip to keep it clean.
  # geom_sf(data = st_as_sf(simple_polys), fill = NA, color = "black", linewidth = 0.5) +
  
  # C. Color Scale (Diverging)
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = limits_f,
    name = "Log-Hazard"
  ) +
  
  # D. Map Decoration
  annotation_scale(location = "bl", width_hint = 0.2, style = "ticks") +
  annotation_north_arrow(
    location = "tr", which_north = "true", 
    pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  
  # E. Theme & Formatting
  theme_minimal() +
  labs(
    title = "Estimated Spatial Log-Hazard",
    subtitle = "San Francisco Structure Fires",
    x = "Longitude", 
    y = "Latitude"
  ) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    legend.position = "right"
  ) +
  
  # F. Coordinate System
  # Ensure geom_sf prints Lat/Lon labels even if data is UTM
  coord_sf(datum = st_crs(4326)) 

# Display Plot
print(p_hat)

ggsave(
  filename = "san_francisco_hat10.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)
################################################################################
######################### PLOTS
################################################################################

# Custom formatters for longitude and latitude
lon_formatter <- function(x) {
  paste0(abs(x), "°", ifelse(x < 0, "W", "E"))
}
lat_formatter <- function(x) {
  paste0(abs(x), "°", ifelse(x < 0, "S", "N"))
}



# 1. Extract Nodes and Segments
# (assuming final_mesh2 is an fdaPDE object)
nodes    <- final_mesh2$nodes     # The coordinates
segments <- final_mesh2$segments  # The indices of points connecting the boundary

# 2. Build the lines
# We loop through the segments and lookup their coordinates
lines_list <- lapply(1:nrow(segments), function(i) {
  # Get the two node indices for this segment
  node_indices <- segments[i, ]
  # Get the coordinates for these nodes
  coords <- nodes[node_indices, ]
  # Create a line string
  sf::st_linestring(coords)
})

# 3. Create the Simple Feature object with the correct CRS
boundary_sf <- sf::st_sfc(lines_list, crs = 32610) %>% 
  sf::st_sf(geometry = .)


plot(boundary_sf)
pts <- sf::st_sample(polys_clean, 100000, type = "regular") 
plot(final_mesh2, pch=".")
points(sf::st_coordinates(pts), col="black", pch='.')
grid <- as.data.frame(sf::st_coordinates(pts))
rm(pts)
names(grid) = c("x", "y")
Agrid = fdaPDE:::CPP_get.psi.Matrix(FEMbasis, as.matrix(grid[,c("x", "y")]))
#Agrid <- eval_FEM_basis_all(FEMbasis, grid[,c("x", "y")])
grid$f_hat = as.vector(Agrid %*% (Z %*% theta_hat))

p_hat <- ggplot(grid, aes(x = x, y = y, color = f_hat)) +
  geom_point(size = 1.2) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    name = "Value"
  ) +
  geom_sf(data = boundary_sf, inherit.aes = FALSE,
          color = "black", linewidth = 1) +
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
  filename = "san_francisco_hat.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

########## NEW PLOT HAT

###############################################
# Use the SAME projected CRS you plot in:
crs_map <- 32610   # <- your UTM

# 1) Get a representative center in a projection (avoid centroid-on-lonlat)
ctr_m <- boundary_sf |>
  st_union() |>
  st_transform(crs_map) |>
  st_centroid()

ctr_ll <- st_transform(ctr_m, 4326)
ctr_xy <- st_coordinates(ctr_ll)[1, ]  # lon, lat

# 2) Go 10 km along TRUE north (bearing 0) in lon/lat
p2_lonlat <- geosphere::destPoint(p = ctr_xy, b = 0, d = 10000)

p1_ll <- st_sfc(st_point(ctr_xy), crs = 4326)
p2_ll <- st_sfc(st_point(p2_lonlat[1, ]), crs = 4326)

# 3) Project both to your map CRS (UTM zone 10N)
p1_m <- st_transform(p1_ll, crs_map)
p2_m <- st_transform(p2_ll, crs_map)

# 4) Angle between true-north vector and +Y axis in projected coords
v <- st_coordinates(p2_m) - st_coordinates(p1_m)  # dx, dy
rot <- atan2(v[1], v[2]) * 180 / pi               # + clockwise from up (your convention)

rot

#####################################################
# 1. Create mesh_edges_m (if not already defined)
# Extract triangles, clean duplicates, and make SF object
tri_matrix <- final_mesh2$triangles
nodes      <- final_mesh2$nodes
all_edges  <- rbind(tri_matrix[, c(1, 2)], tri_matrix[, c(2, 3)], tri_matrix[, c(3, 1)])
all_edges  <- t(apply(all_edges, 1, sort)) # Sort to identify duplicates
unique_edges <- unique(all_edges)

lines_list <- lapply(1:nrow(unique_edges), function(i) {
  sf::st_linestring(nodes[unique_edges[i, ], ])
})
# Important: Mesh nodes are in 26910 (meters), so we define them as such
mesh_edges_m <- sf::st_sfc(lines_list, crs = 26910) %>% 
  sf::st_sf(geometry = .) %>%
  sf::st_transform(32610) # Transform to target CRS (WGS84 UTM 10N)

# 2. Correctly Create pts_sf
# grid coordinates are in 26910 (from polys_clean), NOT 4326
pts_sf <- st_as_sf(grid, coords = c("x", "y"), crs = 26910) 

# Now transform to your target plotting CRS (32610)
pts_m <- st_transform(pts_sf, crs = 32610)

# 3. Create Boundary (Projected)
# Ensure boundary matches the target CRS
boundary_m <- st_transform(boundary_sf, 32610) 

# 4. Plot
range_f <- round(range(pts_m$f_hat, na.rm = TRUE), 1)

p_hat <- ggplot() +
  # Plot mesh first (background)
  geom_sf(data = mesh_edges_m, color = "grey80", linewidth = 0.3) +
  # Plot spatial field
  geom_sf(data = pts_m, aes(color = f_hat), size = 1.5, shape = 15) + 
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = range_f,
    breaks = c(range_f[1], 0, range_f[2]),
    name = "Value"
  ) +
  # Plot boundary on top
  geom_sf(data = boundary_m, color = "black", linewidth = 1, fill = NA) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  ) +
  annotation_scale(location = "br", unit_category = "metric") +
  annotation_north_arrow(location = "tl", rotation = -rot)

p_hat

ggsave(
  filename = "san_francisco_hat2.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

##################
L = max(abs(range_f))+0.05
brks <- round(seq(-L, L, 0.1),2)
length(brks)

brks
#brks <- sort(unique(c(brks, 0)))  # ensure 0 shows up even if seq skipped it

p_hat <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, aes(color = f_hat), size = 2) +
  scale_color_steps2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    #limits = c(-L,L),
    breaks = brks,
    #show.limits = TRUE,
    name = "Value"
  ) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  ) +
  annotation_scale(location = "br", unit_category = "metric") +
  annotation_north_arrow(location = "tl", rotation = -rot)

print(p_hat)
# optional: same legend styling you used before
p_hat <- p_hat + guides(
  color = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(10, "cm"),
    barwidth  = unit(0.6, "cm"),
    label.theme = element_text(size = 9),
  )
)


p_hat

ggsave(
  filename = "san_francisco_hat4.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)




# --- 0) Choose a plotting CRS (match your workflow)
# In your script: simple_polys is 26910; you also use 32610 for plotting.
# Pick ONE and keep everything consistent:
target_crs <- 32610  # WGS84 / UTM zone 10N
# target_crs <- st_crs(simple_polys)  # alternatively: stay in 26910

# --- 1) Mesh edges -> sf lines (from triangles)
tri <- final_mesh2$triangles
nodes <- final_mesh2$nodes

all_edges <- rbind(tri[, c(1,2)], tri[, c(2,3)], tri[, c(3,1)])
all_edges <- t(apply(all_edges, 1, sort))      # canonical direction
edges <- unique(all_edges)

mesh_edges_sf <- st_sfc(
  lapply(seq_len(nrow(edges)), function(i) st_linestring(nodes[edges[i,], ])),
  crs = st_crs(simple_polys)  # mesh nodes are in the same CRS as simple_polys (26910 in your script)
) |>
  st_sf(geometry = _) |>
  st_transform(target_crs)

# --- 2) Boundary sf (clean + reliable): take boundary from simple_polys
boundary_sf <- st_boundary(simple_polys) |>
  st_as_sf() |>
  st_make_valid() |>
  st_transform(target_crs)

# --- 3) Observed points -> sf, then reproject to target_crs
# Option A: if df_analysis exists and has lon/lat columns
pts_sf <- st_as_sf(df_analysis, coords = c("lon", "lat"), crs = 4326, remove = FALSE) |>
  st_transform(target_crs)

# Option B: if you prefer df_first_response / df (same idea)
# pts_sf <- st_as_sf(df_first_response, coords = c("lon","lat"), crs = 4326, remove = FALSE) |>
#   st_transform(target_crs)

# --- 4) Plot
p_mesh <- ggplot() +
  geom_sf(data = mesh_edges_sf, color = "grey70", linewidth = 0.3) +
  geom_sf(data = pts_sf, color = "blue", size = 1.2, alpha = 0.7) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +
  labs(x = "Easting (m)", y = "Northing (m)") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  ) +
  annotation_scale(location = "br", unit_category = "metric") +
  annotation_north_arrow(location = "tl", rotation = -rot)

print(p_mesh)
pts_sf$time <- as.numeric(unlist(pts_sf$time, use.names = FALSE))
p_mesh <- ggplot() +
  geom_sf(data = mesh_edges_sf, color = "grey70", linewidth = 0.3) +
  geom_sf(
    data = pts_sf,
    aes(color = time),
    size = 1.2,
    alpha = 0.8
  ) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  ) +
  scale_color_viridis_c(option = "inferno", name = "Minutes") +
  annotation_scale(location = "br", unit_category = "metric") +
  annotation_north_arrow(location = "tl", rotation = -rot)

print(p_mesh)


ggsave(
  filename = "san_francisco_pmesh.png",
  plot     = p_mesh,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)


p_mesh2 <- ggplot() +
  geom_sf(data = mesh_edges_sf, color = "grey70", linewidth = 0.3) +
  geom_sf(
    data = pts_sf,
    aes(color = time),
    size = 1.2,
    alpha = 0.8
  ) +
  geom_sf(data = boundary_sf, fill = NA, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  ) +
  scale_color_viridis_c(
    option = "inferno",
    name = "Minutes",
    limits = c(0, 25),
    oob = scales::squish,
    breaks = c(0, 5, 10, 15, 20, 25),
    labels = c("0", "5", "10", "15", "20", "25+")
  )+
  annotation_scale(location = "br", unit_category = "metric") +
  annotation_north_arrow(location = "tl", rotation = -rot)

print(p_mesh2)


ggsave(
  filename = "san_francisco_pmesh2.png",
  plot     = p_mesh2,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)


library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)

# --- 1) Data Prep ---

# A) Locator Map Background
ca <- ne_states(country = "united states of america", returnclass = "sf") |>
  subset(name == "California") |>
  st_transform(4326)

# B) Zoom Map Background (High Res)
# Ensure rnaturalearthhires is installed for scale = 10
land_highres <- ne_download(scale = 10, type = "land", category = "physical", returnclass = "sf") |>
  st_make_valid() |>
  st_transform(4326)

# C) Process Mesh (REMOVE ISLANDS)
# Unite all simple polys
mesh_combined <- st_union(simple_polys) |> st_make_valid()

# Cast to separate polygons and keep only the largest one
mesh_poly_native <- st_cast(st_sfc(mesh_combined), "POLYGON") |> st_as_sf()
#mesh_parts$area <- st_area(mesh_parts)
#mesh_poly_native <- mesh_parts[which.max(mesh_parts$area), ]

# Transform result to WGS84
mesh_poly_wgs <- st_transform(mesh_poly_native, 4326)


# --- 2) Calculate Zoom Limits & Locator Box ---
bb <- st_bbox(mesh_poly_wgs)

# Padding (0.6 gives a nice margin around the shape)
pad_x <- 0.6 * (bb$xmax - bb$xmin)
pad_y <- 0.6 * (bb$ymax - bb$ymin)

xlim_zoom <- c(bb$xmin - pad_x, bb$xmax + pad_x)
ylim_zoom <- c(bb$ymin - pad_y, bb$ymax + pad_y)

# Create the black box for the locator map
zoom_box <- st_polygon(list(rbind(
  c(xlim_zoom[1], ylim_zoom[1]),
  c(xlim_zoom[2], ylim_zoom[1]),
  c(xlim_zoom[2], ylim_zoom[2]),
  c(xlim_zoom[1], ylim_zoom[2]),
  c(xlim_zoom[1], ylim_zoom[1])
))) |> st_sfc(crs = 4326)


# --- 3) Plotting ---

# Left Panel: Main Map
p_ca <- ggplot() +
  geom_sf(data = ca, fill = "grey95", color = "grey80", linewidth = 0.3) +
  geom_sf(data = mesh_poly_wgs, fill = "red", color = NA) +
  geom_sf(data = zoom_box, fill = NA, color = "black", linewidth = 0.8) +
  theme_void() +
  theme(plot.margin = margin(r = 10))

# Right Panel: Zoom Map
p_zoom <- ggplot() +
  # High res land background
  geom_sf(data = land_highres, fill = "grey95", color = "grey80", linewidth = 0.3) +
  # The clean mesh (no islands)
  geom_sf(data = mesh_poly_wgs, fill = "red", alpha = 0.9, color = "black", linewidth = 0.2) +
  coord_sf(xlim = xlim_zoom, ylim = ylim_zoom, expand = FALSE) +
  theme_void() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))

# Combine
p_final <- plot_grid(
  p_ca, p_zoom,
  ncol = 2,
  rel_widths = c(1.3, 1),
  labels = NULL
)

print(p_final)
ggsave(
  filename = "san_francisco_inset.png",
  plot     = p_final,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)
# Create a polygon object for the "Locator Box" on the main map
# This rectangle represents exactly what is shown in the zoom panel


colnames(Xcov)
ngrid = nrow(grid)


#emergency        temp      precip        wind  is_holiday     sin_24h     cos_24h     sin_12h     cos_12h
Xgrid = matrix(c(1,10,0,0,0,
              sin(1 * 2 * pi * 12 / 24),
              cos(1 * 2 * pi * 12 / 24),
              sin(1 * 2 * pi * 12 / 24),
              cos(1 * 2 * pi * 12 / 24)), nrow = 1)

grid$risk12 = exp(as.numeric(Xgrid%*% beta_hat) + grid$f_hat)
pts_m$risk12 = grid$risk12
global_min = min(grid$risk12)
global_max = max(grid$risk12)
# Define breaks (you can adjust how many steps you want)
brks_risk <- round(seq(global_min, global_max, length.out = 12), 1)

p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = risk12), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
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
    limits = c(global_min, global_max),
    name = "Risk"
  ) +
  annotation_scale(location = "br", unit_category = "metric") +
  annotation_north_arrow(location = "tl", rotation = -rot)

# Optional: match the colorbar style of your other plots
p_risk <- p_risk + guides(
  color = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(10, "cm"),
    barwidth = unit(0.6, "cm"),
    label.theme = element_text(size = 9)
  )
)

p_risk
ggsave(
  filename = "san_francisco_risk12_level.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)


