PATH = ""
{
  # Package names
  packages <- c(
    "doParallel", "foreach", "sp", "sf", "ggplot2", "survival", "fdaPDE", "scales",
    "ggspatial", "Matrix", "pracma", "data.table", "ggpubr", "spatsurv", "spatstat",
    "leaflet", "fields", "raster", "lubridate", "dplyr", "viridis", "rmapshaper","doRNG","geosphere"
  )
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages],repos='http://cran.us.r-project.org')
  }
  # Packages loadingseq(0.01,0.99,length.out=10)#
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
  
  negPL_pen(par = start,X = X,Phi = Phi,Kred = Kred,time = time, status = status,λ = lambda)
  
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
# helper for NULL-coalescing
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



data("fstimes")

data("fs")
fs <- fs[fs$Description=="Fire Station", ]
fscoords <- cbind(fs$Easting,fs$Northing)
chull <- convexhull.xy(rbind(coordinates(fstimes), fscoords))
win <- expandwinPerfect(chull, TRUE, 1.2)
fsppp <- ppp(x=fscoords[, 1], fscoords[, 2], window=win)
fsintens <- density.ppp(fsppp)
fsintens <- raster(fsintens)
proj4string(fsintens) <- CRS("+init=epsg:27700")
fstimes$fsintens <- raster::extract(fsintens, fstimes) * 100000000

for(i in 1:4) {
  fstimes@data[paste("s", i, sep="")] <- sin(i * 2 * pi * fstimes$timenumeric / 24)
  fstimes@data[paste("c", i, sep="")] <- cos(i * 2 * pi * fstimes$timenumeric / 24)
}

df = data.frame(time = fstimes$S[,1],
                status = fstimes$S[,2],
                Easting_m = fstimes$Easting_m,
                Northing_m = fstimes$Northing_m,
                fsintens = fstimes$fsintens,
                s1 = fstimes$s1, s2 = fstimes$s2, s3 = fstimes$s3, s4 = fstimes$s4,
                c1 = fstimes$c1, c2 = fstimes$c2, c3 = fstimes$c3, c4 = fstimes$c4)
rm(fstimes)
rm(fs)
# Convert df to sf object with BNG coordinates
df_sf <- st_as_sf(df, coords = c("Easting_m", "Northing_m"), crs = 27700)

# Transform to WGS84
df_ll <- st_transform(df_sf, 4326)

# Extract lon/lat as columns
coords <- st_coordinates(df_ll)

# Bind back into original df
df <- df %>%
  mutate(
    lon = coords[,1],
    lat = coords[,2]
  )

head(df)
df$time = df$time/60

#london <- st_read(paste(PATH,"London_Borough_Excluding_MHW.shp", sep = ""))
london <- st_read(paste(PATH,"London_Ward_CityMerged.shp", sep = ""))


latlon_crs <- st_crs(4326)      # WGS84 for latitude/longitude

# 3) Transform London boundary and data points to lat/lon
london_ll <- st_transform(london, latlon_crs)
plot(london_ll, axes = TRUE, graticule = TRUE, main = "London in WGS84 (lon/lat)")
plot(st_boundary(st_union(london_ll)))
london_bd_simp = rmapshaper::ms_simplify(st_boundary(st_union(london_ll)), keep = 0.01) # 0.05
plot(london_bd_simp)
dim(st_coordinates(london_bd_simp))

bd_coords = st_coordinates(london_bd_simp)[,-3]
# avoid repetead nodes
bd_coords = bd_coords[-1,]

# check --- 
plot(bd_coords)
# points(bd_coords[1,1], bd_coords[1,2], col="red", pch=16)
# points(bd_coords[2,1], bd_coords[2,2], col="blue", pch=16)
# points(bd_coords[nrow(bd_coords),1], bd_coords[nrow(bd_coords),2], pch=16)
# text(bd_coords[,1], bd_coords[,2], labels = 1:nrow(bd_coords), pos = 3, cex = 0.7)
# 
# points(bd_coords[seq(2,66,2),1], bd_coords[seq(2,66,2),2], col = "green", pch = 16)
# bd_coords = bd_coords[-c(seq(2,66,2), seq(180,430,2)),]
# 
# plot(bd_coords)
# points(bd_coords[1,1], bd_coords[1,2], col="red", pch=16)
# points(bd_coords[2,1], bd_coords[2,2], col="blue", pch=16)
# points(bd_coords[nrow(bd_coords),1], bd_coords[nrow(bd_coords),2], pch=16)
# text(bd_coords[,1], bd_coords[,2], labels = 1:nrow(bd_coords), pos = 3, cex = 0.7)

segments = cbind(1:(nrow(bd_coords)-1), 2:nrow(bd_coords))
segments =  rbind(segments, c(nrow(bd_coords), 1))
mesh_2 = fdaPDE::create.mesh.2D(nodes = bd_coords, segments = segments)
plot(mesh_2)

# mesh_2 is your mesh created via create.mesh.2D()
triangles <- mesh_2$triangles
coords <- mesh_2$nodes

# Function to compute area of each triangle
triangle_area <- function(tri) {
  p1 <- coords[tri[1], ]
  p2 <- coords[tri[2], ]
  p3 <- coords[tri[3], ]
  0.5 * abs(det(matrix(c(p2 - p1, p3 - p1), ncol = 2)))
}

# Apply over all triangles
areas <- apply(triangles, 1, triangle_area)
(AREA <- sum(areas))

mesh_2 = fdaPDE::refine.mesh.2D(mesh_2, minimum_angle = 30, 
                                maximum_area =  0.1*AREA/(nrow(df)^alpha))
dim(mesh_2$nodes)
plot(mesh_2)

nrow(df)
# mesh_2_inla = inla.mesh.create(
#   loc = mesh_2$nodes, 
#   tv = mesh_2$triangles,
#   crs = 4326
# )
# plot(mesh_2_inla)



###############################################
# 1) Pick a representative point (map center) in lon/lat
ctr_ll <- boundary_m |>
  sf::st_union() |>
  sf::st_centroid() |>
  sf::st_transform(4326)

ctr_xy <- sf::st_coordinates(ctr_ll)[1, ]  # c(lon, lat)

# 2) Go 10 km along *true* north from the center (bearing 0°) in lon/lat
p2_lonlat <- geosphere::destPoint(p = ctr_xy, b = 0, d =50000)  # 10 km

p1_ll <- sf::st_sfc(sf::st_point(ctr_xy), crs = 4326)
p2_ll <- sf::st_sfc(sf::st_point(p2_lonlat[1, ]), crs = 4326)

# 3) Project both to your UTM (grid) coordinates
p1_utm <- sf::st_transform(p1_ll, 32630)
p2_utm <- sf::st_transform(p2_ll, 32630)

# 4) Vector and rotation: angle between true-north vector and +Y axis (grid north)
v   <- sf::st_coordinates(p2_utm) - sf::st_coordinates(p1_utm)  # c(dx, dy)
rot <- atan2(v[1], v[2]) * 180 / pi  # degrees; positive = clockwise from up

#####################################################


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
poly <- boundary_sf |>
  st_union() |>          # dissolve segments into one geometry
  st_line_merge() |>     # connect into longer rings
  st_cast("POLYGON") |>  # cast to polygon(s)
  st_make_valid()

# Convert your df to sf points
pts_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = st_crs(boundary_sf))

# Spatial join or logical check
inside <- st_within(pts_sf, poly, sparse = FALSE)
fit_std = coxph(Surv(time,status)~fsintens+s1+c1+s2+c2+s3+c3+s4+c4, data = df)
bhat_std <- stats::coef(fit_std)
df = df[inside[,1],]

nrow(df)
# --- plot ---
p_mesh <- ggplot() +
  geom_sf(data = mesh_edges_sf, color = "grey60", linewidth = 0.5) +
  geom_point(data = df,
             aes(x = lon, y = lat),
             color = "blue",
             size = 1) +
  geom_sf(data = boundary_sf, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +   # ensures map aspect ratio is respected
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    legend.position = "none"
  )+
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr", rotation = -rot)

p_mesh

ggsave(
  filename = "london_mesh.png",
  plot     = p_mesh,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)


################## NEW PLOT MESH
ctr <- st_coordinates(st_centroid(st_union(boundary_sf)))
utm_zone <- floor((ctr[1] + 180) / 6) + 1
epsg <- if (ctr[2] >= 0) 32600 + utm_zone else 32700 + utm_zone  # 326xx for N, 327xx for S
target_crs <- st_crs(epsg)

st_crs(mesh_edges_sf) <- 4326
st_crs(boundary_sf)   <- 4326

# also convert your points to sf in WGS84
pts_sf <- st_as_sf(df, coords = c("lon","lat"), crs = 4326)

# 3) Reproject everything
mesh_edges_m <- st_transform(mesh_edges_sf, target_crs)
boundary_m   <- st_transform(boundary_sf,   target_crs)
pts_m        <- st_transform(pts_sf,        target_crs)

# 4) Plot in meters (axis units are now meters; aspect is true)
p_mesh <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  #geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, size = 1.2, color = "blue") +  # border
  #geom_sf(data = pts_m, aes(color = "blue"), size = 2) +
  #scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "Censored") +
  geom_sf(data = boundary_sf, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr", rotation = -rot)

p_mesh




ggsave(
  filename = "london_mesh2.png",
  plot     = p_mesh,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

################## plot time

p_time <- ggplot() +
  geom_sf(data = mesh_edges_sf, color = "grey60", linewidth = 0.5) +
  geom_point(data = df,
             aes(x = lon, y = lat, color = time),
             size = 1) +
  geom_sf(data = boundary_sf, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +   # ensures map aspect ratio is respected
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_viridis_c(option = "plasma", name = "Minutes")

p_time
ggsave(
  filename = "london_time.png",
  plot     = p_time,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)


#################### NEW PLOT TIME

p_time <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = time), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_viridis_c(option = "inferno", name = "Minutes") +
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr",rotation = -rot)


p_time

ggsave(
  filename = "london_time2.png",
  plot     = p_time,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)



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

coords = as.matrix(df[,c("lon", "lat")])
## 4. Design matrices (link data to basis) -----------------------
tmp = fdaPDE:::CPP_search_points(mesh_2, coords)
sum(tmp == -1)
A    <- eval_FEM_basis_all(FEMbasis, coords)
#A    <- INLA::inla.spde.make.A(mesh_2_inla, loc = as.matrix(coords))
sum(is.na(A))
Phi  <- A %*% Z                                  # n × (K-1)  constrained basis


Xcov <- as.matrix(df[, c("fsintens","s1", "c1", "s2", "c2", "s3", "c3", "s4", "c4")])         # add your covariates here
fit_std = coxph(Surv(time,status)~fsintens+s1+c1+s2+c2+s3+c3+s4+c4, data = df)
bhat_std <- stats::coef(fit_std)

(lambda_grid <- exp(seq(log(0.05), log(50), length.out = 10)) * (nrow(df)^-gamma) * AREA ) 

# ## 7. Fit with optim() (BFGS) ------------------------------------

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
# 
# print(paste("lambda chosen = ", lambda))
lambda = lambda_grid[10]
# ----- Final fit at selected lambda -----
start.time <- Sys.time()
print(start.time)
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
std_beta1 <- sqrt(diag(solve(H[1:p, 1:p])))
std_beta2 <- sqrt(diag(solve(H)[1:p, 1:p]))

# Name them
names(beta_hat) <- colnames(Xcov)
names(std_beta1) <- colnames(Xcov)
names(std_beta2) <- colnames(Xcov)

# View result
print(beta_hat)
print(std_beta1)
print(std_beta2)

# Also write output to a text file
sink("data_empirica10_results.txt")  # Start writing output to a file
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

# Predict spatial field on full grid


#save.image(file='empirical.RData')
## 8. Post-processing: predict f(s) on a grid --------------------
#load("C:/Users/utente/Documents/empirical.RData")
# 1. Leggi il shapefile

poly <- boundary_sf |>
  st_union() |>
  st_line_merge() |>
  st_polygonize() |>
  st_make_valid()

pts <- sf::st_sample(poly, 100000, type = "regular") 
plot(mesh_2, pch=".")
points(sf::st_coordinates(pts), col="black", pch='.')
grid <- as.data.frame(sf::st_coordinates(pts))
rm(pts)
names(grid) = c("x", "y")
Agrid <- eval_FEM_basis_all(FEMbasis, grid[,c("x", "y")])
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
  filename = "london_hat.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

########## NEW PLOT HAT

# also convert your points to sf in WGS84
pts_sf <- st_as_sf(grid, coords = c("x","y"), crs = 4326)

pts_m        <- st_transform(pts_sf,        target_crs)
range_f <- round(range(pts_m$f_hat, na.rm = TRUE),1)

p_hat <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  #geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = f_hat), size = 2) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = range_f,               # ensure full range is shown
    breaks = c(range_f[1], 0, range_f[2]),  # show min, mid, max on legend
    name = "Value"
  )  +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr",rotation = -rot)


p_hat

ggsave(
  filename = "london_hat2.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

##################
brks <- round(seq(-0.3, 0.3, 0.04),2)
#brks <- sort(unique(c(brks, 0)))  # ensure 0 shows up even if seq skipped it

p_hat <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, aes(color = f_hat), size = 2) +
  scale_color_steps2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0,
    limits = c(max(abs(range_f)), -max(max(abs(range_f)))),
    breaks = brks,
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
  annotation_scale(location = "bl", unit_category = "metric") +
  annotation_north_arrow(location = "tr", rotation = -rot)

# optional: same legend styling you used before
p_hat <- p_hat + guides(
  color = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = unit(10, "cm"),
    barwidth  = unit(0.6, "cm"),
    label.theme = element_text(size = 9)
  )
)

p_hat
ggsave(
  filename = "london_hat3.png",
  plot     = p_hat,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

############### plot RISK by time


pts <- SpatialPoints(as.matrix(grid[,c("x", "y")]), proj4string = CRS("+init=epsg:4326"))
pts_27700 <- spTransform(pts, CRS("+init=epsg:27700"))


grid$fsintens <- raster::extract(fsintens, pts_27700) * 100000000

for(TIME in c(0,6,12,18)){
  for(i in 1:4) {
    grid[,paste("s", i, sep="")] <- sin(i * 2 * pi * TIME / 24)
    grid[,paste("c", i, sep="")] <- cos(i * 2 * pi * TIME / 24)
  }
  Xcov = as.matrix(grid[,colnames(Xcov)])
  grid[,paste("risk", TIME, sep="")] <-  exp( as.numeric(Xcov %*% beta_hat) + grid$f_hat )
}
head(grid)
pts_sf <- st_as_sf(grid, coords = c("x","y"), crs = 4326)
pts_m        <- st_transform(pts_sf,        target_crs)

risk_cols <- grep("^risk", colnames(grid), value = TRUE)
global_min <- min(grid[, risk_cols], na.rm = TRUE)
global_max <- max(grid[, risk_cols], na.rm = TRUE)


p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = risk0), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_viridis_c(
    option = "inferno",
    name   = "Risk",
    limits = c(global_min, global_max)   # 🔑 ensures consistent scale
  ) +
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr",rotation = -rot)

p_risk


ggsave(
  filename = "london_risk0.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)

p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = risk6), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_viridis_c(
    option = "inferno",
    name   = "Risk",
    limits = c(global_min, global_max)   # 🔑 ensures consistent scale
  ) +
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr",rotation = -rot)

p_risk


ggsave(
  filename = "london_risk6.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)
p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = risk12), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_viridis_c(
    option = "inferno",
    name   = "Risk",
    limits = c(global_min, global_max)   # 🔑 ensures consistent scale
  ) +
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr",rotation = -rot)

p_risk


ggsave(
  filename = "london_risk12.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)
p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = risk18), size = 2) +
  geom_sf(data = boundary_m, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +  # keeps 1:1 in projected units
  labs(x = "Longitude", y = "Latitude") +
  scale_x_continuous(labels = lon_formatter) +
  scale_y_continuous(labels = lat_formatter) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10)
  )+
  scale_color_viridis_c(
    option = "inferno",
    name   = "Risk",
    limits = c(global_min, global_max)   # 🔑 ensures consistent scale
  ) +
  annotation_scale(location = "bl", unit_category = "metric") +  # optional: scale bar
  annotation_north_arrow(location = "tr",rotation = -rot)

p_risk


ggsave(
  filename = "london_risk18.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)




# Define breaks (you can adjust how many steps you want)
brks_risk <- round(seq(global_min, global_max, length.out = 20), 2)

p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = risk0), size = 2) +
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
  annotation_scale(location = "bl", unit_category = "metric") +
  annotation_north_arrow(location = "tr")

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
  filename = "london_risk0_level.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)
########### risk 6

p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = risk6), size = 2) +
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
  annotation_scale(location = "bl", unit_category = "metric") +
  annotation_north_arrow(location = "tr")

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
  filename = "london_risk6_level.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)


########### risk 12

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
  annotation_scale(location = "bl", unit_category = "metric") +
  annotation_north_arrow(location = "tr")

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

#p_risk
ggsave(
  filename = "london_risk12_level.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)


########### risk 18

p_risk <- ggplot() +
  geom_sf(data = mesh_edges_m, color = "grey60", linewidth = 0.5) +
  geom_sf(data = pts_m, size = 3, color = "black") +  # border
  geom_sf(data = pts_m, aes(color = risk18), size = 2) +
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
  annotation_scale(location = "bl", unit_category = "metric") +
  annotation_north_arrow(location = "tr")

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
  filename = "london_risk18_level.png",
  plot     = p_risk,
  width    = 10,
  height   = 8,
  dpi      = 300,
  units    = "in"
)
