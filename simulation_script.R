# ================================
# Full Simulation Script (Long/Tidy) with CV-parallel lambda search
# ================================
# 
# Parse command-line arguments
#  This script was run using the following package version 
#  parallel     4.5.1
#  doParallel   1.0.17
#  doRNG        1.8.6.2
#  foreach      1.5.2
#  mgcv         1.9.3
#  sp           2.2.0
#  sf           1.0.21
#  ggplot2      4.0.0
#  survival     3.8.3
#  Matrix       1.7.3
#  pracma       2.4.4
#  fdaPDE       1.1.21
#  tibble       3.3.0
#  purrr        1.1.0
#  dplyr        1.1.4
#  data.table   1.17.8
#  tidyr        1.3.1
#  readr        2.1.5
#  stringr      1.5.2
#  mgcv         1.9.3



# ------- Packages -------
packages <- c(
  "parallel","doParallel", "doRNG" ,"foreach","mgcv","sp","sf",
  "ggplot2","survival","Matrix","pracma", "fdaPDE",
  "tibble","purrr","dplyr","data.table","tidyr","readr","stringr","mgcv"
)




installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages], repos = "http://cran.us.r-project.org")
}
invisible(lapply(packages, library, character.only = TRUE))

# ------- Helpers -------
now_string <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
set.seed(1234)
verbose = FALSE
round_down <- function(x, digits = 1) {
  factor <- 10^digits
  floor(x * factor) / factor
}

# Define the indicator functions 1A, 1B, and 1C
f = function(x,y){fs.test(x,y,r0 = 0.03)-0.05755041}
data(horseshoe2D)
horseshoe_sp_2 <- sp::SpatialPolygons(
  list(
    sp::Polygons(
      list(
        sp::Polygon(horseshoe2D$boundary_nodes)  # Combine x and y into polygon coordinates
      ), 
      ID = "horseshoe2D"  # Assign an ID to the polygon
    )
  )
)


horseshoe_sf <- sf::st_as_sf(horseshoe_sp_2)
horseshoe_boundary <- sf::st_boundary(sf::st_union(sf::st_make_valid(horseshoe_sf)))

mesh_2 <- fdaPDE::create.mesh.2D(
  horseshoe2D$boundary_nodes,
  segments = horseshoe2D$boundary_segments,
  order = 1
)

# Refine mesh based on sample size
mesh_2 <- fdaPDE::refine.mesh.2D(
  mesh_2,
  minimum_angle = 30,
  maximum_area = 0.001
)

plot(mesh_2)

# Compute total area
int_triangles <- mesh_2$triangles
int_coords_mesh <- mesh_2$nodes
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

integrate_over_omega_from_barycenters(f(int_barycenters[,1], int_barycenters[,2]))

(AREA =  integrate_over_omega_from_barycenters(rep(1,nrow(int_barycenters))))

grid <- as.data.frame(int_barycenters)

# Rename coordinate columns
colnames(grid) = c("x", "y")
grid$id = seq_len(nrow(grid))
plot(grid$x, grid$y)



f_true <- f(grid$x,grid$y)
grid$f_true = f_true

# ------- Data generation for one simulation -------
make_spatial_cox_df <- function(n = 500, lambda0 = 0.05,
                                beta1 = 0.03, beta2 = -0.5,
                                censor_rate = 0.02, seed = 1234) {
  set.seed(seed)
  sample_indices <- sample(nrow(grid), n, replace = FALSE)
  x <- grid$x[sample_indices]
  y <- grid$y[sample_indices]
  f_true <- f(x, y)
  X1 <- rnorm(n, mean = 50, sd = 10)
  X2 <- rbinom(n, 1, 0.5)
  eta <- beta1 * X1 + beta2 * X2 + f_true
  event_time  <- rexp(n, rate = lambda0 * exp(eta))
  censor_time <- rexp(n, rate = censor_rate)
  time   <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)
  tibble::tibble(time, status, x, y, X1, X2, f_true,sample_indices)
}

# ------- Cox partial likelihood helpers -------
cumexp <- function(eta) rev(cumsum(rev(exp(eta))))  # ∑_{j ∈ R(t_i)} e^{η_j}

negPL_pen <- function(par, X, Phi, Kred, time, status, λ) {
  p      <- ncol(X)
  beta   <- par[seq_len(p)]
  theta  <- par[-seq_len(p)]
  eta    <- as.vector(X %*% beta + Phi %*% theta)
  
  ord    <- order(time)
  eta    <- eta[ord]; 
  status <- status[ord]
  
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
# (parallelize over lambda values; folds stay serial for warm starts)
logPL_norm <- function(beta, theta, X, Phi, time, status) {
  n   <- nrow(X)
  eta <- as.vector(X %*% beta + Phi %*% theta)
  
  ord <- order(time)
  eta <- eta[ord]; status <- status[ord]
  
  logrisk <- log(cumexp(eta))               # log ∑_{j∈R(t_i)} e^{η_j}
  ll_vec  <- status * (eta - (logrisk - log(n)))
  (1/n) * sum(ll_vec)
}

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
    control = modifyList(list(maxit = 2000, reltol = 1e-7), control),
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
    lambda_grid <- seq(0.1 * lam0, 10 * lam0, length.out = 10)
  }
  
  # random folds
  folds <- sample(rep(1:K, length.out = n))
  
  # results container
  cv_results <- vector("list", length(lambda_grid))
  
  # loop over lambda values
  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    score_g <- 0
    for (k in 1:K) {
      idx_tr <- which(folds != k)
      warm <- NULL
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
    
    cv_results[[i]] <- list(lambda = lam, cvpld = score_g)
  }
  
  cvpld   <- vapply(cv_results, `[[`, numeric(1), "cvpld")
  lambdas <- vapply(cv_results, `[[`, numeric(1), "lambda")
  
  list(
    lambda_opt = lambdas[which.min(cvpld)],
    grid       = lambdas,
    cvpld      = cvpld
  )
}


# ================================
# Fast CV-PLD with parallel folds
# ================================
# Needs: fit_pen_cox(), logPL_norm() already defined in the environment
# Parallel backend: doParallel + foreach (+ doRNG for reproducibility)

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




# Evaluate ALL basis functions of an fdaPDE FEMbasis at given locations
# - FEMbasis: object from create.FEM.basis()
# - locations: N x d matrix of coordinates (d = 2 for 2D meshes, 3 for 3D)
# Returns: N x nbasis numeric matrix, with NA for points outside the mesh.
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




# ------- Simulation settings -------
alpha  <- 0.45  # in (1/8, 1/2)
gamma  <- .55  # > 1/2
lambda0       <- 0.0001
beta_true     <- c(.25, -1)   # (beta1, beta2)
censor_rates  <- c(0.27,2) # ~15%, ~30% expected (realised may vary)
sample_sizes  <- c(250,500,1000, 2000)
num_simulations <- 250

inputs <- as.data.frame(as.matrix(expand.grid(
  censor_rate   = censor_rates,
  sample_size   = sample_sizes,
  id_simulation = seq_len(num_simulations)
)))


# ------- Serial sims: Proposed Spatial PH + Standard PH -------
results_list <- vector("list", length = nrow(inputs))
pb = txtProgressBar(min = 0, max = nrow(inputs), initial = 0, style = 3) 
# Precompute mesh and FEM structures by sample size
unique_sizes <- unique(inputs$sample_size)
mesh_cache <- list()


for (ss in unique_sizes) {
  message("Precomputing mesh and FEM structures for sample size = ", ss)
  
  # Base mesh
  mesh_2 <- fdaPDE::create.mesh.2D(
    horseshoe2D$boundary_nodes,
    segments = horseshoe2D$boundary_segments,
    order = 1
  )
  
  # Compute total area
  triangles <- mesh_2$triangles
  coords_mesh <- mesh_2$nodes
  triangle_area <- function(tri) {
    p1 <- coords_mesh[tri[1], ]
    p2 <- coords_mesh[tri[2], ]
    p3 <- coords_mesh[tri[3], ]
    0.5 * abs(det(matrix(c(p2 - p1, p3 - p1), ncol = 2)))
  }
  areas <- apply(triangles, 1, triangle_area)
  total_area <- sum(areas)
  
  # Refine mesh based on sample size
  mesh_2 <- fdaPDE::refine.mesh.2D(
    mesh_2,
    minimum_angle = 30,
    maximum_area = 0.1*total_area / (ss^alpha)
  )
  plot(mesh_2)
  # Create FEM basis
  FEMbasis <- fdaPDE::create.FEM.basis(mesh = mesh_2)
  nrow(FEMbasis$mesh$nodes)
  # Mass and stiffness matrices
  M0 <- fdaPDE:::CPP_get.FEM.Mass.Matrix(FEMbasis)
  M1 <- fdaPDE:::CPP_get.FEM.Stiff.Matrix(FEMbasis)
  
  # Penalty matrix
  K4 <- M1 %*% solve(M0) %*% M1
  
  # Sum-to-zero constraint
  m_vec <- rowSums(M0)
  Z <- pracma::null(t(m_vec))
  Kred <- t(Z) %*% K4 %*% Z
  Agrid <- eval_FEM_basis_all(FEMbasis, grid[, c("x", "y")])
  #Agrid = NULL
  # Cache everything
  mesh_cache[[as.character(ss)]] <- list(
    mesh = mesh_2,
    FEMbasis = FEMbasis,
    M0 = M0,
    K4 = K4,
    Z = Z,
    Kred = Kred,
    Agrid = Agrid
  )
}
sapply(unique_sizes, function(ss) 
  nrow(mesh_cache[[as.character(ss)]]$FEMbasis$mesh$nodes))
  # ------------------------------------------------------------------------------
# Main simulation loop
# ------------------------------------------------------------------------------

for (i in seq_len(nrow(inputs))) {
  setTxtProgressBar(pb, i)
  
  censor_rate   <- inputs$censor_rate[i]
  sample_size   <- inputs$sample_size[i]
  id_simulation <- inputs$id_simulation[i]
  
  # Retrieve cached mesh & matrices for this sample size
  mc <- mesh_cache[[as.character(sample_size)]]
  FEMbasis <- mc$FEMbasis
  Z        <- mc$Z
  Kred     <- mc$Kred
  Agrid    <- mc$Agrid
  
  # --- Simulate data ---
  df <- make_spatial_cox_df(
    n = sample_size,
    censor_rate = censor_rate,
    seed = id_simulation,
    beta1 = beta_true[1],
    beta2 = beta_true[2],
    lambda0 = lambda0
  )
  #plot(mc$mesh)
  #points(df$x,df$y, col = ifelse(df$status == 1, "blue", "red"))
  censor_lvl <- 1 - mean(df$status)
  coords <- as.matrix(df[, c("x", "y")])
  
  # --- GAM model ---
  fitgam <- gam(
    Surv(time, status) ~ X1 + X2 + s(x, y, bs = "tp", k = 50),
    family = cox.ph(), data = df, method = "REML"
  )
  
  # --- Proposed spatial PH ---
  A <- Agrid[df$sample_indices,]
  Phi <- A %*% Z
  Xcov <- as.matrix(df[, c("X1", "X2")])
  
  # Lambda grid
  lambda_grid <- exp(seq(log(0.05), log(50), length.out = 10)) * (sample_size^-gamma) * AREA#exp(seq(log(10^-5), log(10^5), length.out = 10))
  # start_time = Sys.time()
  # cvsel <- select_lambda_cvpld(
  #   X = Xcov,
  #   Phi = Phi,
  #   Kred = Kred,
  #   time = df$time,
  #   status = df$status,
  #   K = 10,                        # adjust if desired
  #   lambda_grid = lambda_grid,           # or provide numeric vector
  #   seed = 1000 + id_simulation   # CV fold seed
  # )
 
  start_time = Sys.time()
  cvsel <- select_lambda_cvpld_fast(
    X = Xcov, Phi = Phi, Kred = as.matrix(Kred),
    time = df$time, status = df$status,
    K = 5,
    lambda_grid = lambda_grid,
    seed = 1000 + id_simulation,
    workers = min(10,parallel::detectCores() - 1),   # parallel
    criterion = "diff"
  )

  lambda <- cvsel$lambda_opt
  
  # Fit penalized Cox
  start <- rep(0, ncol(Xcov) + ncol(Phi))
  ctime <- Sys.time()
  fit <- fit_pen_cox(
    X = Xcov, Phi = Phi, Kred = Kred,
    time = df$time, status = df$status,
    lambda = lambda, start = start, hessian = TRUE
  )
  Sys.time() - ctime
  
  beta_hat_sp <- fit$beta
  theta_hat   <- fit$theta
  H <- fit$hessian
  
  
  vcov_full <- solve(H)
  vcov_beta <- vcov_full[1:2,1:2]
  std_beta <- sqrt(diag(vcov_beta))
  
  
  # --- Predict spatial field ---
  f_hat <- as.vector(Agrid %*% (Z %*% theta_hat))
  grid$f_hat <- f_hat
  grid$X1 <- 0
  grid$X2 <- 0
  
  pt <- predict(fitgam, newdata = grid, type = "terms", se.fit = TRUE, terms = "s(x,y)")
  f_gam <- drop(pt$fit)
  f_gam <- f_gam - integrate_over_omega_from_barycenters(f_gam) / AREA
  
  grid$f_gam <- f_gam
  grid$diff <- grid$f_hat - grid$f_true
  integrate_over_omega_from_barycenters(grid$diff^2)
  
  # p_diff <- ggplot(grid, aes(x = x, y = y, color = f_hat)) +
  #   geom_point(size = 1.2) +
  #   scale_color_gradient2(
  #     low = "blue", mid = "white", high = "red",
  #     midpoint = 0,
  #     #limits = lims,
  #     name = "Value"
  #   ) +
  #   geom_sf(data = horseshoe_boundary, inherit.aes = FALSE,
  #           color = "black", linewidth = 1) +
  #   coord_sf() +
  #   theme_minimal()
  # p_diff
  # 
  # --- Standard PH model ---
  fit_std <- survival::coxph(Surv(time, status) ~ X1 + X2, data = df)
  bhat_std <- coef(fit_std)
  
  # --- Summaries ---
  summary_df_i <- tibble::tibble(
    id_simulation = id_simulation,
    censor_rate   = censor_rate,
    sample_size   = sample_size,
    censor_level  = censor_lvl,
    lambda        = lambda,
    model         = c("spatial", "standard", "gam"),
    beta1_hat     = c(beta_hat_sp[1], unname(bhat_std["X1"]), fitgam$coefficients[1]),
    beta2_hat     = c(beta_hat_sp[2], unname(bhat_std["X2"]), fitgam$coefficients[2]),
    std_beta1     = std_beta[1],
    std_beta2     = std_beta[2]
  )
  
  fhat_df_i <- tibble::tibble(
    id_simulation = id_simulation,
    sample_size   = sample_size,
    censor_rate   = censor_rate,
    grid_id       = seq_len(nrow(grid)),
    f_hat         = f_hat,
    f_gam         = f_gam
  )
  
  results_list[[i]] <- list(summary = summary_df_i, fhat = fhat_df_i)
}

# =========================
# Bind + Sanity-safe joins
# =========================

length(results_list)

# # Split the combined list into the two long tibbles we need
# is_summary <- purrr::map_lgl(results_list, ~ is.data.frame(.x) && any(c("beta1_hat","beta2_hat") %in% names(.x)))
# summaries  <- dplyr::bind_rows(results_list[is_summary])
# fhats_long <- dplyr::bind_rows(results_list[!is_summary])
# bind all summary tibbles
summaries <- results_list %>%
  map("summary") %>%
  list_rbind()
# bind all fhat tibbles
fhats_long <- results_list %>%
  map("fhat") %>%
  list_rbind()
# Helper: nice label for target censoring level used in the sim setup
censor_label <- function(rate) ifelse(rate == censor_rates[1], "15%", "30%")

# ---- Coefficient metrics (bias & MSE) ----
beta_true <- beta_true  # keep for NSE clarity

coef_stats <- summaries %>%
  mutate(
    bias_beta1 = beta1_hat - beta_true[1],
    bias_beta2 = beta2_hat - beta_true[2],
    mse_beta   = (beta1_hat - beta_true[1])^2 + (beta2_hat - beta_true[2])^2
  ) %>%
  group_by(sample_size, censor_rate, model) %>%
  summarise(
    `Bias beta1` = mean(bias_beta1),
    `Bias beta2` = mean(bias_beta2),
    `MSE beta`   = mean(mse_beta),
    .groups = "drop"
  ) %>%
  mutate(
    SampleSize = sample_size,
    CensorRate = censor_label(censor_rate)
  ) %>%
  select(SampleSize, CensorRate, model,  `Bias beta1`, `Bias beta2`, `MSE beta`)


# --- 95% coverage for the proposed (spatial) method ---
z <- qnorm(0.975)  # 1.959964...


cov_spatial <- summaries %>%
  dplyr::filter(model == "spatial") %>%
  dplyr::mutate(
    lo1 = beta1_hat - z * std_beta1,
    hi1 = beta1_hat + z * std_beta1,
    lo2 = beta2_hat - z * std_beta2,
    hi2 = beta2_hat + z * std_beta2,
    cover1 = (beta_true[1] >= lo1 & beta_true[1] <= hi1),
    cover2 = (beta_true[2] >= lo2 & beta_true[2] <= hi2)
  ) %>%
  dplyr::group_by(sample_size, censor_rate) %>%
  dplyr::summarise(
    `CoverageB beta1 (95%)` = mean(cover1, na.rm = TRUE),
    `CoverageB beta2 (95%)` = mean(cover2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SampleSize = sample_size,
    CensorRate = censor_label(censor_rate)
  ) %>%
  dplyr::select(SampleSize, CensorRate,
                `Coverage beta1 (95%)`, `Coverage beta2 (95%)`)



# ---- Spatial-field metrics (computed only for the spatial model) ----
# Tie f_true to grid_id once
truth_grid <- tibble(grid_id = grid$id, f_true = f_true)

# --- L2 metrics (from averaged field & SD field) ---
field_stats2 <- fhats_long %>%
  dplyr::left_join(truth_grid, by = "grid_id") %>%
  dplyr::group_by(sample_size, censor_rate,id_simulation) %>%
  dplyr::summarise(
    normL2 = integrate_over_omega_from_barycenters( (f_hat-f_true)^2),
    .groups   = "drop"
  ) %>%
  dplyr::group_by(sample_size, censor_rate) %>%
  dplyr::summarise(
    `L2 RMSE` = mean(normL2),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SampleSize = sample_size,
    CensorRate = censor_label(censor_rate)
  ) %>%
  dplyr::select(SampleSize, CensorRate,`L2 RMSE`)


# --- L2 metrics (from averaged field & SD field) ---
field_stats2gam <- fhats_long %>%
  dplyr::left_join(truth_grid, by = "grid_id") %>%
  dplyr::group_by(sample_size, censor_rate,id_simulation) %>%
  dplyr::summarise(
    normL2 = integrate_over_omega_from_barycenters( (f_gam-f_true)^2),
    .groups   = "drop"
  ) %>%
  dplyr::group_by(sample_size, censor_rate) %>%
  dplyr::summarise(
    `L2 RMSE GAM` = mean(normL2),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SampleSize = sample_size,
    CensorRate = censor_label(censor_rate)
  ) %>%
  dplyr::select(SampleSize, CensorRate,`L2 RMSE GAM`)

# # --- Linfinity metrics (per-simulation max over grid, then average across sims) ---
# field_statsInf <- fhats_long %>%
#   dplyr::left_join(truth_grid, by = "grid_id") %>%
#   dplyr::group_by(sample_size, censor_rate, id_simulation) %>%   # <-- important
#   dplyr::summarise(
#     Linf_sim = max(abs(f_true - f_hat)),
#     .groups  = "drop"
#   ) %>%
#   dplyr::group_by(sample_size, censor_rate) %>%
#   dplyr::summarise(
#     `Linf Bias Spatial Field` = mean(Linf_sim),
#     `Linf SD Spatial Field`   = sd(Linf_sim),   # use var(Linf_sim) if you truly want variance
#     .groups = "drop"
#   ) %>%
#   dplyr::mutate(
#     SampleSize = sample_size,
#     CensorRate = censor_label(censor_rate)
#   ) %>%
#   dplyr::select(SampleSize, CensorRate,
#                 `Linf Bias Spatial Field`, `Linf SD Spatial Field`)

# --- Merge ---
field_stats <- field_stats2 %>%
  dplyr::left_join(field_stats2gam, by = c("SampleSize","CensorRate"))

# ---- Combine into the requested wide table ----
coef_spatial  <- coef_stats %>%
  filter(model == "spatial") %>%
  select(SampleSize, CensorRate,
         `Bias beta1 (Spatial)` = `Bias beta1`,
         `Bias beta2 (Spatial)` = `Bias beta2`,
         `MSE beta (Spatial)`   = `MSE beta`)

coef_standard <- coef_stats %>%
  filter(model == "standard") %>%
  select(SampleSize, CensorRate,
         `Bias beta1 (Standard)` = `Bias beta1`,
         `Bias beta2 (Standard)` = `Bias beta2`,
         `MSE beta (Standard)`   = `MSE beta`)

coef_gam <- coef_stats %>%
  filter(model == "gam") %>%
  select(SampleSize, CensorRate,
         `Bias beta1 (GAM)` = `Bias beta1`,
         `Bias beta2 (GAM)` = `Bias beta2`,
         `MSE beta (GAM)`   = `MSE beta`)

final_table <- coef_spatial %>%
  left_join(field_stats, by = c("SampleSize","CensorRate")) %>%
  left_join(coef_standard, by = c("SampleSize","CensorRate")) %>%
  left_join(coef_gam, by = c("SampleSize","CensorRate")) %>%
  arrange(SampleSize, CensorRate) %>%
  select(
    `Sample Size` = SampleSize,
    `Censor rate` = CensorRate,
    `Bias beta1 (Spatial)`,
    `Bias beta2 (Spatial)`,
    `MSE beta (Spatial)`,
    `Bias beta1 (GAM)`,
    `Bias beta2 (GAM)`,
    `MSE beta (GAM)`,
    `L2 RMSE`,
    `L2 RMSE GAM`,
    `Bias beta1 (Standard)`,
    `Bias beta2 (Standard)`,
    `MSE beta (Standard)`
  )


final_table <- final_table %>%
  dplyr::left_join(cov_spatial,
                   by = c("Sample Size" = "SampleSize",
                          "Censor rate" = "CensorRate"))
# View / Save
print(as.data.frame(final_table))
final_table = as.data.frame(final_table)
write.csv(final_table, "final_table.csv")

load("C:/Users/utente/Documents/myEnvironment.RData")
verbose = FALSE
if(verbose){

  avg_fields <- fhats_long %>%
    dplyr::left_join(truth_grid, by = "grid_id") %>%
    dplyr::group_by(sample_size, censor_rate, grid_id) %>%
    dplyr::summarise(
      f_hat_mean  = mean(f_hat, na.rm = TRUE),
      f_true_mean = mean(f_true, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      SampleSize = sample_size,
      CensorRate = censor_label(censor_rate)
    ) %>%
    dplyr::select(SampleSize, CensorRate, grid_id, f_hat_mean, f_true_mean)


  avg_field_list <- avg_fields %>%
    dplyr::group_split(SampleSize, CensorRate)
  head(avg_field_list[[1]])
  names(grid)
  grid = left_join( grid[, c("x", "y", "id")], as.data.frame(avg_field_list[[1]]), by = c("id" = "grid_id"))
  head(grid)

  beta_true = c(.25,-1)
  censor_rate   <- censor_rates[1]
  sample_size   <- 500
  id_simulation <- 1
  df <- make_spatial_cox_df(
    n = sample_size,
    censor_rate = censor_rate,
    seed = id_simulation,
    beta1 = beta_true[1],
    beta2 = beta_true[2],
    lambda0 = 0.0001
  )
  (censor_lvl <- 1 - mean(df$status))
  coords <- as.matrix(df[, c("x", "y")])

  mc <- mesh_cache[[as.character(sample_size)]]
  FEMbasis <- mc$FEMbasis
  Z        <- mc$Z
  Kred     <- mc$Kred
  Agrid    <- mc$Agrid
  mesh_2 = mc$mesh
  plot(mesh_2, pch=".")
  points(df$x, df$y, col = ifelse(df$status, "blue", "red"), pch = 16)
  
  ################ PLOT MESH WITH STATUS POINTS
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
    # boundary emphasized (optional)
    geom_sf(data = boundary_sf, color = "black", linewidth = 2) +
    # your points, blue if status==1 else red
    geom_point(data = df,
               aes(x = x, y = y, color = status == 0),
               size = 2.5) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                        name = "Censored") +
    coord_sf() +
    theme_minimal()
    #theme(legend.position = "none")

  #ggsave("mesh_full_with_points.png", p, width = 7, height = 6, dpi = 300)
  p_mesh
  ggsave(
    filename = "p_mesh.png",
    plot     = p_mesh,
    width    = 10,
    height   = 8,
    dpi      = 300,
    units    = "in"
  )

  
  ################ PLOT F TRUE
  
  pts <- sf::st_sample(horseshoe_sf, 10000, type = "regular") 
  plot(mesh_2, pch=".")
  points(sf::st_coordinates(pts), col="black", pch='.')
  grid2 <- as.data.frame(sf::st_coordinates(pts))
  rm(pts)
  names(grid2) = c("x", "y")
  Agrid2 <- eval_FEM_basis_all(FEMbasis, grid2[,c("x", "y")])
  grid2$f_true = Agrid2%*% qr.solve(Agrid, grid$f_true_mean)
  grid2$f_hat_mean = Agrid2%*% qr.solve(Agrid, grid$f_hat_mean)
  grid2$diff2 =  Agrid2%*% qr.solve(Agrid, grid$f_hat_mean-grid$f_true_mean)
  # find the global min/max across all three variables
  # lims <- range(
  #   c(grid2$f_true, grid$f_hat, grid$f_hat - grid$f_true),
  #   na.rm = TRUE
  # )
  lims = round(range(c(grid2$f_true, grid2$f_hat_mean)),1) + c(-.1,+.1)
  p_true <- ggplot(grid2, aes(x = x, y = y, color = f_true)) +
    geom_point(size = 2) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = lims,
      name = "Value"
    ) +
    geom_sf(data = horseshoe_boundary, inherit.aes = FALSE,
            color = "black", linewidth = 2) +
    coord_sf() +
    theme_minimal()

  p_true
  
  
 
  ggsave(
    filename = "p_true.png",
    plot     = p_true,
    width    = 10,
    height   = 8,
    dpi      = 300,
    units    = "in"
  )
  
  p_hat <- ggplot(grid2, aes(x = x, y = y, color = f_hat_mean)) +
    geom_point(size = 2) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = lims,
      name = "Value"
    ) +
    geom_sf(data = horseshoe_boundary, inherit.aes = FALSE,
            color = "black", linewidth = 2) +
    coord_sf() +
    theme_minimal()
  p_hat
  ggsave(
    filename = "p_hat.png",
    plot     = p_hat,
    width    = 10,
    height   = 8,
    dpi      = 300,
    units    = "in"
  )
  lims = round(range(grid2$diff2),1) + c(-.1,.1)
  p_diff <- ggplot(grid2, aes(x = x, y = y, color = diff2)) +
    geom_point(size = 2) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0,
      limits = lims,
      name = "Value"
    ) +
    geom_sf(data = horseshoe_boundary, inherit.aes = FALSE,
            color = "black", linewidth = 2) +
    coord_sf() +
    theme_minimal()
  p_diff
  ggsave(
    filename = "p_diff.png",
    plot     = p_diff,
    width    = 10,
    height   = 8,
    dpi      = 300,
    units    = "in"
  )

}


