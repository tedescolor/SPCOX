# ==============================================================================
# script_functions.R
# Penalized Cox Model with Optimized Rcpp/Armadillo Backend & Parallel CV
# ==============================================================================

# 1. Load Required Libraries
if (!requireNamespace("Rcpp", quietly = TRUE)) install.packages("Rcpp")
if (!requireNamespace("RcppArmadillo", quietly = TRUE)) install.packages("RcppArmadillo")
if (!requireNamespace("doRNG", quietly = TRUE)) install.packages("doRNG")

library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(foreach)
library(doParallel)
library(doRNG)

# ==============================================================================
# 2. C++ Implementation (Stored as string for parallel export)
# ==============================================================================

cpp_code_source <- '
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------------
// Objective Function: Negative Penalized Partial Likelihood
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double negPL_pen_cpp(arma::vec par, arma::mat X, arma::mat Phi, arma::mat Kred,
                     arma::vec status, double lambda) {

  int p = X.n_cols;
  int n = X.n_rows;

  // Split parameters into beta and theta
  arma::vec beta = par.head(p);
  arma::vec theta = par.tail(par.n_elem - p);

  // Linear predictor: eta = X*beta + Phi*theta
  arma::vec eta = X * beta + Phi * theta;
  arma::vec exp_eta = exp(eta);

  // OPTIMIZATION: Single pass backwards to compute Partial Likelihood
  // Data MUST be pre-sorted by time (descending or ascending handled by loop)
  // Here we assume input is sorted by time ascending (earliest to latest).
  // We iterate backwards (latest to earliest) to build the risk set.
  
  double current_risk_sum = 0.0;
  double log_likelihood = 0.0;

  for (int i = n - 1; i >= 0; --i) {
    current_risk_sum += exp_eta[i]; // Add current individual to risk set

    // If this individual had an event, add their contribution
    if (status[i] == 1) {
      // contribution = eta_i - log(sum_{j in Risk} exp(eta_j))
      log_likelihood += (eta[i] - log(current_risk_sum));
    }
  }

  // Penalty term: -0.5 * lambda * theta^T * K * theta
  double pen = -0.5 * lambda * as_scalar(theta.t() * Kred * theta);

  // Return NEGATIVE objective for minimization: - (LL + Pen)
  return -log_likelihood - pen; 
}

// -----------------------------------------------------------------------------
// Gradient Function
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec negPL_pen_grad_cpp(arma::vec par, arma::mat X, arma::mat Phi, arma::mat Kred,
                             arma::vec status, double lambda) {

  int p = X.n_cols;
  int q = Phi.n_cols;
  int n = X.n_rows;

  arma::vec beta = par.head(p);
  arma::vec theta = par.tail(par.n_elem - p);

  arma::vec eta = X * beta + Phi * theta;
  arma::vec exp_eta = exp(eta);

  // Initialize gradients
  arma::vec grad_beta = arma::zeros(p);
  arma::vec grad_theta = arma::zeros(q);

  // OPTIMIZATION: Single pass backwards
  double risk_sum = 0.0;
  arma::vec x_weighted_sum = arma::zeros(p);
  arma::vec phi_weighted_sum = arma::zeros(q);

  for (int i = n - 1; i >= 0; --i) {
    double w = exp_eta[i];

    // Update Risk Set Sums
    risk_sum += w;
    x_weighted_sum += w * X.row(i).t();      
    phi_weighted_sum += w * Phi.row(i).t(); 

    if (status[i] == 1) {
      // Gradient contribution: covariate_i - (weighted_sum / risk_sum)
      grad_beta += (X.row(i).t() - (x_weighted_sum / risk_sum));
      grad_theta += (Phi.row(i).t() - (phi_weighted_sum / risk_sum));
    }
  }

  // Subtract Penalty Gradient from theta: lambda * K * theta
  grad_theta -= lambda * (Kred * theta);

  // Combine beta and theta gradients
  arma::vec total_grad = join_cols(grad_beta, grad_theta);

  // Return NEGATIVE gradient for minimizer
  return -total_grad;
}
'

# Compile C++ code for the current session
Rcpp::sourceCpp(code = cpp_code_source)

# ==============================================================================
# 3. Optimized R Wrappers
# ==============================================================================

# --- Optimized Fit Function ---
fit_pen_cox_opt <- function(X, Phi, Kred, time, status, lambda,
                            start = NULL, control = list(), hessian = FALSE) {
  
  Kred <- as.matrix(Kred)
  Phi  <- as.matrix(Phi)
  pX   <- ncol(X)
  pPhi <- ncol(Phi)
  
  if (is.null(start)) start <- rep(0, pX + pPhi)
  
  # CRITICAL: Pre-sort data by time once
  ord <- order(time)
  X_ord      <- X[ord, , drop = FALSE]
  Phi_ord    <- Phi[ord, , drop = FALSE]
  status_ord <- status[ord]
  
  # Wrappers calling the C++ functions
  fn_obj <- function(par) {
    negPL_pen_cpp(par, X_ord, Phi_ord, Kred, status_ord, lambda)
  }
  
  fn_grad <- function(par) {
    negPL_pen_grad_cpp(par, X_ord, Phi_ord, Kred, status_ord, lambda)
  }
  
  # Run Optimization
  opt <- optim(
    par     = start,
    fn      = fn_obj,
    gr      = fn_grad,
    method  = "BFGS",
    #control = modifyList(list(maxit = 500, reltol = 1e-7), control),
    hessian = hessian
  )
  
  list(
    beta        = opt$par[seq_len(pX)],
    theta       = opt$par[-seq_len(pX)],
    value       = opt$value,
    counts      = opt$counts,
    convergence = opt$convergence,
    hessian     = opt$hessian
  )
}

# --- Optimized Log-Likelihood Helper (for CV scoring) ---
logPL_cpp <- function(beta, theta, X, Phi, time, status) {
  # Sort data for C++ requirement
  ord <- order(time)
  
  # Call negPL_pen_cpp with lambda = 0
  # It returns: - (LogLikelihood - 0) = -LogLikelihood
  neg_ll <- negPL_pen_cpp(
    par    = c(beta, theta),
    X      = X[ord, , drop = FALSE],
    Phi    = Phi[ord, , drop = FALSE],
    Kred   = matrix(0, nrow=ncol(Phi), ncol=ncol(Phi)), # Dummy Kred
    status = status[ord],
    lambda = 0
  )
  
  # Return Normalized Log-Likelihood (+LL / n)
  return( (-neg_ll) / nrow(X) )
}

# ==============================================================================
# 4. Fast Parallel Cross-Validation
# ==============================================================================

select_lambda_cvpld_opt <- function(
    X, Phi, Kred, time, status,
    K = 5,
    lambda_grid = NULL,
    seed = 123,
    workers = max(1, parallel::detectCores() - 1),
    criterion = c("val","diff") 
) {
  criterion <- match.arg(criterion)
  
  # Sanity checks
  stopifnot(is.matrix(X), is.matrix(Phi), is.matrix(Kred))
  stopifnot(nrow(X) == length(time), length(status) == length(time))
  
  set.seed(seed)
  n <- nrow(X)
  
  if (is.null(lambda_grid)) {
    lam0 <- n^(-0.5)
    lambda_grid <- seq(0.1 * lam0, 10 * lam0, length.out = 10)
  } else {
    lambda_grid <- as.numeric(lambda_grid)
  }
  
  # Folds
  folds <- sample(rep(seq_len(K), length.out = n))
  idx_tr_list  <- lapply(seq_len(K), function(k) which(folds != k))
  idx_val_list <- lapply(seq_len(K), function(k) which(folds == k))
  
  # Parallel Setup
  use_parallel <- workers > 1
  if (use_parallel) {
    cl <- parallel::makeCluster(workers)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    doParallel::registerDoParallel(cl)
    parallel::clusterSetRNGStream(cl, iseed = seed)
    
    # 1. Export the C++ code string to workers
    parallel::clusterExport(cl, "cpp_code_source", envir = environment())
    
    # 2. Compile C++ on every worker
    parallel::clusterEvalQ(cl, {
      library(Rcpp)
      library(RcppArmadillo)
      sourceCpp(code = cpp_code_source)
    })
    
    # 3. Export R helper functions and data
    to_export <- c("X", "Phi", "Kred", "time", "status",
                   "fit_pen_cox_opt", "logPL_cpp", "criterion")
    parallel::clusterExport(cl, varlist = to_export, envir = environment())
  }
  
  warm_list <- vector("list", K)
  cvpld <- numeric(length(lambda_grid))
  
  # Iterate over lambda grid
  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    
    # Define work for one fold
    one_fold <- function(k) {
      idx_tr  <- idx_tr_list[[k]]
      idx_val <- idx_val_list[[k]]
      
      # Fit (Train)
      fit_k <- fit_pen_cox_opt(
        X      = X[idx_tr, , drop = FALSE],
        Phi    = Phi[idx_tr, , drop = FALSE],
        Kred   = Kred,
        time   = time[idx_tr],
        status = status[idx_tr],
        lambda = lam,
        start  = warm_list[[k]]
      )
      
      # Score
      if (criterion == "val") {
        # Val-only: -2 * LL_val
        ll_val <- logPL_cpp(
          fit_k$beta, fit_k$theta,
          X[idx_val, , drop = FALSE],
          Phi[idx_val, , drop = FALSE],
          time[idx_val], status[idx_val]
        )
        crit <- (-2) * ll_val
        
      } else { 
        # Diff: -2 * (LL_full - LL_train)
        ll_full <- logPL_cpp(fit_k$beta, fit_k$theta, X, Phi, time, status)
        ll_tr   <- logPL_cpp(
          fit_k$beta, fit_k$theta,
          X[idx_tr, , drop = FALSE],
          Phi[idx_tr, , drop = FALSE],
          time[idx_tr], status[idx_tr]
        )
        crit <- (-2) * (ll_full - ll_tr)
      }
      
      list(crit = crit, warm = c(fit_k$beta, fit_k$theta))
    }
    
    # Execute folds (Parallel or Serial)
    if (use_parallel) {
      fold_res <- doRNG::`%dorng%`(
        foreach::foreach(k = seq_len(K)),
        one_fold(k)
      )
    } else {
      fold_res <- lapply(seq_len(K), one_fold)
    }
    
    # Aggregate results
    cvpld[i]   <- sum(vapply(fold_res, `[[`, numeric(1), "crit"))
    warm_list  <- lapply(fold_res, `[[`, "warm")
  }
  
  list(
    lambda_opt = lambda_grid[which.min(cvpld)],
    grid       = lambda_grid,
    cvpld      = cvpld
  )
}