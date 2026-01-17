#' ThesisPackage
#'
#' Source functions for data generation, Lavaan model specifications and plotting.
#'
#' @docType _PACKAGE
#' @name ThesisPackage
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline labs
#'   theme_minimal theme element_text xlim ylim
NULL


# =============================================================================
# Data Generating Mechanism
# =============================================================================

#' Data Generating Mechanism
#'
#' Simulates data according to procedure by Bailey et al,2024.
#'
#' @param N Number of participants.
#' @param autocorr_effects Autoregressive effect size.
#' @param timepoints Number of timepoints.
#' @param n_x_confounders Number of X confounders.
#' @param n_y_confounders Number of Y confounders.
#' @param beta_focal Cross-lagged focal effect.
#' @param beta_covfocal Covariance among focal variables.
#' @param beta_covnofocal Covariance among non-focal variables.
#' @param sigma Residual variance.
#' @param seed Random seed.
#'
#' @return A list containing simulated data, parameters, and coefficient matrix.
#' @export
DGM <- function(
    N = 20000,
    autocorr_effects = 0.2,
    timepoints = 100,
    n_x_confounders = 5,
    n_y_confounders = 5,
    beta_focal = 0.1,
    beta_covfocal = 0.1,
    beta_covnofocal = 0.01,
    sigma = 0.5,
    seed = 427
) {
  
  set.seed(seed)
  
  T <- timepoints
  c1 <- n_x_confounders
  c2 <- n_y_confounders
  alpha <- autocorr_effects
  
  irow <- function(x) matrix(x, nrow = sqrt(length(x)), byrow = TRUE)
  row  <- function(x) as.vector(t(x))
  
  v <- c1 + c2 + 2
  
  A <- matrix(NA, v, v)
  diag(A) <- alpha
  A[lower.tri(A)] <- beta_covnofocal
  A[upper.tri(A)] <- beta_covnofocal
  
  c1t <- 3:(2 + c1)
  c2t <- (3 + c1):v
  
  A[c1t, 1] <- beta_covfocal
  A[1, c1t] <- beta_covfocal
  A[c2t, 2] <- beta_covfocal
  A[2, c2t] <- beta_covfocal
  
  A[c1t, c1t][lower.tri(A[c1t, c1t])] <- beta_covfocal
  A[c1t, c1t][upper.tri(A[c1t, c1t])] <- beta_covfocal
  A[c2t, c2t][lower.tri(A[c2t, c2t])] <- beta_covfocal
  A[c2t, c2t][upper.tri(A[c2t, c2t])] <- beta_covfocal
  
  A[1, 2] <- beta_focal
  A[2, 1] <- beta_focal
  
  max_eigenvalue <- max(eigen(A)$values)
  if (max_eigenvalue >= 1) {
    warning("Maximum eigenvalue >= 1; system may not be stationary.")
  }
  
  Sigma <- diag(rep(sigma, v))
  Sigma1 <- irow(solve(diag(v^2) - t(A) %x% t(A)) %*% row(Sigma))
  Mean1 <- matrix(0, v, 1)
  
  df <- matrix(NA, nrow = N, ncol = 2 * T)
  colnames(df) <- c(paste0("x", 1:T), paste0("y", 1:T))
  
  D <- mvtnorm::rmvnorm(N, mean = Mean1, sigma = Sigma1)
  df[, 1] <- D[, 1]
  df[, 1 + T] <- D[, 2]
  
  for (i in 2:T) {
    D <- D %*% t(A) + mvtnorm::rmvnorm(N, sigma = Sigma)
    df[, i] <- D[, 1]
    df[, i + T] <- D[, 2]
  }
  
  list(
    data = df,
    parameters = list(
      N = N,
      timepoints = T,
      n_x_confounders = c1,
      n_y_confounders = c2,
      autocorr_effects = alpha,
      beta_focal = beta_focal,
      beta_covfocal = beta_covfocal,
      beta_covnofocal = beta_covnofocal,
      sigma = sigma,
      StatCOVmatrix = Sigma1,
      max_eigenvalue = max_eigenvalue
    ),
    coefficient_matrix = A
  )
}


# =============================================================================
# Util
# =============================================================================

#' Create Metadata Container
#'
#' Creates a list containing session metadata and an empty results object.
#'
#' @return A list with metadata and results.
#' @export
create_metadata <- function() {
  list(
    metadata = list(
      timestamp = Sys.time(),
      r_version = R.version.string,
      session_info = sessionInfo()
    ),
    results = data.frame()
  )
}


#' Save Results and Clear Memory
#'
#' Saves an R object to disk as an RDS file and triggers garbage collection.
#'
#' @param results Object to save.
#' @param filename_prefix Prefix for the saved file.
#'
#' @return Filename of saved object.
#' @export
save_clear <- function(results, filename_prefix = "sim_results") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0(filename_prefix, "_", timestamp, ".rds")
  saveRDS(results, file = filename)
  gc(FALSE)
  filename
}


# =============================================================================
# Model Specification Generators
# =============================================================================


#' Generate CLPM model specification
#'
#' Generates a Cross-Lagged Panel Model (CLPM) specification in lavaan syntax.
#'
#' @param T Integer. Number of timepoints (>= 2).
#'
#' @return A character string containing the lavaan model specification.
#' @export
generate_clpm <- function(T){
  
  # Input validation
  if (!is.numeric(T) || T <= 1 || T %% 1 != 0) {
    stop("T must be an integer greater or equal to 2.")
  }
  
  # Autoregressive and cross-lagged paths
  if (T > 1) {
    timepoints_reg <- 2:T
    lagged_timepoints <- 1:(T - 1)
    
    paths_x <- paste(
      sprintf("x%d ~ ax*x%d + by*y%d",
              timepoints_reg, lagged_timepoints, lagged_timepoints),
      collapse = "\n"
    )
    
    paths_y <- paste(
      sprintf("y%d ~ bx*x%d + ay*y%d",
              timepoints_reg, lagged_timepoints, lagged_timepoints),
      collapse = "\n"
    )
  } else {
    paths_x <- ""
    paths_y <- ""
  }
  
  # Variances
  var_x <- paste(sprintf("x%d ~~ x%d", 1:T, 1:T), collapse = "\n")
  var_y <- paste(sprintf("y%d ~~ y%d", 1:T, 1:T), collapse = "\n")
  
  # Covariances
  cov_xy <- paste(
    sprintf("x%d ~~ x%d_y%d_cov*y%d", 1:T, 1:T, 1:T, 1:T),
    collapse = "\n"
  )
  
  components <- c(
    "# 1. Autoregressive and cross-lagged paths",
    paths_x,
    paths_y,
    "\n# 2. Variances of observed variables",
    var_x,
    var_y,
    "\n# 3. Covariances between x and y at each timepoint",
    cov_xy
  )
  
  model_string <- paste(components[components != ""], collapse = "\n")
  return(model_string)
}

# ------------------------------------------------------------------------------

#' Generate RI-CLPM model specification
#'
#' Generates a Random-Intercept Cross-Lagged Panel Model (RI-CLPM)
#' specification in lavaan syntax.
#'
#' @param T Integer. Number of timepoints (>= 3).
#'
#' @return A character string containing the lavaan model specification.
#' @export
generate_riclpm <- function(T) {
  
  if (!is.numeric(T) || T <= 2 || T %% 1 != 0) {
    stop("T must be an integer greater or equal to 3.")
  }
  
  ri_x <- sprintf("RIx =~ %s",
                  paste(sprintf("1*x%d", 1:T), collapse = " + "))
  ri_y <- sprintf("RIy =~ %s",
                  paste(sprintf("1*y%d", 1:T), collapse = " + "))
  
  wx_defs <- paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse = "\n")
  wy_defs <- paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse = "\n")
  
  timepoints_reg <- 2:T
  lagged_timepoints <- 1:(T - 1)
  
  paths_wx <- paste(
    sprintf("wx%d ~ ax*wx%d + by*wy%d",
            timepoints_reg, lagged_timepoints, lagged_timepoints),
    collapse = "\n"
  )
  
  paths_wy <- paste(
    sprintf("wy%d ~ bx*wx%d + ay*wy%d",
            timepoints_reg, lagged_timepoints, lagged_timepoints),
    collapse = "\n"
  )
  
  res_covs <- paste(sprintf("wx%d ~~ ur*wy%d", 2:T, 2:T), collapse = "\n")
  
  var_wx <- paste(sprintf("wx%d ~~ wx%d", 1:T, 1:T), collapse = "\n")
  var_wy <- paste(sprintf("wy%d ~~ wy%d", 1:T, 1:T), collapse = "\n")
  
  zero_var_x <- paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse = "\n")
  zero_var_y <- paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse = "\n")
  
  components <- c(
    "# 1. Random Intercepts", ri_x, ri_y,
    "\n# 2. Within-person components (measurement)", wx_defs, wy_defs,
    "\n# 3. Autoregressive and cross-lagged paths", paths_wx, paths_wy,
    "\n# 4. Covariances",
    "wx1 ~~ T1_cov*wy1",
    res_covs,
    "\n# 5. (Co)variances of Random Intercepts",
    "RIx ~~ varRIx*RIx",
    "RIy ~~ varRIy*RIy",
    "RIx ~~ covRI*RIy",
    "\n# 6. (Residual) variances of within-person components",
    var_wx,
    var_wy,
    "\n# 7. Fix observed variances to zero",
    zero_var_x,
    zero_var_y,
    "\n# 8. Fix RI covariances with first state to zero",
    "wx1 ~~ 0*RIx", "wx1 ~~ 0*RIy",
    "wy1 ~~ 0*RIx", "wy1 ~~ 0*RIy"
  )
  
  model_string <- paste(components[components != ""], collapse = "\n")
  return(model_string)
}

# ------------------------------------------------------------------------------

#' Generate segmented RI-CLPM (MRI-CLPM) model specification
#'
#' Generates a segmented (multi–random-intercept) RI-CLPM specification.
#'
#' @param T Integer. Total number of timepoints (>= 3).
#' @param timespan Integer. Length of each segment.
#' @param across_seg_cov Character. One of "equal", "zero", "free",
#'   "toeplitz", or "toeplitz_cross".
#'
#' @return A character string containing the lavaan model specification.
#' @export
generate_segmented_riclpm <- function(T, timespan = 3, across_seg_cov = "equal") {
  
  if (!is.numeric(T) || T <= 2 || T %% 1 != 0) {
    stop("T must be an integer greater or equal to 3.")
  }
  if (!is.numeric(timespan) || timespan <= 0 || timespan %% 1 != 0) {
    stop("timespan must be a positive integer.")
  }
  if (timespan > T) {
    stop("timespan cannot be greater than T.")
  }
  if (!across_seg_cov %in% c("equal", "zero", "free", "toeplitz", "toeplitz_cross")) {
    stop("across_seg_cov must be either 'equal', 'zero', 'free', 'toeplitz', or 'toeplitz_cross'.")
  }
  
  n_segments <- ceiling(T / timespan)
  timepoint_to_segment <- rep(1:n_segments, each = timespan)[1:T]
  
  ri_x_defs <- ri_y_defs <- character()
  
  for (seg in 1:n_segments) {
    tp <- which(timepoint_to_segment == seg)
    ri_x_defs <- c(
      ri_x_defs,
      sprintf("RIx%d =~ %s", seg, paste(sprintf("1*x%d", tp), collapse = " + "))
    )
    ri_y_defs <- c(
      ri_y_defs,
      sprintf("RIy%d =~ %s", seg, paste(sprintf("1*y%d", tp), collapse = " + "))
    )
  }
  
  wx_defs <- paste(sprintf("wx%d =~ 1*x%d", 1:T, 1:T), collapse = "\n")
  wy_defs <- paste(sprintf("wy%d =~ 1*y%d", 1:T, 1:T), collapse = "\n")
  
  timepoints_reg <- 2:T
  lagged_timepoints <- 1:(T - 1)
  
  paths_wx <- paste(
    sprintf("wx%d ~ ax*wx%d + by*wy%d",
            timepoints_reg, lagged_timepoints, lagged_timepoints),
    collapse = "\n"
  )
  paths_wy <- paste(
    sprintf("wy%d ~ bx*wx%d + ay*wy%d",
            timepoints_reg, lagged_timepoints, lagged_timepoints),
    collapse = "\n"
  )
  
  res_covs <- paste(sprintf("wx%d ~~ ur*wy%d", 2:T, 2:T), collapse = "\n")
  
  var_wx <- paste(sprintf("wx%d ~~ wx%d", 1:T, 1:T), collapse = "\n")
  var_wy <- paste(sprintf("wy%d ~~ wy%d", 1:T, 1:T), collapse = "\n")
  
  zero_var_x <- paste(sprintf("x%d ~~ 0*x%d", 1:T, 1:T), collapse = "\n")
  zero_var_y <- paste(sprintf("y%d ~~ 0*y%d", 1:T, 1:T), collapse = "\n")
  
  ri_vars <- c()
  for (seg in 1:n_segments) {
    ri_vars <- c(
      ri_vars,
      sprintf("RIx%d ~~ varRIx%d*RIx%d", seg, seg, seg),
      sprintf("RIy%d ~~ varRIy%d*RIy%d", seg, seg, seg),
      sprintf("RIx%d ~~ covRI%d*RIy%d", seg, seg, seg)
    )
  }
  
  ri_first_state_zero <- c()
  for (seg in 1:n_segments) {
    ri_first_state_zero <- c(
      ri_first_state_zero,
      sprintf("wx1 ~~ 0*RIx%d", seg),
      sprintf("wx1 ~~ 0*RIy%d", seg),
      sprintf("wy1 ~~ 0*RIx%d", seg),
      sprintf("wy1 ~~ 0*RIy%d", seg)
    )
  }
  
  components <- c(
    "# 1. Random Intercepts (Segmented)",
    paste(c(ri_x_defs, ri_y_defs), collapse = "\n"),
    "\n# 2. Within-person components",
    wx_defs,
    wy_defs,
    "\n# 3. Autoregressive and cross-lagged paths",
    paths_wx,
    paths_wy,
    "\n# 4. Covariances",
    "wx1 ~~ wy1",
    res_covs,
    "\n# 5. (Co)variances of Random Intercepts",
    paste(ri_vars, collapse = "\n"),
    "\n# 6. (Residual) variances of within-person components",
    var_wx,
    var_wy,
    "\n# 7. Fix observed variances to zero",
    zero_var_x,
    zero_var_y,
    "\n# 8. Fix RI covariances with first state to zero",
    paste(ri_first_state_zero, collapse = "\n")
  )
  
  model_string <- paste(components[components != ""], collapse = "\n")
  return(model_string)
}


# =============================================================================
# Plotting Functions
# =============================================================================

#' Plot X and Y Time Series
#'
#' Plots longitudinal X and Y trajectories for randomly sampled participants.
#'
#' @param data Matrix or data frame of simulated data.
#' @param timepoints Number of timepoints.
#' @param n_participants Number of participants to plot.
#'
#' @export
plot_xy_sequences <- function(data, timepoints, n_participants = 5) {
  participants <- sample(seq_len(nrow(data)), n_participants)
  
  par(mfrow = c(1, 2))
  matplot(t(data[participants, 1:timepoints]), type = "l", main = "X")
  matplot(t(data[participants, (timepoints + 1):(2 * timepoints)]), type = "l", main = "Y")
  par(mfrow = c(1, 1))
}


#' Plot Random Intercept Variance
#'
#' Line plot of RI variance across conditions.
#'
#' @param data Data frame.
#' @param x_var X-axis variable name.
#' @param group_var Grouping variable name.
#' @param x_label X-axis label.
#' @param group_label Legend label.
#' @param plot_title Plot title.
#'
#' @export
plot_varRIx <- function(data, x_var, group_var, x_label, group_label, plot_title,
                        xlim = NULL, ylim = NULL) {
  
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data[[x_var]],
      y = varRIx,
      color = factor(.data[[group_var]])
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = plot_title,
      x = x_label,
      y = "varRIx",
      color = group_label
    ) +
    ggplot2::theme_minimal()
  
  if (!is.null(xlim)) p <- p + ggplot2::xlim(xlim)
  if (!is.null(ylim)) p <- p + ggplot2::ylim(ylim)
  
  p
}


#' Plot Cross-lagged Bias
#'
#' Plots bias in cross-lagged parameters bx and by.
#'
#' @param data Data frame containing coefficient lists.
#' @param x_var X-axis variable.
#' @param group_var Grouping variable.
#' @param true_bx True bx value.
#' @param true_by True by value.
#' @param x_label X-axis label.
#' @param group_label Legend label.
#' @param title_suffix Title suffix.
#'
#' @export
plot_crosslag_bias <- function(
    data, x_var, group_var, true_bx, true_by,
    x_label, group_label, title_suffix,
    xlim = NULL, ylim = NULL
) {
  
  data$bx <- sapply(data$coefficients, function(x) x[["bx"]])
  data$by <- sapply(data$coefficients, function(x) x[["by"]])
  
  data$bias_bx <- data$bx - true_bx
  data$bias_by <- data$by - true_by
  
  list(
    bx_plot = ggplot2::ggplot(data, ggplot2::aes(.data[[x_var]], bias_bx)) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 0),
    
    by_plot = ggplot2::ggplot(data, ggplot2::aes(.data[[x_var]], bias_by)) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = 0)
  )
}

#' Function template
#'
#' A description of the function
#'
#' @param input Declare input parameters
#'
#' @return the output
#' @export
template <- function(input){
  output <- input * 2
  return(output)
}
