# ============================================================================
# Extreme Value Volatility Estimators Using Null Space Decomposition
# Author: Abhishek Chand
# Date: 2024
# 
# This code implements the methodology from:
# "A Computational Method to Generate Families of Extreme Value Volatility Estimators"
# ============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS)  # For matrix operations
library(quadprog)  # For quadratic programming
library(gridExtra)  # For combining plots

# Set seed for reproducibility
set.seed(42)

# ============================================================================
# SECTION 1: VISUALIZATION FUNCTIONS
# ============================================================================

#' Create Figure 1: Security Price Movement and Extreme Values
#' 
#' This function generates a diagram showing how OHLC values are extracted
#' from continuous price paths during a trading day
#' 
#' @param save_plot Logical, whether to save the plot to file
#' @param filename Character, filename for saved plot
#' @return ggplot object
create_price_movement_diagram <- function(save_plot = TRUE, 
                                          filename = "figure1_price_movement.pdf") {
  
  # Simulate a single day's price path
  n_steps <- 1000
  dt <- 1/n_steps
  
  # Parameters
  mu <- 0.1  # Annual drift (10%)
  sigma <- 0.2  # Annual volatility (20%)
  
  # Time points
  t <- seq(0, 1, length.out = n_steps)
  
  # Generate Brownian motion
  dW <- rnorm(n_steps - 1, 0, sqrt(dt))
  W <- c(0, cumsum(dW))
  
  # Price path (log-prices)
  P <- 100 * exp((mu - 0.5 * sigma^2) * t + sigma * W)
  
  # Trading starts at f = 0.3 (30% of day is closed)
  f <- 0.3
  trading_start <- floor(f * n_steps)
  
  # Extract extreme values
  C0 <- P[1]  # Yesterday's close
  O <- P[trading_start]  # Today's open
  
  # During trading hours
  trading_prices <- P[trading_start:n_steps]
  H <- max(trading_prices)  # High
  L <- min(trading_prices)  # Low
  C <- P[n_steps]  # Close
  
  # Find positions of H and L
  H_idx <- trading_start - 1 + which.max(trading_prices)
  L_idx <- trading_start - 1 + which.min(trading_prices)
  
  # Create data frame for plotting
  price_data <- data.frame(
    time = t,
    price = P,
    period = ifelse(t < f, "Closed", "Trading")
  )
  
  # Create the plot
  p <- ggplot(price_data, aes(x = time, y = price)) +
    # Shade the closed period
    annotate("rect", xmin = 0, xmax = f, ymin = -Inf, ymax = Inf, 
             fill = "gray90", alpha = 0.5) +
    
    # Price path
    geom_line(aes(color = period), size = 1) +
    
    # Mark extreme values
    geom_point(data = data.frame(x = c(0, f, t[H_idx], t[L_idx], 1),
                                 y = c(C0, O, H, L, C),
                                 label = c("C₀", "O", "H", "L", "C")),
               aes(x = x, y = y), size = 3, color = "red") +
    
    # Add labels
    geom_text(data = data.frame(x = c(0, f, t[H_idx], t[L_idx], 1),
                                y = c(C0, O, H, L, C),
                                label = c("C₀", "O", "H", "L", "C")),
              aes(x = x, y = y, label = label), 
              vjust = -1, hjust = 0.5, size = 4) +
    
    # Horizontal lines for normalized values
    geom_hline(yintercept = O, linetype = "dashed", alpha = 0.3) +
    
    # Annotations for normalized values
    annotate("segment", x = 1.05, xend = 1.05, y = O, yend = H,
             arrow = arrow(ends = "both", length = unit(0.1, "inches")),
             color = "blue") +
    annotate("text", x = 1.08, y = (O + H)/2, label = "u = H - O", 
             angle = 90, vjust = -0.5, color = "blue") +
    
    annotate("segment", x = 1.12, xend = 1.12, y = L, yend = O,
             arrow = arrow(ends = "both", length = unit(0.1, "inches")),
             color = "blue") +
    annotate("text", x = 1.15, y = (L + O)/2, label = "d = L - O", 
             angle = 90, vjust = -0.5, color = "blue") +
    
    # Labels and theme
    labs(x = "Time (fraction of day)",
         y = "Price ($)",
         title = "Security Price Movement and Extreme Values",
         subtitle = "Extracting OHLC information from continuous price paths") +
    scale_color_manual(values = c("Closed" = "gray60", "Trading" = "black"),
                       name = "Period") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12)) +
    scale_x_continuous(breaks = c(0, f, 0.5, 0.75, 1),
                       labels = c("0", "f", "0.5", "0.75", "1")) +
    coord_cartesian(xlim = c(0, 1.2))
  
  # Save if requested
  if (save_plot) {
    ggsave(filename, plot = p, width = 10, height = 6)
    message(paste("Plot saved to:", filename))
  }
  
  return(p)
}

# ============================================================================
# SECTION 2: SIMULATION FUNCTIONS
# ============================================================================

#' Simulate one day of price data
#' 
#' @param n_steps Number of time steps per day
#' @param sigma True volatility
#' @param mu Drift parameter
#' @param f Fraction of day market is closed
#' @return List containing OHLC values and normalized statistics
simulate_day <- function(n_steps = 700000, sigma = 0.2, mu = 0, f = 0.3) {
  
  dt <- 1/n_steps
  
  # Generate Brownian motion
  dW <- rnorm(n_steps - 1, 0, sqrt(dt))
  W <- c(0, cumsum(dW))
  
  # Price path (log-prices, starting at 0)
  P <- (mu - 0.5 * sigma^2) * seq(0, 1, length.out = n_steps) + sigma * W
  
  # Trading starts at time f
  trading_start <- floor(f * n_steps)
  
  # Extract values
  C0 <- P[1]  # Yesterday's close (= 0 for first day)
  O <- P[trading_start]  # Open
  
  # During trading hours
  trading_prices <- P[trading_start:n_steps]
  H <- max(trading_prices)
  L <- min(trading_prices)
  C <- P[n_steps]
  
  # Normalized values
  u <- H - O
  d <- L - O  
  c <- C - O
  o <- O - C0
  
  # Quadratic terms
  quadratic_terms <- c(
    u2 = u^2,
    d2 = d^2,
    c2 = c^2,
    ud = u*d,
    uc = u*c,
    dc = d*c,
    o2 = o^2
  )
  
  return(list(
    OHLC = c(O = O, H = H, L = L, C = C, C0 = C0),
    normalized = c(u = u, d = d, c = c, o = o),
    quadratic = quadratic_terms
  ))
}

#' Simulate many days and compute expectations
#' 
#' @param n_days Number of days to simulate
#' @param n_steps Number of time steps per day
#' @param sigma True volatility
#' @param mu Drift parameter
#' @param f Fraction of day market is closed
#' @return Data frame with expectations and covariances
simulate_expectations <- function(n_days = 20000, n_steps = 700000, 
                                  sigma = 0.2, mu = 0, f = 0.3) {
  
  # Pre-allocate matrix for quadratic terms
  quad_terms <- matrix(NA, nrow = n_days, ncol = 7)
  colnames(quad_terms) <- c("u2", "d2", "c2", "ud", "uc", "dc", "o2")
  
  # Simulate days
  message("Simulating ", n_days, " days with ", n_steps, " steps each...")
  pb <- txtProgressBar(min = 0, max = n_days, style = 3)
  
  for (i in 1:n_days) {
    day_data <- simulate_day(n_steps, sigma, mu, f)
    quad_terms[i, ] <- day_data$quadratic
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Compute expectations (normalized by sigma^2)
  expectations <- colMeans(quad_terms) / sigma^2
  
  # Compute variance-covariance matrix
  cov_matrix <- cov(quad_terms) / sigma^4
  
  return(list(
    expectations = expectations,
    covariance = cov_matrix,
    raw_data = quad_terms
  ))
}

# ============================================================================
# SECTION 3: NULL SPACE COMPUTATION
# ============================================================================

#' Construct the constraint matrix K for null space computation
#' 
#' @param terms Character vector of terms to include (e.g., c("u2", "d2", "ud"))
#' @param n_conditions Number of different volatility/drift conditions
#' @return Matrix K for finding null space
construct_K_matrix <- function(terms = c("u2", "d2", "ud"), n_conditions = 6) {
  
  # Different volatility values
  sigmas <- seq(0.1, 0.6, length.out = n_conditions)
  
  K <- matrix(NA, nrow = n_conditions, ncol = length(terms))
  
  for (i in 1:n_conditions) {
    # Get expectations for this volatility
    sim_data <- simulate_expectations(n_days = 1000, n_steps = 10000, 
                                      sigma = sigmas[i])
    K[i, ] <- sim_data$expectations[terms]
  }
  
  colnames(K) <- terms
  return(K)
}

#' Find null space basis vectors
#' 
#' @param K Constraint matrix
#' @return Matrix whose columns are basis vectors for null space
find_null_space <- function(K) {
  
  # Use SVD to find null space
  svd_K <- svd(K)
  
  # Tolerance for determining numerical rank
  tol <- max(dim(K)) * max(svd_K$d) * .Machine$double.eps
  
  # Rank is number of singular values above tolerance
  rank <- sum(svd_K$d > tol)
  
  # Null space dimension
  null_dim <- ncol(K) - rank
  
  if (null_dim == 0) {
    warning("No null space found - system is fully determined")
    return(NULL)
  }
  
  # Null space basis vectors are right singular vectors corresponding to zero singular values
  null_basis <- svd_K$v[, (rank + 1):ncol(K), drop = FALSE]
  
  return(null_basis)
}

#' Find a particular solution to Ka = sigma^2
#' 
#' @param K Constraint matrix
#' @param sigma2 Vector of sigma^2 values
#' @return Particular solution vector
find_particular_solution <- function(K, sigma2) {
  
  # Use least squares to find a particular solution
  particular <- solve(t(K) %*% K) %*% t(K) %*% sigma2
  
  return(as.vector(particular))
}

# ============================================================================
# SECTION 4: ESTIMATOR OPTIMIZATION
# ============================================================================

#' Find minimum variance estimator using null space approach
#' 
#' @param terms Character vector of terms to include
#' @param cov_matrix Variance-covariance matrix of quadratic terms
#' @param verbose Logical, whether to print results
#' @return List containing optimal weights and efficiency
find_optimal_estimator <- function(terms = c("u2", "d2", "ud"), 
                                   cov_matrix = NULL,
                                   verbose = TRUE) {
  
  # If no covariance matrix provided, simulate to get it
  if (is.null(cov_matrix)) {
    message("Simulating to obtain covariance matrix...")
    sim_data <- simulate_expectations(n_days = 20000, n_steps = 100000)
    cov_matrix <- sim_data$covariance[terms, terms]
  }
  
  # Construct constraint matrix
  K <- construct_K_matrix(terms)
  
  # Find null space basis
  null_basis <- find_null_space(K)
  
  # Find particular solution
  sigma2_vec <- rep(1, nrow(K))  # Arbitrary sigma^2 values
  particular <- find_particular_solution(K, sigma2_vec)
  
  if (is.null(null_basis)) {
    # No null space - particular solution is unique
    optimal_weights <- particular
  } else {
    # Form the complete solution: a = particular + null_basis * alpha
    # Need to find optimal alpha
    
    # Quadratic form: (p + B*alpha)' * Cov * (p + B*alpha)
    # Minimizing gives: alpha* = -(B'*Cov*B)^(-1) * B'*Cov*p
    
    B <- null_basis
    p <- particular
    
    BtCovB <- t(B) %*% cov_matrix %*% B
    BtCovp <- t(B) %*% cov_matrix %*% p
    
    alpha_optimal <- -solve(BtCovB) %*% BtCovp
    
    # Optimal weights
    optimal_weights <- p + B %*% alpha_optimal
  }
  
  # Compute variance of optimal estimator
  estimator_variance <- t(optimal_weights) %*% cov_matrix %*% optimal_weights
  
  # Efficiency relative to close-close estimator (variance = 2*sigma^4)
  efficiency <- 2 / estimator_variance
  
  if (verbose) {
    cat("\n=== Optimal Estimator ===\n")
    cat("Terms:", terms, "\n")
    cat("Weights:\n")
    for (i in 1:length(terms)) {
      cat(sprintf("  %s: %.6f\n", terms[i], optimal_weights[i]))
    }
    cat(sprintf("Variance: %.6f\n", estimator_variance))
    cat(sprintf("Efficiency: %.3f\n", efficiency))
  }
  
  return(list(
    weights = as.vector(optimal_weights),
    variance = as.numeric(estimator_variance),
    efficiency = as.numeric(efficiency),
    terms = terms
  ))
}

# ============================================================================
# SECTION 5: KNOWN ESTIMATORS FOR COMPARISON
# ============================================================================

#' Parkinson estimator
#' 
#' @param u Normalized high
#' @param d Normalized low
#' @return Volatility estimate
parkinson_estimator <- function(u, d) {
  return((u - d)^2 / (4 * log(2)))
}

#' Garman-Klass estimator
#' 
#' @param u Normalized high
#' @param d Normalized low  
#' @param c Normalized close
#' @return Volatility estimate
garman_klass_estimator <- function(u, d, c) {
  return(0.511 * (u - d)^2 - 0.019 * (c * (u + d) - 2 * u * d) - 0.383 * c^2)
}

#' Compare known estimators with our optimal estimator
#' 
#' @param n_days Number of simulation days
#' @return Data frame comparing estimators
compare_estimators <- function(n_days = 10000) {
  
  message("Simulating data for comparison...")
  
  # True volatility
  true_sigma <- 0.2
  
  # Simulate many days
  estimates <- data.frame(
    classical = numeric(n_days),
    parkinson = numeric(n_days),
    garman_klass = numeric(n_days),
    optimal_hl = numeric(n_days),
    optimal_ohlc = numeric(n_days)
  )
  
  for (i in 1:n_days) {
    day <- simulate_day(sigma = true_sigma)
    
    # Classical (close-close)
    estimates$classical[i] <- day$normalized["c"]^2
    
    # Parkinson
    estimates$parkinson[i] <- parkinson_estimator(
      day$normalized["u"], 
      day$normalized["d"]
    )
    
    # Garman-Klass
    estimates$garman_klass[i] <- garman_klass_estimator(
      day$normalized["u"],
      day$normalized["d"],
      day$normalized["c"]
    )
  }
  
  # Compute optimal estimators
  opt_hl <- find_optimal_estimator(terms = c("u2", "d2", "ud"), verbose = FALSE)
  opt_ohlc <- find_optimal_estimator(terms = c("u2", "d2", "c2", "ud", "uc", "dc"), 
                                     verbose = FALSE)
  
  # Apply optimal estimators
  for (i in 1:n_days) {
    day <- simulate_day(sigma = true_sigma)
    
    # Optimal HL
    quad_hl <- c(day$quadratic["u2"], day$quadratic["d2"], day$quadratic["ud"])
    estimates$optimal_hl[i] <- sum(opt_hl$weights * quad_hl)
    
    # Optimal OHLC
    quad_ohlc <- c(day$quadratic["u2"], day$quadratic["d2"], day$quadratic["c2"],
                   day$quadratic["ud"], day$quadratic["uc"], day$quadratic["dc"])
    estimates$optimal_ohlc[i] <- sum(opt_ohlc$weights * quad_ohlc)
  }
  
  # Compute variances and efficiencies
  results <- data.frame(
    Estimator = c("Classical", "Parkinson", "Garman-Klass", 
                  "Optimal (HL)", "Optimal (OHLC)"),
    Variance = c(
      var(estimates$classical),
      var(estimates$parkinson),
      var(estimates$garman_klass),
      var(estimates$optimal_hl),
      var(estimates$optimal_ohlc)
    )
  )
  
  results$Efficiency <- results$Variance[1] / results$Variance
  
  return(results)
}

# ============================================================================
# SECTION 6: MAIN ANALYSIS
# ============================================================================

# Example usage and analysis
if (FALSE) {  # Set to TRUE to run examples
  
  # 1. Create price movement diagram
  fig1 <- create_price_movement_diagram(save_plot = TRUE)
  print(fig1)
  
  # 2. Find optimal Parkinson-type estimator (using only H and L)
  opt_parkinson <- find_optimal_estimator(terms = c("u2", "d2", "ud"))
  
  # 3. Find optimal Garman-Klass-type estimator (using OHLC)
  opt_gk <- find_optimal_estimator(terms = c("u2", "d2", "c2", "ud", "uc", "dc"))
  
  # 4. Compare all estimators
  comparison <- compare_estimators(n_days = 5000)
  print(comparison)
  
  # 5. Create comparison plot
  comparison_plot <- ggplot(comparison, aes(x = reorder(Estimator, Efficiency), 
                                            y = Efficiency)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = sprintf("%.2f", Efficiency)), 
              vjust = -0.5, size = 4) +
    labs(x = "Estimator", y = "Efficiency (relative to classical)",
         title = "Efficiency Comparison of Volatility Estimators") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(comparison_plot)
}

# ============================================================================
# END OF MAIN CODE
# ============================================================================