# Simulation: Compare power of Hosmer-Lemeshow vs Ebrahim-Farrington
## Four scenarios: 1) Omitted quadratic term, 2) Omitted binary interaction, 
##                  3) Omitted continuous interaction, 4) Alternative link functions
## Methodology inspired by Hosmer et al. (1997)
##
## INCREMENTAL SIMULATION FEATURE:
## - Results are automatically saved to CSV files for each scenario and sample size
## - When re-running, only missing scenarios will be executed
## - Use show_simulation_status() to check progress
## - Use load_all_results() to load existing results without running simulations
## - Professional A4 summary plot with 4x3 grid (4 scenarios × 3 sample sizes)
## - Customizable test selection: modify 'include_tests' list to choose which tests appear in A4 plot
## - Use set_test_selection() function for quick test configuration
## ---- archive paths (inserted by make_archive.py; the convention is in README.md, "How to run") -------------
EDGE_ARCHIVE_ROOT <- local({
  r <- Sys.getenv("EDGE_ARCHIVE_ROOT")
  if (!nzchar(r)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    r <- if (length(f)) file.path(dirname(normalizePath(f[1], winslash = "/")), "..", "..") else getwd()
  }
  r <- normalizePath(r, winslash = "/", mustWork = FALSE)
  Sys.setenv(EDGE_ARCHIVE_ROOT = r)        # so worker processes started from here resolve the same root
  r
})
edge_path <- function(...) file.path(EDGE_ARCHIVE_ROOT, ...)
edge_battery <- function(...) {           # the author's battery/ folder: analysis/ -> results/analysis, rest -> results/blocks
  if (...length() == 0L) return(edge_path("results", "blocks"))
  p <- file.path(...)
  ifelse(p == "analysis" | startsWith(p, "analysis/"), edge_path("results", p), edge_path("results", "blocks", p))
}
edge_out <- function(...) { d <- edge_path("output", ...); dir.create(d, showWarnings = FALSE, recursive = TRUE); d }
## ---------------------------------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(parallel)
  library(ggplot2)
  # Hard-require these libraries per user request (no manual fallbacks)
  library(ResourceSelection)
  library(ebrahim.gof)
  # Additional libraries for professional plotting
  library(gridExtra)
  library(grid)
  library(cowplot)
})

# Source external tests per user request
try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "pigeonheyse.R"))), silent = TRUE)
try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "Hosmer (H) (equal width interval).R"))), silent = TRUE)

# Declare globals used in ggplot aesthetics and external functions to satisfy linters
utils::globalVariables(c("J", "Power", "Test", "Link", "I", "K", "Test_Label", "N_Label", "x", "y",
                         "pigeon_heyse_test", "hosmer_lemeshow_equal_intervals"))

# --- Configuration ---
sample_sizes <- c(500, 1000, 5000)
num_groups <- 10
alpha <- 0.05
num_trials <- 10000
J_values <- c(0.01, 0.02, 0.03, 0.05, 0.10, 0.20, 0.40)
I_values <- c(0.1, 0.3, 0.5, 0.7)
K_values <- c(0.1, 0.3, 0.5, 0.7)  # For continuous interaction scenario

# --- Test Selection for A4 Professional Plot ---
# Set to TRUE to include the test in the A4 summary plot, FALSE to exclude
include_tests <- list(
  "HL" = TRUE,                    # Hosmer-Lemeshow
  "EF" = TRUE,                    # Ebrahim-Farrington
  "Pigeon-Heyse" = TRUE,          # Pigeon-Heyse
  "HL_equal_width" = FALSE         # HL Equal-Width
)

# Example configurations (uncomment to use):
# Only traditional tests:
# include_tests <- list("HL" = TRUE, "EF" = FALSE, "Pigeon-Heyse" = FALSE, "HL_equal_width" = FALSE)

# Only modern tests:
# include_tests <- list("HL" = FALSE, "EF" = TRUE, "Pigeon-Heyse" = TRUE, "HL_equal_width" = FALSE)

# Compare HL variants:
# include_tests <- list("HL" = TRUE, "EF" = FALSE, "Pigeon-Heyse" = FALSE, "HL_equal_width" = TRUE)

# Function to quickly set test selection
set_test_selection <- function(HL = TRUE, EF = TRUE, PigeonHeyse = TRUE, HL_EqualWidth = TRUE) {
  include_tests <<- list(
    "HL" = HL,
    "EF" = EF,
    "Pigeon-Heyse" = PigeonHeyse,
    "HL_equal_width" = HL_EqualWidth
  )
  message("Test selection updated:")
  for (test in names(include_tests)) {
    status <- if (include_tests[[test]]) "INCLUDED" else "EXCLUDED"
    message("  ", test, ": ", status)
  }
}

# Output path helpers
get_base_dir <- function() {
  base_dir <- dirname(sys.frame(1)$ofile %||% "")
  if (is.null(base_dir) || base_dir == "") base_dir <- "."
  base_dir
}

plot_path_for_n <- function(n) {
  file.path(get_base_dir(), sprintf("Power_OmittedQuadratic_n%d.png", n))
}

plot_path_for_n_interaction <- function(n) {
  file.path(get_base_dir(), sprintf("Power_Interaction_n%d.png", n))
}

plot_path_for_n_altlink <- function(n) {
  file.path(get_base_dir(), sprintf("Power_AltLink_n%d.png", n))
}

plot_path_for_n_continuous_interaction <- function(n) {
  file.path(get_base_dir(), sprintf("Power_ContinuousInteraction_n%d.png", n))
}

# CSV path helpers
csv_path_omitted_quad <- function(n) {
  file.path(get_base_dir(), sprintf("Results_OmittedQuadratic_n%d.csv", n))
}

csv_path_interaction <- function(n) {
  file.path(get_base_dir(), sprintf("Results_Interaction_n%d.csv", n))
}

csv_path_alt_link <- function(n) {
  file.path(get_base_dir(), sprintf("Results_AltLink_n%d.csv", n))
}

csv_path_continuous_interaction <- function(n) {
  file.path(get_base_dir(), sprintf("Results_ContinuousInteraction_n%d.csv", n))
}

# Master results file
csv_path_master <- function() {
  file.path(get_base_dir(), "Master_Simulation_Results.csv")
}

# Null-coalescing helper
`%||%` <- function(a, b) if (is.null(a)) b else a

# Function to check if results exist and load them
load_existing_results <- function(scenario, n) {
  csv_path <- switch(scenario,
    "omitted_quad" = csv_path_omitted_quad(n),
    "interaction" = csv_path_interaction(n),
    "continuous_interaction" = csv_path_continuous_interaction(n),
    "alt_link" = csv_path_alt_link(n),
    stop("Unknown scenario: ", scenario)
  )
  
  if (file.exists(csv_path)) {
    message("Loading existing results from: ", csv_path)
    return(read.csv(csv_path, stringsAsFactors = FALSE))
  } else {
    message("No existing results found for ", scenario, " with n=", n)
    return(NULL)
  }
}

# Function to save results to CSV
save_results <- function(results, scenario, n) {
  csv_path <- switch(scenario,
    "omitted_quad" = csv_path_omitted_quad(n),
    "interaction" = csv_path_interaction(n),
    "continuous_interaction" = csv_path_continuous_interaction(n),
    "alt_link" = csv_path_alt_link(n),
    stop("Unknown scenario: ", scenario)
  )
  
  write.csv(results, csv_path, row.names = FALSE)
  message("Saved results to: ", csv_path)
}

# Function to get missing scenarios
get_missing_scenarios <- function(existing_results, scenario, n) {
  if (is.null(existing_results)) {
    # No existing results, need to run all scenarios
    return(switch(scenario,
      "omitted_quad" = J_values,
      "interaction" = I_values,
      "continuous_interaction" = K_values,
      "alt_link" = alt_links,
      stop("Unknown scenario: ", scenario)
    ))
  }
  
  # Check which scenarios are missing
  all_scenarios <- switch(scenario,
    "omitted_quad" = J_values,
    "interaction" = I_values,
    "continuous_interaction" = K_values,
    "alt_link" = alt_links,
    stop("Unknown scenario: ", scenario)
  )
  
  completed_scenarios <- switch(scenario,
    "omitted_quad" = unique(existing_results$J),
    "interaction" = unique(existing_results$I),
    "continuous_interaction" = unique(existing_results$K),
    "alt_link" = unique(existing_results$Link),
    stop("Unknown scenario: ", scenario)
  )
  
  missing <- setdiff(all_scenarios, completed_scenarios)
  if (length(missing) > 0) {
    message("Missing scenarios for ", scenario, " (n=", n, "): ", paste(missing, collapse = ", "))
  } else {
    message("All scenarios completed for ", scenario, " (n=", n, ")")
  }
  
  return(missing)
}

logit <- function(p) log(p/(1 - p))
inv_logit <- function(eta) 1/(1 + exp(-eta))

# Solve for (beta0, beta1, beta2) given constraints on probabilities at x = -1.5, 3, -3
# Constraints: pi(-1.5)=0.05, pi(3)=0.95, pi(-3)=J
solve_quadratic_betas <- function(J) {
  y1 <- logit(0.05)   # at x=-1.5
  y2 <- logit(0.95)   # at x=3
  y3 <- logit(J)      # at x=-3
  X <- rbind(
    c(1, -1.5, (-1.5)^2),
    c(1, 3.0,  (3.0)^2),
    c(1, -3.0, (-3.0)^2)
  )
  y <- c(y1, y2, y3)
  beta <- solve(X, y)
  names(beta) <- c("beta0", "beta1", "beta2")
  beta
}

# Data generator
generate_data <- function(beta, n) {
  x <- runif(n, min = -3, max = 3)
  eta_true <- beta["beta0"] + beta["beta1"] * x + beta["beta2"] * x^2
  p_true <- inv_logit(eta_true)
  y <- rbinom(n, size = 1, prob = p_true)
  data.frame(x = x, y = y)
}

# --- Interaction omitted scenario ---
# Constraints (Hosmer 1997):
# pi(-3,0)=0.1, pi(-3,1)=0.1, pi(3,0)=0.2, pi(3,1)=0.2+I
# Model: logit pi = b0 + b1 x + b2 d + b3 x d, with b2 = 3 b3 from equal-prob at x=-3
solve_interaction_betas <- function(I) {
  eta_m3_0 <- logit(0.1)        # x=-3, d=0
  eta_p3_0 <- logit(0.2)        # x= 3, d=0
  eta_p3_1 <- logit(0.2 + I)    # x= 3, d=1

  b0 <- (eta_m3_0 + eta_p3_0) / 2
  b1 <- (eta_p3_0 - eta_m3_0) / 6
  # From eta at (3,1): b0 + 3 b1 + 6 b3 = eta_p3_1
  b3 <- (eta_p3_1 - (b0 + 3 * b1)) / 6
  b2 <- 3 * b3
  beta <- c(beta0 = b0, beta1 = b1, beta2 = b2, beta3 = b3)
  beta
}

generate_data_interaction <- function(beta, n) {
  x <- runif(n, min = -3, max = 3)
  d <- rbinom(n, size = 1, prob = 0.5)
  eta <- beta["beta0"] + beta["beta1"] * x + beta["beta2"] * d + beta["beta3"] * x * d
  p <- inv_logit(eta)
  y <- rbinom(n, size = 1, prob = p)
  data.frame(x = x, d = d, y = y)
}

# --- Continuous interaction scenario ---
# Constraints: similar to binary interaction but with z ~ N(0,1)
# pi(-3,-2) = 0.1, pi(-3,0) = 0.1, pi(3,0) = 0.2, pi(3,2) = 0.2+K
solve_continuous_interaction_betas <- function(K) {
  eta_m3_m2 <- logit(0.1)        # x=-3, z=-2
  eta_m3_0 <- logit(0.1)         # x=-3, z=0  
  eta_p3_0 <- logit(0.2)         # x= 3, z=0
  eta_p3_p2 <- logit(0.2 + K)    # x= 3, z=2
  
  # System of equations: eta = b0 + b1*x + b2*z + b3*x*z
  # (-3,-2): b0 - 3*b1 - 2*b2 + 6*b3 = eta_m3_m2
  # (-3, 0): b0 - 3*b1 + 0*b2 + 0*b3 = eta_m3_0
  # ( 3, 0): b0 + 3*b1 + 0*b2 + 0*b3 = eta_p3_0
  # ( 3, 2): b0 + 3*b1 + 2*b2 + 6*b3 = eta_p3_p2
  
  X <- rbind(
    c(1, -3, -2,  6),  # x=-3, z=-2
    c(1, -3,  0,  0),  # x=-3, z=0
    c(1,  3,  0,  0),  # x=3,  z=0
    c(1,  3,  2,  6)   # x=3,  z=2
  )
  y <- c(eta_m3_m2, eta_m3_0, eta_p3_0, eta_p3_p2)
  beta <- solve(X, y)
  names(beta) <- c("beta0", "beta1", "beta2", "beta3")
  beta
}

generate_data_continuous_interaction <- function(beta, n) {
  x <- runif(n, min = -3, max = 3)
  z <- rnorm(n, mean = 0, sd = 1)  # Standard normal
  eta <- beta["beta0"] + beta["beta1"] * x + beta["beta2"] * z + beta["beta3"] * x * z
  p <- inv_logit(eta)
  y <- rbinom(n, size = 1, prob = p)
  data.frame(x = x, z = z, y = y)
}

# --- Alternative link functions scenario ---
alt_links <- c("probit", "cloglog", "stukel_long", "stukel_short", "stukel_asym")

generate_data_alt_link <- function(n, link_name = "probit", beta0 = 0, beta1 = 0.8) {
  x <- runif(n, min = -3, max = 3)
  eta <- beta0 + beta1 * x
  # Stukel inverse link helper
  inv_stukel <- function(eta_vec, a1, a2) {
    z <- numeric(length(eta_vec))
    pos <- eta_vec >= 0
    if (any(pos)) {
      if (abs(a1) < 1e-12) {
        z[pos] <- eta_vec[pos]
      } else {
        disc <- pmax(0, 1 + 2 * a1 * eta_vec[pos])
        z[pos] <- (-1 + sqrt(disc)) / a1
      }
    }
    if (any(!pos)) {
      if (abs(a2) < 1e-12) {
        z[!pos] <- eta_vec[!pos]
      } else {
        disc <- pmax(0, 1 + 2 * a2 * eta_vec[!pos])
        z[!pos] <- (-1 + sqrt(disc)) / a2
      }
    }
    plogis(z)
  }
  p <- switch(
    link_name,
    probit = pnorm(eta),
    cloglog = 1 - exp(-exp(eta)),
    stukel_long = inv_stukel(eta, a1 = -1.0, a2 = -1.0),
    stukel_short = inv_stukel(eta, a1 = 1.0, a2 = 1.0),
    stukel_asym = inv_stukel(eta, a1 = -1.0, a2 = 1.0),
    stop("Unsupported alt link: ", link_name)
  )
  y <- rbinom(n, size = 1, prob = p)
  data.frame(x = x, y = y)
}

# Fit misspecified model (omit quadratic): y ~ x
fit_misspecified <- function(dat) {
  glm(y ~ x, data = dat, family = binomial())
}



# Hosmer-Lemeshow using ResourceSelection only
hl_test <- function(y, p_hat, G = 10) {
  out <- ResourceSelection::hoslem.test(y, p_hat, g = G)
  list(stat = as.numeric(out$statistic), df = as.numeric(out$parameter), pval = as.numeric(out$p.value))
}

# Ebrahim-Farrington using ebrahim.gof only
ef_test <- function(y, p_hat, G = 10) {
  out <- ebrahim.gof::ef.gof(y = y, p = p_hat, G = G)
  pval <- if (!is.null(out$p.value)) as.numeric(out$p.value) else if (!is.null(out$pval)) as.numeric(out$pval) else if (!is.null(out$p)) as.numeric(out$p) else NA_real_
  if (is.na(pval)) stop("ef.gof did not return a recognizable p-value field.")
  list(stat = NA_real_, pval = pval)
}

# One trial given J
run_trial <- function(J, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  beta <- solve_quadratic_betas(J)
  dat <- generate_data(beta, n)
  fit <- fit_misspecified(dat)
  p_hat <- as.numeric(fitted(fit))
  hl <- hl_test(dat$y, p_hat, G = num_groups)
  ef <- ef_test(dat$y, p_hat, G = num_groups)
  # Additional tests: Pigeon-Heyse and HL equal-width
  p_ph <- try(pigeon_heyse_test(data.frame(y = dat$y), fit, g = num_groups)$p_value, silent = TRUE)
  p_hl_eqw <- try(hosmer_lemeshow_equal_intervals(p_hat, dat$y, num_groups = num_groups)$p.value, silent = TRUE)
  c(
    hl_sig = as.integer(hl$pval < alpha),
    ef_sig = as.integer(ef$pval < alpha),
    ph_sig = as.integer(!inherits(p_ph, "try-error") && !is.na(p_ph) && p_ph < alpha),
    hl_eqw_sig = as.integer(!inherits(p_hl_eqw, "try-error") && !is.na(p_hl_eqw) && p_hl_eqw < alpha)
  )
}

# One trial for interaction scenario given I
run_trial_interaction <- function(I, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  beta <- solve_interaction_betas(I)
  dat <- generate_data_interaction(beta, n)
  # Misspecified fit omits interaction term
  fit <- glm(y ~ x + d, data = dat, family = binomial())
  p_hat <- as.numeric(fitted(fit))
  hl <- hl_test(dat$y, p_hat, G = num_groups)
  ef <- ef_test(dat$y, p_hat, G = num_groups)
  p_ph <- try(pigeon_heyse_test(data.frame(y = dat$y), fit, g = num_groups)$p_value, silent = TRUE)
  p_hl_eqw <- try(hosmer_lemeshow_equal_intervals(p_hat, dat$y, num_groups = num_groups)$p.value, silent = TRUE)
  c(
    hl_sig = as.integer(hl$pval < alpha),
    ef_sig = as.integer(ef$pval < alpha),
    ph_sig = as.integer(!inherits(p_ph, "try-error") && !is.na(p_ph) && p_ph < alpha),
    hl_eqw_sig = as.integer(!inherits(p_hl_eqw, "try-error") && !is.na(p_hl_eqw) && p_hl_eqw < alpha)
  )
}

# One trial for continuous interaction scenario given K
run_trial_continuous_interaction <- function(K, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  beta <- solve_continuous_interaction_betas(K)
  dat <- generate_data_continuous_interaction(beta, n)
  # Misspecified fit omits interaction term
  fit <- glm(y ~ x + z, data = dat, family = binomial())
  p_hat <- as.numeric(fitted(fit))
  hl <- hl_test(dat$y, p_hat, G = num_groups)
  ef <- ef_test(dat$y, p_hat, G = num_groups)
  p_ph <- try(pigeon_heyse_test(data.frame(y = dat$y), fit, g = num_groups)$p_value, silent = TRUE)
  p_hl_eqw <- try(hosmer_lemeshow_equal_intervals(p_hat, dat$y, num_groups = num_groups)$p.value, silent = TRUE)
  c(
    hl_sig = as.integer(hl$pval < alpha),
    ef_sig = as.integer(ef$pval < alpha),
    ph_sig = as.integer(!inherits(p_ph, "try-error") && !is.na(p_ph) && p_ph < alpha),
    hl_eqw_sig = as.integer(!inherits(p_hl_eqw, "try-error") && !is.na(p_hl_eqw) && p_hl_eqw < alpha)
  )
}

# Run simulations in parallel over trials for a single J
simulate_for_J <- function(J, n, trials = num_trials) {
  cores <- max(1L, detectCores(logical = TRUE) - 1L)
  cl <- makeCluster(cores)
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  on.exit(try(stopCluster(cl), silent = TRUE))
  clusterSetRNGStream(cl, iseed = 20250911)
  n_current <- n
  clusterEvalQ(cl, {
    try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "pigeonheyse.R"))), silent = TRUE)
    try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "Hosmer (H) (equal width interval).R"))), silent = TRUE)
  })
  clusterExport(cl, varlist = c(
    "n_current", "num_groups", "alpha", "J", "logit", "inv_logit",
    "solve_quadratic_betas", "generate_data", "fit_misspecified"
    , "hl_test", "ef_test", "run_trial",
    "pigeon_heyse_test", "hosmer_lemeshow_equal_intervals"
  ), envir = environment())
  res <- parSapply(cl, 1:trials, function(tk) run_trial(J, n_current))
  data.frame(
    J = J,
    Test = c("HL", "EF", "Pigeon-Heyse", "HL_equal_width"),
    Power = c(
      mean(res["hl_sig", ]),
      mean(res["ef_sig", ]),
      mean(res["ph_sig", ]),
      mean(res["hl_eqw_sig", ])
    )
  )
}

# Parallel sims for interaction I
simulate_for_I <- function(I, n, trials = num_trials) {
  cores <- max(1L, detectCores(logical = TRUE) - 1L)
  cl <- makeCluster(cores)
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  on.exit(try(stopCluster(cl), silent = TRUE))
  clusterSetRNGStream(cl, iseed = 20250911)
  n_current <- n
  clusterEvalQ(cl, {
    try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "pigeonheyse.R"))), silent = TRUE)
    try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "Hosmer (H) (equal width interval).R"))), silent = TRUE)
  })
  clusterExport(cl, varlist = c(
    "n_current", "num_groups", "alpha", "I", "logit", "inv_logit",
    "solve_interaction_betas", "generate_data_interaction",
    "hl_test", "ef_test", "run_trial_interaction",
    "pigeon_heyse_test", "hosmer_lemeshow_equal_intervals"
  ), envir = environment())
  res <- parSapply(cl, 1:trials, function(tk) run_trial_interaction(I, n_current))
  data.frame(
    I = I,
    Test = c("HL", "EF", "Pigeon-Heyse", "HL_equal_width"),
    Power = c(
      mean(res["hl_sig", ]),
      mean(res["ef_sig", ]),
      mean(res["ph_sig", ]),
      mean(res["hl_eqw_sig", ])
    )
  )
}

# Parallel sims for continuous interaction K
simulate_for_K <- function(K, n, trials = num_trials) {
  cores <- max(1L, detectCores(logical = TRUE) - 1L)
  cl <- makeCluster(cores)
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  on.exit(try(stopCluster(cl), silent = TRUE))
  clusterSetRNGStream(cl, iseed = 20250911)
  n_current <- n
  clusterEvalQ(cl, {
    try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "pigeonheyse.R"))), silent = TRUE)
    try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "Hosmer (H) (equal width interval).R"))), silent = TRUE)
  })
  clusterExport(cl, varlist = c(
    "n_current", "num_groups", "alpha", "K", "logit", "inv_logit",
    "solve_continuous_interaction_betas", "generate_data_continuous_interaction",
    "hl_test", "ef_test", "run_trial_continuous_interaction",
    "pigeon_heyse_test", "hosmer_lemeshow_equal_intervals"
  ), envir = environment())
  res <- parSapply(cl, 1:trials, function(tk) run_trial_continuous_interaction(K, n_current))
  data.frame(
    K = K,
    Test = c("HL", "EF", "Pigeon-Heyse", "HL_equal_width"),
    Power = c(
      mean(res["hl_sig", ]),
      mean(res["ef_sig", ]),
      mean(res["ph_sig", ]),
      mean(res["hl_eqw_sig", ])
    )
  )
}

# Parallel sims for alternative link functions
simulate_for_alt_link <- function(link_name, n, trials = num_trials) {
  cores <- max(1L, detectCores(logical = TRUE) - 1L)
  cl <- makeCluster(cores)
  parallel::clusterExport(cl, c("EDGE_ARCHIVE_ROOT", "edge_path", "edge_battery", "edge_out"), envir = environment())  # archive: the path helpers, for the workers
  on.exit(try(stopCluster(cl), silent = TRUE))
  clusterSetRNGStream(cl, iseed = 20250911)
  n_current <- n
  ln <- link_name
  clusterEvalQ(cl, {
    try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "pigeonheyse.R"))), silent = TRUE)
    try(suppressWarnings(source(edge_path("code", "legacy_not_deposited", "Hosmer (H) (equal width interval).R"))), silent = TRUE)
  })
  clusterExport(cl, varlist = c(
    "n_current", "num_groups", "alpha", "ln",
    "generate_data_alt_link", "hl_test", "ef_test",
    "pigeon_heyse_test", "hosmer_lemeshow_equal_intervals"
  ), envir = environment())
  res <- parSapply(cl, 1:trials, function(tk) {
    dat <- generate_data_alt_link(n_current, link_name = ln)
    fit <- glm(y ~ x, data = dat, family = binomial(link = "logit"))
    p_hat <- as.numeric(fitted(fit))
    hl <- hl_test(dat$y, p_hat, G = num_groups)
    ef <- ef_test(dat$y, p_hat, G = num_groups)
    p_ph <- try(pigeon_heyse_test(data.frame(y = dat$y), fit, g = num_groups)$p_value, silent = TRUE)
    p_hl_eqw <- try(hosmer_lemeshow_equal_intervals(p_hat, dat$y, num_groups = num_groups)$p.value, silent = TRUE)
    c(
      hl_sig = as.integer(hl$pval < alpha),
      ef_sig = as.integer(ef$pval < alpha),
      ph_sig = as.integer(!inherits(p_ph, "try-error") && !is.na(p_ph) && p_ph < alpha),
      hl_eqw_sig = as.integer(!inherits(p_hl_eqw, "try-error") && !is.na(p_hl_eqw) && p_hl_eqw < alpha)
    )
  })
  data.frame(
    Link = link_name,
    Test = c("HL", "EF", "Pigeon-Heyse", "HL_equal_width"),
    Power = c(
      mean(res["hl_sig", ]),
      mean(res["ef_sig", ]),
      mean(res["ph_sig", ]),
      mean(res["hl_eqw_sig", ])
    )
  )
}

# Main execution with CSV saving and incremental running
run_all <- function() {
  all_master_results <- list()
  
  for (n in sample_sizes) {
    message("\n=== Processing sample size n = ", n, " ===")
    
    # --- 1. Omitted Quadratic Scenario ---
    message("\n--- Omitted Quadratic Scenario ---")
    existing_quad <- load_existing_results("omitted_quad", n)
    missing_J <- get_missing_scenarios(existing_quad, "omitted_quad", n)
    
    if (length(missing_J) > 0) {
      message("Running simulations for missing J values...")
      new_results_list <- lapply(missing_J, function(J) simulate_for_J(J, n = n))
      new_results <- do.call(rbind, new_results_list)
      new_results$J <- as.numeric(new_results$J)
      new_results$N <- n
      new_results$Scenario <- "omitted_quad"
      
      # Combine with existing results
      if (!is.null(existing_quad)) {
        existing_quad$Scenario <- "omitted_quad"
        results <- rbind(existing_quad, new_results)
      } else {
        results <- new_results
      }
      
      # Save updated results
      save_results(results, "omitted_quad", n)
    } else {
      results <- existing_quad
      results$Scenario <- "omitted_quad"
    }
    
    # Create and save plot
    p <- ggplot(results, aes(x = .data$J, y = .data$Power, color = .data$Test)) +
      geom_line() +
      geom_point() +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      scale_x_continuous(breaks = J_values) +
      labs(
        title = sprintf("Power vs J (Omitted Quadratic, n=%d, %d trials)", n, num_trials),
        x = "J (pi(-3))",
        y = "Power (proportion p < 0.05)",
        color = "Test"
      ) +
      theme_minimal(base_size = 12)
    
    out_plot_path <- plot_path_for_n(n)
    ggsave(out_plot_path, p, width = 8, height = 5, dpi = 150)
    message("Saved plot to: ", normalizePath(out_plot_path))
    
    all_master_results[[paste0("quad_", n)]] <- results
    
    # --- 2. Interaction Omitted Scenario ---
    message("\n--- Interaction Omitted Scenario ---")
    existing_int <- load_existing_results("interaction", n)
    missing_I <- get_missing_scenarios(existing_int, "interaction", n)
    
    if (length(missing_I) > 0) {
      message("Running simulations for missing I values...")
      i_results_list <- lapply(missing_I, function(I) simulate_for_I(I, n = n))
      i_new_results <- do.call(rbind, i_results_list)
      i_new_results$I <- as.numeric(i_new_results$I)
      i_new_results$N <- n
      i_new_results$Scenario <- "interaction"
      
      # Combine with existing results
      if (!is.null(existing_int)) {
        existing_int$Scenario <- "interaction"
        i_results <- rbind(existing_int, i_new_results)
      } else {
        i_results <- i_new_results
      }
      
      # Save updated results
      save_results(i_results, "interaction", n)
    } else {
      i_results <- existing_int
      i_results$Scenario <- "interaction"
    }
    
    # Create and save plot
    p2 <- ggplot(i_results, aes(x = .data$I, y = .data$Power, color = .data$Test)) +
      geom_line() +
      geom_point() +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      scale_x_continuous(breaks = I_values) +
      labs(
        title = sprintf("Power vs I (Interaction Omitted, n=%d, %d trials)", n, num_trials),
        x = "I (increment at x=3, d=1)",
        y = "Power (proportion p < 0.05)",
        color = "Test"
      ) +
      theme_minimal(base_size = 12)
    
    out_plot_path2 <- plot_path_for_n_interaction(n)
    ggsave(out_plot_path2, p2, width = 8, height = 5, dpi = 150)
    message("Saved plot to: ", normalizePath(out_plot_path2))
    
    all_master_results[[paste0("int_", n)]] <- i_results
    
    # --- 3. Continuous Interaction Omitted Scenario ---
    message("\n--- Continuous Interaction Omitted Scenario ---")
    existing_cont_int <- load_existing_results("continuous_interaction", n)
    missing_K <- get_missing_scenarios(existing_cont_int, "continuous_interaction", n)
    
    if (length(missing_K) > 0) {
      message("Running simulations for missing K values...")
      k_results_list <- lapply(missing_K, function(K) simulate_for_K(K, n = n))
      k_new_results <- do.call(rbind, k_results_list)
      k_new_results$K <- as.numeric(k_new_results$K)
      k_new_results$N <- n
      k_new_results$Scenario <- "continuous_interaction"
      
      # Combine with existing results
      if (!is.null(existing_cont_int)) {
        existing_cont_int$Scenario <- "continuous_interaction"
        k_results <- rbind(existing_cont_int, k_new_results)
      } else {
        k_results <- k_new_results
      }
      
      # Save updated results
      save_results(k_results, "continuous_interaction", n)
    } else {
      k_results <- existing_cont_int
      k_results$Scenario <- "continuous_interaction"
    }
    
    # Create and save plot
    p3 <- ggplot(k_results, aes(x = .data$K, y = .data$Power, color = .data$Test)) +
      geom_line() +
      geom_point() +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      scale_x_continuous(breaks = K_values) +
      labs(
        title = sprintf("Power vs K (Continuous Interaction Omitted, n=%d, %d trials)", n, num_trials),
        x = "K (continuous interaction effect)",
        y = "Power (proportion p < 0.05)",
        color = "Test"
      ) +
      theme_minimal(base_size = 12)
    
    out_plot_path3 <- plot_path_for_n_continuous_interaction(n)
    ggsave(out_plot_path3, p3, width = 8, height = 5, dpi = 150)
    message("Saved plot to: ", normalizePath(out_plot_path3))
    
    all_master_results[[paste0("cont_int_", n)]] <- k_results
    
    # --- 4. Alternative Link Functions Scenario ---
    message("\n--- Alternative Link Functions Scenario ---")
    existing_alt <- load_existing_results("alt_link", n)
    missing_links <- get_missing_scenarios(existing_alt, "alt_link", n)
    
    if (length(missing_links) > 0) {
      message("Running simulations for missing link functions...")
      l_results_list <- lapply(missing_links, function(lnk) simulate_for_alt_link(lnk, n = n))
      l_new_results <- do.call(rbind, l_results_list)
      l_new_results$N <- n
      l_new_results$Scenario <- "alt_link"
      
      # Combine with existing results
      if (!is.null(existing_alt)) {
        existing_alt$Scenario <- "alt_link"
        l_results <- rbind(existing_alt, l_new_results)
      } else {
        l_results <- l_new_results
      }
      
      # Save updated results
      save_results(l_results, "alt_link", n)
    } else {
      l_results <- existing_alt
      l_results$Scenario <- "alt_link"
    }
    
    # Create and save plot
    p3 <- ggplot(l_results, aes(x = .data$Link, y = .data$Power, group = .data$Test, color = .data$Test)) +
      geom_line(aes(group = .data$Test)) +
      geom_point(position = position_dodge(width = 0.2)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
      labs(
        title = sprintf("Power under Alternative Links (n=%d, %d trials)", n, num_trials),
        x = "Link Function",
        y = "Power (proportion p < 0.05)",
        color = "Test"
      ) +
      theme_minimal(base_size = 12)
    
    out_plot_path3 <- plot_path_for_n_altlink(n)
    ggsave(out_plot_path3, p3, width = 8, height = 5, dpi = 150)
    message("Saved plot to: ", normalizePath(out_plot_path3))
    
    all_master_results[[paste0("alt_", n)]] <- l_results
  }
  
  # Combine all results and save master file
  if (length(all_master_results) > 0) {
    # Ensure all data frames have the same column structure before combining
    standardize_columns <- function(df, scenario_type) {
      if (scenario_type == "omitted_quad" && !"J" %in% names(df)) df$J <- NA
      if (scenario_type == "interaction" && !"I" %in% names(df)) df$I <- NA
      if (scenario_type == "continuous_interaction" && !"K" %in% names(df)) df$K <- NA
      if (scenario_type == "alt_link" && !"Link" %in% names(df)) df$Link <- NA
      
      # Ensure all expected columns exist
      expected_cols <- c("Test", "Power", "N", "Scenario", "J", "I", "K", "Link")
      for (col in expected_cols) {
        if (!col %in% names(df)) df[[col]] <- NA
      }
      
      # Reorder columns consistently
      df[, expected_cols, drop = FALSE]
    }
    
    # Standardize each result set
    for (i in seq_along(all_master_results)) {
      scenario_type <- unique(all_master_results[[i]]$Scenario)[1]
      all_master_results[[i]] <- standardize_columns(all_master_results[[i]], scenario_type)
    }
    
    master_results <- do.call(rbind, all_master_results)
    master_csv_path <- csv_path_master()
    write.csv(master_results, master_csv_path, row.names = FALSE)
    message("\n=== Saved master results to: ", master_csv_path, " ===")
    
    return(master_results)
  } else {
    message("No results to save")
    return(NULL)
  }
}

# Convenience function to load all existing results without running new simulations
load_all_results <- function() {
  all_results <- list()
  
  for (n in sample_sizes) {
    # Load omitted quadratic results
    quad_results <- load_existing_results("omitted_quad", n)
    if (!is.null(quad_results)) {
      quad_results$Scenario <- "omitted_quad"
      all_results[[paste0("quad_", n)]] <- quad_results
    }
    
    # Load interaction results
    int_results <- load_existing_results("interaction", n)
    if (!is.null(int_results)) {
      int_results$Scenario <- "interaction"
      all_results[[paste0("int_", n)]] <- int_results
    }
    
    # Load continuous interaction results
    cont_int_results <- load_existing_results("continuous_interaction", n)
    if (!is.null(cont_int_results)) {
      cont_int_results$Scenario <- "continuous_interaction"
      all_results[[paste0("cont_int_", n)]] <- cont_int_results
    }
    
    # Load alternative link results
    alt_results <- load_existing_results("alt_link", n)
    if (!is.null(alt_results)) {
      alt_results$Scenario <- "alt_link"
      all_results[[paste0("alt_", n)]] <- alt_results
    }
  }
  
  if (length(all_results) > 0) {
    combined_results <- do.call(rbind, all_results)
    message("Loaded ", nrow(combined_results), " existing result rows from CSV files")
    return(combined_results)
  } else {
    message("No existing results found")
    return(NULL)
  }
}

# Function to show simulation status
show_simulation_status <- function() {
  cat("=== Simulation Status ===\n")
  
  for (n in sample_sizes) {
    cat("\nSample size n =", n, ":\n")
    
    # Check omitted quadratic
    existing_quad <- load_existing_results("omitted_quad", n)
    missing_J <- get_missing_scenarios(existing_quad, "omitted_quad", n)
    cat("  Omitted Quadratic:", 
        length(J_values) - length(missing_J), "/", length(J_values), "completed\n")
    
    # Check interaction
    existing_int <- load_existing_results("interaction", n)
    missing_I <- get_missing_scenarios(existing_int, "interaction", n)
    cat("  Interaction:", 
        length(I_values) - length(missing_I), "/", length(I_values), "completed\n")
    
    # Check continuous interaction
    existing_cont_int <- load_existing_results("continuous_interaction", n)
    missing_K <- get_missing_scenarios(existing_cont_int, "continuous_interaction", n)
    cat("  Continuous Interaction:", 
        length(K_values) - length(missing_K), "/", length(K_values), "completed\n")
    
    # Check alternative links
    existing_alt <- load_existing_results("alt_link", n)
    missing_links <- get_missing_scenarios(existing_alt, "alt_link", n)
    cat("  Alternative Links:", 
        length(alt_links) - length(missing_links), "/", length(alt_links), "completed\n")
  }
  
  # Check master file
  master_path <- csv_path_master()
  if (file.exists(master_path)) {
    cat("\nMaster results file exists:", master_path, "\n")
  } else {
    cat("\nMaster results file not found\n")
  }
}

# Professional A4 Summary Plot Function
create_summary_plot <- function(master_results = NULL) {
  # Load results if not provided
  if (is.null(master_results)) {
    master_results <- load_all_results()
    if (is.null(master_results)) {
      stop("No results available. Please run simulations first.")
    }
  }
  
  # Filter data based on selected tests
  selected_tests <- names(include_tests)[unlist(include_tests)]
  plot_data <- master_results[master_results$Test %in% selected_tests, ]
  
  if (nrow(plot_data) == 0) {
    stop("No tests selected! Please set at least one test to TRUE in include_tests.")
  }
  
  # Create scenario labels for better presentation
  plot_data$Scenario_Label <- factor(plot_data$Scenario, 
    levels = c("omitted_quad", "interaction", "continuous_interaction", "alt_link"),
    labels = c("Omitted Quadratic\n(Nonlinearity)", 
               "Omitted Interaction\n(Binary Effect Modification)",
               "Omitted Continuous Interaction\n(Continuous Effect Modification)",
               "Alternative Links\n(Link Misspecification)"))
  
  # Create sample size labels
  plot_data$N_Label <- factor(paste0("n = ", plot_data$N), 
    levels = paste0("n = ", sort(unique(plot_data$N))))
  
  # Create test labels with better formatting (only for selected tests)
  test_label_mapping <- c("HL" = "Hosmer-Lemeshow", 
                         "EF" = "Ebrahim-Farrington", 
                         "Pigeon-Heyse" = "Pigeon-Heyse", 
                         "HL_equal_width" = "HL Equal-Width")
  
  plot_data$Test_Label <- factor(plot_data$Test,
    levels = selected_tests,
    labels = test_label_mapping[selected_tests])
  
  # Reorder data to ensure Ebrahim-Farrington is plotted last (on top)
  plot_data <- plot_data[order(plot_data$Test != "EF"), ]
  
  # Create custom color palette with Ebrahim-Farrington as red
  create_custom_colors <- function(selected_tests) {
    # Base colors for tests (excluding red for EF)
    base_colors <- c("blue", "darkgreen", "black", "purple", "orange", "brown", "pink")
    
    colors <- character(length(selected_tests))
    names(colors) <- test_label_mapping[selected_tests]
    
    # Assign red to Ebrahim-Farrington
    if ("EF" %in% selected_tests) {
      colors["Ebrahim-Farrington"] <- "red"
    }
    
    # Assign other colors to remaining tests
    other_tests <- names(colors)[names(colors) != "Ebrahim-Farrington"]
    if (length(other_tests) > 0) {
      colors[other_tests] <- base_colors[1:length(other_tests)]
    }
    
    return(colors)
  }
  
  # Get the custom color palette
  custom_colors <- create_custom_colors(selected_tests)
  
  # Prepare scenario-specific data for plotting
  create_scenario_plot <- function(scenario_name, scenario_label) {
    data_subset <- plot_data[plot_data$Scenario == scenario_name, ]
    
    if (scenario_name == "omitted_quad") {
      # For omitted quadratic, plot against J values
      # Split data for layering (EF on top)
      data_other <- data_subset[data_subset$Test != "EF", ]
      data_ef <- data_subset[data_subset$Test == "EF", ]
      
      p <- ggplot(data_subset, aes(x = J, y = Power, color = Test_Label)) +
        # Plot other tests first
        geom_line(data = data_other, size = 0.8) +
        geom_point(data = data_other, size = 1.5) +
        # Plot EF on top with thicker line
        geom_line(data = data_ef, size = 1.2) +
        geom_point(data = data_ef, size = 2) +
        scale_x_continuous(breaks = J_values, labels = sprintf("%.2f", J_values)) +
        scale_color_manual(values = custom_colors) +
        labs(x = "J (π(-3))", title = scenario_label) +
        theme_minimal(base_size = 8) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 7),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "gray80")
        )
    } else if (scenario_name == "interaction") {
      # For interaction, plot against I values
      # Split data for layering (EF on top)
      data_other <- data_subset[data_subset$Test != "EF", ]
      data_ef <- data_subset[data_subset$Test == "EF", ]
      
      p <- ggplot(data_subset, aes(x = I, y = Power, color = Test_Label)) +
        # Plot other tests first
        geom_line(data = data_other, size = 0.8) +
        geom_point(data = data_other, size = 1.5) +
        # Plot EF on top with thicker line
        geom_line(data = data_ef, size = 1.2) +
        geom_point(data = data_ef, size = 2) +
        scale_x_continuous(breaks = I_values) +
        scale_color_manual(values = custom_colors) +
        labs(x = "I (Binary Interaction Effect)", title = scenario_label) +
        theme_minimal(base_size = 8) +
        theme(
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 7),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "gray80")
        )
    } else if (scenario_name == "continuous_interaction") {
      # For continuous interaction, plot against K values
      # Split data for layering (EF on top)
      data_other <- data_subset[data_subset$Test != "EF", ]
      data_ef <- data_subset[data_subset$Test == "EF", ]
      
      p <- ggplot(data_subset, aes(x = K, y = Power, color = Test_Label)) +
        # Plot other tests first
        geom_line(data = data_other, size = 0.8) +
        geom_point(data = data_other, size = 1.5) +
        # Plot EF on top with thicker line
        geom_line(data = data_ef, size = 1.2) +
        geom_point(data = data_ef, size = 2) +
        scale_x_continuous(breaks = K_values) +
        scale_color_manual(values = custom_colors) +
        labs(x = "K (Continuous Interaction Effect)", title = scenario_label) +
        theme_minimal(base_size = 8) +
        theme(
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 7),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "gray80")
        )
    } else if (scenario_name == "alt_link") {
      # For alternative links, plot as bar chart
      p <- ggplot(data_subset, aes(x = Link, y = Power, fill = Test_Label)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        scale_fill_manual(values = custom_colors) +
        labs(x = "Link Function", title = scenario_label) +
        theme_minimal(base_size = 8) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          axis.title = element_text(size = 7),
          plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "gray80")
        )
    }
    
    # Common formatting
    p <- p +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), 
                        labels = sprintf("%.1f", seq(0, 1, 0.2))) +
      labs(y = "Power") +
      facet_wrap(~ N_Label, nrow = 1, scales = "free_x")
    
    return(p)
  }
  
  # Create individual plots for each scenario
  p1 <- create_scenario_plot("omitted_quad", "Omitted Quadratic (Nonlinearity)")
  p2 <- create_scenario_plot("interaction", "Omitted Binary Interaction (Binary Effect Modification)")
  p3 <- create_scenario_plot("continuous_interaction", "Omitted Continuous Interaction (Continuous Effect Modification)")
  p4 <- create_scenario_plot("alt_link", "Alternative Links (Link Misspecification)")
  
  # Create a shared legend based on selected tests
  num_selected_tests <- length(selected_tests)
  if (num_selected_tests > 0) {
    # Create dummy data for legend with only selected tests
    legend_data <- data.frame(
      Test_Label = factor(test_label_mapping[selected_tests], levels = test_label_mapping[selected_tests]),
      x = rep(1, num_selected_tests),
      y = rep(1, num_selected_tests)
    )
  } else {
    stop("No tests selected for legend creation!")
  }
  
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Test_Label)) +
    geom_point(size = 3) +
    scale_color_manual(values = custom_colors) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.5, "cm")
    ) +
    guides(color = guide_legend(nrow = 1))
  
  legend <- get_legend(legend_plot)
  
  # Combine all plots
  
  # Create title
  title <- textGrob("Goodness-of-Fit Test Power Comparison Across Scenarios", 
                   gp = gpar(fontsize = 16, fontface = "bold"))
  
  # Create subtitle with selected tests information
  selected_test_names <- test_label_mapping[selected_tests]
  test_info <- paste(selected_test_names, collapse = ", ")
  subtitle_text <- sprintf("Monte Carlo Simulation Results (%s trials per scenario)\nIncluded Tests: %s", 
                          format(num_trials, big.mark = ","), test_info)
  subtitle <- textGrob(subtitle_text, 
                      gp = gpar(fontsize = 10, fontface = "italic"))
  
  # Arrange plots in grid
  combined_plot <- plot_grid(
    plot_grid(title, subtitle, ncol = 1, rel_heights = c(1, 0.5)),
    p1, p2, p3, p4,
    legend,
    ncol = 1, 
    rel_heights = c(0.8, 3, 3, 3, 3, 0.8)
  )
  
  # Save as high-quality A4 PDF
  output_path <- file.path(get_base_dir(), "Professional_Summary_A4.pdf")
  ggsave(output_path, combined_plot, 
         width = 8.27, height = 11.69,  # A4 dimensions in inches
         units = "in", dpi = 300, device = "pdf")
  
  # Also save as high-quality PNG
  output_path_png <- file.path(get_base_dir(), "Professional_Summary_A4.png")
  ggsave(output_path_png, combined_plot, 
         width = 8.27, height = 11.69,  # A4 dimensions in inches
         units = "in", dpi = 300, device = "png")
  
  message("Professional A4 summary plot saved to:")
  message("  PDF: ", output_path)
  message("  PNG: ", output_path_png)
  
  return(combined_plot)
}

# Create custom A4 plot with specific test selection (without changing global settings)
create_custom_summary_plot <- function(master_results = NULL, 
                                      custom_tests = list("HL" = TRUE, "EF" = TRUE, "Pigeon-Heyse" = TRUE, "HL_equal_width" = TRUE),
                                      output_filename = "Custom_Summary_A4") {
  # Temporarily store original settings
  original_tests <- include_tests
  
  # Set custom test selection
  include_tests <<- custom_tests
  
  # Create the plot
  tryCatch({
    plot <- create_summary_plot(master_results)
    
    # Save with custom filename
    output_path_pdf <- file.path(get_base_dir(), paste0(output_filename, ".pdf"))
    output_path_png <- file.path(get_base_dir(), paste0(output_filename, ".png"))
    
    ggsave(output_path_pdf, plot, width = 8.27, height = 11.69, units = "in", dpi = 300, device = "pdf")
    ggsave(output_path_png, plot, width = 8.27, height = 11.69, units = "in", dpi = 300, device = "png")
    
    message("Custom A4 summary plot saved to:")
    message("  PDF: ", output_path_pdf)
    message("  PNG: ", output_path_png)
    
    return(plot)
  }, finally = {
    # Restore original settings
    include_tests <<- original_tests
  })
}

# Helper function to extract legend from ggplot
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

if (identical(environment(), globalenv())) {
  # Show current status first
  show_simulation_status()
  
  # Run simulations (only missing ones will be executed)
  cat("\n=== Starting simulation run ===\n")
  out <- run_all()
  
  # Show summary
  cat("\n=== Final Results Summary ===\n")
  if (!is.null(out)) {
    cat("Total result rows:", nrow(out), "\n")
    cat("Scenarios:", paste(unique(out$Scenario), collapse = ", "), "\n")
    cat("Sample sizes:", paste(unique(out$N), collapse = ", "), "\n")
    cat("Tests:", paste(unique(out$Test), collapse = ", "), "\n")
    
    # Create professional summary plot
    cat("\n=== Creating Professional A4 Summary Plot ===\n")
    tryCatch({
      create_summary_plot(out)
      cat("Professional summary plot created successfully!\n")
    }, error = function(e) {
      cat("Error creating summary plot:", e$message, "\n")
      cat("Individual plots are still available in the working directory.\n")
    })
  }
}


