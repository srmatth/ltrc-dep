# 2024-03-26_beta_1.rds

```
N <- 1000
true_beta <- c(1)
n <- 300
final_res <- list()

for (i in 1:N) {
  cat("Working on iteration ", i, "\n")
  x <- matrix(runif(length(true_beta)*n, -3, 3), ncol = length(true_beta))
  eps <- rnorm(n)
  y <- x %*% true_beta + eps

  ## Induce left-truncation into our data
  # t_b <- runif(n, -1, 0) # -abs(rnorm(n, -0.5))
  # t <- x %*% true_beta + t_b
  t <- runif(n, -6, 0)
  obs_indx <- which(y > t)

  # c_b <- runif(n, 2, 5) # abs(rnorm(n, 1.5))
  # c <- x %*% true_beta + c_b
  c <- runif(n, 0, 6)
  delta <- as.numeric(y < c)
  y_cens <- ifelse(y < c, y, c)

  y_obs <- y_cens[obs_indx]
  x_obs <- x[obs_indx,]
  t_obs <- t[obs_indx]
  delta_obs <- delta[obs_indx]

  lin_mod <- lm(y_obs ~ x_obs)
  beta <- lin_mod$coefficients[-1]
  gamma <- rnorm(6)

  dat <- data.frame(
    y_obs, delta_obs, x_obs
  )

  tictoc::tic()
  res <- ltrc(formula = survival::Surv(y_obs, delta_obs) ~ x_obs, t_obs, data = dat, n_start = 20, int_knots = 10)
  time <- tictoc::toc()

  tmp_dat <- list(
    n_obs_start = n,
    n_obs = length(y_obs),
    pct_truncated = 1 - length(y_obs) / n,
    pct_censored = 1 - mean(delta_obs),
    true_beta = true_beta,
    naive_beta = lin_mod$coefficients[2:(length(true_beta) + 1)],
    time = time$toc - time$tic,
    lnlklhd_start = res$lnlklhd0,
    beta_start = res$theta0[1:length(true_beta)],
    converge = res$converge,
    errored = res$err_occurred,
    n_iter = res$nbriter,
    lnlklhd = res$lnlklhd,
    beta = res$theta[1:length(true_beta)],
    gamma = res$theta[-(1:length(true_beta))],
    score = res$score,
    basis_boundary_knots = res$boundary_knots,
    basis_internal_knots = res$knots
  )

  final_res[[i]] <- tmp_dat
}
```

# 2024-03-27_beta_2_neg1.rds

```
N <- 1000
true_beta <- c(2, -1)
n <- 300
final_res <- list()

for (i in 1:N) {
  cat("Working on iteration ", i, "\n")
  x <- matrix(runif(length(true_beta)*n, -3, 3), ncol = length(true_beta))
  eps <- rnorm(n)
  y <- x %*% true_beta + eps
  
  ## Induce left-truncation into our data
  # t_b <- runif(n, -1, 0) # -abs(rnorm(n, -0.5))
  # t <- x %*% true_beta + t_b
  t <- runif(n, -6, 0)
  obs_indx <- which(y > t)
  
  # c_b <- runif(n, 2, 5) # abs(rnorm(n, 1.5))
  # c <- x %*% true_beta + c_b
  c <- runif(n, 0, 6)
  delta <- as.numeric(y < c)
  y_cens <- ifelse(y < c, y, c)
  
  y_obs <- y_cens[obs_indx]
  x_obs <- x[obs_indx,]
  t_obs <- t[obs_indx]
  delta_obs <- delta[obs_indx]
  
  lin_mod <- lm(y_obs ~ x_obs)
  beta <- lin_mod$coefficients[-1]
  gamma <- rnorm(6)
  
  dat <- data.frame(
    y_obs, delta_obs, x_obs
  )
  
  tictoc::tic()
  res <- ltrc(formula = survival::Surv(y_obs, delta_obs) ~ x_obs, t_obs, data = dat, n_start = 20, int_knots = 10)
  time <- tictoc::toc()
  
  tmp_dat <- list(
    n_obs_start = n,
    n_obs = length(y_obs),
    pct_truncated = 1 - length(y_obs) / n,
    pct_censored = 1 - mean(delta_obs),
    true_beta = true_beta,
    naive_beta = lin_mod$coefficients[2:(length(true_beta) + 1)],
    time = time$toc - time$tic,
    lnlklhd_start = res$lnlklhd0,
    beta_start = res$theta0[1:length(true_beta)],
    converge = res$converge,
    errored = res$err_occurred,
    n_iter = res$nbriter,
    lnlklhd = res$lnlklhd,
    beta = res$theta[1:length(true_beta)],
    gamma = res$theta[-(1:length(true_beta))],
    score = res$score,
    inform = res$inform,
    basis_boundary_knots = res$boundary_knots,
    basis_internal_knots = res$knots
  )
  
  final_res[[i]] <- tmp_dat
}
```

# 2024-03-27_beta_1_ev_error.rds

```
N <- 1000
true_beta <- c(1)
n <- 300
final_res <- list()
error_dist <- function(n) {
  coin_flip <- runif(n)
  ep1 <- rnorm(n)
  ep2 <- rnorm(n, 0, 3)
  pmax(ifelse(coin_flip > 0.5, ep1, ep2), -5)
}

for (i in 1:N) {
  cat("Working on iteration ", i, "\n")
  x <- matrix(runif(length(true_beta)*n, -3, 3), ncol = length(true_beta))
  eps <- error_dist(n)
  y <- x %*% true_beta + eps
  
  ## Induce left-truncation into our data
  # t_b <- runif(n, -1, 0) # -abs(rnorm(n, -0.5))
  # t <- x %*% true_beta + t_b
  t <- runif(n, -8, 0)
  obs_indx <- which(y > t)
  
  # c_b <- runif(n, 2, 5) # abs(rnorm(n, 1.5))
  # c <- x %*% true_beta + c_b
  c <- runif(n, 0, 6)
  delta <- as.numeric(y < c)
  y_cens <- ifelse(y < c, y, c)
  
  y_obs <- y_cens[obs_indx]
  x_obs <- x[obs_indx,]
  t_obs <- t[obs_indx]
  delta_obs <- delta[obs_indx]
  
  lin_mod <- lm(y_obs ~ x_obs)
  beta <- lin_mod$coefficients[-1]
  gamma <- rnorm(6)
  
  dat <- data.frame(
    y_obs, delta_obs, x_obs
  )
  
  tictoc::tic()
  res <- ltrc(formula = survival::Surv(y_obs, delta_obs) ~ x_obs, t_obs, data = dat, n_start = 20, int_knots = 10)
  time <- tictoc::toc()
  
  tmp_dat <- list(
    n_obs_start = n,
    n_obs = length(y_obs),
    pct_truncated = 1 - length(y_obs) / n,
    pct_censored = 1 - mean(delta_obs),
    true_beta = true_beta,
    naive_beta = lin_mod$coefficients[2:(length(true_beta) + 1)],
    time = time$toc - time$tic,
    lnlklhd_start = res$lnlklhd0,
    beta_start = res$theta0[1:length(true_beta)],
    converge = res$converge,
    errored = res$err_occurred,
    n_iter = res$nbriter,
    lnlklhd = res$lnlklhd,
    beta = res$theta[1:length(true_beta)],
    gamma = res$theta[-(1:length(true_beta))],
    score = res$score,
    inform = res$inform,
    basis_boundary_knots = res$boundary_knots,
    basis_internal_knots = res$knots
  )
  
  final_res[[i]] <- tmp_dat
}
```

# 2024-03-28_beta_1.rds

```
N <- 500
true_beta <- c(1)
n <- 500
final_res <- list()
error_dist <- function(n) {
  # coin_flip <- runif(n)
  # ep1 <- rnorm(n)
  # ep2 <- rnorm(n, 0, 3)
  # ifelse(coin_flip > 0.5, ep1, ep2)
  rnorm(n)
}

for (i in 1:N) {
  cat("Working on iteration ", i, "\n")
  x <- matrix(runif(length(true_beta)*n, -3, 3), ncol = length(true_beta))
  eps <- error_dist(n)
  y <- x %*% true_beta + eps

  ## Induce left-truncation into our data
  # t_b <- runif(n, -1, 0) # -abs(rnorm(n, -0.5))
  # t <- x %*% true_beta + t_b
  t <- runif(n, -6, 0)
  obs_indx <- which(y > t)

  # c_b <- runif(n, 2, 5) # abs(rnorm(n, 1.5))
  # c <- x %*% true_beta + c_b
  c <- runif(n, 0, 6)
  delta <- as.numeric(y < c)
  y_cens <- ifelse(y < c, y, c)

  y_obs <- y_cens[obs_indx]
  x_obs <- x[obs_indx,]
  t_obs <- t[obs_indx]
  delta_obs <- delta[obs_indx]

  lin_mod <- lm(y_obs ~ x_obs)
  beta <- lin_mod$coefficients[-1]
  gamma <- rnorm(6)

  dat <- data.frame(
    y_obs, delta_obs, x_obs
  )

  tictoc::tic()
  res <- ltrc(formula = survival::Surv(y_obs, delta_obs) ~ x_obs, t_obs, data = dat, n_start = 20, int_knots = 5)
  time <- tictoc::toc()

  rc_res <- ltrc(formula = survival::Surv(y_obs, delta_obs) ~ x_obs, rep(-10, length(delta_obs)), data = dat, n_start = 20, int_knots = 5)

  tmp_dat <- list(
    n_obs_start = n,
    n_obs = length(y_obs),
    pct_truncated = 1 - length(y_obs) / n,
    pct_censored = 1 - mean(delta_obs),
    true_beta = true_beta,
    naive_beta = lin_mod$coefficients[2:(length(true_beta) + 1)],
    rc_beta = rc_res$theta[1:length(true_beta)],
    rc_score = rc_res$score,
    time = time$toc - time$tic,
    lnlklhd_start = res$lnlklhd0,
    beta_start = res$theta0[1:length(true_beta)],
    converge = res$converge,
    errored = res$err_occurred,
    n_iter = res$nbriter,
    lnlklhd = res$lnlklhd,
    beta = res$theta[1:length(true_beta)],
    gamma = res$theta[-(1:length(true_beta))],
    score = res$score,
    inform = res$inform,
    basis_boundary_knots = res$boundary_knots,
    basis_internal_knots = res$knots
  )

  final_res[[i]] <- tmp_dat
}

readr::write_rds(final_res, "inst/sim_res/2024-03-28_beta_1.rds")
```

