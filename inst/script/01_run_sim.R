devtools::load_all()
library(extRemes)

N_CAL <- 10
N <- 510
true_beta <- c(1)
n <- 500
final_res <- list()
error_dist <- function(n) {
  # coin_flip <- runif(n)
  # ep1 <- rnorm(n)
  # ep2 <- rnorm(n, 0, 3)
  # err_vec <- ifelse(coin_flip > 0.5, ep1, ep2)
  # err_vec <- runif(n, -3, 3)
  # err_vec <- rnorm(n)
  err_vec <- revd(n)
  return(err_vec)
}
cal_knots <- numeric(N_CAL)
FINAL_KNOTS <- 5

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

  if (i <= N_CAL) {
    cat("Calibrating the number of knots using Cross Validation: ", i, " of ", N_CAL, "\n")
    best_knots <- 5
    best_lklhd <- -Inf
    tictoc::tic()
    res <- ltrc(formula = survival::Surv(y_obs, delta_obs) ~ x_obs, t_obs, data = dat, n_start = 5, int_knots = 0)
    time <- tictoc::toc()
    cal_knots[i] <- length(res$knots)
    if (i == N_CAL) {
      FINAL_KNOTS <- round(mean(cal_knots))
      cat("Finished Cross-Validation: FINAL_KNOTS = ", FINAL_KNOTS, "\n")
    }
  } else {
    tictoc::tic()
    res <- ltrc(formula = survival::Surv(y_obs, delta_obs) ~ x_obs, t_obs, data = dat, n_start = 20, int_knots = FINAL_KNOTS)
    time <- tictoc::toc()

    rc_res <- ltrc(formula = survival::Surv(y_obs, delta_obs) ~ x_obs, rep(-10, length(delta_obs)), data = dat, n_start = 20, int_knots = FINAL_KNOTS)

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

    final_res[[i-N_CAL]] <- tmp_dat
  }
}

readr::write_rds(final_res, "inst/sim_res/2024-04-03_beta_1_ev_error.rds")

test_res <- readRDS("inst/sim_res/2024-03-25_beta_1.rds")
N <- length(test_res)

betas <- numeric(N)
naive_betas <- numeric(N)
gammas <- matrix(NA, nrow = N, ncol = 9)
for (i in 1:(N-10)) {
  naive_betas[i] <- final_res[[i]]$naive_beta
  betas[i] <- final_res[[i]]$beta
  gammas[i,] <- final_res[[i]]$gamma
  # print(i$score)
  # print(i$pct_truncated)
  # print(i$pct_censored)
}

plot(density(betas))
plot(density(gammas[,1]))
plot(density(naive_betas))

betas <- matrix(NA, nrow = N, ncol = 2)
naive_betas <- matrix(NA, nrow = N, ncol = 2)
gammas <- matrix(NA, nrow = N, ncol = 13)
for (i in 1:N) {
  naive_betas[i,] <- final_res[[i]]$naive_beta
  betas[i,] <- final_res[[i]]$beta
  gammas[i,] <- final_res[[i]]$gamma
  # print(i$score)
  # print(i$pct_truncated)
  # print(i$pct_censored)
}

plot(density(betas[,1]))
plot(density(betas[,2]))
plot(density(gammas[,1]))
plot(density(naive_betas))

mean(betas[,1])
mean(betas[,2])
