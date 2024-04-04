sim_res_to_df <- function(sim_res) {
  res <- data.frame()
  p <- length(sim_res[[1]]$beta)
  if (p == 1) {
    for (i in sim_res) {

      tmp <- data.frame(
        n_obs_start = i$n_obs_start,
        n_obs = i$n_obs,
        pct_truncated = i$pct_truncated,
        pct_censored = i$pct_censored,
        time = i$time,
        lnlklhd_start = i$lnlklhd_start,
        lnlklhd_end = i$lnlklhd,
        n_iter = i$n_iter,
        p_beta = length(i$beta),
        p_gamma = length(i$gamma),
        norm_score = sqrt(sum(i$score * i$score)),
        norm_rc_score = sqrt(sum(i$rc_score * i$rc_score)),
        true_beta = i$true_beta,
        naive_beta = i$naive_beta,
        rc_beta = i$rc_beta,
        beta_start = i$beta_start,
        beta = i$beta,
        beta_score = i$score[1],
        beta_var = solve(i$inform)[1,1]
      )
      res <- rbind(res, tmp)
    }
  } else if (p == 2) {
    for (i in sim_res) {

      tmp <- data.frame(
        n_obs_start = i$n_obs_start,
        n_obs = i$n_obs,
        pct_truncated = i$pct_truncated,
        pct_censored = i$pct_censored,
        time = i$time,
        lnlklhd_start = i$lnlklhd_start,
        lnlklhd_end = i$lnlklhd,
        n_iter = i$n_iter,
        p_beta = length(i$beta),
        p_gamma = length(i$gamma),
        norm_score = sqrt(sum(i$score * i$score)),
        true_beta_1 = i$true_beta[1],
        true_beta_2 = i$true_beta[2],
        naive_beta_1 = i$naive_beta[1],
        naive_beta_2 = i$naive_beta[2],
        beta_start_1 = i$beta_start[1],
        beta_start_2 = i$beta_start[2],
        beta_1 = i$beta[1],
        beta_2 = i$beta[2]
      )
      res <- rbind(res, tmp)
    }
  } else {
    print("Not yet implemented for multi-dimensional beta")
    return(NULL)
  }
  return(res)
}

predict_hazard <- function(ltrc_obj, vals, within_bounds = FALSE) {
  if (within_bounds) {
    x <- vals[vals < ltrc_obj$basis_boundary_knots[2]]
    res <- exp(predict_basis(x, knots = ltrc_obj$basis_internal_knots,
                      boundary_knots = ltrc_obj$basis_boundary_knots) %*% ltrc_obj$gamma)
    return(c(res, rep(NA, length(vals) - length(x))))
  }
  exp(predict_basis(vals, knots = ltrc_obj$basis_internal_knots,
                    boundary_knots = ltrc_obj$basis_boundary_knots) %*% ltrc_obj$gamma)
}

true_hazard <- function(x) {
  # PDF of the standard normal distribution
  pdf <- dnorm(x)

  # Survival function of the standard normal distribution
  survival_function <- 1 - pnorm(x)

  # Hazard function
  hazard <- pdf / survival_function

  return(hazard)
}
