newtraph <- function(y, X, t, delta, beta, gamma, lklhd, verbose = FALSE, EPS = 1e-4, iter_max = 50,
                     TAU = 3, K = 3) {
  ## Evaluate likelihood for starting values of theta
  # errs <- y - X * beta
  init_res <- lklhd(beta, gamma, X, y, t, delta, K = K)

  ## initialize values
  converge <- TRUE
  delta_vec <- rep(1, length(c(beta, gamma)))
  iters <- 0
  theta <- c(beta, gamma)
  theta_new <- theta
  new_res <- init_res
  err_occurred <- FALSE
  step_halving_start <- 0
  should_end <- FALSE
  p <- length(beta)
  p_gamma <- length(gamma)

  while(sum(abs(delta_vec) > EPS) >= 1){
    iters <- iters + 1
    # print(iters)
    theta_old <- theta_new
    lnlklhd_old <- new_res$lnlklhd
    inform_reg <- new_res$inform # + lambda * diag(nrow(new_res$inform))

    delta_vec <- tryCatch({
      solve(inform_reg) %*% new_res$score
    },
    error = function(e) {
      print(e)
      err_occurred <<- TRUE
      rep(0, length(c(beta, gamma)))
    })
    # print(delta_vec)
    theta_new <- theta_old + delta_vec
    # errs <- y - X * theta_new[1:length(beta)]
    new_res <- tryCatch({
      lklhd(theta_new[1:p], theta_new[(p+1):(p+p_gamma)], X, y, t, delta, K = K)
    },
    error = function(e) {
      err_occurred <<- TRUE
      print(e)
      new_res
    })
    lnlklhd_new <- new_res$lnlklhd
    step_size <- 1
    while (is.nan(lnlklhd_new) | lnlklhd_new - lnlklhd_old <= -EPS) {
      # print("Step Halving!")
      if (verbose) cat("Performing Step Halving (new log-likelihood", round(lnlklhd_new, 6), ")\n")
      step_size <- step_size / 2
      theta_new <- theta_old + delta_vec * step_size
      ## if (abs(step_halving_start - lnlklhd_new) <= EPS) break
      # errs <- y - X * theta_new[1:length(beta)]
      new_res <- lklhd(theta_new[1:p], theta_new[(p+1):(p+p_gamma)], X, y, t, delta, K = K)
      lnlklhd_new <- new_res$lnlklhd
    }

    if (abs(step_halving_start - lnlklhd_new) <= EPS) {
      warning("Stopping after step halving failed!")
      converge <- FALSE
      break
    }

    step_halving_start <- lnlklhd_new

    if (verbose == TRUE) {
      print(paste0("Iteration: ", iters, "; Theta Vector: ", stringr::str_c(theta_new, collapse = ", ")))
      print("Log Likelihood:")
      print(new_res$lnlklhd)
      print("Score Vector: ")
      print(new_res$score)
      print("Fisher's Information Matrix:")
      print(new_res$inform)
    }
    if (iters >= iter_max) {
      warning("Stopping after reaching the maximum number of iterations!")
      converge <- FALSE
      break
    }
  }
  if (err_occurred) converge <- FALSE

  return(
    list(
      theta0 = theta,
      lnlklhd0 = init_res$lnlklhd,
      converge = converge,
      nbriter = iters,
      theta = as.vector(theta_new),
      lnlklhd = new_res$lnlklhd,
      score = new_res$score,
      inform = new_res$inform,
      err_occurred = err_occurred,
      knots = new_res$basis$knots,
      boundary_knots = new_res$basis$boundary_knots
    )
  )
}
