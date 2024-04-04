lnlklhd <- function(beta, gamma, X, y, t, delta, K = 3, knots = NULL, boundary_knots = NULL) {
  ## Error checking and argument massaging
  if (!is.matrix(X)) {
    X <- matrix(X, ncol = 1)
  }

  ## Getting some important values
  N <- nrow(X)
  p <- ncol(X)
  p_gamma <- length(gamma)
  e_b <- y -  X %*% beta
  t_b <- t - X %*% beta
  tracker <- data.frame(
    is_e = c(rep(0, N), rep(1, N)),
    order = c(1:N, 1:N),
    val = c(t_b, e_b)
  ) %>%
    dplyr::arrange(val) %>%
    dplyr::mutate(
      to_next = dplyr::lead(val, default = max(val)) - val,
      to_prior = val - dplyr::lag(val, default = min(val)),
      width = (to_next + to_prior) / 2
    )
  if ((!is.null(knots)) & (!is.null(boundary_knots))) {
    K <- length(knots)
    basis <- list(
      basis = predict_basis(tracker$val, knots, boundary_knots),
      d_1 = predict_basis(tracker$val, knots, boundary_knots, deriv = 1),
      d_2 = predict_basis(tracker$val, knots, boundary_knots, deriv = 2),
      knots = knots,
      boundary_knots = boundary_knots
    )
  } else {
    basis <- get_basis(c(min(t_b), e_b), K = K)
    new_basis <- list(
      basis = predict_basis(tracker$val, basis$knots, basis$boundary_knots),
      d_1 = predict_basis(tracker$val, basis$knots, basis$boundary_knots, deriv = 1),
      d_2 = predict_basis(tracker$val, basis$knots, basis$boundary_knots, deriv = 2),
      knots = basis$knots,
      boundary_knots = basis$boundary_knots
    )
    basis <- new_basis
  }
  e_indx <- which(tracker$is_e == 1)
  e_order <- tracker$order[which(tracker$is_e == 1)]
  e_indx <- data.frame(e_indx, e_order) %>%
    dplyr::arrange(e_order) %>%
    dplyr::pull(e_indx)
  t_indx <- which(tracker$is_e == 0)
  t_order <- tracker$order[which(tracker$is_e == 0)]
  t_indx <- data.frame(t_indx, t_order) %>%
    dplyr::arrange(t_order) %>%
    dplyr::pull(t_indx)

  ## define function we want to integrate over
  to_int <- function(s) {
    exp(predict_basis(s, basis$knots, basis$boundary_knots) %*% gamma)
  }
  to_int_2 <- function(s) {
    matrix(predict_basis(s, basis$knots, basis$boundary_knots), ncol = 1) * as.numeric(exp(predict_basis(s, basis$knots, basis$boundary_knots) %*% gamma))
  }

  ## Set up computation of integrals using midpoint method
  midpoints <- 0.5*(tracker$val[2:nrow(tracker)] + tracker$val[1:(nrow(tracker)-1)])
  f_midpoints <- to_int(midpoints)
  f_midpoints_2 <- matrix(to_int_2(midpoints), ncol = K + 3, byrow = FALSE)
  widths <- tracker$val[2:(2*N)] - tracker$val[1:(2*N-1)]
  int_1 <- matrix(f_midpoints * widths, ncol = 1)
  int_2 <- f_midpoints_2 * widths

  ## Check for testing purposes
  # integrate(to_int, lower = -14, 2.35)

  ## Get indicator matrix for the observations
  ind_mat <- sapply(1:N, FUN = function(.x) {
    as.numeric(midpoints >= t_b[.x] & midpoints < e_b[.x])
  })
  new_basis <- predict_basis(midpoints, basis$knots, basis$boundary_knots)

  lklhd <- sum(delta * (basis$basis %*% gamma)[e_indx]) - sum(t(int_1) %*% ind_mat)

  score_1 <- t(X) %*% matrix(-delta * (basis$d_1 %*% gamma)[e_indx] +
                               exp(basis$basis %*% gamma)[e_indx] -
                               exp(basis$basis %*% gamma)[t_indx], ncol = 1)
  score_2 <- rowSums(t(delta * basis$basis[e_indx,]) - t(int_2) %*% ind_mat)
  score <- c(score_1, score_2)

  tst <- exp(basis$basis %*% gamma)[e_indx]
  tst2 <- (basis$d_1 %*% gamma)[e_indx]
  tst3 <- exp(basis$basis %*% gamma)[t_indx]
  tst4 <- (basis$d_1 %*% gamma)[t_indx]
  tst5 <- (basis$d_2 %*% gamma)[e_indx]
  tst6 <- basis$basis[e_indx,]
  tst7 <- basis$basis[t_indx,]
  tst8 <- basis$d_1[e_indx,]
  inform_1 <- matrix(0, nrow = p, ncol = p)
  inform_2 <- matrix(0, nrow = p, ncol = p_gamma)
  inform_4 <- matrix(0, nrow = p_gamma, ncol = p_gamma)

  for (i in 1:N) {
    inform_1 <- inform_1 + X[i,] %*% t(X[i,]) * (delta[i] * tst5[i] - tst[i] * tst2[i] + tst3[i] * tst4[i])
    inform_2 <- inform_2 + X[i,] %*% matrix(-tst8[i,] + tst6[i,] * tst[i] + tst7[i,] * tst3[i], nrow = 1)
  }
  for (j in 1:p_gamma) {
    inform_4[j,] <- -rowSums(t(new_basis[,j] * as.numeric(int_1) * new_basis) %*% ind_mat)
  }
  inform_3 <- t(inform_2)

  inform <- -rbind(cbind(inform_1, inform_2), cbind(inform_3, inform_4))

  list(
    lnlklhd = lklhd,
    score = score,
    inform = inform,
    basis = basis
  )

}
