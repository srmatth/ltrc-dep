lnlklhd_integration <- function(beta, gamma, X, y, t, delta, K = 3) {
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
  basis <- get_basis(tracker$val, K = K)

  ## define function we want to integrate over
  to_int <- function(s, k, b_k, g) {
    exp(predict_basis(s, k, b_k) %*% g)
  }
  to_int_2 <- function(s, k, b_k, g, K) {
    predict_basis(s, basis$knots, basis$boundary_knots)[,K] * as.numeric(exp(predict_basis(s, basis$knots, basis$boundary_knots) %*% gamma))
  }
  to_int_3 <- function(s, k, b_k, g, K, L) {
    predict_basis(s, basis$knots, basis$boundary_knots)[,L] * predict_basis(s, basis$knots, basis$boundary_knots)[,K] * as.numeric(exp(predict_basis(s, basis$knots, basis$boundary_knots) %*% gamma))
  }

  lnlklhd <- 0
  score <- rep(0, p + p_gamma)
  inform <- matrix(0, nrow = p + p_gamma, ncol = p + p_gamma)
  e_basis <- predict_basis(e_b, basis$knots, basis$boundary_knots)
  t_basis <- predict_basis(t_b, basis$knots, basis$boundary_knots)
  e_basis_d <- predict_basis(e_b, basis$knots, basis$boundary_knots, deriv = 1)
  t_basis_d <- predict_basis(t_b, basis$knots, basis$boundary_knots, deriv = 1)
  e_basis_d_2 <- predict_basis(e_b, basis$knots, basis$boundary_knots, deriv = 2)
  for (i in 1:N) {
    int_1 <- integrate(
      f = to_int,
      lower = t_b[i],
      upper = e_b[i],
      k = basis$knots,
      b_k = basis$boundary_knots,
      g = gamma
    )[[1]]
    lnlklhd <- lnlklhd + delta[i] * e_basis[i,] %*% gamma - int_1
    for (j in 1:p) {
      score[j] <- score[j] + X[i,j] * (-delta[i] * e_basis_d[i,] %*% gamma +
                                         exp(e_basis[i,] %*% gamma) - exp(t_basis[i,] %*% gamma))
      for (k in 1:p) {
        inform[j,k] <- inform[j,k] + X[i,j] * X[i,k] *
          (delta[i] * e_basis_d_2[i,] %*% gamma -
             exp(e_basis[i,] %*% gamma) * e_basis_d[i,] %*% gamma +
             exp(t_basis[i,] %*% gamma) * t_basis_d[i,] %*% gamma)
        if (j != k) {
          inform[k,j] <- inform[j,k]
        }
      }
      for (k in 1:p_gamma) {
        inform[j, p+k] <- inform[j, p+k] + X[i,j] * (-e_basis_d[i,k] +
                                                       e_basis[i,k] * exp(e_basis[i,] %*% gamma) +
                                                       t_basis[i,k] * exp(t_basis[i,] %*% gamma))
        inform[p+k, j] <- inform[j, p+k]
      }
    }

    for (j in 1:p_gamma) {
      int_2 <- integrate(
        f = to_int_2,
        lower = t_b[i],
        upper = e_b[i],
        k = basis$knots,
        b_k = basis$boundary_knots,
        g = gamma,
        K = j
      )[[1]]
      score[p+j] <- score[p+j] + e_basis[i, j] - int_2

      for (m in 1:p_gamma) {
        int_3 <- integrate(
          f = to_int_3,
          lower = t_b[i],
          upper = e_b[i],
          k = basis$knots,
          b_k = basis$boundary_knots,
          g = gamma,
          K = j,
          L = m
        )[[1]]
        inform[p + j, p+m] = inform[p + j, p+m] - int_3
        inform[p + m, p+j] = inform[p + j, p+m] - int_3
      }
    }
  }
  inform = -inform

  list(
    lnlklhd = lklhd,
    score = score,
    inform = inform,
    basis = basis
  )

}
