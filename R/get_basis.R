get_basis <- function(eps, K = 3) {
  knots <- get_knots(eps, K)
  basis <- splines2::bSpline(eps, knots = knots)
  d_1 <- splines2::dbs(eps, knots = knots, derivs = 1L)
  d_2 <- splines2::dbs(eps, knots = knots, derivs = 2L)

  return(
    list(
      basis = basis,
      d_1 = d_1,
      d_2 = d_2,
      knots = knots,
      boundary_knots = c(min(eps), max(eps))
    )
  )
}

predict_basis <- function(x, knots, boundary_knots, deriv = 0) {
  if (deriv == 0) {
    val <- splines2::bSpline(x, knots = knots, Boundary.knots = boundary_knots)
  } else if (deriv == 1) {
    val <- splines2::dbs(x, knots = knots, Boundary.knots = boundary_knots, derivs = 1L)
  } else if (deriv == 2) {
    val <- splines2::dbs(x, knots = knots, Boundary.knots = boundary_knots, derivs = 2L)
  }
  return(val)
}


get_knots <- function(eps, K = 3) {
  # Ensure the vector is sorted
  eps <- sort(eps)

  # Number of elements in each group, at least approximately
  n <- length(eps)
  group_size <- n / (K + 1)

  # Calculate cut points (not including the minimum and maximum)
  cut_points <- numeric(K)
  for (i in 1:K) {
    cut_point_index <- round(i * group_size)
    cut_points[i] <- eps[cut_point_index]
  }

  return(cut_points)
}
