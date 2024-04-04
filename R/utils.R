# Function to extract necessary components
extract_surv_components <- function(formula, data) {
  # Validate input
  if (!inherits(formula, "formula")) {
    stop("First argument must be a formula.")
  }
  if (!is.data.frame(data)) {
    stop("Second argument must be a data frame.")
  }

  # Extract the response part of the formula (left-hand side)
  response <- formula[[2]]

  # Ensure the response is a Surv object
  if (!inherits(eval(response, envir = data), "Surv")) {
    stop("The response side of the formula must be a Surv object.")
  }

  # Extract response (Survival times and event indicator)
  surv_object <- eval(response, envir = data)
  survival_times <- surv_object[,1]
  event_indicator <- surv_object[,2]

  # Extract covariate names (right-hand side of the formula)
  covariates <- all.vars(formula[[3]])

  # Prepare the covariate matrix
  covariate_matrix <- model.matrix(formula, data)[, -1] # Removing the intercept

  # Return a list containing the components
  return(list(survival_times = survival_times,
              event_indicator = event_indicator,
              covariate_matrix = covariate_matrix))
}

# Example usage:
# Assuming 'data' is your data frame containing the survival times, event indicators,
# and covariates (e.g., age, treatment).
# formula <- Surv(time, status) ~ age + treatment
# components <- extract_surv_components(formula, data)


rect_integration = function(x, f)
{
  # check if the variable of integration is numeric
  if (!is.numeric(x))
  {
    stop('The first argument is not numeric.')
  }
  # check if f is a function
  if (!is.function(f))
  {
    stop('The second argument is not a function.')
  }
  ### finish checks
  # obtain length of variable of integration and integrand
  n.points = length(x)
  # midpoints
  midpoints = 0.5*(x[2:n.points] + x[1:(n.points-1)])
  # function evaluated at midpoints
  f.midpoints = f(midpoints)
  # calculate the widths of the intervals between adjacent pairs of points along the variable of integration
  interval.widths = x[2:n.points] - x[1:(n.points-1)]
  # implement rectangular integration
  # calculate the sum of all areas of rectangles that are used to approximate the integral
  rectangular.integral = sum(interval.widths * f.midpoints)
  # print the definite integral
  return(rectangular.integral)
}
