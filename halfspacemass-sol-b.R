## only updated functions for evaluate_depth()

# evaluate depth for given test data, halfspaces, and metric of choice
# (package "foreach" used for parallelization, register backends to use)
# input: data: test data (matrix of number-of-test-data x space-dimension)
#        halfspaces: output of train_depth-function
#        metric: metric of choice (halfspace mass following Chen et al. or
#                Monte-Carlo approximate Tukey halfspace depth)
# output: vector of halfspace masses for each point in data
evaluate_depth <- function(data, halfspaces, metric = c("mass", "depth")) {
  require("foreach")
  ## input check ##
  if (!is.matrix(data)) {
    data <- try(as.matrix(data))
    if (inherits(data, "try-error")) {
      stop("<data> not convertable to a matrix")
    } else {
      message("<data> converted to matrix")
    }
  }
  checkmate::assert_matrix(
    data, mode = "numeric", any.missing = FALSE, min.rows = 1, min.cols = 1
  )
  if (!inherits(halfspaces, "trained_halfspaces")) {
    stop("<halfspaces> is not output of the function 'train_depth()'")
  }
  metric <- match.arg(metric)
  #################
  n_halfspaces <- length(halfspaces)
  n_data <- nrow(data)

  masses_of_points <- foreach(
    i = seq_len(n_data),
    .combine = c,
    .export = c("calculate_mass", "project_point_onto_directions")
  ) %dopar%
    calculate_mass(point = data[i, , drop = FALSE], halfspaces = halfspaces,
                   metric = metric)

  masses_of_points
}

# calculate mass for one point with respect to several halfspaces
# input: point: data point (matrix of 1 x space-dimension)
#        halfspaces: output of train_depth-function
#        metric: metric of choice (halfspace mass following Chen et al. or
#                Monte-Carlo approximate Tukey halfspace depth)
# output: vector of halfspace masses of input point
calculate_mass <- function(point, halfspaces, metric) {
  n_halfspaces <- length(halfspaces)

  directions <- t(sapply(halfspaces, "[[", "direction"))
  projected_points <- project_point_onto_directions(point = point, directions = directions)

  halfspaces_metrics <- t(sapply(halfspaces, function(x) {
    c(x[["split_point"]], x[["count_lower"]], x[["count_higher"]])
  }))
  colnames(halfspaces_metrics) <- c("split_point", "count_lower", "count_higher")

  if (metric == "mass") {
    masses_in_halfspaces <- ifelse(
      projected_points < halfspaces_metrics[, "split_point", drop = FALSE],
      halfspaces_metrics[, "count_lower", drop = FALSE],
      halfspaces_metrics[, "count_higher", drop = FALSE]
    )
    return(sum(masses_in_halfspaces) / n_halfspaces)
  } else {
    masses_in_halfspaces <- ifelse(
      projected_points < halfspaces_metrics[, "split_point", drop = FALSE],
      halfspaces_metrics[, "count_lower", drop = FALSE],
      NA
    )
    return(min(masses_in_halfspaces, na.rm = TRUE))
  }
}

# scalar projection of one point onto several vectors
# input: point: vector of space dimension
#        directions: matrix of vectors (number-of-vectors x space-dimension)
# output: vector of projections
project_point_onto_directions <- function(point, directions) {
  norms_of_directions <- sqrt(rowSums(directions^2))
  projected_points <- as.vector(point %*% t(directions)) / norms_of_directions
  projected_points
}

