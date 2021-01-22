# train halfspace depths following algorithm 1 in
# "Chen, B., Ting, K.M., Washio, T. et al. (2015); Half-space mass: a maximally
# robust and efficient data depth method; Machine Learning, 100(2):677–699"
# (package "doRNG" used for parallelization, register backends to use)
# input: data: training data (matrix of number-points x space-dimension),
#        n_halfspace: number of halfspaces to sample,
#        subsample: proportion of data used to calculate each halfspace
#                   (t in paper),
#        scope: parameter determining size of convex region covering a data
#               (lambda in paper),
#        seed: for random numbers
# output: matrix with random direction in data space (direction) split point
#         (split_point), number of points lower (count_lower) and higher
#         (count_higher) than the split point
train_depth <- function(data, n_halfspace, subsample = 1, scope = 1, seed = 123456) {
  require("doRNG")
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
  checkmate::assert_count(n_halfspace, positive = TRUE)
  checkmate::assert_number(subsample, lower = 0, upper = 1)
  checkmate::assert_number(scope, lower = 1)
  checkmate::assert_integerish(seed)
  #################
  set.seed(seed)

  trained_halfspaces <- foreach(
    seq_len(n_halfspace),
    .export = c("train_one_halfspace", "select_split_point",
                "project_points_onto_direction")
  ) %dorng%
    train_one_halfspace(
      data = data, n_halfspace = n_halfspace,
      subsample = subsample, scope = scope
    )

  attr(trained_halfspaces, "rng") <- NULL
  class(trained_halfspaces) <- "trained_halfspaces"
  trained_halfspaces
}

# train one halfspace with a sample from data
# input: data, n_halfspace, subsample, scope as for train_depth-function
# output: list with random direction in data space (random_direction) split point
#         (split_point), number of points lower (count_lower) and higher
#         (count_higher) than the split point for one halfspace
train_one_halfspace <- function(data, n_halfspace, subsample, scope) {
  space_dimension <- ncol(data)
  random_direction <- runif(space_dimension, min = -100, 100)
  subsample_size <- subsample * nrow(data)
  random_sample <- data[sample(nrow(data), subsample_size), , drop = FALSE]

  projected_points <- project_points_onto_direction(
    data = random_sample,
    direction = random_direction
  )

  split_point <- select_split_point(projected_points, scope)

  count_lower <- max(c(0, sum(projected_points < split_point)), na.rm = TRUE) /
    subsample_size
  count_higher <- (length(projected_points) - count_lower) /
    subsample_size

  list(direction = random_direction, split_point = split_point,
       count_lower = count_lower, count_higher = count_higher)
}

# scalar projection of data points onto one vector
# input: data: data points (matrix of number-points x space-dimension),
#        direction: vector
# output: vector of projected points
project_points_onto_direction <- function(data, direction) {
  as.vector(data %*% direction) / as.numeric(sqrt(crossprod(direction)))
}

# calculate the split point
# input: projected_points: projected points in sample as vector,
#        scope: parameter determining size of convex region covering a data
#               (lambda in paper)
# output: split point (scalar)
select_split_point <- function(projected_points, scope) {
  max_of_projections <- max(projected_points)
  min_of_projections <- min(projected_points)
  mid_of_projections <- 0.5 * (max_of_projections + min_of_projections)

  limit_distance <- 0.5 * scope * (max_of_projections - min_of_projections)
  lower_limit <- mid_of_projections - limit_distance
  upper_limit <- mid_of_projections + limit_distance

  runif(1, min = lower_limit, max = upper_limit)
}


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

