# Testing approximate vs exact Halfspace Depth in 3D
library(scatterplot3d)

# make depth::depth accept matrices containing points u from test as rows:
hsdepth <- function(test, train) {
  apply(test, 1, depth::depth, x = train)
}

# wrap train and test steps into a single function:
my_hsdepth <- function(test, train, ...) {
  evaluate_depth(
    data = test,
    halfspaces = train_depth(train, ...),
    metric = "depth"
  )
}


# generate list of test cases:
set.seed(121133)
n <- 100
theta <- runif(n, 0, 2 * pi)
u <- runif(n, -1, 1)
x <- sqrt(1 - u^2) * cos(theta)
y <- sqrt(1 - u^2) * sin(theta)
z <- u

data_list <- list(
  sphere = cbind(x, y, z),
  gaussian = cbind(rnorm(100), rnorm(100), rnorm(100)),
  clustered = cbind(
    rnorm(100, mean = rep(c(0, 2, 0, 0), l = 100), sd = .2),
    rnorm(100, mean = rep(c(0, 0, 2, 0), l = 100), sd = .2),
    rnorm(100, mean = rep(c(0, 0, 0, 2), l = 100), sd = .2)
  ),
  grid = expand.grid(seq(-2, 2, l = 5), seq(-2, 2, l = 5), seq(-2, 2, l = 5))
)
data_list <- lapply(data_list, as.data.frame)

layout(matrix(1:4, 2, 2))
lapply(data_list, scatterplot3d)

# test approximate equivalence for the 4 datasets above:
# (values are scaled differently, so compare ranks and/or check correlation)
for (train in data_list) {
  approx_depth <- my_hsdepth(
    test = data_list$grid, train = train,
    n_halfspace = 1e4, scope = 1, seed = 23
  )
  exact_depth <- hsdepth(test = data_list$grid, train = train)
  testthat::expect_true(
    cor(
      approx_depth,
      exact_depth
    ) > 0.99
  )
  testthat::expect_equivalent(
    rank(approx_depth),
    rank(exact_depth),
    tol = nrow(data_list$grid) / 20
  )
}
