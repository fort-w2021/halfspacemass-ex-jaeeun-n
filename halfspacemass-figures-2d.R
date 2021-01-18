# copy-pasted code from "topdown-halfspacemass-ex.Rmd"

### preparation ####
library(ggplot2)
theme_set(theme_minimal())

# visualize half-space depth/mass values for 2D-data on a grid of points
# if no grid is provided, a grid over min/max +/- .2 * range of data is created.
#   data: a 2d data.frame with columns z1 and z2
#   points: add points in data to plot?
#   metric: argument for evaluate_depth
plot_depth <- function(halfspaces, data, grid = NULL, points = TRUE,
                       gridlength = 70, metric = "mass") {
  # not doing input checks here because this is just tooling/diagnostics...
  if (is.null(grid)) {
    range_1 <- range(data$z1)
    range_2 <- range(data$z2)
    grid <- expand.grid(
      z1 = seq(range_1[1] - .2 * diff(range_1),
               range_1[2] + .2 * diff(range_1),
               length = gridlength
      ),
      z2 = seq(range_2[1] - .2 * diff(range_2),
               range_2[2] + .2 * diff(range_2),
               length = gridlength
      )
    )
  }
  grid_depth <- evaluate_depth(
    data = as.matrix(grid),
    halfspaces = halfspaces,
    metric = metric)
  grid_halfspaces <- cbind(grid, depth = grid_depth)
  # use colors as in Chen et al.:
  spectralcolors <- c(
    "darkblue", "blue", "cyan", "lightgreen",
    "yellow", "orange", "red", "darkred"
  )
  p <- ggplot(grid_halfspaces, aes(x = z1, y = z2)) +
    geom_tile(aes(fill = depth, colour = depth)) +
    scale_fill_gradientn(metric, colors = spectralcolors) +
    scale_colour_gradientn(metric, colors = spectralcolors)

  if (points & !is.null(data)) {
    p <- p +
      geom_point(data = data,
                 aes(x = z1, y = z2),
                 colour = rgb(1, 1, 1, .8))
  }
  p
}


### figure 3 ####
library(gridExtra)
data_fig3 <- data.frame(z1 = c(-2, -.5, .5, 2), z2 = 0)
grid_fig3 <- expand.grid(z1 = seq(-3, 3, l = 51),
                         z2 = seq(-3, 3, l = 51))

depth_fig3 <- train_depth(data_fig3, n_halfspace = 1e4,
                          scope = 1, seed = 4163)
# need scope > 1 for reliable halfspace _depth_ approximation:
depth_fig3_scope15 <- train_depth(data_fig3,
                                  n_halfspace = 1e3, scope = 1.5,
                                  seed = 4163)
gridExtra::grid.arrange(
  plot_depth(depth_fig3_scope15,
             data = data_fig3, grid = grid_fig3,
             metric = "depth") +
    ggtitle("Tukey Halfspace Depth"),
  plot_depth(depth_fig3, data = data_fig3, grid = grid_fig3) +
    ggtitle("Halfspace Mass (Chen et al.)"),
  nrow = 1
)


### figure 5 ####
set.seed(187471431)
# 2D standard Normal:
cluster <- data.frame(
  z1 = rnorm(50) / 2,
  z2 = rnorm(50) / 2,
  group = "cluster"
)
# polar coordinates: points with distance 3 to 5 from the origin, at 90° - 270°:
left_anomalies <- data.frame(
  angle = runif(10, pi / 2, 3 * pi / 2),
  length = runif(10, 3, 5)
)
# convert to cartesian coords
left_anomalies <- with(left_anomalies, data.frame(
  z1 = length * cos(angle),
  z2 = length * sin(angle),
  group = "anomaly"
))
# ~ N_2(\mu = (6,0), \Sigma = I_2)
right_anomalies <- data.frame(
  z1 = rnorm(20) / 5 + 6,
  z2 = rnorm(20) / 5,
  group = "anomaly"
)
data_fig5 <- rbind(cluster,
                   left_anomalies,
                   right_anomalies)

hs_fig5 <- train_depth(data_fig5[, 1:2],
                       n_halfspace = 1e4, subsample = .5,
                       seed = 4165
)
fig5 <- plot_depth(hs_fig5, data = data_fig5[, 1:2], points = FALSE)
# can't assign two colour scales to one plot, so plot 2 groups separately:
fig5 +
  geom_point(
    data = subset(data_fig5, group == "cluster"),
    aes(x = z1, y = z2), color = rgb(0, 0, 1, .5)
  ) +
  geom_point(
    data = subset(data_fig5, group == "anomaly"),
    aes(x = z1, y = z2), color = rgb(1, 0, 0, .5)
  )



