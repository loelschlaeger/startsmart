### Illustration of the type-1 and type-3 initialization strategies

library("ggplot2")

set.seed(1)

N <- 300
Sigma <- matrix(c(2, 1, 1, 2), nrow = 2)
L <- t(chol(Sigma))

colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")

### two type-1 regressors of the first alternative, the systematic utility of
### the second alternative is fixed at zero; the same values are also read as
### the values of one type-3 regressor for the two alternatives, for which
### beta_3 = (beta_1[1], -beta_1[2]) gives the same utility difference

beta_1 <- c(1, -2)
beta_3 <- c(1, 2)

data <- data.frame(x_1 = numeric(N), x_2 = numeric(N), y = numeric(N))
for (n in seq_len(N)) {
  x_n <- stats::rnorm(2)
  U_n <- c(sum(x_n * beta_1), 0) + L %*% stats::rnorm(2)
  data[n, ] <- c(x_n, which.max(U_n))
}
y <- data$y

X <- matrix(0, nrow = 2 * N, ncol = 2)
X[c(TRUE, FALSE), ] <- as.matrix(data[, 1:2])
X_c <- scale(X, scale = FALSE)
y_c <- as.vector(t(diag(2)[y, ]))
y_c <- y_c - mean(y_c)
alpha_hat <- as.vector(solve(crossprod(X_c), crossprod(X_c, y_c)))
alpha_hat <- alpha_hat / sqrt(sum(alpha_hat^2))

M <- stats::cov(diag(2)[y, , drop = FALSE], as.matrix(data[, 1:2]))
s_hat <- svd(M)$v[, nrow(M)]
s_hat <- s_hat * sign(s_hat[1])

cat(sprintf(
  "type-1: beta direction = (%.2f, %.2f), estimate = (%.2f, %.2f)\n",
  beta_1[1] / sqrt(sum(beta_1^2)), beta_1[2] / sqrt(sum(beta_1^2)),
  alpha_hat[1], alpha_hat[2]
))
cat(sprintf(
  "type-3: beta = (%.2f, %.2f), estimate = (%.2f, %.2f)\n",
  beta_3[1], beta_3[2], 1, s_hat[1] / s_hat[2]
))

### panels

delta_v <- seq(-4, 4, by = 2)
cov_space <- ggplot(data, aes(x = x_1, y = x_2)) +
  geom_abline(
    slope = beta_3[1] / beta_3[2], intercept = -delta_v / beta_3[2], linetype = "dashed",
    linewidth = ifelse(delta_v == 0, 0.8, 0.3)
  ) +
  labs(x = expression(X[1]), y = expression(X[2])) +
  scale_x_continuous(breaks = seq(-2, 2, 2), limits = c(-3.5, 3.5)) +
  scale_y_continuous(breaks = seq(-2, 2, 2), limits = c(-3.5, 3.5)) +
  coord_fixed() +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.key.height = unit(3, "mm"),
    legend.margin = margin(t = -4),
    panel.grid.minor = element_blank(),
    plot.margin = margin(2, 5.5, 2, 5.5),
    plot.tag = element_text(face = "bold", size = 9),
    plot.tag.position = c(0.02, 0.98)
  )

points_a <- data.frame(
  x_1 = c(0, 2, 2), x_2 = c(0, 0, 1), label = c(0, 2, 0),
  nudge_x = c(-0.1, 0.5, 0.55), nudge_y = c(-0.9, -0.9, 0.9)
)
arrows_a <- data.frame(x = 0, y = 0, xend = 1.9, yend = c(0, 0.95))
panel_a <- cov_space +
  geom_segment(
    data = arrows_a, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
    colour = "white", linewidth = 1.6
  ) +
  geom_segment(
    data = arrows_a, aes(x = x, y = y, xend = xend, yend = yend), inherit.aes = FALSE,
    arrow = arrow(length = unit(0.06, "inches")), linewidth = 0.4
  ) +
  geom_point(data = points_a, aes(x = x_1, y = x_2), size = 1.4) +
  geom_label(
    data = points_a,
    aes(x = x_1 + nudge_x, y = x_2 + nudge_y, label = paste0("Delta[2]*V==", label)),
    parse = TRUE, size = 2.3, linewidth = 0, label.padding = unit(0.6, "mm")
  ) +
  labs(tag = "A")

panel_b <- cov_space +
  geom_point(aes(color = x_1 * beta_3[1] - x_2 * beta_3[2]), size = 0.8) +
  scale_color_gradientn(
    colours = c(colors[2], colors[2], "grey80", colors[1], colors[1]),
    values = scales::rescale(c(-9.5, -1.5, 0, 1.5, 9.5)),
    limits = c(-9.5, 9.5), breaks = seq(-6, 6, by = 3)
  ) +
  labs(color = expression(Delta[2] * V), tag = "B")

unit_length <- function(v, length = 3) v / sqrt(sum(v^2)) * length
ends <- rbind(
  unit_length(beta_1), unit_length(alpha_hat),
  unit_length(1 / beta_3), unit_length(s_hat)
)
arrows <- data.frame(xend = ends[, 1], yend = ends[, 2], color = colors[c(3, 4, 3, 4)])
arrow_labels <- data.frame(
  x = c(1.95, 3.15), y = c(-3.15, 1.95), label = c("beta^'(1)'", "s")
)

panel_c <- cov_space +
  geom_point(aes(color = factor(y, levels = c(1, 2))), size = 0.8, alpha = 0.35) +
  scale_colour_manual(
    values = colors[1:2], guide = guide_legend(override.aes = list(alpha = 1, size = 1.5))
  ) +
  geom_segment(
    data = arrows, aes(x = 0, y = 0, xend = xend, yend = yend), inherit.aes = FALSE,
    colour = "white", linewidth = 1.6
  ) +
  geom_segment(
    data = arrows, aes(x = 0, y = 0, xend = xend, yend = yend), inherit.aes = FALSE,
    colour = arrows$color, linewidth = 0.7,
    arrow = arrow(length = unit(0.08, "inches"), type = "closed")
  ) +
  geom_label(
    data = arrow_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE,
    parse = TRUE, size = 2.6, linewidth = 0, label.padding = unit(0.4, "mm")
  ) +
  labs(color = expression(y), tag = "C")

### layout: one row of three panels

width <- 5.4
grDevices::pdf(NULL, width = width, height = 5)
grobs <- lapply(list(panel_a, panel_b, panel_c), ggplot2::ggplotGrob)
g <- do.call(cbind, c(grobs, size = "max"))
is_null <- function(u) grid::unitType(u) == "null"
fixed_w <- sum(grid::convertWidth(g$widths[!is_null(g$widths)], "in", valueOnly = TRUE))
fixed_h <- sum(grid::convertHeight(g$heights[!is_null(g$heights)], "in", valueOnly = TRUE))
invisible(dev.off())
height <- fixed_h + (width - fixed_w) / 3

grDevices::cairo_pdf("figures/illustration.pdf", width = width, height = height)
grid::grid.draw(g)
invisible(dev.off())
