### Evaluate the application results and create the figure and the table

library(ggplot2)

set.seed(1)

app <- "cars_us"
ll_tol <- 2

### read results ---------------------------------------------------------------

### the runs are stored as one archive per result directory
results_path <- function(dir) {
  if (dir.exists(dir)) return(dir)
  archive <- paste0(dir, ".tar.gz")
  if (!file.exists(archive)) stop("Results '", dir, "' not found.")
  utils::untar(archive, exdir = tempdir())
  file.path(tempdir(), dir)
}

### the two initializations are stored in separate backups
outputs <- list()
for (init in c("proposed", "random")) {
  Sim <- oeli::Simulator$new(
    use_backup = results_path(file.path("application", paste0(app, "_", init)))
  )
  incomplete <- sum(Sim$cases$.error, na.rm = TRUE) +
    sum(Sim$cases$.pending, na.rm = TRUE)
  if (incomplete > 0) {
    warning(
      incomplete, " ", init, " case(s) failed or are missing.", call. = FALSE
    )
  }
  outputs <- c(outputs, Sim$results$.output)
}
theta_names <- outputs[[1L]]$parameter_names
q <- length(theta_names)

attribute_labels <- c(
  hev = "HEV", phev10 = "PHEV 10", phev20 = "PHEV 20", phev40 = "PHEV 40",
  bev75 = "BEV 75", bev100 = "BEV 100", bev150 = "BEV 150",
  phevFastcharge = "PHEV fast charging", bevFastcharge = "BEV fast charging",
  american = "American", japanese = "Japanese", chinese = "Chinese",
  skorean = "South Korean", price = "Price", opCost = "Operating cost",
  accelTime = "Acceleration time"
)
is_cov <- grepl("^[OL][0-9]+$", theta_names)
label_tex <- ifelse(
  is_cov, sub("^([OL])([0-9]+)$", "$\\1_{\\2}$", theta_names),
  attribute_labels[theta_names]
)
label_plotmath <- ifelse(
  is_cov, sub("^([OL])([0-9]+)$", "\\1[\\2]", theta_names),
  paste0("\"", attribute_labels[theta_names], "\"")
)
if (anyNA(label_tex)) stop("Parameter name(s) without label.")

res <- do.call(rbind, lapply(seq_along(outputs), function(k) {
  o <- outputs[[k]]
  data.frame(
    case = k,
    init = if (o$init_type == "proposed") "Proposed" else "Random",
    ll_init = o$ll_init,
    ll_fit = o$ll_fit,
    sec = o$sec_init + o$sec_fit,
    sec_init = o$sec_init,
    iterations = o$nlm_iterations,
    code = o$nlm_code,
    stringsAsFactors = FALSE
  )
}))
res$init <- factor(res$init, levels = c("Random", "Proposed"))
ll_best <- max(res$ll_fit)
theta_hat <- outputs[[which.max(res$ll_fit)]]$theta_fit
res$success <- res$ll_fit >= ll_best - ll_tol

### console summary ------------------------------------------------------------

cat(sprintf("\n%s: reference value l(theta_hat) = %.2f, q = %d\n", app, ll_best, q))
for (level in levels(res$init)) {
  d <- res[res$init == level, ]
  s <- d[d$success, ]
  cat(sprintf(
    paste(
      "%-9s | n = %3d | successful %3.0f%% | median iterations %4.0f |",
      "median time %6.1f min | initialization %4.1f min |",
      "per successful estimate %5.0f iterations, %6.1f min\n"
    ),
    level, nrow(d), 100 * mean(d$success),
    if (nrow(s)) median(s$iterations) else NA_real_,
    if (nrow(s)) median(s$sec) / 60 else NA_real_,
    median(d$sec_init) / 60, sum(d$iterations) / sum(d$success),
    sum(d$sec) / sum(d$success) / 60
  ))
}
cat(sprintf(
  "objective gap at the initial values: proposed %.1f, median random %.1f\n",
  ll_best - res$ll_init[res$init == "Proposed"][1],
  median(ll_best - res$ll_init[res$init == "Random"])
))
cat(sprintf(
  "failed random runs: smallest gap at termination %.1f, median iterations %.1f\n",
  min(ll_best - res$ll_fit[res$init == "Random" & !res$success]),
  median(res$iterations[res$init == "Random" & !res$success])
))

### figure ---------------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)

init_colors <- c("Random" = "black", "Proposed" = "#08519c")
outcome_shapes <- c("successful" = 16, "failed" = 4)
outcome_colors <- c("successful" = "grey20", "failed" = "firebrick")
theme_paper <- theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -6),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
  )

dev <- do.call(rbind, lapply(seq_along(outputs), function(k) {
  o <- outputs[[k]]
  data.frame(
    parameter = factor(theta_names, levels = theta_names),
    deviation = o$theta_init - theta_hat,
    init = res$init[k],
    stringsAsFactors = FALSE
  )
}))
p_dev <- ggplot(dev, aes(x = parameter, y = deviation, color = init)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_jitter(
    data = dev[dev$init == "Random", ],
    width = 0.2, height = 0, size = 0.6, alpha = 0.3
  ) +
  geom_point(
    data = dev[dev$init == "Proposed", ], size = 1.6, shape = 18
  ) +
  scale_color_manual(
    values = init_colors, name = NULL,
    guide = guide_legend(override.aes = list(size = 1.6, alpha = 1, shape = 16))
  ) +
  scale_x_discrete(labels = parse(text = label_plotmath)) +
  labs(x = "Parameter", y = expression(theta^{(0)} - hat(theta))) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7))

perf <- rbind(
  transform(res, metric = "Computation time in minutes", value = sec / 60),
  transform(res, metric = "Iterations", value = iterations)
)
perf$metric <- factor(
  perf$metric, levels = c("Computation time in minutes", "Iterations")
)
perf$outcome <- factor(
  ifelse(perf$success, "successful", "failed"),
  levels = c("successful", "failed")
)

enough <- names(which(table(res$init) >= 5L))
p_perf <- ggplot(perf, aes(x = init, y = value)) +
  geom_boxplot(
    data = perf[perf$init %in% enough, ],
    fill = NA, outlier.shape = NA, linewidth = 0.3
  ) +
  geom_jitter(
    data = perf[perf$init %in% enough, ],
    aes(shape = outcome, color = outcome),
    width = 0.2, height = 0, size = 0.8, alpha = 0.7, stroke = 0.5
  ) +
  geom_point(
    data = perf[!perf$init %in% enough, ],
    size = 1.8, color = init_colors[["Proposed"]]
  ) +
  facet_wrap(~ metric, scales = "free_y", strip.position = "left") +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000),
    labels = scales::label_number(drop0trailing = TRUE, big.mark = ",")
  ) +
  scale_shape_manual(values = outcome_shapes, name = NULL) +
  scale_color_manual(values = outcome_colors, name = NULL) +
  labs(x = "Initialization", y = NULL) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

cairo_pdf(file.path("figures", paste0("app_", app, ".pdf")), width = 5.4, height = 3.2)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(
  1, 2, widths = grid::unit(c(0.58, 0.42), "npc")
)))
print(p_dev, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_perf, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
invisible(dev.off())

### table: estimates and proposed initial values --------------------------------

dir.create("tables", showWarnings = FALSE)
fmt <- function(x, digits = 2) {
  formatC(x, format = "f", digits = digits, big.mark = ",")
}
num <- function(x) paste0("$", fmt(x), "$")
theta_0 <- outputs[[which(res$init == "Proposed")[1L]]]$theta_init
n_blocks <- 4L
blocks <- split(seq_len(q), ceiling(seq_len(q) / ceiling(q / n_blocks)))
row_block <- function(idx) c(
  paste0("Parameter & ", paste(label_tex[idx], collapse = " & "), " \\\\"),
  "\\midrule",
  paste0("$\\hat{\\theta}$ & ", paste(num(theta_hat[idx]), collapse = " & "), " \\\\"),
  paste0("$\\theta^{(0)}$ & ", paste(num(theta_0[idx]), collapse = " & "), " \\\\")
)
writeLines(c(
  paste0("\\begin{tabular}{l", strrep("r", max(lengths(blocks))), "}"),
  "\\toprule",
  unlist(lapply(seq_along(blocks), function(b) {
    c(if (b > 1L) "\\midrule", row_block(blocks[[b]]))
  })),
  "\\bottomrule",
  "\\end{tabular}"
), file.path("tables", paste0("app_", app, "_estimates.tex")))

invisible(res)
