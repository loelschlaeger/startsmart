### Evaluate the simulation results and create figures and tables for the paper

library(ggplot2)

set.seed(1)

ll_tol <- 5
dist_tol <- 5

scenarios <- data.frame(
  id = c("S1", "S2", "S3", "S4", "S5", "S6"),
  label = c(
    "type_1", "type_2", "type_3", "mixed_fixed", "random_type_1",
    "mixed_random"
  ),
  effects = c(
    "10 type-I",
    "5 type-II",
    "5 type-III",
    "3 type-I, 3 type-II, 3 type-III",
    "5 type-I, all random",
    "2 type-I, 2 type-II, 2 type-III, one of each random"
  ),
  J = c(5, 5, 5, 5, 5, 3),
  N = c(5000, 5000, 5000, 5000, 200, 200),
  Tp = c(1, 1, 1, 1, 10, 10),
  objective = c(
    "log-likelihood", "log-likelihood", "log-likelihood", "log-likelihood",
    "pairwise CML", "pairwise CML"
  ),
  stringsAsFactors = FALSE
)

init_levels <- c("Random", "Proposed", "DGP")

### read results ---------------------------------------------------------------

### the runs are stored as one archive per result directory
results_path <- function(dir) {
  if (dir.exists(dir)) return(dir)
  archive <- paste0(dir, ".tar.gz")
  if (!file.exists(archive)) stop("Results '", dir, "' not found.")
  utils::untar(archive, exdir = tempdir())
  file.path(tempdir(), dir)
}

read_scenario <- function(label) {
  Sim <- oeli::Simulator$new(use_backup = results_path(file.path("sim", label)))
  results <- Sim$results
  incomplete <- sum(Sim$cases$.error, na.rm = TRUE) +
    sum(Sim$cases$.pending, na.rm = TRUE)
  if (incomplete > 0) {
    warning(
      incomplete, " case(s) failed or are missing in scenario '", label, "'.",
      call. = FALSE
    )
  }
  rows <- lapply(seq_len(nrow(results)), function(k) {
    o <- results$.output[[k]]
    dgp <- o$theta_dgp
    rel_dist <- function(theta) sqrt(sum((theta - dgp)^2)) / sqrt(sum(dgp^2))
    best <- max(o$ll_fit_random, o$ll_fit_smart, o$ll_fit_true)
    one <- function(init, suffix, sec_init) {
      data.frame(
        label = label,
        case = results$.case[k],
        init = init,
        q = length(dgp),
        sec_init = sec_init,
        sec_fit = o[[paste0("sec_fit_", suffix)]],
        iterations = o[[paste0("nlm_iterations_", suffix)]],
        code = o[[paste0("nlm_code_", suffix)]],
        ll_init = o[[paste0("ll_init_", suffix)]],
        ll_fit = o[[paste0("ll_fit_", suffix)]],
        ll_best = best,
        dist_init = rel_dist(o[[paste0("theta_init_", suffix)]]),
        dist_fit = rel_dist(o[[paste0("theta_fit_", suffix)]]),
        stringsAsFactors = FALSE
      )
    }
    rbind(
      one("Random", "random", 0),
      one("Proposed", "smart", o$sec_init_smart),
      one("DGP", "true", 0)
    )
  })
  do.call(rbind, rows)
}

res <- do.call(rbind, lapply(scenarios$label, read_scenario))
res$init <- factor(res$init, levels = init_levels)
res$id <- factor(
  scenarios$id[match(res$label, scenarios$label)], levels = scenarios$id
)
res$sec_total <- res$sec_init + res$sec_fit
res$success <- res$ll_fit >= res$ll_best - ll_tol & res$dist_fit <= dist_tol

### classify the outcome of every run
###   diverged         the parameter norm ran away, see condition (ii) above
###   iteration limit  nlm code 4, the optimizer was still moving when it stopped
###   stalled          nlm codes 3 and 5, the last step failed to improve or the
###                    maximum step size was exceeded repeatedly
###   inferior optimum nlm codes 1 and 2, the optimizer converged to a
###                    stationary point below the best objective value
res$outcome <- ifelse(
  res$success, "successful",
  ifelse(
    res$dist_fit > dist_tol, "diverged",
    ifelse(
      res$code == 4, "iteration limit",
      ifelse(res$code %in% c(3, 5), "stalled", "inferior optimum")
    )
  )
)
outcome_levels <- c(
  "successful", "diverged", "iteration limit", "stalled", "inferior optimum"
)
res$outcome <- factor(res$outcome, levels = outcome_levels)

### figure: optimization outcomes ---------------------------------------------

dir.create("figures", showWarnings = FALSE)

long <- rbind(
  transform(res, metric = "Total computation time in seconds", value = sec_total),
  transform(res, metric = "Iterations until termination", value = iterations)
)
long$metric <- factor(
  long$metric,
  levels = c("Total computation time in seconds", "Iterations until termination")
)
long$outcome <- factor(
  ifelse(long$success, "successful", "failed"),
  levels = c("successful", "failed")
)

theme_paper <- theme_bw(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 8),
    legend.position = "bottom",
    legend.margin = margin(t = -6),
    panel.grid.minor = element_blank()
  )

p_opt <- ggplot(long, aes(x = init, y = value)) +
  geom_boxplot(fill = NA, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(
    aes(shape = outcome, color = outcome),
    width = 0.2, height = 0, size = 0.9, alpha = 0.6, stroke = 0.5
  ) +
  facet_grid(
    metric ~ id, scales = "free_y", switch = "y",
    labeller = labeller(metric = label_wrap_gen(20))
  ) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000), labels = scales::comma) +
  scale_shape_manual(values = c(16, 4), name = NULL) +
  scale_color_manual(values = c("grey20", "firebrick"), name = NULL) +
  labs(x = "Initialization", y = NULL) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("figures/sim_optimization.pdf", p_opt, width = 5.4, height = 3.9, device = cairo_pdf)

### tables -----------------------------------------------------------------------

dir.create("tables", showWarnings = FALSE)

fmt <- function(x, digits = 1) formatC(x, format = "f", digits = digits, big.mark = ",")
pct <- function(x) paste0(round(100 * x), "\\%")
agg <- function(v, f, subset = rep(TRUE, nrow(res))) {
  tapply(v[subset], list(res$id[subset], res$init[subset]), f)
}

## scenario overview
q_by_label <- tapply(res$q, res$label, unique)
sc <- scenarios
sc$q <- q_by_label[sc$label]
scen_rows <- sprintf(
  "%s & %s & %d & %s & %d & %d & %s \\\\",
  sc$id, sc$effects, sc$J, fmt(sc$N, 0), sc$Tp, sc$q, sc$objective
)
writeLines(c(
  "\\begin{tabular}{llrrrrl}",
  "\\toprule",
  "Scenario & Coefficients & $J$ & $N$ & $T$ & $q$ & Objective \\\\",
  "\\midrule",
  scen_rows,
  "\\bottomrule",
  "\\end{tabular}"
), "tables/sim_scenarios.tex")

## initial values
is_prop <- res$init == "Proposed"
med_init_time <- tapply(res$sec_init[is_prop], res$id[is_prop], median)
med_dist <- agg(res$dist_init, median)
med_gap <- agg(res$ll_best - res$ll_init, median)
ini_rows <- sprintf(
  "%s & %s & %s & %s & %s & %s \\\\",
  scenarios$id,
  fmt(med_init_time[scenarios$id]),
  fmt(med_dist[scenarios$id, "Random"], 2), fmt(med_dist[scenarios$id, "Proposed"], 2),
  fmt(med_gap[scenarios$id, "Random"], 0), fmt(med_gap[scenarios$id, "Proposed"], 0)
)
writeLines(c(
  "\\begin{tabular}{lrrrrr}",
  "\\toprule",
  " & Computation time & \\multicolumn{2}{c}{Relative distance} & \\multicolumn{2}{c}{$l(\\hat{\\theta}) - l(\\theta^{(0)})$} \\\\",
  "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6}",
  "Scenario & Proposed & Random & Proposed & Random & Proposed \\\\",
  "\\midrule",
  ini_rows,
  "\\bottomrule",
  "\\end{tabular}"
), "tables/sim_initial.tex")

## optimization outcomes
succ <- agg(res$success, mean)
med_time_s <- agg(res$sec_total, median, res$success)
med_iter_s <- agg(res$iterations, median, res$success)
res_rows <- sprintf(
  "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  scenarios$id,
  pct(succ[scenarios$id, "Random"]), pct(succ[scenarios$id, "Proposed"]), pct(succ[scenarios$id, "DGP"]),
  fmt(med_time_s[scenarios$id, "Random"]), fmt(med_time_s[scenarios$id, "Proposed"]), fmt(med_time_s[scenarios$id, "DGP"]),
  fmt(med_iter_s[scenarios$id, "Random"], 0), fmt(med_iter_s[scenarios$id, "Proposed"], 0), fmt(med_iter_s[scenarios$id, "DGP"], 0)
)
writeLines(c(
  "\\begin{tabular}{lrrrrrrrrr}",
  "\\toprule",
  " & \\multicolumn{3}{c}{(a) Success rate} & \\multicolumn{3}{c}{(b) Total time in seconds} & \\multicolumn{3}{c}{(c) Iterations} \\\\",
  "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7} \\cmidrule(lr){8-10}",
  "Scenario & Random & Proposed & DGP & Random & Proposed & DGP & Random & Proposed & DGP \\\\",
  "\\midrule",
  res_rows,
  "\\bottomrule",
  "\\end{tabular}"
), "tables/sim_results.tex")

## outcomes of the optimization runs
outcome_rows <- character()
for (level in init_levels) {
  outcome_rows <- c(outcome_rows, sprintf(
    "\\multicolumn{%d}{l}{\\textit{%s initialization}} \\\\",
    length(outcome_levels) + 1L, level
  ))
  for (id in scenarios$id) {
    counts <- table(res$outcome[res$id == id & res$init == level])
    outcome_rows <- c(outcome_rows, paste0(
      id, " & ", paste(as.integer(counts), collapse = " & "), " \\\\"
    ))
  }
  if (level != init_levels[length(init_levels)]) {
    outcome_rows <- c(outcome_rows, "\\midrule")
  }
}
writeLines(c(
  paste0("\\begin{tabular}{l", strrep("r", length(outcome_levels)), "}"),
  "\\toprule",
  paste0(
    "Scenario & Successful & Diverged & Iteration limit & Stalled &",
    " Inferior optimum \\\\"
  ),
  "\\midrule",
  outcome_rows,
  "\\bottomrule",
  "\\end{tabular}"
), "tables/sim_outcomes.tex")

### console summary ---------------------------------------------------------------

cat("Success rates (objective gap <=", ll_tol, "and relative distance <=", dist_tol, "):\n")
print(round(succ, 2))
cat("\nMedian total time, all runs:\n"); print(round(agg(res$sec_total, median), 1))
cat("\nMedian total time, successful runs:\n"); print(round(med_time_s, 1))
cat("\nMedian iterations, all runs:\n"); print(round(agg(res$iterations, median)))
cat("\nMedian iterations, successful runs:\n"); print(round(med_iter_s))
cat("\nMedian initial distance:\n"); print(round(med_dist, 2))
cat("\nMedian objective gap at initialization:\n"); print(round(med_gap))
cat("\nMedian initialization time (Proposed):\n"); print(round(med_init_time, 1))
cat("\nnlm termination codes by init:\n"); print(table(res$init, res$code, res$id))
cat("\noutcomes by init and scenario:\n"); print(table(res$outcome, res$init, res$id))
cat("\nFailed runs (init != Random):\n")
print(res[!res$success & res$init != "Random",
          c("id", "case", "init", "ll_fit", "ll_best", "code", "dist_init", "dist_fit")])

invisible(res)
