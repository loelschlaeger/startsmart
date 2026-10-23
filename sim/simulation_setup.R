### Simulation setup
###
### Run simulation settings:
###
###   Rscript -e 'source("sim/simulation_setup.R")'
###
### Optional variables can be defined before sourcing:
###
###   install_dep <- TRUE        # install package dependencies (default TRUE)
###   settings_to_run <- 1:6     # subset of settings (default: all)
###   runs <- 100L               # replications per setting (default 100)
###   workers <- NULL            # parallel workers (default: cores - 1)
###
### The working directory must be the project root. Results are saved as
### sim/<label>/case_*.rds. A setting whose output directory already exists is
### skipped.

if (!exists("install_dep")) install_dep <- TRUE
if (!exists("settings_to_run")) settings_to_run <- NULL
if (!exists("runs")) runs <- 100L
if (!exists("workers")) workers <- NULL

### install dependencies: Rprobit is installed from the tarball in the project
### root; all other packages come from CRAN
if (isTRUE(install_dep)) {
  repos <- "https://cloud.r-project.org"
  install_missing <- function(pkgs) {
    pkgs <- trimws(sub("[(].*", "", pkgs))
    pkgs <- setdiff(pkgs, c("R", "Rprobit", rownames(installed.packages(priority = "base"))))
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0) install.packages(missing, repos = repos)
  }
  tarball <- list.files(pattern = "^Rprobit_.*\\.tar\\.gz$")
  if (length(tarball) != 1L) {
    stop("Expected exactly one Rprobit tarball in the project root.")
  }
  tmp <- tempfile()
  untar(tarball, files = "Rprobit/DESCRIPTION", exdir = tmp)
  desc <- read.dcf(file.path(tmp, "Rprobit", "DESCRIPTION"))
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(desc))
  deps <- unlist(strsplit(gsub("\n", " ", paste(desc[, fields], collapse = ",")), ","))
  install_missing(deps[nzchar(trimws(deps))])
  install_missing(c("checkmate", "future", "progressr"))
  install.packages(tarball, repos = NULL, type = "source")
}

### define simulation setup
### P[1,2,3][f,r]: number of fixed/random type-1 (type-2, type-3) effects
### ASC: either "no", "f" (fixed), or "r" (random) ASCs (in addition to P2f/r)
### J: number of choice alternatives
### N: number of deciders
### Tp: number of choice occasions per decider
### cml_type: "probit" for full probit likelihood, "FP" for full pairwise CML
probit_mle <- function(P1f, P1r, P2f, P2r, P3f, P3r, ASC, J, N, Tp, cml_type) {

  ### input checks
  checkmate::assert_count(P1f)
  checkmate::assert_count(P1r)
  checkmate::assert_count(P2f)
  checkmate::assert_count(P2r)
  checkmate::assert_count(P3f)
  checkmate::assert_count(P3r)
  checkmate::assert_choice(ASC, c("no", "f", "r"))
  checkmate::assert_int(J, lower = 2)
  checkmate::assert_count(N, positive = TRUE)
  checkmate::assert_count(Tp, positive = TRUE)
  checkmate::assert_choice(cml_type, c("probit", "FP"))

  ### build formula, re, and alt_names
  x <- function(i, f, r) {
    if (f + r) paste0(letters[seq_len(f + r)], i) else character()
  }
  p <- function(x) if (length(x)) paste(x, collapse = "+") else "0"
  x1 <- x(1, P1f, P1r)
  x2 <- x(2, P2f, P2r)
  x3 <- x(3, P3f, P3r)
  b <- if (length(x2)) {
    paste0(p(x2), ifelse(ASC == "no", "+0", ""))
  } else {
    if (ASC == "no") "0" else "1"
  }
  formula <- stats::as.formula(paste("choice ~", p(x1), "|", b, "|", p(x3)))
  re <- c(tail(x1, P1r), tail(x2, P2r), tail(x3, P3r), if(ASC == "r") "ASC")
  alt_names <- LETTERS[seq_len(J)]

  ### setup Rprobit object
  r_sim <- Rprobit::setup(
    formula = formula,
    df = NULL,
    re = re,
    alt_names = alt_names,
    corr_re = TRUE,
    error_struc = "full",
    N = N,
    Tp = Tp
  )
  theta_dgp <- r_sim$choice_data$theta_dgp

  ### nlm setting
  iterlim <- 1000L
  gradtol <- 1e-8
  steptol <- 1e-8

  ### estimation with random initialization
  r_sim_random <- r_sim$clone(deep = TRUE)
  theta_init_random <- rnorm(r_sim_random$mod$length_theta)
  r_sim_random$fit(
    theta_init = theta_init_random, cml_type = cml_type, verbose = FALSE,
    iterlim = iterlim, gradtol = gradtol, steptol = steptol
  )
  theta_init_random <- r_sim_random$theta_init
  theta_fit_random <- r_sim_random$theta_fit
  ll_init_random <- r_sim_random$optimization_details$ll_init
  ll_fit_random <- r_sim_random$ll_fit
  sec_fit_random <- r_sim_random$optimization_details$sec_fit
  nlm_code_random <- r_sim_random$optimization_details$code
  nlm_iterations_random <- r_sim_random$optimization_details$iterations

  ### estimation with smart initialization
  r_sim_smart <- r_sim$clone(deep = TRUE)
  r_sim_smart$fit(
    theta_init = NULL, cml_type = cml_type, verbose = FALSE,
    iterlim = iterlim, gradtol = gradtol, steptol = steptol
  )
  theta_init_smart <- r_sim_smart$theta_init
  theta_fit_smart <- r_sim_smart$theta_fit
  ll_init_smart <- r_sim_smart$optimization_details$ll_init
  ll_fit_smart <- r_sim_smart$ll_fit
  sec_init_smart <- r_sim_smart$optimization_details$sec_init
  sec_fit_smart <- r_sim_smart$optimization_details$sec_fit
  nlm_code_smart <- r_sim_smart$optimization_details$code
  nlm_iterations_smart <- r_sim_smart$optimization_details$iterations
  init_type_smart <- r_sim_smart$optimization_details$init_type

  ### estimation initialized at the true DGP parameters
  r_sim_true <- r_sim$clone(deep = TRUE)
  r_sim_true$fit(
    theta_init = theta_dgp, cml_type = cml_type, verbose = FALSE,
    iterlim = iterlim, gradtol = gradtol, steptol = steptol
  )
  theta_init_true <- r_sim_true$theta_init
  theta_fit_true <- r_sim_true$theta_fit
  ll_init_true <- r_sim_true$optimization_details$ll_init
  ll_fit_true <- r_sim_true$ll_fit
  sec_fit_true <- r_sim_true$optimization_details$sec_fit
  nlm_code_true <- r_sim_true$optimization_details$code
  nlm_iterations_true <- r_sim_true$optimization_details$iterations

  ### return results
  list(
    "theta_dgp" = theta_dgp,
    "theta_init_random" = theta_init_random,
    "theta_init_smart" = theta_init_smart,
    "theta_init_true" = theta_init_true,
    "theta_fit_random" = theta_fit_random,
    "theta_fit_smart" = theta_fit_smart,
    "theta_fit_true" = theta_fit_true,
    "ll_init_random" = ll_init_random,
    "ll_init_smart" = ll_init_smart,
    "ll_init_true" = ll_init_true,
    "ll_fit_random" = ll_fit_random,
    "ll_fit_smart" = ll_fit_smart,
    "ll_fit_true" = ll_fit_true,
    "sec_init_smart" = sec_init_smart,
    "sec_fit_random" = sec_fit_random,
    "sec_fit_smart" = sec_fit_smart,
    "sec_fit_true" = sec_fit_true,
    "nlm_code_random" = nlm_code_random,
    "nlm_code_smart" = nlm_code_smart,
    "nlm_code_true" = nlm_code_true,
    "nlm_iterations_random" = nlm_iterations_random,
    "nlm_iterations_smart" = nlm_iterations_smart,
    "nlm_iterations_true" = nlm_iterations_true,
    "init_type_smart" = init_type_smart
  )
}

### define simulation settings
settings <- list(
  list(
    "label" = "type_1",
    "args" = list(
      "P1f" = list(10), "P1r" = list(0),
      "P2f" = list(0), "P2r" = list(0),
      "P3f" = list(0), "P3r" = list(0),
      "ASC" = list("no"), "J" = list(5),
      "N" = list(5000), "Tp" = list(1),
      "cml_type" = list("probit")
    )
  ),
  list(
    "label" = "type_2",
    "args" = list(
      "P1f" = list(0), "P1r" = list(0),
      "P2f" = list(5), "P2r" = list(0),
      "P3f" = list(0), "P3r" = list(0),
      "ASC" = list("no"), "J" = list(5),
      "N" = list(5000), "Tp" = list(1),
      "cml_type" = list("probit")
    )
  ),
  list(
    "label" = "type_3",
    "args" = list(
      "P1f" = list(0), "P1r" = list(0),
      "P2f" = list(0), "P2r" = list(0),
      "P3f" = list(5), "P3r" = list(0),
      "ASC" = list("no"), "J" = list(5),
      "N" = list(5000), "Tp" = list(1),
      "cml_type" = list("probit")
    )
  ),
  list(
    "label" = "mixed_fixed",
    "args" = list(
      "P1f" = list(3), "P1r" = list(0),
      "P2f" = list(3), "P2r" = list(0),
      "P3f" = list(3), "P3r" = list(0),
      "ASC" = list("no"), "J" = list(5),
      "N" = list(5000), "Tp" = list(1),
      "cml_type" = list("probit")
    )
  ),
  list(
    "label" = "random_type_1",
    "args" = list(
      "P1f" = list(0), "P1r" = list(5),
      "P2f" = list(0), "P2r" = list(0),
      "P3f" = list(0), "P3r" = list(0),
      "ASC" = list("no"), "J" = list(5),
      "N" = list(200), "Tp" = list(10),
      "cml_type" = list("FP")
    )
  ),
  list(
    "label" = "mixed_random",
    "args" = list(
      "P1f" = list(1), "P1r" = list(1),
      "P2f" = list(1), "P2r" = list(1),
      "P3f" = list(1), "P3r" = list(1),
      "ASC" = list("no"), "J" = list(3),
      "N" = list(200), "Tp" = list(10),
      "cml_type" = list("FP")
    )
  )
)
if (is.null(settings_to_run)) settings_to_run <- seq_along(settings)
checkmate::assert_integerish(
  settings_to_run, lower = 1L, upper = length(settings), unique = TRUE
)
checkmate::assert_count(runs, positive = TRUE)

### setup parallel framework
if (is.null(workers)) workers <- max(1, future::availableCores() - 1)
future::plan(future::multicore, workers = workers)

### setup progress updates
progressr::handlers(global = TRUE)
progressr::handlers(
  progressr::handler_progress(format = ">> :percent, :eta to go :message")
)

### save system information (identical for all settings, since one machine)
path <- file.path("sim", "system_information.txt")
sink(path, split = FALSE)
cat("date/time:")
print(Sys.time())
cat("\nsettings:", settings_to_run, "\n")
cat("\nR session information:\n")
print(sessionInfo())
cat("\npackage versions:\n")
version_of <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(packageVersion(pkg))
  } else {
    NA_character_
  }
}
print(vapply(c("Rprobit", "RprobitB", "oeli"), version_of, character(1)))
cat("\nlibrary paths:\n")
print(.libPaths())
cat("\nparallel computation:\n")
print(future::plan())
sink()

### run simulations, one setting after the other
for (setting in settings_to_run) {
  setting_label <- settings[[setting]]$label
  setting_args <- settings[[setting]]$args
  setting_path <- file.path("sim", setting_label)
  if (dir.exists(setting_path)) {
    message("Skipping setting '", setting_label, "', output directory exists.")
    next
  }
  message("Starting setting '", setting_label, "' at ", format(Sys.time()))
  time_start <- Sys.time()
  Sim <- oeli::Simulator$new(verbose = TRUE)
  do.call(Sim$define, c(list("f" = probit_mle), setting_args))
  Sim$go(runs = runs, backup = TRUE, path = setting_path)
  cat(
    sprintf(
      "\nsetting '%s': %d runs finished at %s after %.1f minutes\n",
      setting_label, runs, format(Sys.time()),
      as.numeric(difftime(Sys.time(), time_start, units = "mins"))
    ),
    file = path, append = TRUE
  )
}
