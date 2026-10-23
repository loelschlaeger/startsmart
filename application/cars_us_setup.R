### Application 2: stated vehicle choices of US car buyers
###
### Fits a normally mixed multinomial probit model with correlated random
### coefficients to the US conjoint data of Helveston et al. (2015) from
### proposed and from random initial values:
###
###   source("application/cars_us_setup.R")
###
### The following variables can be defined before sourcing:
###
###   install_dep <- TRUE   # install package dependencies (default TRUE)
###   runs <- 100L          # random initializations (default 100)
###   workers <- NULL       # parallel workers (default: cores - 1)
###
### The data are the `cars_us` data set of the logitr package. The working
### directory must be the project root, which must contain the Rprobit tarball.
### Results are saved as case_*.rds in application/cars_us_proposed and
### application/cars_us_random, which must not already exist, because backups
### are not overwritten. System information and the run time are written to
### application/system_information_cars_us.txt.

if (!exists("install_dep")) install_dep <- TRUE
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
  install_missing(c("checkmate", "future", "ggplot2", "logitr", "oeli", "progressr"))
  install.packages(tarball, repos = NULL, type = "source")
}

checkmate::assert_count(runs, positive = TRUE)

### setup parallel framework
if (is.null(workers)) workers <- max(1, future::availableCores() - 1)
future::plan(future::multicore, workers = workers)

### setup progress updates
progressr::handlers(global = TRUE)
progressr::handlers(
  progressr::handler_progress(format = ">> :percent, :eta to go :message")
)

### define the application
### init_type: "proposed" for the initialization procedure, "random" for
###   initial values drawn from a standard normal distribution
cars_us_mle <- function(init_type) {

  ### input checks
  checkmate::assert_choice(init_type, c("proposed", "random"))

  ### choice data: 384 US respondents each chose 15 times between three vehicle
  ### profiles
  utils::data("cars_us", package = "logitr", envir = environment())
  raw <- as.data.frame(cars_us)
  raw <- raw[order(raw$id, raw$obsnum), ]
  raw$alt <- stats::ave(seq_len(nrow(raw)), raw$id, raw$obsnum, FUN = seq_along)
  attributes <- c(
    "hev", "phev10", "phev20", "phev40", "bev75", "bev100", "bev150",
    "phevFastcharge", "bevFastcharge", "price", "opCost", "accelTime",
    "american", "japanese", "chinese", "skorean"
  )
  df <- unique(raw[c("id", "obsnum")])
  df$occasion <- stats::ave(seq_len(nrow(df)), df$id, FUN = seq_along)
  key <- paste(df$id, df$obsnum)
  for (a in 1:3) {
    rows <- raw[raw$alt == a, ]
    k <- match(key, paste(rows$id, rows$obsnum))
    for (variable in attributes) {
      df[[paste0(variable, "_A", a)]] <- as.numeric(rows[[variable]][k])
    }
  }
  chosen_rows <- raw[raw$choice == 1, ]
  df$choice <- paste0("A", chosen_rows$alt[match(key, paste(chosen_rows$id, chosen_rows$obsnum))])

  ### scale the regressors
  ### the price (in thousands of dollars), the operating cost (in cents per
  ### mile), and the acceleration time (in seconds) are divided by their
  ### standard deviations; the powertrain, charging, and brand indicators are
  ### left unscaled
  standardize <- function(columns) {
    values <- unlist(df[columns], use.names = FALSE)
    deviation <- stats::sd(values)
    if (is.finite(deviation) && deviation > 0) deviation else 1
  }
  for (variable in c("price", "opCost", "accelTime")) {
    columns <- paste0(variable, "_A", 1:3)
    df[columns] <- df[columns] / standardize(columns)
  }

  ### setup Rprobit object
  ### the profiles are unlabelled, so the model has no alternative-specific
  ### constants; the coefficients of price, operating cost, and acceleration
  ### time are correlated Gaussian random effects
  r <- Rprobit::setup(
    formula = choice ~ hev + phev10 + phev20 + phev40 + bev75 + bev100 + bev150 +
      phevFastcharge + bevFastcharge + price + opCost + accelTime +
      american + japanese + chinese + skorean | 0,
    df = df,
    re = c("price", "opCost", "accelTime"),
    column_choice = "choice",
    column_decider_id = "id",
    column_occasion_id = "occasion",
    corr_re = TRUE,
    error_struc = "full"
  )

  ### nlm setting
  iterlim <- 1000L
  gradtol <- 1e-8
  steptol <- 1e-8

  ### initial values
  theta_init <- if (init_type == "random") {
    stats::rnorm(r$mod$length_theta)
  } else {
    NULL
  }

  ### estimation by the full pairwise composite marginal likelihood, because
  ### the panel dimension makes the full likelihood expensive
  r$fit(
    theta_init = theta_init, cml_type = "FP", verbose = FALSE,
    iterlim = iterlim, gradtol = gradtol, steptol = steptol
  )
  details <- r$optimization_details

  ### return results
  list(
    "init_type" = init_type,
    "parameter_names" = c(
      sub("^b\\[[0-9]+\\]:", "", rownames(r$mod$H_b)),
      if (r$mod$length_theta_O > 0) paste0("O", seq_len(r$mod$length_theta_O)),
      if (r$mod$length_theta_L > 0) paste0("L", seq_len(r$mod$length_theta_L))
    ),
    "theta_init" = r$theta_init,
    "theta_fit" = r$theta_fit,
    "ll_init" = details$ll_init,
    "ll_fit" = r$ll_fit,
    "sec_init" = details$sec_init,
    "sec_fit" = details$sec_fit,
    "nlm_code" = details$code,
    "nlm_iterations" = details$iterations
  )
}

### save system information
path <- file.path("application", "system_information_cars_us.txt")
sink(path, split = FALSE)
cat("application: cars_us\n")
cat("date/time:")
print(Sys.time())
cat("\nrandom initializations:", runs, "\n")
cat("\nR session information:\n")
print(sessionInfo())
cat("\npackage versions:\n")
version_of <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(packageVersion(pkg)) else NA_character_
}
print(vapply(c("Rprobit", "RprobitB", "oeli", "logitr"), version_of, character(1)))
cat("\nlibrary paths:\n")
print(.libPaths())
cat("\nparallel computation:\n")
print(future::plan())
sink()

### run the application
### the proposed initialization is computed once, the random initialization
### `runs` times; the two arms are stored in separate backups
time_start <- Sys.time()
for (init in c("proposed", "random")) {
  Sim <- oeli::Simulator$new(verbose = TRUE)
  Sim$define(f = cars_us_mle, init_type = list(init))
  Sim$go(
    runs = if (init == "proposed") 1L else runs, backup = TRUE,
    path = file.path("application", paste0("cars_us_", init))
  )
}
cat(
  sprintf(
    "\n1 proposed and %d random initializations finished at %s after %.1f minutes\n",
    runs, format(Sys.time()),
    as.numeric(difftime(Sys.time(), time_start, units = "mins"))
  ),
  file = path, append = TRUE
)
