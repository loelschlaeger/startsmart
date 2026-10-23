# Application setup

This directory compares initialization strategies for probit model estimation
on an empirical data set, using the **Rprobit** package.

Two initializations are compared:

1. **Random initialization:** the initial values are drawn with `rnorm()`.
   This is repeated `runs` times.
2. **Proposed initialization:** the initialization procedure is used. This is
   run once.

## Application: cars_us

The `cars_us` data are the US part of the conjoint survey of Helveston et al.
(2015), *Transportation Research Part A* 73, 96-112: 384 respondents each chose
15 times between three vehicle options described by powertrain (conventional,
hybrid, plug-in hybrid with 10, 20, or 40 miles of electric range, battery
electric with 75, 100, or 150 miles of range), fast-charging capability, price
in thousands of dollars, operating cost in cents per mile, acceleration time
in seconds, and brand origin (German, American, Japanese, Chinese, South
Korean). The data are distributed with the **logitr** package.

The profiles are unlabelled, so the position of a profile in a choice task
serves as the alternative identifier and the model has no alternative-specific
constants. All attributes enter with generic coefficients (type 1). The
coefficients of price, operating cost, and acceleration time carry correlated
Gaussian random effects (`re = c("price", "opCost", "accelTime")`,
`corr_re = TRUE`). The error covariance matrix is flexible and the parameter 
vector has 24 entries. Because the panel dimension makes the full likelihood 
expensive, the full pairwise composite marginal likelihood is maximized 
(`cml_type = "FP"`).

## Scaling

Dividing by the standard deviation puts the parameters in the range of a draw
from the standard normal distribution and is therefore the numerically most
favourable case for random initialization. This is deliberate: the proposed
initialization is compared against a random draw on the terms of the latter.

## Start the application

From the project root, which must contain the Rprobit tarball:

```r
source("application/cars_us_setup.R")
```

The following variables can be defined before sourcing:

```r
install_dep <- TRUE   # install package dependencies (default TRUE)
runs <- 100L          # random initializations (default 100)
workers <- NULL       # parallel workers (default: cores - 1)
```

Package dependencies are installed from CRAN, except **Rprobit**, which is
installed from the tarball `Rprobit_<version>.tar.gz` in the project root.

The output directories `application/cars_us_proposed` and
`application/cars_us_random` must not already exist, because backups are not
overwritten. Parallel execution uses `future::multicore` and therefore requires
a Unix-like operating system.

## Output

Successful runs are saved as `case_*.rds` in `application/cars_us_proposed`
and `application/cars_us_random`. Each result contains the initialization
type, the parameter names, the initial and the fitted parameter vector, the
log-likelihoods at initialization and after fitting, the initialization and
optimization times, and the `nlm` convergence code and iteration count. System
information and the total run time are written to
`application/system_information_cars_us.txt`.

## Start simulation on AWS instance

Below, replace `IP` by the public address of your instance and `KEY` by the path 
to your private key file.

### 1. Upload the application files

Open the **AWS Management Console** and go to **EC2**. Click **Instances** in
the left-hand menu. Select your running EC2 instance. In the **Details**
section, look for **Public IPv4 address**. Copy this IP address. The remote
paths below are relative to the home directory of the `ubuntu` user. Run:

```
ssh -i KEY ubuntu@IP "mkdir -p startsmart/application"
scp -i KEY application/cars_us_setup.R ubuntu@IP:startsmart/application/
scp -i KEY Rprobit_1.0.0.tar.gz ubuntu@IP:startsmart/
```

### 2. Install R

Connect to the instance:

```
ssh -i KEY ubuntu@IP
```

On the instance, run:

```bash
sudo apt update
sudo apt install -y r-base r-base-dev build-essential gfortran cmake pkg-config libnlopt-dev libcurl4-openssl-dev libssl-dev libxml2-dev libuv1-dev imagemagick libmagick++-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev libcairo2-dev libjpeg-dev libpng-dev libtiff-dev
```

Create a library for the packages:

```bash
mkdir -p /home/ubuntu/R/library
echo 'R_LIBS_USER=/home/ubuntu/R/library' >> ~/.Renviron
```

### 3. Run the application

```bash
cd ~/startsmart
nohup Rscript -e 'source("application/cars_us_setup.R")' > cars_us.log 2>&1 &
```

It will continue running even if the SSH connection is closed. To follow the
progress:

```bash
tail -f cars_us.log
```

### 4. Download the results

Exit the instance:

```bash
exit
```

Run:

```
scp -r -i KEY ubuntu@IP:startsmart/application application_new
```
