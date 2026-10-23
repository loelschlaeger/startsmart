# Simulation setup

This simulation compares initialization strategies for probit model estimation 
using the **Rprobit** package:

1. **Random initialization:** the initial values are drawn with `rnorm()`.
2. **Smart initialization:** the initialization procedure is used.
3. **True initialization:** the true parameter vector `theta_dgp` is used.

## Parameters

- `P1f/P1r`, `P2f/P2r`, `P3f/P3r`: numbers of fixed/random type-1/2/3 effects,
- `ASC`: no, fixed, or random alternative-specific constants;
- `J`: alternatives;
- `N`: decision makers;
- `Tp`: choice occasions per decision maker;
- `cml_type`: full probit likelihood (`"probit"`) or full pairwise composite
  likelihood (`"FP"`).
  
## Model details

Random effects are Gaussian and correlated (`corr_re = TRUE`), and the error
covariance $\Sigma$ is flexible (`error_struc = "full"`). The first alternative 
is the reference alternative, so its type-2 coefficients, including the 
ASC, are fixed to zero. For $\Sigma$, the first row and column are fixed to 
zero, which expresses utilities relative to the reference alternative, 
and $\Sigma_{2,2}=1$ fixes the utility scale.

## Scenarios

Each scenario runs 100 times:

| `setting` | Label | Effects | `length(theta)` | `J` | `N` | `Tp` | CML type |
|---:|---|---|---:|---:|---:|---:|---|
| 1 | `type_1` | 10 fixed type-1 | 19 | 5 | 5000 | 1 | `probit` |
| 2 | `type_2` | 5 fixed type-2 | 29 | 5 | 5000 | 1 | `probit` |
| 3 | `type_3` | 5 fixed type-3 | 34 | 5 | 5000 | 1 | `probit` |
| 4 | `mixed_fixed` | 3 fixed effects of each type | 39 | 5 | 5000 | 1 | `probit` |
| 5 | `random_type_1` | 5 random type-1 | 29 | 5 | 200 | 10 | `FP` |
| 6 | `mixed_random` | 1 fixed + 1 random effect of each type | 35 | 3 | 200 | 10 | `FP` |

## Start simulation

All scenarios are run one after the other on the same machine with a single
call from the project root:

```r
source("sim/simulation_setup.R")
```

The following variables can be defined before sourcing:

```r
install_dep <- TRUE      # install package dependencies (default TRUE)
settings_to_run <- 1:6   # subset of settings (default: all)
runs <- 100L             # replications per setting (default 100)
workers <- NULL          # parallel workers (default: cores - 1)
```

Package dependencies are installed from CRAN, except **Rprobit**, which is
installed from the tarball `Rprobit_<version>.tar.gz` in the project root. 

A scenario whose output directory `sim/<label>` already exists is skipped, so
that an interrupted run can be resumed. Existing directories from an earlier run
must therefore be removed before recomputing. Parallel execution uses
`future::multicore` and therefore requires a Unix-like operating system.

## Output

Successful runs are saved as `sim/<scenario-label>/case_*.rds`, for example:

```r
case <- readRDS("sim/type_1/case_001.rds")
result <- case$.output
```

Each result contains:

- DGP, initial, and fitted parameter vectors for all three initializations;
- log-likelihoods at initialization and after fitting;
- initialization and optimization times;
- `nlm` convergence codes and iteration counts;
- `init_type_smart`, indicating automatic initialization or unit fallback;
- fitted `Rprobit` objects with suffixes `_random`, `_smart`, and `_true`.

System information, package versions, and the run time of each scenario are
written to `sim/system_information.txt`.

## Start simulation on AWS instance

Below, replace `IP` by the public address of your instance and `KEY` by the path 
to your private key file.

### 1. Upload simulation files

Open the **AWS Management Console** and go to **EC2**. Click **Instances** in
the left-hand menu. Select your running EC2 instance. In the **Details**
section, look for **Public IPv4 address**. Copy this IP address.

The simulation script expects the project layout `startsmart/sim/simulation_setup.R`
and `startsmart/Rprobit_<version>.tar.gz`. The remote paths below are relative
to the home directory of the `ubuntu` user. Run:

```
ssh -i KEY ubuntu@IP "mkdir -p startsmart/sim"
scp -i KEY sim/simulation_setup.R ubuntu@IP:startsmart/sim/
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

### 3. Run simulation script

Start the simulation:

```bash
cd ~/startsmart
nohup Rscript -e 'source("sim/simulation_setup.R")' > simulation.log 2>&1 &
```

The simulation will continue running even if the SSH connection is closed. To
follow the simulation progress:

```bash
tail -f simulation.log
```

### 4. Download the results

Exit the instance:

```bash
exit
```

Run:

```
scp -r -i KEY ubuntu@IP:startsmart/sim sim_new
```
