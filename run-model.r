library(rstan)
library(data.table)
library(optparse)

# Commandline options and parsing
parser <- OptionParser()
parser <- add_option(parser, c("-D", "--debug"), action="store_true",
                     help="Perform a debug run of the model")
parser <- add_option(parser, c("-F", "--full"), action="store_true",
                     help="Perform a full run of the model")
parser <- add_option(parser, c("-S", "--fixseed"), action="store_true",
                     help="Fix the RNG seed during sampling to 1")
parser <- add_option(parser, c("-I", "--scaleifr"), type="double", default=1,
                     help="Scale the ifr values")
parser <- add_option(parser, c("-G", "--gradeffect"), type="double", default=Inf,
                     help="Gradual full effectiveness of interventions")
parser <- add_option(parser, c("-P", "--forecast"), type="integer", default=0,
                     help="Number of days to forecast")
cmdoptions <- parse_args(parser, args=commandArgs(trailingOnly = TRUE),
                         positional_arguments=TRUE)

# Default run parameters for the model
DEBUG <- !is.null(cmdoptions$options$debug)
FULL <- !is.null(cmdoptions$options$full)
if (DEBUG && FULL) {
  stop("Setting both debug and full run modes at once is invalid")
}
seed <- if (is.null(cmdoptions$options$fixseed))
            sample.int(.Machine$integer.max, 1) else 1

if (length(cmdoptions$args) == 0) {
  stop("Model name must be specified")
}
StanModel <- cmdoptions$args[1]

print(sprintf("Running %s", StanModel))
if (DEBUG) {
  print("Running in DEBUG mode")
} else if (FULL) {
  print("Running in FULL mode")
}
print(sprintf("Setting random seed to %d", seed))

source("utils/process-covariates.r")
source("utils/covariate-size-effects.r")
source("utils/plot-2-panel.r")
source("utils/plot-3-panel.r")
source("utils/utils.r")

# Read which countires to use
countries <- readRDS('data/regions.rds')

# Read deaths data for regions
d <- readRDS('data/COVID-19-up-to-date.rds')

# Read IFR and pop by country
ifr.by.country <- readRDS('data/popt-ifr.rds')
if (cmdoptions$options$scaleifr != 1) {
  print(sprintf("Scaling IFR by %g", cmdoptions$options$scaleifr))
  ifr.by.country$ifr <- ifr.by.country$ifr * cmdoptions$options$scaleifr
}

# Read interventions
interventions <- readRDS('data/interventions.rds')

# Maximum number of days to simulate
N2 <- (max(d$DateRep) - min(d$DateRep) + 1)[[1]]

processed_data <- process_covariates(countries=countries,
                                     interventions=interventions,
                                     d=d, ifr.by.country=ifr.by.country, N2=N2)

forecast <- cmdoptions$options$forecast
if (forecast > 0) {
  ## when forecasting, the stan data shouldn't contain the last `forecast` days,
  ## but processed_data contains all dates, so that we can compare observed
  ## and predicted deaths
  data_forecast <- process_covariates(countries=countries,
                                      interventions=interventions,
                                      d=d[d$DateRep <= tail(d$DateRep, 1) - forecast, ],
                                      ifr.by.country=ifr.by.country, N2=N2)
  processed_data$stan_data <- data_forecast$stan_data
}
stan_data <- processed_data$stan_data

if (grepl("nolastintervention", StanModel)) {
  ## data without the last intervention column for Sweden
  stan_data$X <- stan_data$X[, , -7]
  stan_data$P <- 6
}

## set the name of the output results
output.name <- StanModel
if (cmdoptions$options$scaleifr != 1)
  output.name <- paste0(output.name, "_scaleifr", cmdoptions$options$scaleifr)
if (!is.infinite(cmdoptions$options$gradeffect)) {
  print(sprintf("Setting gradual delay of effects to %g", cmdoptions$options$gradeffect))
  stan_data <- gradual.npi.effect(stan_data, cmdoptions$options$gradeffect)
  output.name <- paste0(output.name, "_gradeffect", cmdoptions$options$gradeffect)
}
if (forecast > 0) {
  print(sprintf("Forecasting last %d days", forecast))
  output.name <- paste0(output.name, "_forecast", forecast)
}

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

print("Compiling Stan model")
model <- stan_model(paste0("stan-models/", StanModel, ".stan"))

if (DEBUG) {
  iter <- 40
  warmup <- 20
  chains <- 2
  control <- list(adapt_delta=0.95, max_treedepth=11)
} else if (FULL) {
  iter <- 1800
  warmup <- 1000
  chains <- 5
  control <- list(adapt_delta=0.99, max_treedepth=20)
} else {
  iter <- 600
  warmup <- 300
  chains <- 4
  control <- list(adapt_delta=0.95, max_treedepth=11)
}
fit <- sampling(model, data=stan_data, seed=seed, refresh=iter / 20,
                iter=iter, warmup=warmup, chains=chains, control=control)

save.results(fit, processed_data, output.name)
write.csv(fit.summary(fit, processed_data, stan_data, forecast),
          file=paste0("results/", output.name, "-summary.csv"))
write.csv(sampler.stats(fit),
          file=paste0("results/", output.name, "-sampler.csv"))
plot.mu.rt.rhat(output.name)
plot_covars(output.name)
make_two_panel_plot(output.name)
make_three_panel_plot(output.name)
