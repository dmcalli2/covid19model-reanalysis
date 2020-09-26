library(bayesplot)
library(ggplot2)

## allow the effect of each NPI to come in gradually
gradual.npi.effect <- function(stan_data, lambda=0.7) {
  grad.effect <- round(1 - exp(-lambda * 1:20), 2)

  ## loop over countries
  for (m in 1:nrow(stan_data$X)) {
    first.day <- apply(stan_data$X[m, , ], 2, function(z) which(z == 1)[1])
    for (p in 1:length(first.day)) {
      if (is.na(first.day[p])) next
      stan_data$X[m, first.day[p] + 1:length(grad.effect) - 1, p] <- grad.effect
    }
  }

  return(stan_data)
}

## save results data
save.results <- function(fit, processed_data, modelname) {
  out <- rstan::extract(fit)
  prediction <- out$prediction
  estimated.deaths <- out$E_deaths
  estimated.deaths.cf <- out$E_deaths0

  stan_data <- processed_data$stan_data
  dates <- processed_data$dates
  deaths_by_country <- processed_data$deaths_by_country
  reported_cases <- processed_data$reported_cases
  countries <- countries$Regions

  save(fit, prediction, dates, reported_cases, deaths_by_country, countries,
       estimated.deaths, estimated.deaths.cf, stan_data, processed_data,
       file=paste0("results/", modelname, "-stanfit.Rdata"))
}

## deaths/counterfactual and MSE
get.deaths.mse <- function(fit, processed_data, forecast=0) {
  cases.pr <- rstan::extract(fit, "prediction")[[1]]
  deaths.pr <- rstan::extract(fit, "E_deaths")[[1]]
  deaths.cf <- rstan::extract(fit, "E_deaths0")[[1]]
  deaths.ob <- processed_data$deaths_by_country
  dates <- processed_data$dates
  S <- nrow(deaths.pr) # number of samples

  ## loop over countries
  cum.deaths.mse <- NULL
  for (m in 1:length(dates)) {
    ## consider only deaths within the modelled period
    N <- length(dates[[m]])
    if (forecast > 0)
      per <- tail(1:N, forecast)
    else
      per <- 1:N
    mse <- apply(deaths.pr[, , m], 1,
                 function(z) mean((z[per] - deaths.ob[[m]][per])^2))
    cum.deaths.mse <- rbind(cum.deaths.mse,
                            c(sum(colMeans(deaths.pr[, 1:N, m])),
                              sum(colMeans(deaths.cf[, 1:N, m])),
                              sum(colMeans(cases.pr[, 1:N, m])),
                              mean(mse)))

  }

  cum.deaths.mse <- as.data.frame(cum.deaths.mse)
  rownames(cum.deaths.mse) <- names(dates)
  colnames(cum.deaths.mse) <- c("deaths", "counterfactual", "infections", "mse")
  cum.deaths.mse[, 1:3] <- round(cum.deaths.mse[, 1:3])
  cum.deaths.mse$mse <- round(cum.deaths.mse$mse, 1)
  return(cum.deaths.mse)
}

## compute ranges of herd immunity threshold from posterior means of parameters
hit <- function(fit) {
  R0 <- rstan::extract(fit, "mu")[[1]]
  alpha <- tryCatch(matrix(rstan::extract(fit, "a_het")[[1]],
                           nrow(R0), ncol(R0)),
                    error=function(e) Inf)
  hit <- 1 - (1 / R0)^(1 - 1 / (1 + alpha))
  list(hit=colMeans(hit), R0=colMeans(R0),
       range.hit=range(colMeans(hit)),
       range.R0=range(colMeans(R0)),
       alpha=mean(alpha))
}

## create some plots
plot.mu.rt.rhat <- function(modelname) {
  filename <- paste0("results/", modelname, "-stanfit.Rdata")
  print(sprintf("loading: %s", filename))
  load(filename)
  print('Generating mu, rt plots')
  out <- rstan::extract(fit)
  mu <- as.matrix(out$mu)
  colnames(mu) <- countries
  g <- mcmc_intervals(mu, prob=0.9)
  ggsave(sprintf("figures/%s-mu.png", modelname), g, width=4, height=6)
  tmp <- lapply(1:length(countries), function(i) (out$Rt_adj[, stan_data$N[i], i]))
  Rt_adj <- do.call(cbind, tmp)
  colnames(Rt_adj) <- countries
  g <- mcmc_intervals(Rt_adj, prob=0.9)
  ggsave(sprintf("figures/%s-final-rt.png", modelname), g, width=4, height=6)

  print("Generating rhat plot")
  g <- mcmc_rhat_hist(rhat(fit))
  ggsave(sprintf("figures/%s-rhat.png", modelname), g, width=4, height=6)
}

## export Stan's neg_binomial_2_lpmf function as negbim2
cat("Exporting Stan's neg_binomial_2_lpmf\n")
rstan::expose_stan_functions(rstan::stanc(model_code='
functions {
  real negbim2(int[] n, real[] mu, real phi) {
    return neg_binomial_2_lpmf(n | mu, phi);
  }
}
model{}
'))

dic <- function(fit, stan_data) {

  ## total likelihood
  lik <- rowSums(rstan::extract(fit, "log_lik")[[1]])

  ## mean deviance
  d1 <- -2 * mean(lik)

  ## deviance of the mean
  d2 <- 0
  for (m in 1:stan_data$M) {
    per <- stan_data$EpidemicStart[m]:stan_data$N[m]
    deaths.ob <- stan_data$deaths[per, m]
    deaths.pr <- rstan::extract(fit, "E_deaths")[[1]][, per, m]
    phi <- rstan::extract(fit, "phi")[[1]]
    d2 <- d2 + negbim2(deaths.ob, colMeans(deaths.pr), mean(phi))
  }
  d2 <- -2 * d2

  return(round(c(pD=d1 - d2, pV=mean(var(-2 * lik)) / 2, dic=2 * d1 - d2, deviance=d1), 1))
}

waic <- function(fit) {
  loo::waic(rstan::extract(fit, "log_lik")[[1]])
}

## collect summary measures of fit
fit.summary <- function(fit, processed_data, stan_data, forecast) {
  dic.model <- dic(fit, stan_data)
  hit.model <- hit(fit)
  waic.model <- waic(fit)
  loo.model <- loo(fit)
  mm <- get.deaths.mse(fit, processed_data, forecast)
  mm.sums <- colSums(mm)
  data.frame(
    Deviance=dic.model["deviance"],
    pD=dic.model["pD"],
    DIC=dic.model["dic"],
    R0=round(mean(hit.model$R0), 2),
    alpha=gsub("Inf", "\U221E", round(mean(hit.model$alpha), 2)),
    H=round(mean(hit.model$hit), 2),
    Infections=mm.sums["infections"],
    UK.Infections=mm["United_Kingdom", "infections"],
    Deaths=mm.sums["deaths"],
    Counterfactual=mm.sums["counterfactual"],
    row.names=deparse(substitute(fit)), stringsAsFactors=FALSE)
}

## collect statistics on the sampler behaviour
sampler.stats <- function(fit) {
    sp <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
    accept.stat <- sapply(sp, function(x) mean(x[, "accept_stat__"]))
    stepsize <- sapply(sp, function(x) mean(x[, "stepsize__"]))
    divergences <- sapply(sp, function(x) sum(x[, "divergent__"]))
    treedepth <- sapply(sp, function(x) max(x[, "treedepth__"]))
    gradients <- sapply(sp, function(z) sum(z[, "n_leapfrog__"]))
    et <- round(rstan::get_elapsed_time(fit), 2)
    res <- cbind(accept.stat, stepsize, divergences, treedepth, gradients, et)
    avg <- colMeans(res)
    tot <- colSums(res)
    round(rbind(res, all=c(avg[1:2], tot[3], max(res[, 4]), tot[5:7])), 4)
}
