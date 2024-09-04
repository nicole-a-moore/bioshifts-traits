## trying a Bayesian approach

#### building the generative model 

## Sample:
## range shift 
rs = rnorm(n = 100, mean = 30, sd = 10)
hist(rs)

## Fixed effects:
## climate velocity 
cv = 
## dispersal potential 
dp 

## Parameters to estimate:
## the y value of the breakpoint in the relationship (this is also the slope after the breakpoint)
B
## the slope before the breakpoint (dp < cv)
s
## the intercept of the relationship before the breakpoint (dp < cv)
int






# Load rstan package
library(rstan)
# Avoid unnecessary recompiling
rstan_options(auto_write = TRUE)
# optional: Distribute work over multiple CPU cores
options(mc.cores = parallel::detectCores())

## set a seed
set.seed(696969)

## specify sample size to simulate 
n <- 300

## generate climate velocity data (exponential distribution, from ~0 to ~7)
#CV <- rnorm(n = n, mean = 5, sd = 1)
CV <- rnorm(n = n, mean = 15, sd = 1)
hist(CV)

## generate dispersal distance data (from 1 - 1000, with most species having small values)
DP <- rnorm(n = n, mean = 15, sd = 10)
hist(DP)

## set NA to 0
DP[is.na(DP)] <- 0

## generate average values on range expansion rate outcome
## when DP < CV, RS = DP
## when DP > CV, RS = CV
meanRS <- ifelse(DP < CV,
                  0 + 1 * DP,
                  CV)
plot(DP, meanRS, main = "Underlying function")


## generate observed values
error <- 0.5
RS <- rnorm(n = n, mean = meanRS, sd = error)
plot(DP, RS, main = "Observed data")

## remove negative obs range shifts 
obsRS <- RS[RS > 0]
DP <- DP[RS > 0]
CV <- CV[RS > 0]
plot(DP, obsRS, main = "Observed data")


model_bp <- '
// You need to specify the kind of input data, incl. number of observations.
data { 
  int<lower=1> N;  // total number of observations (integer); at least 1
  real RS[N];     // outcome variable with N elements (real-valued)
  real DP[N];     // predictor variable with N elements (real-valued)
  real CV[N];     // predictor variable with N elements (real-valued)
}

// the parameters to be estimated from the data
parameters { 
  real intercept;                 // = predicted RS when DP = 0
  real slope_before;              // slope when DP < CV
  real slope_after;               // slope when DP > CV
  real<lower = 0> error;          // standard deviation of residuals
                                  //  (always positive, hence <lower = 0>)
} 

// Functions of estimated parameters.
transformed parameters{
  vector[N] conditional_mean; // the estimated average RS for each DP observation

  // conditional_mean depends on whether DP is less than or greater than CV
  for (i in 1:N) {
    if (DP[i] < CV[i]) {
      conditional_mean[i] = slope_before * DP[i] + intercept;
    } else {
      conditional_mean[i] = slope_after * DP[i] + slope_before * DP[bp] + intercept;
    }
  }
}

// The model itself specifies how the data are expected to have
// been generated and what the prior expectations for the model parameters are.
model {
  // Set priors
  intercept ~ normal(0, 1);  // Average RS at breakpoint
  slope_before ~ normal(1, 1);  // Slope before breakpoint
  slope_after ~ normal(0, 1);   // Slope after breakpoint
  error ~ normal(0, 0.5);        // Residual error, likely between 0 and 2*20
  
  // How the data are expected to have been generated:
  // normal distribution with mu = conditional_mean and 
  // std = error, estimated from data.
  for (i in 1:N) {
    RS[i] ~ normal(conditional_mean[i], error);
  }
}

generated quantities {
  vector[N] sim_RS;               // Simulate new data using estimated parameters.
  vector[N] log_lik;               // Useful for model comparisons; not done here.
  vector[50] sim_conditional_mean; // Useful for plotting.

  // Compute conditional means for DPs between 1 and 50.
  for (i in 1:50) {
    if (i < CV[i]) {
      sim_conditional_mean[i] = intercept + slope_before * i;
    } else {
      sim_conditional_mean[i] = intercept + slope_after * i;
    }
  }

  for (i in 1:N) {
    sim_RS[i] = normal_rng(conditional_mean[i], error);
    log_lik[i] = normal_lpdf(RS[i] | conditional_mean[i], error);
  }
}
'


data_list <- list(
  DP = DP,
  RS = obsRS,
  CV = CV,
  N = length(DP)
)
fit_bp_sim <- stan(model_code = model_bp, 
                   data = data_list)

print(fit_bp_sim,
      par = c("intercept", "slope_before", "slope_after", "error"))


shinystan::launch_shinystan(fit_bp_sim)

# rstan's 'extract' is likely to conflict with another function
# called 'extract', so specify the package, too.
simulated_data <- rstan::extract(fit_bp_sim)$sim_RS
# simulated_data is a matrix with 4000 rows and 80 columns.
# For the plot, I select 8 rows at random:
simulated_data <- simulated_data[sample(1:4000, 8), ]

par(mfrow = c(3, 3))

# Plot the observed data
plot(data_list$DP, data_list$RS,
     xlab = "DP", ylab = "RS",
     main = "observed")

# Plot the simulated data
for (i in 1:8) {
  plot(data_list$DP, simulated_data[i, ],
       xlab = "DP", ylab = "RS",
       main = "simulated")
}




##########################
library(rethinking)

## make fake dispersal distance data
dp = rlnorm(200, 1, 3)
hist(dp)
## make climate velocity data
cv = runif(200, 0, 8)
hist(cv)
## make fake range shift
rs = ifelse(dp > cv, cv, dp)
hist(rs)

## plot 
data <- data.frame(dp = dp, cv = cv, rs = rs) 
data %>%
  ggplot(aes(x = dp, y = rs, colour = cv)) +
  geom_point() +
  coord_equal()


set.seed(2971)

N <- 100 # 100 lines

a <- rnorm(N, 0, 1) # intercept - mean 0
b <- rlnorm(N, 1, 0.3) # first slope - constrain so always positive using log normal distribution
b2 <- rnorm(N, 0, 0.05) # second slope - mean 0
bp <- rlnorm(N, 1, 1) # breakpoint - mean = mean climate velocity, constrain so always positive using log normal distribution

bp <- cv

## plot
plot(NULL, xlim=c(0,70), ylim=c(-20,30),
     xlab= "potential dispersal rate", ylab="range expansion rate")

abline(h=0 , lty=2)
abline(h=272, lty=1, lwd=0.5)

for (i in 1:N) {
  curve(ifelse(x < first(which(x > bp[i])),
                (a[i] + b[i]*(x)),
                (a[i] + b[i]*(first(which(x > bp[i]))) + b2[i]*(x))),
  from = min(dp), 
  to = max(dp), 
  add = TRUE,
  col = col.alpha("black",0.2),
  n = 70)
}


## fit the model
m4.3 <- quap(alist(
  rs ~ dnorm(mu, sigma),
  mu <- if_else(dp < cv, 
               a + b*(dp),
               0*dp),
  a ~ dnorm(0, 1),
  b ~ dlnorm(1, 0.3),
  sigma ~ dunif(0, 3)), 
  data = data)





model_bp <- '
// You need to specify the kind of input data, incl. number of observations.
data { 
  int<lower=1> N;  // total number of observations (integer); at least 1
  real RS[N];     // outcome variable with N elements (real-valued)
  real DP[N];     // predictor variable with N elements (real-valued)
  real CV[N];     // predictor variable with N elements (real-valued)
}

// the parameters to be estimated from the data
parameters { 
  real intercept;                 // predicted RS when DP = 0
  real slope_before;              // slope when DP < CV
  real slope_after;               // slope when DP > CV
  real<lower = 0> error;          // standard deviation of residuals
                                  //  (always positive, hence <lower = 0>)
} 

// Functions of estimated parameters.
transformed parameters{
  vector[N] conditional_mean; // the estimated average RS for each DP observation

  // conditional_mean depends on whether DP is less than or greater than CV
  for (i in 1:N) {
    if (DP[i] < CV[i]) {
      conditional_mean[i] = intercept + slope_before * (DP[i] - CV[i]);
    } else {
      conditional_mean[i] = intercept + slope_after * (DP[i] - CV[i]);
    }
  }
}

// The model itself specifies how the data are expected to have
// been generated and what the prior expectations for the model parameters are.
model {
  // Set priors
  intercept ~ normal(0, 1);  // Average RS at breakpoint
  slope_before ~ normal(1, 1);  // Slope before breakpoint
  slope_after ~ normal(0, 1);   // Slope after breakpoint
  error ~ normal(0, 0.5);        // Residual error, likely between 0 and 2*20
  
  // How the data are expected to have been generated:
  // normal distribution with mu = conditional_mean and 
  // std = error, estimated from data.
  for (i in 1:N) {
    RS[i] ~ normal(conditional_mean[i], error);
  }
}

generated quantities {
  vector[N] sim_RS;               // Simulate new data using estimated parameters.
  vector[N] log_lik;               // Useful for model comparisons; not done here.
  vector[50] sim_conditional_mean; // Useful for plotting.

  // Compute conditional means for DPs between 1 and 50.
  for (i in 1:50) {
    if (i < CV[i]) {
      sim_conditional_mean[i] = intercept + slope_before * i;
    } else {
      sim_conditional_mean[i] = intercept + slope_after * i;
    }
  }

  for (i in 1:N) {
    sim_RS[i] = normal_rng(conditional_mean[i], error);
    log_lik[i] = normal_lpdf(RS[i] | conditional_mean[i], error);
  }
}
'

data_list <- list(
  DP = dp,
  RS = rs,
  CV = cv,
  N = length(dp)
)
fit_bp_sim <- stan(model_code = model_bp, 
                   data = data_list)
print(fit_bp_sim,
      par = c("intercept", "slope_before", "slope_after", "error"))

# rstan's 'extract' is likely to conflict with another function
# called 'extract', so specify the package, too.
simulated_data <- rstan::extract(fit_bp_sim)$sim_RS
# simulated_data is a matrix with 4000 rows and 80 columns.
# For the plot, I select 8 rows at random:
simulated_data <- simulated_data[sample(1:4000, 8), ]

par(mfrow = c(3, 3))

# Plot the observed data
plot(log(data_list$DP), data_list$RS,
     xlab = "DP", ylab = "RS",
     main = "observed")

# Plot the simulated data
for (i in 1:8) {
  plot(log(data_list$DP), simulated_data[i, ],
       xlab = "DP", ylab = "RS",
       main = "simulated")
}

