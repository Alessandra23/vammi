model_code <- "
model
{
  # Likelihood
   for (i in 1:N) {
    # Model for phenotype
    y[i] ~ dnorm(mu[i], sigma^-2)
    mu[i] = muall + g[gen[i]] + e[env[i]] + blin[i]
    blin[i] = sum(lambda[1:Q] * gamma[gen[i],1:Q] * delta[env[i],1:Q])

   }

   muall ~ dnorm(90, 10^-2) # grand mean

   # Prior on genotype effect
  for(i in 1:I) {
    g[i] ~ dnorm(0, sigma_g^-2) # Prior on genotype effect
  }

  for(i in 1:J) {
    e[i] ~ dnorm(0, sigma_e^-2) # Prior on genotype effect
  }

  # Priors on gamma
  for(q in 1:Q){
    for(i in 1:I){
      thetaG[i,q] ~ dnorm(0,1)
    }
    mG[q] = sum(thetaG[1:I,q])/I
    for(i in 1:I){
    thetaGNew[i,q] = thetaG[i,q] - mG[q]
    }
    sqrtThetaG[q] = sqrt(1/(sum(thetaGNew[1:I,q]^2 + 0.000001)))
    for(i in 1:I){
      gamma[i,q] = thetaGNew[i,q]*sqrtThetaG[q]
    }
  }

   # Priors on delta
   for(q in 1:Q){
    for(j in 1:J){
      thetaD[j,q] ~ dnorm(0,1)
    }
    mD[q] = sum(thetaD[1:J,q])/J
    for(j in 1:J){
    thetaDNew[j,q] = thetaD[j,q] - mD[q]
    }
    sqrtThetaD[q] = sqrt(1/(sum(thetaDNew[1:J,q]^2+ 0.000001)))
    for(j in 1:J){
      delta[j,q] = thetaDNew[j,q]*sqrtThetaD[q]
    }
  }

  # Prior on eigenvalues
  for(q in 1:Q) {
    lambda_raw[q] ~ dnorm(0, 100^-2)T(0,)
  }
  lambda = sort(lambda_raw)

  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma_e ~ dt(0, 10^-2, 1)T(0,)
  sigma_g ~ dt(0, 10^-2, 1)T(0,)
}
"


# Set up the data
model_data <- list(N = length(dat$y), y = dat$y, J = dat$n_e, I = dat$n_g, Q = dat$Q,
                   gen = dat$x$g, env = dat$x$e)

# Choose the parameters to watch
model_parameters <- c("g", "e", "mu", "sigma", "blin")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)


qplot(model_run$BUGSoutput$mean$g, dat$g) + theme_bw() + geom_abline()



# OLD CODE ----------------------------------------------------------------


# Run some VB on a simple AMMI-type model with the interaction only

# Clear workspace and load packages
rm(list = ls())
library(tidyverse)
library(truncnorm)
library(gridExtra)
library(R2jags)

# Simulate data -----------------------------------------------------------

# Simulate some data

set.seed(123)
n_g <- 8 # Number of genotypes
n_e <- 5 # Number of environment
n <- n_g * n_e # Number of obs altogether
m_g <- 2
s_g <- 1
m_e <- 2
s_e <- 1
alpha <- 5
beta <- 1
tau <- rgamma(1, alpha, beta) # Residual precision
omega <- 1
lambda <- 1 # Only eigenvalue
gamma <- c(1, rnorm(n_g - 1)) # Genotype effects - first has to be positive
delta <- rnorm(n_e) # Environment effects
gen_env <- expand.grid(1:n_g, 1:n_e) # Create grid to look up values
genotype <- gen_env[, 1]
environment <- gen_env[, 2]

# Finally create y
g <- rnorm(n_g, m_g, s_g)
e <- rnorm(n_e, m_e, s_e)
blin <- y <- rep(NA, n)
for (i in 1:n) {
  blin[i] <- lambda * gamma[genotype[i]] * delta[environment[i]]
}

mu_ij <- g[genotype] + e[environment] + blin
y <- rnorm(n, mean = mu_ij, sd = 1 / sqrt(tau))
df <- data.frame(genotype, environment, y)


# Do some plots if required
# p1 <- qplot(x = factor(genotype), y = y, data = df, geom = "boxplot")
# p2 <- qplot(x = factor(environment), y = y, data = df, geom = "boxplot")



# VB version --------------------------------------------------------------

# Set up a place to store the values
n_iter <- 100

# Starting values
mu_g <- g
tau_g <- 1
mu_g_sq <- mu_g^2 + 1/tau_g
mu_e <- e
tau_e <- 1
mu_e_sq <- mu_e^2 + 1/tau_e
alpha_tau <- alpha + n / 2 - 1
tau_lambda <- 100
mu_lambda <- 1
mu_lambda_sq <- mu_lambda^2 + 1 / tau_lambda
tau_gamma <- 100
mu_gamma <- rep(1, n_g)
mu_gamma_sq <- mu_gamma^2 + 1 / tau_gamma
tau_delta <- 100
mu_delta <- delta
mu_delta_sq <- mu_delta^2 + 1 / tau_delta


for (t in 1:n_iter) {

  # Update tau
  beta_tau <- beta + (sum(y^2) + mu_lambda_sq * sum(mu_gamma_sq[genotype] * mu_delta_sq[environment]) +
                        n_e*sum(mu_g_sq) + 2*sum(mu_g)*sum(mu_e) + n_g*sum(mu_e_sq) +
                        2*mu_lambda*sum(mu_g*mu_gamma[genotype] * mu_delta[environment]) +
                        2*mu_lambda*sum(mu_e*mu_gamma[genotype] * mu_delta[environment])-
                        2* sum(y*mu_g) - 2*sum(y*mu_e) -
                        2 * mu_lambda * sum(y * mu_gamma[genotype] * mu_delta[environment])) / 2
  mu_tau <- alpha_tau / beta_tau
  #mu_tau <- 0.5

  # Update g
  tau_g <- n_e*mu_tau + 1/tau_g
  for(i in 1:n_g){
    mu_g[i] <- (mu_tau*(sum(y[genotype == i]) - sum(mu_e) -
                          mu_lambda*sum(mu_gamma[genotype] * mu_delta)))/2*tau_g
  }
  #mu_g <- mu_g - mean(mu_g)
  mu_g_sq <- mu_g^2 + 1/tau_g

  # Update e
  tau_e <- n_g*mu_tau + 1/tau_e
  for(j in 1:n_e){
    mu_e[j] <- (mu_tau*(sum(y[environment == j]) - sum(mu_g) -
                          mu_lambda*sum(mu_gamma[genotype]*mu_delta[environment])))/2*tau_e
  }
  #mu_e <- mu_e - mean(mu_e)
  mu_e_sq <- mu_e^2 + 1/tau_e

  # Update lambda
  tau_lambda <- mu_tau * sum(mu_gamma_sq[genotype] * mu_delta_sq[environment]) + omega
  mu_lambda <- mu_tau*(sum(y * mu_gamma[genotype] * mu_delta[environment]) -
                          sum(g*mu_gamma[genotype] * mu_delta[environment]) -
                          sum(e*mu_gamma[genotype] * mu_delta[environment])) / 2*tau_lambda
  mu_lambda_sq <- mu_lambda^2 + 1 / tau_lambda

  # Update gamma
  tau_gamma <- mu_tau * mu_lambda_sq * sum(mu_delta_sq) + 1
  for (i in 1:n_g) {
    mu_gamma[i] <- -mu_tau * (-mu_lambda * sum(y[genotype == i] * mu_delta) +
                             mu_lambda*mu_g[i]*sum(mu_delta) +
                             mu_lambda*sum(mu_e*mu_delta)) / 2*tau_gamma
  }
  mu_gamma_sq <- mu_gamma^2 + 1 / tau_gamma

  # Update delta
  tau_delta <- mu_tau * mu_lambda_sq * sum(mu_gamma_sq) + 1
  for (j in 1:n_e) {
    mu_delta[j] <- -mu_tau *(-mu_lambda * sum(y[environment == j] * mu_gamma) +
                              mu_lambda*sum(g*gamma) +
                               mu_lambda*mu_e[j]*sum(gamma)) / 2*tau_delta
  }
  mu_delta_sq <- mu_delta^2 + 1 / tau_delta
}

# Check plots -------------------------------------------------------------

# Check the bi-linear mean
vb_blin <- mu_lambda * mu_gamma[genotype] * mu_delta[environment]
qplot(blin, vb_blin) + geom_abline() + theme_bw()

qplot(g, mu_g) + geom_abline() + theme_bw()
qplot(e, mu_e) + geom_abline() + theme_bw()





