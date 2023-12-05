# Run some VB on a simple AMMI-type model with the interaction only

# Clear workspace and load packages
rm(list = ls())
library(tidyverse)
library(truncnorm)
library(gridExtra)
library(R2jags)

# Simulate data -----------------------------------------------------------

set.seed(123)
n_g <- 8 # Number of genotypes
n_e <- 5 # Number of environment
n <- n_g * n_e # Number of obs altogether
m_g <- 10
s_g <- 1
m_e <- 10
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

# VB version --------------------------------------------------------------

# Set up a place to store the values
n_iter <- 1000

# Starting values
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
  beta_tau <- beta + (sum(y^2) + mu_lambda_sq * sum(mu_gamma_sq[genotype] * mu_delta_sq[environment]) - 2 * mu_lambda * sum(y * mu_gamma[genotype] * mu_delta[environment])) / 2
  mu_tau <- alpha_tau / beta_tau

  # Update lambda
  tau_lambda <- mu_tau * sum(mu_gamma_sq[genotype] * mu_delta_sq[environment]) + omega
  mu_lambda <- mu_tau * sum(y * mu_gamma[genotype] * mu_delta[environment]) / tau_lambda
  mu_lambda_sq <- mu_lambda^2 + 1 / tau_lambda

  # Update gamma
  tau_gamma <- mu_tau * mu_lambda_sq * sum(mu_delta_sq) + 1
  for (i in 1:n_g) {
    mu_gamma[i] <- mu_tau * mu_lambda * sum(y[genotype == i] * mu_delta) / tau_gamma
  }
  mu_gamma_sq <- mu_gamma^2 + 1 / tau_gamma

  # Update delta
  tau_delta <- mu_tau * mu_lambda_sq * sum(mu_gamma_sq) + 1
  for (j in 1:n_e) {
    mu_delta[j] <- mu_tau * mu_lambda * sum(y[environment == j] * mu_gamma) / tau_delta
  }
  mu_delta_sq <- mu_delta^2 + 1 / tau_delta
}

# Check plots -------------------------------------------------------------

# Check the bi-linear mean
vb_blin <- mu_lambda * mu_gamma[genotype] * mu_delta[environment]
qplot(blin, vb_blin) + geom_abline()
