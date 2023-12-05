# VB code

# Run some VB on a simple AMMI-type model with the interaction only

# Clear workspace and load packages
rm(list = ls())
library(tidyverse)
library(truncnorm)
library(gridExtra)
library(R2jags)
library(ggplot2)


vbAMMI <- function(n_iter = 1000,
                   data,
                   hyperparameters = list(),
                   initial_values = list()) {

  y <- data$y
  n <- length(data$y)
  n_g <- data$n_g
  n_e <- data$n_e
  Q <- data$Q
  genotype <- data$x[, 1]
  environment <- data$x[, 2]

  # reading the hyperparameters
  s_mu <- hyperparameters$s_mu
  m_mu <- hyperparameters$m_mu
  s_g  <- hyperparameters$s_g
  m_g  <- hyperparameters$m_g
  s_e  <- hyperparameters$s_e
  m_e  <- hyperparameters$m_e
  s_lambda <- hyperparameters$s_lambda
  beta <- hyperparameters$beta

  # reading the initial values
  mu_mu <-  initial_values$mu_mu
  tau_mu <- initial_values$tau_mu
  mu_g <- initial_values$mu_g
  tau_g <- initial_values$tau_g
  mu_e <- initial_values$mu_e
  tau_e <- initial_values$tau_e
  alpha_tau <-initial_values$alpha_tau
  tau_lambda <- initial_values$tau_lambda
  mu_lambda <- initial_values$mu_lambda
  mu_lambda_sq <- initial_values$mu_lambda_sq
  tau_gamma <-initial_values$tau_gamma
  mu_gamma <- initial_values$mu_gamma
  mu_gamma_sq <- initial_values$mu_gamma_sq
  tau_delta <-initial_values$tau_delta
  mu_delta <- initial_values$mu_delta
  mu_delta_sq <- initial_values$mu_delta_sq
  mu_mu_sq <- initial_values$mu_mu_sq
  mu_g_sq <- initial_values$mu_g_sq
  mu_e_sq <- initial_values$mu_e_sq

  # running
  for (iter in 1:n_iter) {
    # Update tau
    beta_tau <- beta + (sum(y^2) - 2*mu_mu*sum(y) - 2*sum(y*mu_g) - 2*sum(y*mu_e) -
                          mu_lambda*sum(y*mu_gamma[genotype]*mu_delta[environment]) +
                          n*mu_mu_sq + 2*n_e*mu_mu*sum(mu_g) + 2*n_g*mu_mu*sum(mu_e) +
                          n_e*sum(mu_g_sq) + 2*sum(mu_g)*sum(mu_e) + n_g*sum(mu_e_sq)+
                          2*mu_mu*mu_lambda*sum(mu_gamma[genotype]*mu_delta[environment]) +
                          2*mu_lambda*sum(mu_g*mu_gamma[genotype]*mu_delta[environment]) +
                          2*mu_lambda*sum(mu_e*mu_gamma[genotype]*mu_delta[environment]) +
                          mu_lambda_sq*sum(mu_gamma_sq[genotype]*mu_delta_sq[environment]))/2
    mu_tau <- alpha_tau / beta_tau

    # Update mu
    tau_mu <- n*mu_tau + (1/s_mu^2)
    mu_mu <- ((m_mu/s_mu)- (mu_tau * (sum(y) + n_e*sum(mu_g) + n_g*sum(mu_e) +
                            mu_lambda*sum(mu_gamma[genotype] * mu_delta[environment]))))/(2*tau_mu)
    mu_mu_sq <- mu_mu^2 + 1/tau_mu

    # Update g
    tau_g <- n_e*mu_tau + 1/s_g^2
    for(i in 1:n_g){
      mu_g[i] <- (mu_tau*(sum(y[genotype == i]) - n_e*mu_mu - sum(mu_e) -
                            mu_lambda*sum(mu_gamma[genotype] * mu_delta)))/(2*tau_g)
    }
    #mu_g <- mu_g - mean(mu_g)
    mu_g_sq <- mu_g^2 + 1/tau_g

    # Update e
    tau_e <- n_g*mu_tau + 1/s_e^2
    for(j in 1:n_e){
      mu_e[j] <- (mu_tau*(sum(y[environment == j]) - n_g*mu_mu - sum(mu_g) -
                            mu_lambda * sum(mu_gamma * mu_delta[environment])))/(2*tau_e)
    }
    #mu_e <- mu_e - mean(mu_e)
    mu_e_sq <- mu_e^2 + 1/tau_e

    # Update lambda
    tau_lambda <- mu_tau * sum(mu_gamma_sq[genotype] * mu_delta_sq[environment]) + 1/(s_lambda^2)
    mu_lambda <- (mu_tau * (sum(y*mu_gamma[genotype] * mu_delta[environment]) -
                            mu_mu*sum(mu_gamma[genotype] * mu_delta[environment]) -
                            sum(mu_g*mu_gamma[genotype] * mu_delta[environment]) -
                            sum(mu_e*mu_gamma[genotype] * mu_delta[environment])))/(2*tau_lambda)
    mu_lambda_sq <- mu_lambda^2 + 1 / tau_lambda

    # Update gamma
    tau_gamma <- mu_tau * mu_lambda_sq * sum(mu_delta_sq) + 1
    for (i in 1:n_g) {
      mu_gamma[i] <- (- mu_tau * (-mu_lambda * sum(y[genotype == i] * mu_delta) +
                        mu_mu*mu_lambda*sum(mu_delta) + mu_lambda*sum(mu_delta*mu_e) +
                        mu_lambda*mu_g[i]*sum(mu_delta))) / (2*tau_gamma)
    }
    mu_gamma <- (mu_gamma - mean(mu_gamma))/sd(mu_gamma)
    mu_gamma_sq <- mu_gamma^2 + 1 / tau_gamma

    # Update delta
    tau_delta <- mu_tau * mu_lambda_sq * sum(mu_gamma_sq) + 1
    for (j in 1:n_e) {
      mu_delta[j] <- (-mu_tau *(-mu_lambda * sum(y[environment == j] * mu_gamma) +
                        mu_mu*mu_lambda*sum(mu_gamma) + mu_lambda*sum(mu_gamma*mu_g) +
                        mu_lambda*mu_e[j]*sum(mu_gamma))) / (2*tau_delta)
    }
    mu_delta <- (mu_delta - mean(mu_delta))/sd(mu_delta)
    mu_delta_sq <- mu_delta^2 + 1 / tau_delta
  }

  return(list(mu_mu = mu_mu,
              tau_mu = tau_mu,
              mu_g = mu_g,
              tau_g = tau_g,
              mu_e = mu_e,
              tau_e = tau_e,
              mu_lambda = mu_lambda,
              tau_lambda = tau_lambda,
              mu_gamma = mu_gamma,
              tau_lambda = tau_lambda,
              mu_delta = mu_delta,
              tau_delta = tau_delta,
              mu_tau = mu_tau
              ))

}



# VB version --------------------------------------------------------------

dat <- generate_data_AMMI(n_g = 25,
                          n_e = 10,
                          m_g = 10,
                          s_g = 1,
                          m_e = 10,
                          s_e = 1,
                          s_y = 1,
                          m_lambda = 12,
                          s_lambda = 1)

# hyperparameters
hyperparameters <- list(s_mu = 10 ,
                        m_mu = 100 ,
                        s_g = 1 ,
                        m_g = 10 ,
                        s_e = 1 ,
                        m_e = 10 ,
                        m_lambda = 12,
                        s_lambda = 1,
                        alpha = 1,
                        beta = 5)

# Starting values
initial_values <- list(mu_mu = 100,
                       tau_mu = 1/100,
                       mu_mu_sq = 100 + 100,
                       mu_g = dat$g,
                       tau_g = 1/100,
                       mu_g_sq = 100 + 100,
                       mu_e = dat$e,
                       tau_e = 1/100,
                       mu_e_sq = 100 + 100,
                       alpha_tau = hyperparameters$alpha + length(dat$y),
                       tau_lambda = 1/100,
                       mu_lambda = 12,
                       mu_lambda_sq = 12^2 + 100, #mu_lambda^2 + 1 / tau_lambda ,
                       tau_gamma = 1/100,
                       mu_gamma = dat$gamma[,1],
                       mu_gamma_sq = dat$gamma[,1]^2 + 100, #mu_gamma^2 + 1 / tau_gamma ,
                       tau_delta = 1/100,
                       mu_delta = dat$delta[,1],
                       mu_delta_sq = dat$delta[,1]^2 + 100) #mu_delta^2 + 1 / tau_delta)


model_run <- vbAMMI(data = dat,
                    hyperparameters = hyperparameters,
                    initial_values = initial_values)

model_run

#values_e <- matrix(rnorm(1000, model_run$mu_e, sqrt(1/model_run$tau_e)), ncol = length(model_run$mu_e))
#apply(values_e, 2, mean)

qplot(dat$e, model_run$mu_e) + geom_abline()
qplot(dat$g, model_run$mu_g) + geom_abline()



vb_blin <- model_run$mu_lambda * model_run$mu_gamma[dat$x[,1]] * model_run$mu_delta[dat$x[,2]]
qplot(dat$blin, vb_blin) + geom_abline()


# Run vammi using different  initial values

# random values

# using classical ammi

# using jags


## setting jags

# Set up the data
model_data <- list(N = length(dat$y), y = dat$y, J = dat$n_e, I = dat$n_g, Q = dat$Q,
                   gen = dat$x$g, env = dat$x$e)

# Choose the parameters to watch
model_parameters <- c("g", "e", "mu", "sigma", "blin")

# setting gibbs

# Set the data
I <- dat$n_g
J <- dat$n_e
set.seed(02)
y <- matrix(rnorm(I * J), nrow = I, ncol = J)
# dat$y
#
# matrix(dat$y, nrow = dat$n_g, ncol = dat$n_e)
#
# dat$y[dat$x$g]
# Set the priors
mu_prior <- list(mean = 0, sigma2 = 1)
g_prior <- list(mean = 10, sigma2 = 1)
e_prior <- list(mean = 10, sigma2 = 1)
lambda_prior <- list(mean = 0, sigma2 = 1)
inv_sigma2_prior <- list(shape = 2, rate = 2)

# Set the number of iterations and burn-in
n_iter <- 6000
burn_in <- 1000

# Run the Gibbs sampler
# posterior_samples <- gibbs_sampler(data = y,
#                                    mu_prior = mu_prior,
#                                    g_prior = g_prior,
#                                    e_prior = e_prior,
#                                    lambda_prior = lambda_prior,
#                                    inv_sigma2_prior = inv_sigma2_prior,
#                                    n_iter = n_iter,
#                                    burn_in = burn_in)


# Comp time

times <- microbenchmark::microbenchmark(vammi = vbAMMI(data = dat,
                                                       hyperparameters = hyperparameters,
                                                       initial_values = initial_values),
                                        gibbs_ammi =  gibbs_sampler(data = y,
                                                                   mu_prior = mu_prior,
                                                                   g_prior = g_prior,
                                                                   e_prior = e_prior,
                                                                   lambda_prior = lambda_prior,
                                                                   inv_sigma2_prior = inv_sigma2_prior,
                                                                   n_iter = n_iter,
                                                                   burn_in = burn_in),
                                        jags_ammi = jags(
                                          data = model_data,
                                          parameters.to.save = model_parameters,
                                          model.file = textConnection(model_code)
                                        ),
                                        times = 1,
                                        unit = 's')


times

times_df <- data.frame(n = c(50, 250, 500, 750, 1000, 5000, 10000, 15000, 20000),
                    vammi = c(0.7600064,),
                    gibbs = c(0.5683700,))




# fit classical ammi ------------------------------------------------------

dat <- generate_data_AMMI(n_g = 25,
                          n_e = 10,
                          m_g = 10,
                          s_g = 1,
                          m_e = 10,
                          s_e = 1,
                          s_y = 2,
                          m_lambda = 12,
                          s_lambda = 1)

classical <- run_classical_AMMI(dat, Q = 1)

# hyperparameters
hyperparameters <- list(s_mu = 10,
                        m_mu = 100 ,
                        s_g = 1 ,
                        m_g = 100 ,
                        s_e = 1 ,
                        m_e = 100 ,
                        m_lambda = 12,
                        s_lambda = 1,
                        alpha = 5,
                        beta = 1)

# Starting values



initial_values_classical <- list(mu_mu = classical$mu_hat ,
                                 tau_mu = 10,
                                 mu_mu_sq = classical$mu_hat^2 + 1/1,
                                 mu_g = classical$g_hat ,
                                 tau_g = var(classical$g_hat) ,
                                 mu_g_sq = classical$g_hat^2 + 1/var(classical$g_hat),
                                 mu_e = classical$e_hat,
                                 tau_e = var(classical$e_hat),
                                 mu_e_sq = classical$e_hat^2 + 1/var(classical$e_hat),
                                 alpha_tau = hyperparameters$alpha + (length(dat$y)) / 2 - 1,
                                 tau_lambda = 100 ,
                                 mu_lambda = classical$lambda_hat ,
                                 mu_lambda_sq = classical$lambda_hat^2 + 1/100, #mu_lambda^2 + 1 / tau_lambda ,
                                 tau_gamma = as.vector(var(classical$gamma_hat)),
                                 mu_gamma = as.vector(classical$gamma_hat),
                                 mu_gamma_sq = as.vector(classical$gamma_hat)^2 + 1/as.vector(var(classical$gamma_hat)), #mu_gamma^2 + 1 / tau_gamma ,
                                 tau_delta = as.vector(var(classical$delta_hat)),
                                 mu_delta = as.vector(classical$delta_hat),
                                 mu_delta_sq = as.vector(classical$delta_hat)^2 + 1/as.vector(var(classical$delta_hat))) #mu_delta^2 + 1 / tau_delta)


vb_classical <- vbAMMI(data = dat,
       hyperparameters = hyperparameters,
       initial_values = initial_values_classical, n_iter = 100)


qplot(vb_classical$mu_g, dat$g) + theme_bw()
qplot(vb_classical$mu_e, dat$e) + theme_bw()

qplot(vb_classical$mu_g - mean(vb_classical$mu_g), dat$g) + theme_bw()

qplot(classical$g_hat, dat$g) + theme_bw()

