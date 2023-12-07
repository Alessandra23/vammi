# Clear workspace and load packages
rm(list = ls())
library(tidyverse)
library(truncnorm)
library(gridExtra)
library(R2jags)

## Generate data

generate_data_AMMI <- function(n_g, # Number of genotypes
                               n_e, # Number of environments
                               m_g,
                               s_g, # standard deviation of g
                               m_e,
                               s_e, # standard deviation of e
                               s_y, # standard deviation of y
                               m_lambda,
                               s_lambda
) {
  # Total number of observations
  N <- n_g * n_e

  # Generate g (genotypes)
  g <- rnorm(n_g, m_g, s_g)
  g <- g - mean(g) # impose the sum-to-zero restriction
  # Generate e (environments)
  e <- rnorm(n_e, m_e, s_e)
  e <- e - mean(e) # impose the sum-to-zero restriction
  # Set the grand mean
  mu <- rnorm(1, 90, 10)

  lambda <- truncnorm::rtruncnorm(n = 1, a = 0, mean = m_lambda, sd = s_lambda)
  # Number of components in the bilinear part
  Q <- length(lambda)

  # Generate gamma
  # gamma <- matrix(NA, nrow = I ,ncol = Q)
  # gamma[1,] <- truncnorm::rtruncnorm(Q, a=0)
  # gamma[-1,] <- rnorm((I-1)*Q)
  gamma <- generate_gamma_delta(n_g, Q)

  # Generate delta
  # delta <- matrix(rnorm(J*Q), nrow = J ,ncol = Q)
  delta <- generate_gamma_delta(n_e, Q)
  # Generate the "design matrix"
  x <- expand.grid(1:n_g, 1:n_e)
  names(x) <- c("g", "e") # g = genotype and e = envorinment
  x$g <- as.factor(x$g)
  x$e <- as.factor(x$e)

  # Generate the interaction/bilinear part
  blin <- rep(0, n_g*n_e)
  for (k in 1:length(lambda)) {
    blin <- blin + lambda[k] * gamma[x[, "g"], k] * delta[x[, "e"], k]
  }

  # Now simulate the response
  mu_ij <- mu + g[x[, "g"]] + e[x[, "e"]] + blin

  # Compute the response for the TRAINING set
  y <- rnorm(N, mu_ij, s_y)

  # Compute the response for the TEST set
  y_test <- rnorm(N, mu_ij, s_y)

  # This is a matrix representation from Alessandra (it works fine)
  # mu_ij <- mu*I1%*%t(J1) + kronecker(g,t(J1)) + kronecker(t(e), (I1)) + gamma%*%diag(lambda)%*%t(delta)
  # y <- rnorm(N, c(mu.Y), s_y)

  return(list(
    y = y,
    y_test = y_test,
    x = x,
    n_g = n_g,
    n_e = n_e,
    Q = Q,
    s_g = s_g,
    s_e = s_e,
    s_y = s_y,
    lambda = lambda,
    g = g,
    e = e,
    gamma = gamma,
    delta = delta,
    blinear = blin
  ))
}

square_root_matrix <- function(x) {
  # When Q = 1, x will be a scalar
  if (nrow(x) == 1) {
    return(sqrt(x))
  }

  # When Q > 1, then x will be a matrix
  if (nrow(x) > 1) {
    # Jordan normal form
    X <- eigen(x)
    P <- X$vectors
    A <- diag(X$values)

    A_sqrt <- diag(sqrt(X$values))
    P_inv <- solve(P)
    x_sqrt <- P %*% A_sqrt %*% P_inv
    return(x_sqrt)
  }
}

generate_gamma_delta <- function(INDEX, Q) {
  first_row <- TRUE

  while (first_row) {
    raw_par <- matrix(rnorm(INDEX * Q), ncol = Q)
    par_mean <- matrix(rep(apply(raw_par, 2, mean), each = nrow(raw_par)), ncol = Q)
    par_aux <- raw_par - par_mean

    # Constraints ----
    # apply(par_aux,2,sum)
    parTpar <- solve(t(par_aux) %*% (par_aux))
    A <- square_root_matrix(parTpar)
    samples <- par_aux %*% A

    # Force the first to be positive
    for (i in 1:nrow(samples)) {
      row1 <- samples[1, ]
      if (all(samples[i, ] > 0)) {
        aux <- samples[i, ]
        samples[1, ] <- aux
        samples[i, ] <- row1
        return(samples)
      }
    }
    # t(samples)%*%samples == 0
    # apply(samples,2,sum) == diag(Q)
  }
}


# vammi

vbAMMI <- function(n_iter = 2000,
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
  for (t in 1:n_iter) {

    # Update tau
    beta_tau <- beta + (sum(y^2) - 2*mu_mu*sum(y) - 2*sum(y*mu_g) - 2*sum(y*mu_e) -
                          2*mu_lambda*sum(y*mu_gamma[genotype]*mu_delta[environment]) +
                          n*mu_mu_sq +  2*n_e*mu_mu*sum(mu_g) +  2*n_g*mu_mu*sum(mu_e) +
                          2*mu_mu*mu_lambda*sum(mu_gamma[genotype]*mu_delta[environment]) +
                          n*sum(mu_g_sq) + 2*sum(mu_g)*sum(mu_e) +
                          2*mu_lambda*sum(mu_g*mu_gamma[genotype]*mu_delta[environment]) +
                          n*sum(mu_e_sq) + 2*mu_lambda*sum(mu_e*mu_gamma[genotype]*mu_delta[environment]) +
                          mu_lambda_sq * sum(mu_gamma_sq[genotype] * mu_delta_sq[environment])) / 2
    mu_tau <- alpha_tau / beta_tau

    # Update mu
    tau_mu <- n*mu_tau + (1/s_mu^2)
    mu_mu <- ((m_mu/s_mu^2) + (mu_tau * (sum(y) - n_e*sum(mu_g) - n_g*sum(mu_e) -
                                           mu_lambda*sum(mu_gamma[genotype] * mu_delta[environment]))))/tau_mu
    mu_mu_sq <- mu_mu^2 + 1/tau_mu

    # Update g
    tau_g <- n_e*mu_tau + 1/s_g^2
    for(i in 1:n_g){
      mu_g[i] <- ((m_g/s_g^2) + mu_tau*(sum(y[genotype == i]) - n_e*mu_mu - sum(mu_e) -
                                          mu_lambda*sum(mu_gamma[genotype] * mu_delta)))/(tau_g)
    }
    mu_g <- mu_g - mean(mu_g)
    mu_g_sq <- mu_g^2 + 1/tau_g

    # Update e
    tau_e <- n_g*mu_tau + 1/s_e^2
    for(j in 1:n_e){
      mu_e[j] <- ((m_e/s_e^2) + mu_tau*(sum(y[environment == j]) - n_g*mu_mu - sum(mu_g) -
                                          mu_lambda*sum(mu_gamma*mu_delta[environment])))/(tau_e)
    }
    mu_e <- mu_e - mean(mu_e)
    mu_e_sq <- mu_e^2 + 1/tau_e

    # Update lambda
    tau_lambda <- mu_tau * sum(mu_gamma_sq[genotype] * mu_delta_sq[environment]) + 1/(s_lambda^2)
    mu_lambda <- mu_tau * sum(y * mu_gamma[genotype] * mu_delta[environment]) / tau_lambda
    mu_lambda_sq <- mu_lambda^2 + 1 / tau_lambda

    # Update gamma
    tau_gamma <- mu_tau * mu_lambda_sq * sum(mu_delta_sq) + 1
    for (i in 1:n_g) {
      mu_gamma[i] <- mu_tau * mu_lambda * sum(y[genotype == i] * mu_delta) / tau_gamma
    }
    mu_gamma <- (mu_gamma - mean(mu_gamma))/sd(mu_gamma)
    mu_gamma_sq <- mu_gamma^2 + 1 / tau_gamma

    # Update delta
    tau_delta <- mu_tau * mu_lambda_sq * sum(mu_gamma_sq) + 1
    for (j in 1:n_e) {
      mu_delta[j] <- mu_tau * mu_lambda * sum(y[environment == j] * mu_gamma) / tau_delta
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


# gibbs

gibbs_sampler <- function(data, mu_prior, g_prior, e_prior,
                          lambda_prior, gamma_prior, delta_prior,
                          inv_sigma2_prior, n_iter, burn_in) {
  # I <- data$n_g
  # J <- data$n_e
  # data <- dat$y

  I <- nrow(data)
  J <- ncol(data)

  # Allocate memory for posterior samples
  mu_samples <- numeric(n_iter)
  g_samples <- matrix(0, n_iter, I)
  e_samples <- matrix(0, n_iter, J)
  lambda_samples <- numeric(n_iter)
  gamma_samples <- matrix(0, n_iter, I)
  delta_samples <- matrix(0, n_iter, J)
  inv_sigma2_samples <- numeric(n_iter)

  # Initialize parameters
  mu <- mu_prior$mean
  g <- rnorm(I, g_prior$mean, sqrt(g_prior$sigma2))
  e <- rnorm(J, e_prior$mean, sqrt(e_prior$sigma2))
  lambda <- rtruncnorm(1, a = 0, mean = lambda_prior$mean, sd = sqrt(lambda_prior$sigma2))
  gamma <- c(rtruncnorm(1, a = 0, mean = 0, sd = 1), rnorm(I - 1, 0, 1))
  delta <- rnorm(J, 0, 1)
  inv_sigma2 <- rgamma(1, inv_sigma2_prior$shape, inv_sigma2_prior$rate)

  # Iterate over the number of iterations
  for (iter in 1:n_iter) {
    # Update mu
    mu_mean <- (sum(data) - sum(g) * I - sum(e) * J - sum(lambda * gamma %*% t(delta))) / (I * J)
    mu <- rnorm(1, mu_mean, sqrt(mu_prior$sigma2))

    # Update g
    for (i in 1:I) {
      g_mean <- (sum(data[i, ]) - I * mu - sum(e) - sum(lambda * gamma[i] * delta)) / J
      g[i] <- rnorm(1, g_mean, sqrt(g_prior$sigma2))
    }

    # Update e
    for (j in 1:J) {
      e_mean <- (sum(data[, j]) - J * mu - sum(g) - sum(lambda * gamma * delta[j])) / I
      e[j] <- rnorm(1, e_mean, sqrt(e_prior$sigma2))
    }

    # Update lambda
    lambda_mean <- sum((data - mu - outer(g, rep(1, J), "+") - outer(rep(1, I), e, "+")) * gamma %*% t(delta)) / sum(gamma^2 %*% t(delta^2))
    lambda <- rtruncnorm(1, a = 0, mean = lambda_mean, sd = sqrt(lambda_prior$sigma2))

    # Update gamma
    gamma[1] <- rtruncnorm(1, a = 0, mean = 0, sd = 1)
    gamma[-1] <- rnorm(I - 1, 0, 1)

    # Update delta
    delta <- rnorm(J, 0, 1)

    # Update inv_sigma2
    epsilon <- data - mu - outer(g, rep(1, J), "+") - outer(rep(1, I), e, "+") - lambda * gamma %*% t(delta)
    inv_sigma2_shape <- inv_sigma2_prior$shape + (I * J) / 2
    inv_sigma2_rate <- inv_sigma2_prior$rate + sum(epsilon^2) / 2
    inv_sigma2 <- rgamma(1, inv_sigma2_shape, inv_sigma2_rate)

    # Store samples
    mu_samples[iter] <- mu
    g_samples[iter, ] <- g
    e_samples[iter, ] <- e
    lambda_samples[iter] <- lambda
    gamma_samples[iter, ] <- gamma
    delta_samples[iter, ] <- delta
    inv_sigma2_samples[iter] <- inv_sigma2
  }

  # Discard burn-in samples
  mu_samples <- mu_samples[-(1:burn_in)]
  g_samples <- g_samples[-(1:burn_in), ]
  e_samples <- e_samples[-(1:burn_in), ]
  lambda_samples <- lambda_samples[-(1:burn_in)]
  gamma_samples <- gamma_samples[-(1:burn_in), ]
  delta_samples <- delta_samples[-(1:burn_in), ]
  inv_sigma2_samples <- inv_sigma2_samples[-(1:burn_in)]

  # Return posterior samples
  list(
    mu = mu_samples,
    g = g_samples,
    e = e_samples,
    lambda = lambda_samples,
    gamma = gamma_samples,
    delta = delta_samples,
    inv_sigma2 = inv_sigma2_samples
  )
}


# jags

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



# VB version --------------------------------------------------------------

n_g = 100
n_e = 10
n <- n_g * n_e

dat <- generate_data_AMMI(n_g = n_g,
                          n_e = n_e,
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
                        alpha = 0.1,
                        beta = 0.1)

# Starting values
initial_values <- list(mu_mu = 10 ,
                       tau_mu = 1 ,
                       mu_mu_sq = 100 + 1/1,
                       mu_g = rep(10, dat$n_g) ,
                       tau_g = 1 ,
                       mu_g_sq = 100 + 1/1,
                       mu_e = rep(10, dat$n_e) ,
                       tau_e = 1 ,
                       mu_e_sq = 100 + 1/1,
                       alpha_tau = hyperparameters$alpha + (length(dat$y)) / 2 - 1,
                       tau_lambda = 100 ,
                       mu_lambda = 1 ,
                       mu_lambda_sq = 1 + 1/100, #mu_lambda^2 + 1 / tau_lambda ,
                       tau_gamma = 100 ,
                       mu_gamma = rep(1, dat$n_g) ,
                       mu_gamma_sq = rep(1, dat$n_g)^2 + 1/100, #mu_gamma^2 + 1 / tau_gamma ,
                       tau_delta = 100 ,
                       mu_delta = rep(1, dat$n_e),
                       mu_delta_sq = rep(1, dat$n_e)^2 + 1 / 100) #mu_delta^2 + 1 / tau_delta)




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
#set.seed(02)
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
                                                       n_iter = 1000,
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

# Including in the times_df the values from times
times_df <- data.frame(n = c(50, 250, 500, 750, 1000, 5000, 10000, 15000, 20000),
                       gibbs = c(3.612318, 8.988203, 20.425634, 45.024673, 43.002456,280.234658,608.452637,1060.537389,1500.948734),
                       jaggs = c(4.734365, 52.588370, 162.456398, 358.973845, 503.77655,829.028374,1485.029374,2657.837465,3786.638466),
                       vammi = c(1.483440, 2.996492, 5.129946, 16.560467, 18.675354, 150.283747,270.638394,480.625384,758.018263))

times_df_small <- times_df[, -3] |> filter(n %in% c(50, 250, 500, 750, 1000))
times_df_small <- times_df_small |> reshape2::melt(id.vars = 'n')
levels(times_df_small$variable) <- c('MCMC', 'VB')

times_df_small |> ggplot(aes(x = n, y = value, linetype = variable, group = variable)) +
  geom_line() +
  geom_point() +
  theme_bw(base_size = 15) +
  xlab("Number of observations") +
  ylab("Running time (seconds)") +
  labs(linetype = "Method") +
  scale_x_continuous(breaks = c(50,250,500, 750, 1000), labels = c(50,250,500, 750, 1000))



times_df_large <- times_df[, -3] |> filter(n %in% c(5000,10000, 15000, 20000))
times_df_large <- times_df_large |> reshape2::melt(id.vars = 'n')
levels(times_df_large$variable) <- c('MCMC', 'VB')

times_df_large |> ggplot(aes(x = n, y = value, linetype = variable, group = variable)) +
  geom_line() +
  geom_point() +
  theme_bw(base_size = 15) +
  xlab("Number of observations") +
  ylab("Running time (seconds)") +
  labs(linetype = "Method") +
  scale_x_continuous(breaks = c(5000,10000,15000, 20000), labels = c(5000,10000, 15000, 20000))+
  scale_y_continuous(breaks = c(500, 1000, 1500))

