library(dplyr)
library(mvtnorm)
library(truncnorm)


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
n_iter <- 1000
burn_in <- 500

# Run the Gibbs sampler
posterior_samples <- gibbs_sampler(data = y,
                                   mu_prior = mu_prior,
                                   g_prior = g_prior,
                                   e_prior = e_prior,
                                   lambda_prior = lambda_prior,
                                   inv_sigma2_prior = inv_sigma2_prior,
                                   n_iter = n_iter,
                                   burn_in = burn_in)

# Extract posterior samples
mu_posterior <- posterior_samples$mu
g_posterior <- posterior_samples$g
e_posterior <- posterior_samples$e
lambda_posterior <- posterior_samples$lambda
gamma_posterior <- posterior_samples$gamma
delta_posterior <- posterior_samples$delta
inv_sigma2_posterior <- posterior_samples$inv_sigma2


# mu_posterior |> hist()
# g_posterior  |> hist()
# e_posterior  |> hist()
# lambda_posterior |> hist()
# gamma_posterior  |> hist()
# delta_posterior  |> hist()
# inv_sigma2_posterior |> density() |> plot()

