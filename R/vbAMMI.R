rm(list = ls())

# Functions to meet the restrictions --------------------------------------

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

# Model G + E + BLIN -------------------------------------------------------------

# Simulate the data

n_g = 25
n_e = 10
m_g = 100
s_g = 1
m_e = 100
s_e = 1
s_y = 1
m_lambda = 60
s_lambda = 2
s_mu = 10
m_mu = 90

n <- n_g * n_e

set.seed(2022)
# Generate g (genotypes)
g <- rnorm(n_g, m_g, s_g)
g <- g - mean(g) # impose the sum-to-zero restriction
# Generate e (environments)
e <- rnorm(n_e, m_e, s_e)
e <- e - mean(e) # impose the sum-to-zero restriction
# Set the grand mean
mu <- rnorm(1, m_mu, s_mu)

lambda <- truncnorm::rtruncnorm(n = 1, a = 0, mean = m_lambda, sd = s_lambda)
# Number of components in the bilinear part
Q <- length(lambda)

# Generate gamma
gamma <- generate_gamma_delta(n_g, Q)

# Generate delta
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
y <- rnorm(n, mu_ij, s_y)

# Compute the response for the TEST set
y_test <- rnorm(n, mu_ij, s_y)

# hist(y)
# hist(blin)
# hist(g)
# hist(e)

genotype <- x[, "g"]
environment <- x[, "e"]

# Fitting the vb

n_iter <- 1000

# Initial values

# Starting values
alpha <- 100#10
beta <- 2#50

alpha_tau <- alpha + n / 2 - 1

tau_lambda <- 5
mu_lambda <- 40
mu_lambda_sq <- mu_lambda^2 + 1 / tau_lambda

tau_gamma <- 100
mu_gamma <- gamma[,1]
mu_gamma_sq <- mu_gamma^2 + 1 / tau_gamma

tau_delta <- 100
mu_delta <- delta[,1]
mu_delta_sq <- mu_delta^2 + 1 / tau_delta

tau_mu <- 100
mu_mu <- 90
mu_mu_sq <- mu_mu^2 + 1/tau_mu

tau_g <- 100
mu_g <-  rep(100,n_g)
mu_g_sq <- mu_g^2 + 1/tau_g

tau_e <- 100
mu_e <-  rep(10,n_e)
mu_e_sq <- mu_e^2 + 1/tau_e

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



qplot(g, mu_g) + geom_abline() + theme_bw()
qplot(e, mu_e) + geom_abline() + theme_bw()

# qplot(g, mu_g - mean(mu_g)) + geom_abline() + theme_bw()
# qplot(e, mu_e - mean(mu_e)) + geom_abline() + theme_bw()

vb_blin <- mu_lambda * mu_gamma[x[,1]] * mu_delta[x[,2]]
y_hat <-  mu_mu + mu_g[x[, "g"]] + mu_e[x[, "e"]] + vb_blin

qplot(blin, vb_blin) + geom_abline() + theme_bw()
qplot(y,y_hat) + geom_abline() + theme_bw()
qplot(blin, vb_blin) + geom_abline() + theme_bw()


caret::RMSE(y,y_hat)
caret::RMSE(blin, vb_blin)


qplot(e, mu_e) + geom_abline() + theme_bw()
qplot(g, mu_g) + geom_abline() + theme_bw()


# intervals

list_g <- list()
for(i in 1:n_g){
  values <- rnorm(100, mu_g[i], sqrt(1/tau_g))
  list_g$lower[[i]] <- quantile(values, 0.05)
  list_g$upper[[i]] <- quantile(values, 0.95)
}

ci_g <- data.frame(true = g, vb = mu_g, lower = unlist(list_g$lower), upper = unlist(list_g$upper))

ci_g %>% ggplot(aes(x = true, y = vb, ymax = upper, ymin = lower)) +
  geom_pointrange(colour = "steelblue", size = 0.25) +
  geom_abline() + theme_bw(base_size = 16) +
  xlab("Genotype effects") +
  ylab("VI")



list_e <- list()
for(i in 1:n_e){
  values <- rnorm(100, mu_e[i], sqrt(1/tau_e))
  list_e$lower[[i]] <- quantile(values, 0.05)
  list_e$upper[[i]] <- quantile(values, 0.95)
}

ci_e <- data.frame(true = e, vb = mu_e, lower = unlist(list_e$lower), upper = unlist(list_e$upper))

ci_e %>% ggplot(aes(x = true, y = vb, ymax = upper, ymin = lower)) +
  geom_pointrange(colour = "steelblue", size = 0.25) +
  geom_abline() + theme_bw(base_size = 16) +
  xlab("Environment effects") +
  ylab("VI")



values_lambda <- rnorm(10000, mu_lambda, sqrt(1/tau_lambda))
lower_lambda <- quantile(values_lambda, 0.05)
upper_lambda <- quantile(values_lambda, 0.95)

list_gamma <- list()
for(i in 1:n_g){
  values <- rnorm(10000, mu_gamma[i], sqrt(1/tau_gamma))
  list_gamma$lower[[i]] <- quantile(values, 0.05)
  list_gamma$upper[[i]] <- quantile(values, 0.95)
}

list_delta <- list()
for(i in 1:n_e){
  values <- rnorm(10000, mu_delta[i], sqrt(1/tau_delta))
  list_delta$lower[[i]] <- quantile(values, 0.05)
  list_delta$upper[[i]] <- quantile(values, 0.95)
}



lower_blin <- lower_lambda * unlist(list_gamma$lower[x[,1]]) * unlist(list_delta$lower[x[,2]])
upper_blin <- upper_lambda * unlist(list_gamma$upper[x[,1]]) * unlist(list_delta$upper[x[,2]])


ci_blin <- data.frame(true = blin, vb = vb_blin, lower = lower_blin, upper = upper_blin)

ci_blin %>% ggplot(aes(x = true, y = vb, ymax = upper, ymin = lower)) +
  geom_pointrange(colour = "steelblue", size = 0.25) +
  geom_abline() + theme_bw(base_size = 16) +
  xlab("Interaction effects") +
  ylab("VI")


values_mu <- rnorm(100, mu_mu, sqrt(1/tau_mu))
lower_mu <- quantile(values_mu, 0.05)
upper_mu <- quantile(values_mu, 0.95)

lower_y_hat <-  lower_mu + unlist(list_g$lower[x[, "g"]]) + unlist(list_e$lower[x[, "e"]]) + lower_blin
upper_y_hat <-  upper_mu + unlist(list_g$upper[x[, "g"]]) + unlist(list_e$upper[x[, "e"]]) + upper_blin

ci_y <- data.frame(true = y, vb = y_hat, lower = lower_y_hat, upper = upper_y_hat)

ci_y %>% ggplot(aes(x = true, y = vb, ymax = upper, ymin = lower)) +
  geom_pointrange(colour = "steelblue", size = 0.25) +
  geom_abline() + theme_bw(base_size = 16) +
  xlab("y") +
  ylab(expression(hat(y)))

