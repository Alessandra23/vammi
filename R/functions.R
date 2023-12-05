# Simulate data

# These functions are from ambarti package (https://github.com/ebprado/AMBARTI)
#' @export
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

#' @description # function from the package AMBARTI
#' @export
run_classical_AMMI <- function(data, Q){

  x_train <- data$x
  y_train <- data$y

  g = as.factor(x_train[,"g"])
  e = as.factor(x_train[,"e"])

  # Fit the linear model
  linear_mod = aov(y_train ~ g + e + g:e)
  # linear_mod = lm(y_train ~ g + e)

  # Get the residuals for the interaction g:e
  #interaction_tab = matrix(residuals(linear_mod), ncol = length(unique(g)), length(unique(env)))
  interaction_tab = model.tables(linear_mod, type='effects', cterms = 'g:e')
  interaction_tab = interaction_tab$tables$`g:e`

  # Get the number of PCs
  if (is.null(data$Q) == FALSE) {Q = data$Q}

  # Run the Singular Value Decomposition (SVD) to compute lambda, gamma, and delta
  sv_dec <- svd(interaction_tab, nu = Q, nv = Q)

  # Get parameter estimates
  # mu_hat     = linear_mod$coefficients[1] # slightly biased compared to mean(y_train)
  mu_hat     = mean(y_train)
  g_hat      = aggregate(x = y_train - mu_hat, by = list(g), FUN = "mean")[,2]
  e_hat      = aggregate(x = y_train - mu_hat, by = list(e), FUN = "mean")[,2]
  lambda_hat = sv_dec$d[1:Q]
  gamma_hat  = -1*sv_dec$u
  delta_hat  = -1*sv_dec$v

  return(list(mu_hat     = mu_hat,
              g_hat      = g_hat,
              e_hat      = e_hat,
              lambda_hat = lambda_hat,
              gamma_hat  = gamma_hat,
              delta_hat  = delta_hat))
}

