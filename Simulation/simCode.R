# Simulation study
library(vammi)
rm(list = ls())

# Simulate data -----------------------------------------------------------
n_g <- c(6, 10, 12, 25, 50, 100, 200) # Number of genotypes
n_e <- c(10, 12, 20, 30, 50, 100) # Number of environments
s_g <- 10 # standard deviation of alpha
s_e <- 10 # standard deviation of alpha
s_y <- 1 # standard deviation of y
lambda <- c(12, 20, 25, c(12,20), c(12,25), c(20,25)) # values for lambda

# Set a seed to make it reproducible
set.seed(02)

# Generate data from the AMMI model used in the paper (data_ng_ne_lambda)
data_25_12_12 <- generate_data_AMMI(n_g[4], n_e[2], s_g, s_e, s_y, lambda[1])
data_25_12_20 <- generate_data_AMMI(n_g[4], n_e[2], s_g, s_e, s_y, lambda[2])




classical <- run_classical_AMMI(data_25_12_12, Q = 1)




set.seed(123)
n_g <- 8 # Number of genotypes
n_e <- 5 # Number of environment
n <- n_g * n_e # Number of obs altogether
alpha <- 0.1
beta <- 0.1
tau <- rgamma(1, alpha, beta) # Residual precision
omega <- 1
lambda <- 1 # Only eigenvalue
gamma <- c(1, rnorm(n_g - 1)) # Genotype effects - first has to be positive
delta <- rnorm(n_e) # Environment effects
x <- expand.grid(1:n_g, 1:n_e) # Create grid to look up values
genotype <- x[, 1]
environment <- x[, 2]

# Finally create y
blin <- y <- rep(NA, n)
for (i in 1:n) {
  blin[i] <- lambda * gamma[genotype[i]] * delta[environment[i]]
  y[i] <- rnorm(1,
                mean = blin[i],
                sd = 1 / sqrt(tau)
  )
}
df <- data.frame(genotype, environment, y)

