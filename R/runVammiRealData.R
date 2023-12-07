rm(list = ls())

# loadRData <- function(fileName){
#   load(fileName)
#   get(ls()[ls() != "fileName"])
# }


# year2010 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2010_data.RData")
# year2011 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2011_data.RData")
# year2012 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2012_data.RData")
# year2013 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2013_data.RData")
# year2014 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2014_data.RData")
# year2015 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2015_data.RData")
# year2016 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2016_data.RData")
# year2017 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2017_data.RData")
# year2018 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2018_data.RData")
# year2019 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2019_data.RData")
#
# # everything (data proc.)
#
load("~/Documents/GitHub/vammi/vammi/Simulation/datIrl_procc.RData")
load("~/Documents/GitHub/vammi/vammi/Simulation/initial_values.RData")



y <- datIrl$y
genotype <- datIrl$g
environment <- datIrl$e
n_g = 85
n_e = 17
m_g = 100
s_g = 1
m_e = 100
s_e = 1
m_lambda = 60
s_lambda = 2
s_mu = 1
m_mu = 10
gamma <- gammas_deltas$gamma
delta <- gammas_deltas$delta
Q <- 1
n <- n_g * n_e
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
                        mu_lambda_sq * sum(mu_gamma_sq[genotype] * mu_delta_sq[environment]))/10
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
vb_blin <- mu_lambda * mu_gamma[genotype] * mu_delta[environment]
y_hat <-  mu_mu + mu_g[genotype] + mu_e[environment] + vb_blin
caret::RMSE(y,y_hat)




