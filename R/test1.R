I <- 3
J <- 2
set.seed(02)
x <- matrix(runif(I*J), nrow=I, ncol=J)
g <- runif(I)

# Calculate the sum using the outer() function
sum(outer(g, rowSums(x)))

g[1]*x[1,1] + g[1]* x[1,2] +
g[2]*x[2,1] + g[2]* x[2,2] +
g[3]*x[3,1] + g[3]*x[3,2]



g <- c(2,4,6)
e <- c(1,2)
gamma <- c(1,5,10)
delta <- c(2,4)

x <- expand.grid(1:I, 1:J)
names(x) <- c("g", "e") # g = genotype and e = envorinment
x$g <- as.factor(x$g)
x$e <- as.factor(x$e)
y <- g[x[, "g"]] + e[x[, "e"]]

gamma
delta

sum(y)
sum(g*y)
sum(delta*sum(g*gamma))

sum(g*gamma[x[, "g"]]*delta[x[, "e"]])

sum(y*gamma[x[, "g"]]*delta[x[, "e"]]) ## this is the right one

sum(delta*y*gamma)

dat <- matrix(y, nrow=I, ncol=J)

g[1]*dat[1,1] + g[1]* dat[1,2] +
g[2]*dat[2,1] + g[2]* dat[2,2] +
g[3]*dat[3,1] + g[3]* dat[3,2]


sum(g[1]*dat[1,] + g[2]*dat[2,] + g[3]*dat[3,])

v <- vector()
for (i in 1:I) {
  v[i] <- sum(g[i]*y[x[, "g"] == i])
}
v |> sum()


sum(y^2)

sum(y*g)

sum(y[x[, "g"] == 1])
dat[1,] |> sum()

ng <- I
ne <- J



# write code --------------------------------------------------------------

a <- 0.1

n <- length(y)
alpha_tau <- n + a


# update tau
beta_tau <- b + (sum(y^2) - 2*mu_mu*sum(y) - 2*sum(y*mu_g) - 2*sum(y*mu_e) -
                   mu_lambda*sum(y*mu_gamma*mu_delta) + n*mu_mu_sq + 2*mu_mu*ne*sum(mu_g) +
                   2*mu_mu*ng*sum(mu_e) + ne*sum(mu_g_sq) + 2*sum(mu_g*mu_e) + ng*sum(mu_e_sq) +
                   2*mu_mu*mu_lambda*sum(mu_gamma*mu_delta) + 2*mu_lambda*sum(mu_g*mu_gamma*mu_delta) +
                   2*mu_lambda*sum(mu_e*mu_gamma*mu_delta) + mu_delta_sq*sum(mu_gamma_sq*mu_delta_sq)
                 )/2
mu_tau <- alpha_tau/beta_tau

# update mu
tau_mu <- n*mu_tau + 1/s_mu_sq
mu_mu <- (m*(1/s_mu_sq) - mu_tau(sum(y) + ne*sum(g) + ng*sum(e) + mu_lambda*sum(mu_gamma*mu_delta)))/(2*(n*mu_tau + 1/s_mu_sq))
mu_mu_sq <- mu_mu^2 + 1/tau_mu


# update g
tau_g <- 2*(mu_tau*ne + 1/(s_g^2))
for(i in 1:ng){
  mu_g[i] <- mu_tau*(sum(y[genotype == i]) - mu_mu*ne - sum(mu_e) - mu_lambda*mu_gamma[i]*sum(mu_delta))/tau_g
}
mu_g_sq <- mu_g^2 + 1/tau_g

# update e
tau_e <- 2*(mu_tau*ng + 1/(s_e^2))
for(i in 1:ng){
  mu_g[i] <- mu_tau*(sum(y[genotype == i]) - mu_mu*ne - sum(mu_e) - mu_lambda*mu_gamma[i]*sum(mu_delta))/tau_g
}
mu_g_sq <- mu_g^2 + 1/tau_g







