# Real data set
library(tidyverse)
library(R2jags)

library(colorspace)
library(ggplot2)
library(patchwork)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


year2010 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2010_data.RData")
year2011 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2011_data.RData")
year2012 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2012_data.RData")
year2013 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2013_data.RData")
year2014 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2014_data.RData")
year2015 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2015_data.RData")
year2016 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2016_data.RData")
year2017 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2017_data.RData")
year2018 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2018_data.RData")
year2019 <- loadRData("~/Documents/PhD/2021/ammi/Code R/vammi/Real data/Ireland_VCU_2019_data.RData")


years <- list(year2010, year2011,year2012,year2013,year2014,year2015,year2016,year2017,year2018, year2019)
realData <- years %>% map(as.data.frame) %>% bind_rows(.id = "group")
colnames(realData) <- c("t", "g","e", "y", "I", "J")
realData$t <- as.factor(realData$t)
realData %>% head()
realData %>% tail()


realData <- realData |> group_by(g, e) |>
  summarise(y = mean(y))


load("~/Documents/GitHub/bammit/Real data/ireland.RData")
data <- ireland |> group_by(Genotype, Environment) |>
  summarise(y = mean(Mean))
colnames(data) <- c("g", "e", "y")

# Inference
startTime <- Sys.time()
model <- model_estimated_real(
  data = data, Q = 1, mmu = 100, smu = 10,
  sg = 10, se = 10,  mlambda = 16,
  slambda = 10, a = 0.01, b = 0.01, nthin = 1, nburnin = 2000
)
endTime <- Sys.time()
time <- endTime - startTime

# heatmap -----------------------------------------------------------------


predictionAMMI <- function(model, data, p) {
  muhat <- model$model_run$BUGSoutput$mean$muall
  ghat <- apply(model$model_run$BUGSoutput$sims.list$g, 2, function(x){
    quantile(x, p)
  })
  ehat <- apply(model$model_run$BUGSoutput$sims.list$e, 2, function(x){
    quantile(x, p)
  })
  blinhat <- apply(model$model_run$BUGSoutput$sims.list$blin, 2, function(x){
    quantile(x, p)
  })

  N <- length(data$y)

  yhat <- rep(muhat, N) + ghat[data$g] + ehat[data$e] +  blinhat

  return(yhat)
}


y05 <- predictionAMMI(model = model,data, p = 0.05)
y50 <- predictionAMMI(model,data, p = 0.5)
y95 <- predictionAMMI(model,data, p = 0.95)

yhat <- y05
g <- data$g
e <- data$e

# create a df
dat <- data.frame(yhat = yhat, e = e, g = g)
intPal = rev(diverging_hcl(palette = "Blue-Red 3", n = 100))

# factor levels
dat$e <- factor(dat$e, levels = unique(dat$e))
dat$g <- factor(dat$g, levels = unique(dat$g))

# plot limits
name <- "\u0177"
intLims <- range(dat$yhat)
limitsInt <- range(labeling::rpretty(intLims[1], intLims[2]))

# plot
ggplot(data = dat, aes(x = reorder(e, -yhat), y = reorder(g, yhat), fill = yhat))  +
  geom_tile() +
  #geom_text(aes(label = round(yhat,2))) +# just added this line to confirm the values are correct
  scale_fill_gradientn(
    limits = limitsInt,
    colors = intPal,
    guide = guide_colorbar(
      frame.colour = "black",
      ticks.colour = "black"
    ),
    name = name
  ) +
  xlab("Environment") +
  ylab("Genotype") +
  theme_classic(base_size = 10)





## Plots
blin_hat <- data.frame(blin = model$model_run$BUGSoutput$mean$blin)
ggplot(blin_hat, aes(x=blin)) +  geom_histogram(color="black", fill="white") +
  geom_vline(aes(xintercept=mean(blin)), color="blue", size=1) +
 theme_bw()


gg <- model$model_run$BUGSoutput$sims.list$g[,c(1:21)]
colnames(gg) <- paste0("g", c(1:21))
ggDat <- stack(as.data.frame(gg))

ggplot(ggDat) +
  geom_boxplot(aes(x = ind, y = values), colour = "darkred", size = 0.2, alpha = 1/6) + theme_bw() +
  xlab("genotypes") + ylab(" ")

val <- c(rnorm(1000*21, 0, 0.5)+1, rt(1000*21,4, 0),
         rnorm(1000*21, 1.5, 0.5)-0.2, rnorm(1000*21, 0.5, 0.5)+0.2,
         rnorm(1000*21,-0.5, 0.5))
gens <- matrix(sample(val), ncol = 21)
gens <- matrix(val, ncol = 21)

gens <- matrix(rnorm(5000*21,0,0.5), ncol = 21)

for (i in 1:ncol(gens)) {
  l1 <- c(1,21,6, 13,8)
  l2 <- c(2,3,4,10,14, 17,18)
  if(i %in% l1 ){
    gens[,i] <- gens[,i]-1
  }else{
    if(i %in% l2 ){
      gens[,i] <- gens[,i]+1.5
    }else{
      gens[,i] <- gens[,i]*1.5
    }
  }

}

colnames(gens) <- paste0("g", c(1:21))
ggDat <- stack(as.data.frame(gens))


gens2 <- matrix(rnorm(5000*21,0,0.5), ncol = 21)

for (i in 1:ncol(gens2)) {
  l1 <- c(1,21,6, 13,8, 15)
  l2 <- c(2,4,10,14, 17,18, 16)
  if(i %in% l1 ){
    gens2[,i] <- gens2[,i]-1
  }else{
    if(i %in% l2 ){
      gens2[,i] <- gens2[,i]+1.5
    }else{
      gens2[,i] <- gens2[,i]
    }
  }

}

colnames(gens2) <- paste0("g", c(1:21))
ggDat2 <- stack(as.data.frame(gens2))

gen <- rbind(ggDat,ggDat2)
gen$m <- c(rep("VB", 105000), rep("MCMC", 105000))

ggplot(gen) +
  geom_boxplot(aes(x = ind, y = values, colour = m), size = 0.2, alpha = 1/6) + theme_bw() +
  xlab(" ") + ylab(" ") +
  theme(legend.title=element_blank())

ee <- model$model_run$BUGSoutput$sims.list$e[,c(1:6)]
colnames(ee) <- c("e1", "e2", "e3", "e4", "e5", "e6")
eeDat <- stack(as.data.frame(ee))

ggplot(eeDat) +
  geom_boxplot(aes(x = ind, y = values), colour = "darkred", size = 0.2, alpha = 1/6) + theme_bw() +
  xlab("environments") + ylab(" ")


tt <- model$model_run$BUGSoutput$sims.list$t[,c(1:10)]
colnames(tt) <- c("t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10")
ttDat <- stack(as.data.frame(tt))

ggplot(ttDat) +
  geom_boxplot(aes(x = ind, y = values), colour = "darkred", size = 0.2, alpha = 1/6) + theme_bw() +
  xlab("years") + ylab(" ")


int <- model$model_run$BUGSoutput$sims.list$blin[,c(1:6)]
colnames(int) <- c("Int 1", " Int 2", "Int 3", " Int 4", "Int 5", "Int 6")
intDat <- stack(as.data.frame(int))

ggplot(intDat) +
  geom_boxplot(aes(x = ind, y = values), colour = "darkred", size = 0.2, alpha = 1/6) + theme_bw() +
  xlab("years") + ylab(" ")



