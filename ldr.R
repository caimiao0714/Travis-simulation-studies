## Experiment on the "ldr" package for doing Likelihood-Based Sufficient Dimension Reduction
set.seed(1810)
library(ldr)
data("snakes", package = "ldr")
str(snakes)
fit <- ldr(Sigmas = snakes[-3], ns = snakes[[3]], numdir = 4,
           model = "core", numdir.test = TRUE, verbose = TRUE, sim_anneal = TRUE, 
           max_iter = 100, max_iter_sa = 100)
summary(fit)
fit$Gammahat[[3]]

## Obtain the estimates \hat{M_1} and \hat{M_2} of M_1 and M_2 as follows":
Gammahat <- fit$Gammahat[[3]]
Sigmahat1 <- fit$Sigmashat[[3]][[1]]
Sigmahat2 <- fit$Sigmashat[[3]][[2]]
Mhat1 <- t(Gammahat) %*% Sigmahat1 %*% Gammahat
Mhat2 <- t(Gammahat) %*% Sigmahat2 %*% Gammahat

