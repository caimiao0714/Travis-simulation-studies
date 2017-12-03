
# load packages ===========================================

library(MASS)
library(Matrix)
library(GMCM)
library(MatchIt)

# generate data ===========================================


# correlations

r1 = 0.3
r2 = 0.5
r3 = 0.8


# block diagnoal correlation matrix

m1 = matrix(r1, nrow=3, ncol=3)
diag(m1) = 1
m2 = matrix(r2, nrow=3, ncol=3)
diag(m2) = 1
m3 = matrix(r3, nrow=3, ncol=3)
diag(m3) = 1

cmat = bdiag(m1, m2, m3)


# covariates

x = mvrnorm(n=1000, mu=rep(0,9), Sigma=cmat)


# pt: the probability to draw the binary treatment
pt = GMCM:::inv.logit(rowSums(x))


# tr: treatment
tr = rbinom(n = 1000, size = 1, prob = pt)

# y
y = rnorm(n = 1000, mean = tr * 3, sd = 1)


# constructing the data.frame
dat <- data.frame(x, tr, y)




#------------------------------------------------------
# regression analysis
reg1t <- lm(y ~ tr + X1 + X2 + X3, data = dat)
summary(reg1t) # treatment, X1, X2, and X3

reg1 <- lm(y ~ X1 + X2 + X3 , data = dat)
summary(reg1) # X1, X2, and X3


reg2t <- lm(y ~ tr + X4 + X5 + X6 , data = dat)
summary(reg2t) # treatment, X4, X5, and X6

reg2 <- lm(y ~ X4 + X5 + X6 , data = dat)
summary(reg2) # X4, X5, and X6


reg3t <- lm(y ~ tr + X7 + X8 + X9 , data = dat)
summary(reg3t) # treatment, X4, X5, and X6

reg3 <- lm(y ~ X7 + X8 + X9 , data = dat)
summary(reg3) # X4, X5, and X6


regtotal <- lm(y ~ tr + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 , data = dat)
summary(regtotal) # all controlling variables



# METHOD 1.1 Principal Component Analysis: Stage 1 - integrate 9Xs in 1 PCA
pca1 <- prcomp(dat[1:9], scale = FALSE)

dat_pca1 <- data.frame(dat$y, dat$tr, pca$x[,1:3])
names(dat_pca1) <- c("y", "tr", paste(rep("pc1", 3), 1:3, sep = ""))

reg_pca1 <- lm(y ~ tr + pc11 + pc12 + pc13, data = dat_pca1)
summary(reg_pca1)

# METHOD 1.2 Principal Component Analysis: Stage 2 - integrate 9Xs in 3 PCAs
pca3_1 <- prcomp(dat[1:3], scale = FALSE)
pca3_2 <- prcomp(dat[4:6], scale = FALSE)
pca3_3 <- prcomp(dat[7:9], scale = FALSE)

dat_pca3 <- data.frame(dat$y, dat$tr, pca3_1$x[,1], pca3_2$x[,1], pca3_3$x[,1])
names(dat_pca3) <- c("y", "tr", paste(rep("pc3", 3), 1:3, sep = ""))

reg_pca3 <- lm(y ~ tr + pc31 + pc32 + pc33, data = dat_pca3)
summary(reg_pca3)

# METHOD 2 Propensity Score 
matchdata1.1 = matchit(tr ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = dat, method = "nearest", ratio = 1)

matchdata1.2 = matchit(tr ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = dat, method = "nearest", ratio = 2)

matchdata1.3 = matchit(tr ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = dat, method = "nearest", ratio = 3)

# question: why do we do the inv.logit of the RowSums of 9Xs?
# problem: PCA, I cannot distinguish the order of the three major components generated in stage 1. I don't know what PC11, PC12, PC13 correspond to. The results seem weired between METHOD1 and METHOD2.
# question: propensity score matching and weighting, do I use all 9 Xs or use the ones after the PCA.
# problem: the control observations are less than the treatment observations
