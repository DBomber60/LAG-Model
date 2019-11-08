a = c(1)

rm(list=ls())
set.seed(1)

# let's simulate some data - two coefficients
nobs = 200
X = matrix(rnorm(nobs * 2), nrow = nobs, ncol = 2)
beta_true = c(1,2)
Y = rbinom(nobs, 1, prob = plogis(X %*% beta_true) )

# prior parameters
mu_0 = c(0,0)
c = 1


# test by fitting traditional logistic regression
m0 = glm(Y ~ X, family = binomial)
summary(m0)

logpost = function(beta_curr, dat) {
  t(X %*% beta_curr) %*% dat - sum(log(1 + exp(X %*% beta_curr))) - .5/c * crossprod(beta_curr)
}

# simulate from multivariate normal


# random walk model
nIter = 10000
lsig = 1

# initialize parameters
beta_curr = c(0,0)
logpi_curr = logpost(beta_curr, dat = Y)

# keep track of results here
Res = array(NA, dim = c(nIter, 2))

for (k in 1:nIter) {
  beta_prop = beta_curr + exp(lsig) * rnorm(2)
  logpi_prop = logpost(beta_prop, dat = Y)
  acc = min(1, exp(logpi_prop - logpi_curr) )
  if (runif(1) < acc) {
    # set beta to the proposed beta
    beta_curr = beta_prop
    logpi_curr = logpi_prop
  }
  Res[k,] = beta_curr
}
Output = Res[2001:nIter,]
hist(Output[,2])


# sim from MVN
mu=c(1,2)
Sigma=matrix(c(1,0.9,0.9,1),2,2,byrow=T)
d=length(mu)
L=t(chol(Sigma)); Z=rnorm(d,0,1)
X=mu+L%*%Z