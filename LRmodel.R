# bayesian logistic regression with multivariate normal prior or beta parameters
# Beta ~ N(mu_0, sigma); mu_0 set to 0 vector here; sigma is a scaled identity matrix


# logpost computes the log of the posterior density of beta vector
# input: current beta value/ binary response vector (Y)/ real data matrix (X)/ scaling factor on 
# prior Sigma matrix (c)
# output: log of the posterior density of beta given data

logpost = function(beta_curr, Y, X, c) {
  t(X %*% beta_curr) %*% Y - sum(log(1 + exp(X %*% beta_curr))) - .5/c * crossprod(beta_curr)
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


# test above with simulated data - two coefficients
set.seed(1)
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




# sim from MVN
mu=c(1,2)
Sigma=matrix(c(1,0.9,0.9,1),2,2,byrow=T)
d=length(mu)
L=t(chol(Sigma)); Z=rnorm(d,0,1)
X=mu+L%*%Z