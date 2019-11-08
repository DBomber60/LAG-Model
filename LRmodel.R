# bayesian logistic regression with multivariate normal prior on beta parameters
# fit using a random walk MH algorithm (normal proposal distribution)
# Beta ~ N(mu_0, sigma); mu_0 set to 0 vector here; sigma is a scaled identity matrix


# logpost computes the log of the posterior density of beta vector
# input: current beta value/ binary response vector (Y)/ real data matrix (X)/ scaling factor on 
# prior Sigma matrix (c)
# output: log of the posterior density of beta given data

logpost = function(beta_curr, Y, X, c) {
  t(X %*% beta_curr) %*% Y - sum(log(1 + exp(X %*% beta_curr))) - .5/c * crossprod(beta_curr)
}

# obtain posterior draws from beta

bayesLR = function(Y, X, nIter, c, lsig) {
  
  p = dim(X)[2] # number of beta parameters
  
  # initialize parameters
  beta_curr = rep(0, p) 
  logpi_curr = logpost(beta_curr, Y, X, c)
  
  # keep track of results here
  Res = array(NA, dim = c(nIter, p))
  
  for (k in 1:nIter) {
    beta_prop = beta_curr + exp(lsig) * rnorm(p) # propose new beta vector from normal distribution
    logpi_prop = logpost(beta_prop, Y, X, c)
    acc = min(1, exp(logpi_prop - logpi_curr) )
    if (runif(1) < acc) {
      # set beta to the proposed beta
      beta_curr = beta_prop
      logpi_curr = logpi_prop
    }
    Res[k,] = beta_curr
  }
  return(Res)
}

# test above with simulated data - two coefficients
set.seed(1)
nobs = 200
X = matrix(rnorm(nobs * 2), nrow = nobs, ncol = 2)
beta_true = c(1,2)
Y = rbinom(nobs, 1, prob = plogis(X %*% beta_true) )

draws = bayesLR(Y, X, nIter = 10000, c = 1, lsig = 1)

# visualize posterior draws
hist(draws[1000:nIter,2])


# test by fitting traditional logistic regression
m0 = glm(Y ~ X, family = binomial)
summary(m0)