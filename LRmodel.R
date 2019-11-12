# bayesian logistic regression with multivariate normal prior on beta parameters
# fit using a random walk MH algorithm (normal proposal distribution)
# Beta ~ N(mu_0, sigma); mu_0 set to 0 vector here; sigma is a scaled identity matrix


# logpost computes the log of the posterior density of beta vector
# input: current beta value/ binary response vector (Y)/ real data matrix (X)/ scaling factor on 
# prior sigma matrix (c) - bigger c implies larger steps through the posterior parameter space
# output: log of the posterior density of beta given data

logpostLR = function(beta_curr, Y, X, c) {
  t(X %*% beta_curr) %*% Y - sum(log(1 + exp(X %*% beta_curr))) - .5/c * crossprod(beta_curr)
}


# logposterior for cox prop hazards model
# data: S, C (clique set), k, current beta value, cn (clique number of current graph)

# beta_curr = rnorm(n)

logpostCOX = function(beta_curr, S, k, cn, Design) {
  denom = log (sum ( exp(Design %*% c(1:cn, beta_curr) ) ) )
  return ( sum (S %*% c(1:cn, beta_curr) - k * denom) )
}



# random walk metropolis hastings algorithm - multivariate normal proposal distribution
# input: logposterior - function that computes the log posterior density

rw_mh = function(Y, X, nIter, c, lsig, logposterior) {
  
  p = dim(X)[2] # number of beta parameters
  
  # initialize parameters
  beta_curr = rep(0, p) 
  logpi_curr = logposterior(beta_curr, Y, X, c)
  
  # keep track of results here
  Res = array(NA, dim = c(nIter, p))
  
  for (k in 1:nIter) {
    beta_prop = beta_curr + exp(lsig) * rnorm(p) # propose new beta vector from normal distribution
    logpi_prop = logposterior(beta_prop, Y, X, c)
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

# random walk metropolis hastings algorithm - multivariate normal proposal distribution
# input: logposterior - function that computes the log posterior density

rw_mh_COX = function(S, k, cn, Design, nIter, lsig) {
  
  p = n #dim(X)[2] # number of beta parameters
  
  # initialize parameters
  beta_curr = rep(0, p) 
  logpi_curr = logpostCOX(beta_curr, S, k, cn, Design)
  
  # keep track of results here
  Res = array(NA, dim = c(nIter, p))
  
  for (j in 1:nIter) {
    beta_prop = beta_curr + exp(lsig) * rnorm(p) # propose new beta vector from normal distribution
    logpi_prop = logpostCOX(beta_curr, S, k, cn, Design)
    acc = min(1, exp(logpi_prop - logpi_curr) )
    if (runif(1) < acc) {
      # set beta to the proposed beta
      beta_curr = beta_prop
      logpi_curr = logpi_prop
    }
    Res[j,] = beta_curr
  }
  return(Res)
}



# test above with simulated data - two coefficients
# set.seed(1)
# nobs = 200
# nIter = 10000
# X = matrix(rnorm(nobs * 2), nrow = nobs, ncol = 2)
# beta_true = c(1,2)
# Y = rbinom(nobs, 1, prob = plogis(X %*% beta_true) )
# 
# draws = rw_mh(Y, X, nIter, c = 1, lsig = 1, logposterior = logpostLR)

# proportion for burn in
#burn_prop = 0.15

# visualize posterior draws
#hist(draws[(nIter * burn_prop):nIter,2])


# test by fitting traditional logistic regression
m0 = glm(Y ~ X, family = binomial)
#summary(m0)
