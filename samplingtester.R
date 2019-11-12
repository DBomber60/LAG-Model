source('LRmodel.R')
source('drivers.R')
source('rewens.R')

set.seed(1) # for reproducibility
library(igraph)
library(tidyverse)
library(rlist)

# set parameters

# PARAMETERS
n = 20 # items
nt = 500 # transactions
gamma <- rnorm(n, -1.5) # graph vertex coefficients
theta=.2

#################### SAMPLE FROM MODEL/ ESTIMATE GRAPH/ GAMMA #########################

g <- g_sample(n, gamma)
# generate design matrix for sampling cliques
C <- lapply(cliques(g), as.vector) # quick check: `table(sapply(C, length))`
cn <- max(sapply(C, length)) - 1 # clique number of `g` (minus 1)

# parameters
alpha <- (1:cn) # cardinality coefficients 
beta_true <- rnorm(n) # item (vertex) coefficients
sampled = lag_sample(G=g, C=C, cn=cn, nt=nt, theta=theta, gamma = gamma, beta = beta_true, alpha = alpha) 


######### BETA SAMPLES ############
beta_curr = rnorm(n)

logpostCOX(beta_curr, S = sampled$S, k = sampled$k, cn, sampled$Design)

# MCMC sampler for gamma samples
nIter = 100000
reg = rw_mh_COX(S = sampled$S, k = sampled$k, cn = cn, Design = sampled$Design, nIter=nIter, lsig = -4)
output = reg[10000:nIter,]
hist(output[,which.max(beta_true)])
beta_true[which.max(beta_true)]




######### GAMMA SAMPLES ############

# build design matrix/ response vector based on graph
design = array(0, dim = c(choose(n,2),n))
Y = array(0, dim = c(choose(n,2),1))

# build design matrix & response vector
row_num = 1
for (i in 1:(n-1)) {
  for(j in (i+1):n) {
    design[row_num, c(i,j)] = 1
    if (get.edge.ids(g, c(i,j)) > 0) {Y[row_num] = 1}
    row_num = row_num + 1
  }
}

# MCMC sampler for gamma samples
nIter = 10000
reg = rw_mh(Y, X = design, nIter, c = 1, lsig = -3, logposterior = logpostLR)
output = reg[2000:nIter,]
hist(output[,which.max(gamma)])
# summary(glm(Y ~ design - 1, family = binomial)) # frequentist check




