setwd("~/Documents/LAG Model")
source('functions3.R')
source('~/Documents/LAG Model/sample/graphsample.R')
source('rewens.R')
library(igraph)
library(tidyverse)
library(rlist)

# toy example - generate a sample from a small set of items
# compute likelihood values of neighboring graphs
set.seed(1)
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

# estimate initial graph
M = t(sampled$D) %*% sampled$D
ghat = g_estimate(sampled$D, M=M, nt = nt)
kshat = kestimate(ghat$graph, sampled$D, cn = ghat$cn)


# graph sampling procedure
currentg = ghat$graph
cn = ghat$cn


# keep track of which edges were added/ removed
changed = list()


for (j in 1:10) {
  # sample graph from chordal neighbors
  nbd = graph_sample(currentg, kshat = kshat, n, cn, beta = beta_true)
  sampledg = sample.int(length(nbd$samplep), 1, prob = exp(nbd$samplep - median(nbd$samplep)))
  a = nbd$acceptance[sampledg] # acceptance function
  # if we accept the new sampled graph ...
  if(runif(1) < min(1,a)) {
    print("accept")
    if(nbd$added[sampledg]==1) {
      newg=add_edges(currentg, nbd$edge_change[[sampledg]])
    } else {newg = currentg - E(currentg,nbd$edge_change[[sampledg]])}
    
    # keep track of the changes made
    changed = list.append(changed, c(nbd$edge_change[[sampledg]], nbd$added[sampledg]))
    
    # update parameters
    ksupdate = sK_update(D = sampled$D, S=kshat$S, k=kshat$k, nbd$edge_change[[sampledg]], 
                         g = currentg, nbd$added[sampledg], cn = cn)
    
    currentg = newg
    Cl <- lapply(cliques(currentg), as.vector) # quick check: `table(sapply(C, length))`
    cn <- max(sapply(Cl, length)) - 1 # clique number of `g` (minus 1)
    kshat = ksupdate
  }
}

