source('drivers.R')
source('rewens.R')
library(igraph)
library(tidyverse)
library(rlist)

# input: graph/ kshat (list with k and S estimates)
# output: sampling probability of each possible edge change (neighboring chordal graphs)

graph_sample = function(gr_current, kshat, n, cn, beta) {

  # get clique number of the current graph
  #Cl <- lapply(cliques(gr_current), as.vector) # quick check: `table(sapply(C, length))`
  #cn <- max(sapply(Cl, length)) - 1 # clique number of `g` (minus 1)
  
  # get chordal neighbors of current graph
  neighbors = lazy_proposal(gr_current, n=n) 
  
  # initialize vectors to store sampling proobabiity (J(G^new | G^ current)) and  
  # acceptance function for each chordal neighbor
  samplep = array(0, dim = length(neighbors$edge_perturb))
  acceptance = array(0, dim = length(neighbors$edge_perturb))
  kdiff = array(NA, dim = length(neighbors$edge_perturb))
  
  # for each chordal neighbor, compute sampling probability and acceptance probability
  for (m in seq_along(neighbors$edge_perturb)) {
    #print(m)
    added = neighbors$added[m] # binary variable
    edge_change = neighbors$edge_perturb[[m]] # which edge is added/ deleted
    ks_temp = sK_update(D = sampled$D, S=kshat$S, k=kshat$k, edge_change, g = gr_current, added, cn = cn)
    kdiff[m] = sum(ks_temp$k - kshat$k)
    if(added==1) {gnew=add_edges(gr_current, edge_change)} else {gnew = gr_current - E(gr_current,edge_change)}
    
    # new design matrix
    newX = getDesign(gnew)
    f = log(sum(exp(newX %*% c(1:(dim(newX)[2]-n),beta))))
    
    # new cn value
    Cl <- lapply(cliques(gnew), as.vector) # quick check: `table(sapply(C, length))`
    cntemp <- max(sapply(Cl, length)) - 1 # clique number of `g` (minus 1)
    
    samplep[m] = -sum(ks_temp$k)*f # log (J(G^{t+1} | J(G^{t})))
    acceptance[m] = min(1e300,exp(sum( (ks_temp$k - kshat$k) * log(theta)  + ks_temp$S %*% c(1:cntemp,beta_true) - kshat$S %*% c(1:cn,beta_true)) ))
  }
  return(list(samplep=samplep, kdiff = kdiff, acceptance=acceptance, edge_change = neighbors$edge_perturb, added = neighbors$added))
}



# testing
a = c(-10, -20, -30)
k = 2
b = exp(a) /sum( exp(a) )
c = exp(a/k) / sum( exp(a/k) )
df = data.frame(a, b, c)

j = a + 10
df$j = exp(j) / sum( exp(j) )

# compute exp(-1000)
expm1(1e-10)

