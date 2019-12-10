source('drivers.R')
source('graphsample.R')
source('rewens.R')
source('graphplot.R')
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
changed2 = list()


for (j in 1:10) {
  print(j)
  # sample graph from chordal neighbors
  nbd = graph_sample(currentg, kshat = kshat, n, cn, beta = beta_true)
  
  # data frame for checking - can comment out
  datcheck = data.frame(nbd$kdiff, nbd$samplep, nbd$acceptance)
  
  # subtract the median sampling probability to prevent computational underflow
  probvec = exp(nbd$samplep - median(nbd$samplep))
  probvec = ifelse(probvec == Inf, 1e300, probvec)
  sampledg = sample.int(length(nbd$samplep), 1, prob = probvec)
  a = nbd$acceptance[sampledg] # acceptance function
  # if we accept the new sampled graph ...
  if(runif(1) < min(1,a)) {
    print("accept")
    if(nbd$added[sampledg]==1) {
      newg=add_edges(currentg, nbd$edge_change[[sampledg]])
    } else {newg = currentg - E(currentg,nbd$edge_change[[sampledg]])}
    
    # keep track of the changes made
    changed2 = list.append(changed2, c(nbd$edge_change[[sampledg]], nbd$added[sampledg], j))
    
    # update parameters
    ksupdate = sK_update(D = sampled$D, S=kshat$S, k=kshat$k, nbd$edge_change[[sampledg]], 
                         g = currentg, nbd$added[sampledg], cn = cn)
    
    currentg = newg
    Cl <- lapply(cliques(currentg), as.vector) # quick check: `table(sapply(C, length))`
    cn <- max(sapply(Cl, length)) - 1 # clique number of `g` (minus 1)
    kshat = ksupdate
  }
}


############ PLOTTING RESULTS (graphplot.R) #############


# read in a list from the graph sampling algorithm that contains: 1) edges that are changed
# 2) whether added or removed; 3) which iteration the change occured on

# traceG.RData contains the edges changed from a true graph simulated with set.seed(1) 
# in gsample_testing.R. gamma <- rnorm(n, -2)

#saveRDS(changed, file = "traceG.RData")
chglist = readRDS("data/traceG.12.9.8pm.RData")

# key inputs: true graph/ initial graph estimate
trueg = g
initg = ghat$graph

# number of initial false positive/ false negative edges
fng = difference(trueg, initg) # G - hat(G)
fpg = difference(initg, trueg) # hat(G) - G

fpfa = gen_fpfa(chgl = chgl.1, fpg, fng)
matdat = do.call(rbind, chgl.1)

x = matdat[,4]
y1 = fpfa$fpa # false positive array
y2 = fpfa$fna # false negative array

df = data.frame(iter = c(0,x), fpa = y1, fna = y2)

write.csv(df, "gsample.csv")





