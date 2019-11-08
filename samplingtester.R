source('LRmodel.R')

library(igraph)
library(tidyverse)
library(rlist)

n = 50 # 30 items
gamma = rnorm(n, -1)

g_sample = function(n, gamma) {
  edges <- c()
  for (i in 1:(n - 1)) { #    for each
    for (j in (i + 1):n) { # vertex pair
      if (rbinom(1, 1, plogis(gamma[i] + gamma[j])) == 1)
        edges <- c(edges, i, j)
    }
  }
  g1 <- make_empty_graph(n, directed = FALSE) %>% add_edges(edges)
  #chordalg1 = is_chordal(g1, fillin = T, newgraph = T)$newgraph
  #g3 = chordalg2$newgraph # ghat
  return(g1)
}

# sample a graph
g = g_sample(n, gamma)

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

nIter = 10000
reg = bayesLR(Y, X = design, nIter, c = 1, lsig = -3)

output = reg[2000:nIter,]

hist(output[,which.max(gamma)])

summary(glm(Y ~ design - 1, family = binomial))

