# [ Notes on the hypergraph model for transactions ]

# Given a collection of `n` items, the model has two main levels: one that
# generates a latent *association* graph `G` on the items, for example,
#
#   logit(P((i, j) in G)) = gamma_i + gamma_j, for each pair of items i, j
#
# and, conditional on G, another level that samples the transactions as an
# Ewens mixture of hyperedges, that is, we first sample the number `K` of
# hyperedges in the transaction, with K ~ Ewens(n, theta), and then sample the
# hyperedges. The only possible hyperedges are *maximal* cliques in `G`, and
# for `c in C(G)`, each (maximal) clique `c` in the clique set of `G`,
# can be sampled with weight `pi_c` proportional to
#
#   exp(alpha_{|c|} + sum_{v in c} beta_v)
#
# that is, the sampled hyperedges correspond to non-zero entries in
# `H ~ MN(K; {pi_c}_{c in C(G)})`.

# inverse logit
ilogit <- function (x) 1 / (1 + exp(-x))

# log (unsigned) Stirling numbers of the first kind:
# s1(n + 1, k) = n * s1(n, k) + s1(n, k - 1) and s1(n, 1) = (n - 1)!
lstirling1 <- local({
  S <- list(c(0), c(0, 0))
  function (n) {
    ns <- length(S)
    if (ns < n) { # allocate?
      for (ni in (ns + 1):n) {
        Si <- numeric(ni); Si[1] <- lfactorial(ni - 1)
        Si1 <- S[[ni - 1]]
        for (ki in 2:(ni - 1)) {
          # Si[ki] <- lse2(Si1[ki - 1], log(ni - 1) + Si1[ki]) # TODO
          u <- log(ni - 1) + Si1[ki] - Si1[ki - 1]
          Si[ki] <- Si1[ki - 1] + log1p(exp(u))
        }
        S[[ni]] <<- Si
      }
    }
    S[[n]]
  }
})

# sample from Ewens distribution with parameter `theta`
rewens <- function (n, theta) {
  s <- integer(length(theta))
  if (length(n) == 1) {
    ls <- lstirling1(n)
    for (i in 1:length(theta)) {
      ti <- theta[i]
      p <- exp(ls + (1:n) * log(ti) - (lgamma(ti + n) - lgamma(ti)))
      s[i] <- sample.int(n, size = 1, prob = p)
    }
  }
  else {
    if (length(theta) != length(n))
      stop("inconsistent argument lengths")
    for (i in 1:length(theta)) {
      ls <- lstirling1(n[i]); ti <- theta[i]
      p <- exp(ls + (1:n[i]) * log(ti) - (lgamma(ti + n[i]) - lgamma(ti)))
      s[i] <- sample.int(n[i], size = 1, prob = p)
    }
  }
  s
}


# sample from model
model1 = function(gamma, n=100, alpha1, beta, theta, nt) {

  n <- 100 # #items
  
  # sample association graph `g`:
  library(igraph)
  #gamma <- rnorm(n, -2) # roughly logit(degree / n)
  edges <- c()
  for (i in 1:(n - 1)) { #    for each
    for (j in (i + 1):n) { # vertex pair
      if (rbinom(1, 1, ilogit(gamma[i] + gamma[j])) == 1)
        edges <- c(edges, i, j)
    }
  }
  g <- make_empty_graph(n, directed = FALSE) %>% add_edges(edges)
  
  # sample number of hyperedges `k`:
  C <- lapply(cliques(g), as.vector) # quick check: `table(sapply(C, length))`
  cn <- max(sapply(C, length)) - 1 # clique number of `g` (minus 1)
  X <- matrix(0, nrow = length(C), ncol = cn + n) # design
  for (i in seq_along(C)) {
    kc <- length(C[[i]]) - 1
    if (kc > 0) X[i, kc] <- 1
    X[i, cn + C[[i]]] <- 1 / (kc + 1)
  }
  
  alpha <- alpha1 * (1:cn) / cn # cardinality coefficients
  #beta <- rnorm(n) # item (vertex) coefficients
  #theta <- .25 # E[K] = sum(theta / (theta + 1:n - 1)) ~ 2.2
  #nt <- 100
  k <- rewens(n, rep(theta, nt)) # #cliques
  
  # sample single transaction
  sample_transaction = function(X,k,alpha, beta) { 
    transaction = array(0, dim = c(1,cn + n)) # transaction vector - row vector with element for each item
    rsk = as.numeric(exp(matrix(X[, 1:cn], ncol = cn) %*% alpha + X[, cn + 1:n] %*% beta))
    draw = sample.int(nrow(X), 1, prob = rsk)
    transaction = transaction + X[draw,]
    if (k>1) {
      for (i in 2:k) {
        # X_new contains rows with emtpy intersection with union of already sampled cliques
        X_new = X[which(apply(X, 1, function(x) length(intersect(which(X[draw,cn + 1:n]>0),which(x[(cn + 1):n]>0))))==0),] 
        rsk = as.numeric(exp(matrix(X_new[, 1:cn], ncol = cn) %*% alpha + X_new[, cn + 1:n] %*% beta))
        draw = sample.int(nrow(X_new), 1, prob = rsk)
        transaction = transaction + X_new[draw,]
        X = X_new
      }
    }  
    return(transaction)
  }
  
  # S holds binary transaction matrix 
  S = array(0, dim = c(nt,n+cn))
  
  for (i in 1:nt) {
    S[i,] = sample_transaction(X, k[i], alpha, beta)
  }
  
  S = S[,-c(1:cn)]
  S1 = apply(S, 1, function(x) ifelse(x>0,1,0))
  
  return(S1)
}

#S = model1(gamma = rnorm(n, -2), n=100, alpha1 = -1, beta = rnorm(100), theta = .25, nt = 100)
#hist(apply(S,1,sum))
#mean(apply(S,1,sum))




