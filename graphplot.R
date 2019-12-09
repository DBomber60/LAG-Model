library(latex2exp)


set.seed(1)
# read in a list from the graph sampling algorithm that contains: 1) edges that are changed
# 2) whether added or removed; 3) which iteration the change occured on

saveRDS(changed, file = "traceG.RData")
chgl = readRDS("traceG.RData")

trueg = g
initg = ghat$graph

# number of initial false positive/ false negative edges
fng = difference(trueg, initg) # G - hat(G)
fpg = difference(initg, trueg) # hat(G) - G

get.edge.ids(fng, c(2,19))

# for each row in mat1, if the edge is added, see if the addition reduces the
# number of false negatives; for each that is removed, see if the deletion
# reduces the number of false positives
testd = mat1[1,]

# input: edge (vector), added (binary), falsepg (graph of false positives), falseng (graph of false negs)
# output: 
echg = function(edge, added, falsepg, falseng) {
  
  result = ""
  
  # if the edge was added and in the false negative (fn) graph,
  # then false negatives -= 1; else, false positives is incremented
  
  if(added == 1) {
    if(get.edge.ids(fng, edge) != 0) {
        result = "fnm1"
        falseng = delete_edges(falseng, edge)
      } else {
        result = "fpp1"
        falsepg = add_edges(falsepg, edge)
      }
    
  # if the edge was removed and is in the false positive (fp) graph,
  # then false positives -= 1; else, false negatives is incremented
    
    } else if(get.edge.ids(fpg, edge) != 0) {
        result = "fpm1" 
        faslepg = delete_edges(faslepg, edge)
    } else {
        result = "fnp1"
        falseng = add_edges(falseng, edge)
    }
  return(list(result = result, falsepg = falsepg, falseng = falseng))
}






echg(edge = c(testd[1], testd[2]), added = testd[3])



# fake data to generate graph 
colors = c("cornflowerblue","darkgreen")

mat1 = do.call(rbind, chgl)

chg = apply(mat1, 1, function(x) echg(edge = c(x[1], x[2]), added = x[3]))


fn_init = gsize(fng)
fp_init = gsize(fpg)

fns = c(fn_init)
fps = c(fp_init)

# need to change the actual graph as well 
gendat = function(fns, fps, chg) {
  for (k in seq_along(chg)) {
    if (chg[k] == "fnm1") {fns = c(fns, -1)
    } else if(chg[k] == "fpm1") {fps = c(fps, -1)
    } else if(chg[k] == "fnp1") {fns = c(fns, 1)
    } else {fps = c(fps,1)}
  }
  return(list(fns = fns, fps = fps))
}

changes = gendat(fns, fps, chg)

x = mat1[,4]
y = rpois(dim(mat1)[1], lambda = 3)

fn = stepfun(x, c(y,1))

plot(fn, xlim = c(1,200), frame = F, do.points = F,
     lty = 1, lwd = 4, col = colors[1], 
     main = "Graph Trace Plot",
     xlab = "Iteration", ylab = "Edge Count")


