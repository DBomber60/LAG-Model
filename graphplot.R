library(latex2exp)
set.seed(1)


# for each row in mat1, if the edge is added, see if the addition reduces the
# number of false negatives; for each that is removed, see if the deletion
# reduces the number of false positives

# input: edge (vector), added (binary), falsepg (graph of false positives), falseng (graph of false negs)
# output: impact on graph estimate relative to true graph (false positives/ false negatvies)
echg = function(edge, added, falsepg, falseng) {
  
  result = ""
  
  # if the edge was added and in the false negative (fn) graph,
  # then false negatives -= 1; else, false positives is incremented
  
  if(added == 1) {
    ind = get.edge.ids(falseng, edge) # index of the added edge in the false negative graph (0 if not in it)
    if(ind != 0) { # if the edge is in the false negative graph (that is, adding it is good since its in gtrue)
        result = "fnm1"
        falseng = falseng %>% delete_edges(ind)
      } else {
        result = "fpp1"
        falsepg = add_edges(falsepg, edge)
      }
    
  # if the edge was removed and is in the false positive (fp) graph,
  # then false positives -= 1; else, false negatives is incremented
    
    } else {
     
        ind =  get.edge.ids(falsepg, edge) # index of the removed edge in the false positive graph
        if(ind != 0) { # if it is in the false pos graph
          result = "fpm1" 
          falsepg = falsepg %>% delete_edges(ind)
        } else {
            result = "fnp1"
            falseng = add_edges(falseng, edge)
        }
      }
  return(list(result = result, falsepg = falsepg, falseng = falseng))
}



# input: list of changes (each item in the list has the edge change/ added/ iteration)
# output: false positive array/ false negative array - length of the list of changes + 1 (for initial values)

gen_fpfa = function(chgl, fpg, fng) {
  
  # arrays to hold data on false positives/ negatives
  fn_init = gsize(fng)
  fp_init = gsize(fpg)

  fna = array(c(fn_init,rep(NA, length(chgl))), dim = c(length(chgl)+1,1))
  fpa = array(c(fp_init, rep(NA, length(chgl))), dim = c(length(chgl)+1,1))
  
  # for each change made to the graph 
  for (k in seq_along(chgl)) {
    dat = chgl[[k]]
    newg = echg(edge = c(dat[1], dat[2]), added = dat[3], fpg, fng)
    if (newg$result == "fpm1") {
      fpa[k+1] = fpa[k] - 1
      fna[k+1] = fna[k]
      fpg = newg$falsepg
    } else if(newg$result == "fpp1") {
      fpa[k+1] = fpa[k] + 1
      fna[k+1] = fna[k]
      fpg = newg$falsepg
    } else if(newg$result == "fnm1") {
      fpa[k+1] = fpa[k]
      fna[k+1] = fna[k] - 1 
      fng = newg$falseng
    } else {
      fpa[k+1] = fpa[k]
      fna[k+1] = fna[k] + 1 
      fng = newg$falseng
    }
  }
  
  return(list(fna=fna, fpa=fpa))

}
  


# colors to generate graph 
colors = c("cornflowerblue","darkgreen")

graphplot = function(iters, fpa, fna) {
  
  fn1 = stepfun(iters, fpa)
  fn2 = stepfun(iters, fna)

  plot(fn1, xlim = c(1,200), ylim = c(0, 15), frame = F, do.points = F,
       lty = 1, lwd = 4, col = colors[1],
       main = "Graph Trace Plot",
       xlab = "Iteration", ylab = "Edge Count")
  plot(fn2, add = T, do.points = F, lty = 3, lwd = 4, col = colors[2])

}

