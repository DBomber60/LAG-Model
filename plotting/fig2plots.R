library(ggplot2)
library(tidyr)
library(latex2exp)
library(dplyr)
library(gridExtra)
source('drivers.R')



##########################################################
################ GRAPH TRACE PLOTS #######################
##########################################################


# generate "true" graph and graph estimates
set.seed(1)
# PARAMETERS
n = 30 # items
nt = 500 # transactions
gamma <- rnorm(n, -2) # graph vertex coefficients
theta=.2

#################### SAMPLE FROM MODEL #########################

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
ghat1 = g_estimate(sampled$D, M=M, nt = nt) # lo pval threshold
ghat2 = g_estimate(sampled$D, M=M, nt = nt, pval = .001) # hi pval threshold

# read data
lowp = readRDS('data/traceG.12.11.RData')
hip = readRDS('data/traceG.12.12(hip).RData')

# generate FP/ FN data arrays
# number of initial false positive/ false negative edges
fng_low = difference(trueg, ghat1$graph) # G - hat(G)
fpg_low = difference(ghat1$graph, trueg) # hat(G) - G

fng_hi = difference(trueg, ghat2$graph) # G - hat(G)
fpg_hi = difference(ghat2$graph, trueg) # hat(G) - G


fpfa_low = gen_fpfa(lowp, fpg_low, fng_low)
fpfa_hi = gen_fpfa(hip, fpg_hi, fng_hi)


matdat_lo = do.call(rbind, lowp)
matdat_hi = do.call(rbind, hip)

x_lo = matdat_lo[,4]
x_hi = matdat_hi[,4]
y1.1 = fpfa_low$fpa # false positive array
y2.1 = fpfa_low$fna # false negative array

y1.2 = fpfa_hi$fpa # false positive array
y2.2 = fpfa_hi$fna # false negative array


##### g1 has false positives over 200 iterations ########

cols = c("cornflowerblue","darkgreen")

g1 = ggplot() + geom_step(mapping = aes(x = x_lo, y = y1.1[1:length(x_lo)]), color = cols[1], size=1 ) + 
  ylim(0,20) + ylab("Edge Count") + xlab("Iteration") + 
  geom_step(mapping = aes(x = x_hi, y = y1.2[1:length(x_hi)]), color = cols[2], size=1 ) + theme_bw() + 
  ggtitle(TeX("Graph Trace Plot: |E($\\hat{G}$) - E(G)|"))

g2 = ggplot() + geom_step(mapping = aes(x = x_lo, y = y2.1[1:length(x_lo)]), color = cols[1], size=1 ) + 
  ylim(0,20) + ylab("") + xlab("Iteration") +
  geom_step(mapping = aes(x = x_hi, y = y2.2[1:length(x_hi)]), color = cols[2], size=1 ) + theme_bw() + 
  ggtitle(TeX("Graph Trace Plot: |E(G) - E($\\hat{G}$)|"))

g=arrangeGrob(g1, g2, nrow=1, ncol=2)

ggsave("graphsample.pdf",g, dpi = 320, height = 4)




####### Boxplots comparing model accuracy

graphest = read.csv('./data/graphest.csv')
otherest = read.csv('./data/consdat2.csv')
newdist = read.csv('./data/newdist2.csv')
gsample = read.csv('./data/gsample.csv')

otherest = otherest[,-1]
graphest = graphest[,-1]
newdist = newdist[,-1]
names(newdist) = c("dense","sparse")
newdist$dense = 1 - newdist$dense
newdist$sparse = 1 - newdist$sparse

long_nd = gather(newdist, type, estimate)
long_nd$parameter = "Match Similarity"
long_ge = gather(graphest, parameter, estimate, c(jaccard, spear))
longdat = rbind(long_nd, long_ge)

p11 = ggplot(filter(longdat, type=="sparse"), aes(x=parameter, y=estimate)) + 
  geom_boxplot() + ylim(0,1) + theme_bw() + ylab("Similarity/ Correlation") + xlab("") + 
  scale_x_discrete(labels=c('jaccard'="Jaccard", "Match Similarity"="Match", "spear" = "Spearman")) + 
  ggtitle("Graph Similarity - Sparse") + theme(plot.title = element_text(size = 11))
p12 = ggplot(filter(longdat, type=="dense"), aes(x=parameter, y=estimate)) + 
  geom_boxplot() + ylim(0,1) + theme_bw() + ylab("") + xlab("") + 
  scale_x_discrete(labels=c('jaccard'="Jaccard", "Match Similarity"="Match", "spear" = "Spearman")) + 
  ggtitle("Graph Similarity - Dense")+ theme(plot.title = element_text(size = 11))


long_oe = gather(otherest, parameter, estimate, c(beta,gamma,theta))

p21 = ggplot(filter(long_oe,type=="sparse"), aes(x=parameter, y=estimate)) + geom_boxplot() + 
  ylim(0,2.5) + ylab("Relative Error") + xlab("") + theme_bw() + 
  ggtitle("Estimation Error - Sparse")+ theme(plot.title = element_text(size = 11))

p22 = ggplot(filter(long_oe,type=="dense"), aes(x=parameter, y=estimate)) + geom_boxplot() + 
  ylim(0,2.5) + ylab("") + theme_bw() + xlab("") + 
  ggtitle("Estimation Error - Dense")+ theme(plot.title = element_text(size = 11))


g=arrangeGrob(p11,p12,g1, p21, p22, g1, nrow=2, ncol=3)

ggsave("this2.pdf",g)







############################# Plots to compare FIM with the LAG ######################## 

library(dplyr)
library(ggplot2)
library(latex2exp)

type = rep('match',9)
comparison = rep(1:3,3)
estimate = rnorm(9)
df1 = data.frame(type=type, comparison=comparison, estimate=estimate)
type = rep('jaccard',9)
comparison = rep(1:3,3)
estimate = rnorm(9,1)
df2 = data.frame(type=type, comparison=comparison, estimate=estimate)
dat = rbind(df1, df2)
dat$comparison = as.factor(dat$comparison)

ggplot(filter(dat, type=="match"), aes(x=comparison, y=estimate)) +
  geom_boxplot() + 
  scale_x_discrete(labels=c('1'= parse(text = TeX('$M(G,\\hat{G})$')),
                            '2'= parse(text = TeX('$M(G,fim)$')),
                            '3'= parse(text = TeX('$M(\\hat{G},fim)')) )) + 
  theme_bw() + ylab("Similarity/ Correlation") + xlab("Comparison") + 
  ggtitle("Set Similarity Comparisons")



