#setwd("~/Documents/LAG Model")
library(ggplot2)
library(tidyr)
library(latex2exp)
library(dplyr)
library(gridExtra)

####### Boxplots comparing model accuracy

graphest = read.csv('graphest.csv')
otherest = read.csv('consdat2.csv')
newdist = read.csv('newdist2.csv')
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
  scale_x_discrete(labels=c('jaccard'="Jaccard index", "spear" = "Spearman Corr.")) + 
  ggtitle("Graph Similarity - Sparse") + theme(plot.title = element_text(size = 11))
p12 = ggplot(filter(longdat, type=="dense"), aes(x=parameter, y=estimate)) + 
  geom_boxplot() + ylim(0,1) + theme_bw() + ylab("") + xlab("") + 
  scale_x_discrete(labels=c('jaccard'="Jaccard Index", "spear" = "Spearman Corr.")) + 
  ggtitle("Graph Similarity - Dense")+ theme(plot.title = element_text(size = 11))


long_oe = gather(otherest, parameter, estimate, c(beta,gamma,theta))

p21 = ggplot(filter(long_oe,type=="sparse"), aes(x=parameter, y=estimate)) + geom_boxplot() + 
  ylim(0,2.5) + ylab("Relative Error") + xlab("") + theme_bw() + 
  ggtitle("Estimation Error - Sparse")+ theme(plot.title = element_text(size = 11))

p22 = ggplot(filter(long_oe,type=="dense"), aes(x=parameter, y=estimate)) + geom_boxplot() + 
  ylim(0,2.5) + ylab("") + theme_bw() + xlab("") + 
  ggtitle("Estimation Error - Dense")+ theme(plot.title = element_text(size = 11))

g=arrangeGrob(p11,p12, p21, p22, nrow=2, ncol=2)

ggsave("this.pdf",g)



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



