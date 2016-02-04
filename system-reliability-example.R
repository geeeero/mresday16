# system reliability examples

#install.packages("devtools")
#library("devtools")
#install_github("louisaslett/ReliabilityTheory")
library("ReliabilityTheory")
library(ggplot2)
library(reshape2)
#source("discrete-reliability-functions.R)

bottomlegend <- theme(legend.position = 'bottom', legend.direction = 'horizontal', legend.title = element_blank())
rightlegend <- theme(legend.title = element_blank())

# produces survival signature matrix for one component of type "name",
# for use in nonParBayesSystemInference()
oneCompSurvSign <- function(name){
  res <- data.frame(name=c(0,1), Probability=c(0,1))
  names(res)[1] <- name
  res
}

# produces data frame with prior and posterior component survival function
# for component of type "name" based on nonParBayesSystemInference() inputs
# for all components except survival signature; alpha, beta must be data frames
# where each column corresponds to the component type, so there must be a match;
# test.data must be a list where each element is named by component
# and contains the corresponding vector of failure times
oneCompPriorPost <- function(name, at.times, test.data, alpha, beta){
  sig <- oneCompSurvSign(name)
  nodata <- list(name=NULL)
  names(nodata) <- name
  a <- alpha[, match(name, names(alpha))]
  b <- beta[, match(name, names(beta))]
  data <- test.data[match(name, names(test.data))]
  prio <- nonParBayesSystemInference(at.times, sig, nodata, a, b)
  post <- nonParBayesSystemInference(at.times, sig,   data, a, b)
  data.frame(Time=rep(at.times,2), SysRel=c(prio, post),
             Item=rep(c("Prior", "Posterior"), each=length(at.times)))
}

ny2a <- function(n,y){
  n*y
}

ny2b <- function(n,y){
  n*(1-y)
}

# -----------------------------------------------------------------

set.seed(424242)
b0 <- graph.formula(s -- 2:3 -- 4 -- 5:6 -- 1 -- t, 2 -- 5, 3 -- 6)
b0 <- setCompTypes(b0, list("T1"=c(2,3,5,6), "T2"=c(4), "T3"=c(1)))
b0nulldata <- list("T1"=NULL, "T2"=NULL, "T3"=NULL)
b0testdata <- list("T1"=runif(20, min=0, max=10), #c(2:6),
                   "T2"=runif(10, min=0, max=10), #c(4.0, 4.2, 4.4, 4.6, 4.8, 5.0),
                   "T3"=runif(10, min=8, max=10)) #(14:19)/2) # T3 late failures
b0dat <- melt(b0testdata); names(b0dat) <- c("x", "Part")
b0dat$Part <- ordered(b0dat$Part, levels=c("T1", "T2", "T3", "System"))
b0sig <- computeSystemSurvivalSignature(b0)
b0t <- seq(0, 10, length.out=101)

b0a <- data.frame(T1=ny2a(n=rep(1,101), y=relbt1), T2=ny2a(n=rep(5,101), y=rel0), T3=ny2a(n=rep(5,101), y=rel3))
b0b <- data.frame(T1=ny2b(n=rep(1,101), y=relbt1), T2=ny2b(n=rep(5,101), y=rel0), T3=ny2b(n=rep(5,101), y=rel3))
b0a[b0a == 1] <- 1-1e-6
b0b[b0b == 0] <-   1e-6

b0T1 <- oneCompPriorPost("T1", b0t, b0testdata, b0a, b0b)
b0T2 <- oneCompPriorPost("T2", b0t, b0testdata, b0a, b0b)
b0T3 <- oneCompPriorPost("T3", b0t, b0testdata, b0a, b0b)
b0prio <- nonParBayesSystemInference(b0t, b0sig, b0nulldata, b0a, b0b)
b0post <- nonParBayesSystemInference(b0t, b0sig, b0testdata, b0a, b0b)

b0df <- rbind(data.frame(b0T1, Part="T1"), data.frame(b0T2, Part="T2"), data.frame(b0T3, Part="T3"),
              data.frame(Time=rep(b0t,2), SysRel=c(b0prio, b0post),
                         Item=rep(c("Prior", "Posterior"), each=length(b0t)), Part="System"))
b0df$Item <- ordered(b0df$Item, levels=c("Prior", "Posterior"))
b0df$Part <- ordered(b0df$Part, levels=c("T1", "T2", "T3", "System"))

tuered <- rgb(0.839,0.000,0.290)
tueblue <- rgb(0.000,0.400,0.800)
tueyellow <- rgb(1.000,0.867,0.000)

sr0 <- ggplot(b0df, aes(x=Time)) + scale_colour_manual(values = c(tuered, tueblue))
sr0 <- sr0 + geom_line(aes(y=SysRel, group=Item, colour=Item)) + facet_wrap(~Part, nrow=2)
sr0 <- sr0 + geom_rug(aes(x=x), data=b0dat) + xlab("Time") + ylab("Reliability")
pdf("figs/sysrelex.pdf", width=6, height=3.75)
sr0 + rightlegend
dev.off()

#