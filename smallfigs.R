# small figures for Beta-Binomial Model
library(ggplot2)
library(reshape2)

updateLuckY <- function (n0, y0, tau, n){ (n0*y0+tau)/(n0+n) }
updateLuckN <- function (n0, n){ n0+n }

nyupdate <- function (pr, data){
  nn <- updateLuckN(pr[1], data[2])
  yn <- updateLuckY(pr[1], pr[2], data[1], data[2])
  c(nn,yn)
}

dbetany <- function(x, ny, ...){
  dbeta(x, shape1=ny[1]*ny[2], shape2=ny[1]*(1-ny[2]), ...)
}

sfn <- 16
sfpr <- c(8, 0.75) # n0, y0
sfdata <- c(10, sfn)  #  s,  n
sfpos <- nyupdate(sfpr, sfdata)

sfbetavec <- seq(0,1, length.out=101)
bdsfpr <- dbetany(sfbetavec, sfpr)
bdsfpo <- dbetany(sfbetavec, sfpos)

sfbinomvec <- 0:sfn
sfbinom <- dbinom(sfbinomvec, size=sfn, prob=0.75)

pdf("figs/smallfig-prior.pdf", width=2, height=2)
qplot(x=sfbetavec, y=bdsfpr, geom="line", ylim=c(0,4.5), xlab="p", ylab="f(p)")
dev.off()
pdf("figs/smallfig-posterior.pdf", width=2, height=2)
qplot(x=sfbetavec, y=bdsfpo, geom="line", ylim=c(0,4.5), xlab="p", ylab="f(p|s)")
dev.off()
pdf("figs/smallfig-binom.pdf", width=2, height=2)
qplot(x=sfbinomvec, y=sfbinom, geom="bar", stat="identity", xlab="s", ylab="f(s|p)")
dev.off()

#