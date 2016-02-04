# discrete reliability function examples
library(ggplot2)
library(reshape2)

rel2haz <- function(relvec){
  ptj   <- relvec[-length(relvec)]
  ptjp1 <- relvec[-1]
  (ptj - ptjp1)/ptj
}

haz2rel <- function(hazvec){
  ptj <- cumprod(1-hazvec)
  c(1, ptj)
}

relplotmaker <- function(xvec, relvec, name){
  pdf(file=paste("figs/", name, ".pdf", sep=""), width=4, height=2.5)
  #pl <- qplot(x=xvec, y=relvec, geom="line", ylim=c(0,1), xlab="t", ylab=expression(p[t]^k))
  df <- data.frame(x=xvec, y=relvec)
  pl <- ggplot(df, aes(x, y)) + geom_line(size=0.3) + geom_point(size=0.2) + 
    ylim(0,1) + xlab("t") + ylab(expression(p[t]^k))
  print(pl)
  dev.off()
}

relxvec <- seq(0, 10, length.out=101)
#hazxvec <- relxvec[-length(relxvec)]
hazxvec <- relxvec[-1]
rel0 <- c(rep(1, 5), rep(0.9, 5), rep(0.8,10), rep(0.75, 15), rep(0.65, 10), rep(0.6, 5),
          rep(0.35, 5), rep(0.3, 10), rep(0.2, 15), rep(0.1, 10), rep(0.05,11))
haz0 <- rel2haz(rel0)
qplot(x=relxvec, y=rel0, geom="line", ylim=c(0,1), xlab="t", ylab=expression(p[t]^k))
qplot(x=hazxvec, y=haz0, geom="line", ylim=c(0,0.2), xlab="t", ylab=expression(h[t]^k))
relplotmaker(relxvec, rel0, "discr-rel-0")

rel1 <- c(1.00-0.2*relxvec[1:11], 0.80-0.08*relxvec[2:31], 0.56-0.25*relxvec[2:11],
          0.31-0.1*relxvec[2:21], 0.11-0.03*relxvec[2:31])
plot(relxvec,rel1, type="l")
haz1 <- rel2haz(rel1)
qplot(x=relxvec, y=rel1, geom="line", ylim=c(0,1), xlab="t", ylab=expression(p[t]^k))
qplot(x=hazxvec, y=haz1, geom="line", ylim=c(0,0.1), xlab="t", ylab=expression(h[t]^k))
relplotmaker(relxvec, rel1, "discr-rel-1")

rel2 <- 1-pweibull(relxvec, shape=1, scale=2)
haz2 <- rel2haz(rel2)
qplot(x=relxvec, y=rel2, geom="line", ylim=c(0,1), xlab="t", ylab=expression(p[t]^k))
qplot(x=hazxvec, y=haz2, geom="line", ylim=c(0,0.1), xlab="t", ylab=expression(h[t]^k))
relplotmaker(relxvec, rel2, "discr-rel-2")

rel2h <- haz2rel(rep(0.05, 100))

rel3 <- 1-pweibull(relxvec, shape=2.5, scale=5)
haz3 <- rel2haz(rel3)
qplot(x=relxvec, y=rel3, geom="line", ylim=c(0,1), xlab="t", ylab=expression(p[t]^k))
qplot(x=hazxvec, y=haz3, geom="line", ylim=c(0,0.1), xlab="t", ylab=expression(h[t]^k))
relplotmaker(relxvec, rel3, "discr-rel-3")

rel4 <- 1-plnorm(relxvec, meanlog=1.5, sdlog=0.5)
haz4 <- rel2haz(rel4)
qplot(x=relxvec, y=rel4, geom="line", ylim=c(0,1), xlab="t", ylab=expression(p[t]^k))
qplot(x=hazxvec, y=haz4, geom="line", ylim=c(0,0.1), xlab="t", ylab=expression(h[t]^k))
relplotmaker(relxvec, rel4, "discr-rel-4")

rel5 <- 1-pgamma(relxvec, shape=2.5, rate=0.75)
haz5 <- rel2haz(rel5)
qplot(x=relxvec, y=rel5, geom="line", ylim=c(0,1), xlab="t", ylab=expression(p[t]^k))
qplot(x=hazxvec, y=haz5, geom="line", ylim=c(0,0.1), xlab="t", ylab=expression(h[t]^k))
relplotmaker(relxvec, rel5, "discr-rel-5")

hazbt1 <- c(0.05-0.02*(1:20)/10, rep(0.01,40), 0.01+0.03*(1:40)/10)
plot(hazxvec, hazbt1, type="l")
relbt1 <- haz2rel(hazbt1)
plot(relxvec, relbt1, type="l")
relplotmaker(relxvec, relbt1, "discr-rel-bt1")

hazbt2 <- c(0.08-0.03*(1:20)/10, rep(0.02,40), 0.02+0.02*(1:40)/10)
plot(hazxvec, hazbt1, type="l")
lines(hazxvec, hazbt2)
relbt2 <- haz2rel(hazbt2)
plot(relxvec, relbt1, type="l")
lines(relxvec, relbt2)


#