x <- rnorm(32)
x.mra <- mra(x)
sum(x - apply(matrix(unlist(x.mra), nrow=32), 1, sum))^2

## Figure 4.19 in Gencay, Selcuk and Whitcher (2001)
data(ibm)     
ibm.returns <- diff(log(ibm))
ibm.volatility <- abs(ibm.returns)
## Haar
ibmv.haar <- mra(ibm.volatility, "haar", 1, "dwt")
names(ibmv.haar) <- c("d1", "d2", "d3", "d4", "s4")
## LA(8)
ibmv.la8 <- mra(ibm.volatility, "la8", 4, "dwt")
names(ibmv.la8) <- c("d1", "d2", "d3", "d4", "s4")
lines(ibmv.haar[[2]],col='red')
plot(ibm.volatility)
dev.off()
par(mfrow=c(1,1))
