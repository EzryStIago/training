

a <- rnorm(1)
b <- runif(1, 0, 4) 
c <- a+b

cat(paste(signif(a,3), "+", signif(b,3), "=", signif(c,3), sep=""))
