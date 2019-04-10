x <- seq(-2, 5, by=0.1)
y <- 9*x + rnorm(length(x))
plot(y ~ x)
