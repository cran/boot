isMatrix <- function(x) length(dim(x)) == 2
arima.sim <- function(n, model, rand.gen=rnorm, innov=rand.gen(n, ...), ...)
{
    x <- ts(c(rnorm(100), innov[1:n]), start = -99)
    if(length(model$ma)) x <- filter(x, c(1, model$ma), sides=1)
    if(length(model$ar)) x <- filter(x, model$ar, method="recursive")
    as.ts(x[-(1:100)])
}
