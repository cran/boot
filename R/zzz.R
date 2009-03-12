isMatrix <- function(x) length(dim(x)) == 2

ssample <- function(x) if(length(x) > 1) sample(x) else x
