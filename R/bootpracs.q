# part of R package boot
# copyright (C) 1997-2001 Angelo J. Canty
# corrections (C) 1997-2007 B. D. Ripley
#
# Unlimited distribution is permitted

# empirical log likelihood ---------------------------------------------------------

EL.profile <- function(y, tmin = min(y) + 0.1, tmax = max(y) - 0.1, n.t=25,
 		u = function(y, t) y - t ) {
#  Calculate the profile empirical log likelihood function
	EL.loglik <- function(lambda) {
		temp <- 1 + lambda * EL.stuff$u
		if (any(temp <= 0))
			out <- NA
		else	out <- - sum(log(1 + lambda * EL.stuff$u))
		out
		}
	EL.paras <- matrix(NA, n.t, 3)
	lam <- 0.001
	for(it in 0:(n.t-1)) {
        	t <- tmin + ((tmax - tmin) * it)/(n.t-1)
        	EL.stuff <- list(u = u(y, t))
#		assign("EL.stuff",EL.stuff, frame=1)
        	EL.out <- nlm(EL.loglik, lam)
		i <- 1
		while (EL.out$code > 2 && (i < 20)) {
			i <- i+1
			lam <- lam/5
			EL.out <- nlm(EL.loglik, lam)
		}
        	EL.paras[1 + it,  ] <- c(t, EL.loglik(EL.out$x), EL.out$x)
		lam <- EL.out$x
	}
	EL.paras[,2] <- EL.paras[,2]-max(EL.paras[,2])
	EL.paras
}




EEF.profile <- function( y, tmin=min(y)+0.1, tmax=max(y)-0.1, n.t=25,
		u=function(y,t) { y-t } ) {
	EEF.paras <- matrix( NA, n.t+1, 4)
	for (it in 0:n.t) {
		t <- tmin + (tmax-tmin)*it/n.t;
		psi <- as.vector(u( y, t ));
		fit <- glm(zero~psi -1,poisson(log));
		f <- fitted(fit);
		EEF.paras[1+it,] <- c( t, sum(log(f)-log(sum(f))), sum(f-1),
			coefficients(fit) );
	}
	EEF.paras[,2] <- EEF.paras[,2] - max(EEF.paras[,2]);
	EEF.paras[,3] <- EEF.paras[,3] - max(EEF.paras[,3]);
	EEF.paras
}

lik.CI <- function(like, lim ) {
#
#  Calculate an interval based on the likelihood of a parameter.
#  The likelihood is input as a matrix of theta values and the
#  likelihood at those points.  Also a limit is input.  Values of
#  theta for which the likelihood is over the limit are then used
#  to estimate the end-points.
#
#  Not that the estimate only works for unimodal likelihoods.
#
	L <- like[,2]
	theta <- like[,1]
	n <- length(L)
	i <- min(c(1:n)[L>lim])
	if (is.na(i))
		stop("likelihood never exceeds ",lim)
	j <- max(c(1:n)[L>lim])
	if (i==j)
		stop("likelihood exceeds ", lim, " at only one point")
	if (i==1) bot <- -Inf
	else {	i <- i+c(-1,0,1)
		x <- theta[i]
		y <- L[i]-lim;
		co <- coefficients( lm(y~x+x^2) );
		bot <- (-co[2]+sqrt( co[2]^2-4*co[1]*co[3] ) )/(2*co[3] )
	}
	if (j==n) top <- Inf
	else {	j <- j+c(-1,0,1)
		x <- theta[j]
		y <- L[j]-lim;
		co <- coefficients( lm(y~x+x^2) );
		top <- (-co[2]-sqrt( co[2]^2-4*co[1]*co[3] ) )/(2*co[3] );
	}
	out <- c( bot, top )
	names(out) <- NULL
	out
}

nested.corr <- function(data,w,t0,M) {
# Statistic for the example nested bootstrap on the cd4 data.
	corr.fun <- function(d, w=rep(1,nrow(d))/nrow(d)) {
		w <- w/sum(w)
		n <- nrow(d)
		m1 <- sum(d[,1]*w)
		m2 <- sum(d[,2]*w)
		v1 <- sum(d[,1]^2*w)-m1^2
		v2 <- sum(d[,2]^2*w)-m2^2
		rho <- (sum(d[,1]*d[,2]*w)-m1*m2)/sqrt(v1*v2)
		i <- rep(1:n,round(n*w))
		us <- (d[i,1]-m1)/sqrt(v1)
		xs <- (d[i,2]-m2)/sqrt(v2)
		L <- us*xs-0.5*rho*(us^2+xs^2)
		c(rho, sum(L^2)/nrow(d)^2)
	}
	n <- nrow(data)
	i <- rep(1:n,round(n*w))
	t <- corr.fun(data,w)
	z <- (t[1]-t0)/sqrt(t[2])
	nested.boot <- boot(data[i,],corr.fun,R=M,stype="w")
	z.nested <- (nested.boot$t[,1]-t[1])/sqrt(nested.boot$t[,2])
	c(z,sum(z.nested<z)/(M+1))
}
