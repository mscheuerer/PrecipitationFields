
###
#   A collection of R-functions used by the CSGD method or for verification plots
#
#      written by Michael Scheuerer, Oct 2016
#


library(Hmisc)


pctg <- function(x,mu,sigma,shift)  {
	return(pgamma(x-shift, scale=(sigma^2)/mu, shape=(mu/sigma)^2))
}


qctg <- function(x,mu,sigma,shift)  {
	pmax(0, shift + qgamma( x, scale=(sigma^2)/mu, shape=(mu/sigma)^2) )
}



gini.md <- function(x,na.rm=FALSE)  {
   if(na.rm & any(is.na(x)))  x <- x[!is.na(x)]
   n <-length(x)
   return(4*sum((1:n)*sort(x,na.last=TRUE))/(n^2)-2*mean(x)*(n+1)/n)
}



wgt.md <- function(x,w,na.rm=FALSE)  {
	if(na.rm) {
		miss <- is.na(x) | w==0
		x <- x[!miss]
		w <- w[!miss]
	}
	x.ord <- order(as.vector(x))
	x.sort <- x[x.ord]
	W <- cumsum(w[x.ord])
	N <- length(W)
	2*sum(W[-N]*(1-W[-N])*diff(x.sort))
}


crps.climo <- function(par,obs)  {
	obs <- sort(obs[!is.na(obs)])
	n <- length(obs)
	k0 <- sum(obs==0)

	shape <- (par[1]/par[2])^2
	scale <- par[1]/shape
	shift <- par[3]

	crps <- numeric(n)

	c.std <- -shift/scale
	y.std <- (obs[(k0+1):n]-shift)/scale

	F.k.c <- pgamma(c.std, shape=shape)
	F.kp1.c <- pgamma(c.std, shape=shape+1)
	F.2k.2c <- pgamma(2*c.std, shape=2*shape)
	B.05.kp05 <- beta(0.5,shape+0.5)

	F.k.y <- pgamma(y.std, shape=shape)
	F.kp1.y <- pgamma(y.std, shape=shape+1)

	crps[1:k0]     <- c.std*(2*F.k.c-1) - shape*(2*F.kp1.c-1+F.k.c^2-2*F.kp1.c*F.k.c) - c.std*F.k.c^2 - (shape/pi)*B.05.kp05*(1-F.2k.2c)
	crps[(k0+1):n] <- y.std*(2*F.k.y-1) - shape*(2*F.kp1.y-1+F.k.c^2-2*F.kp1.c*F.k.c) - c.std*F.k.c^2 - (shape/pi)*B.05.kp05*(1-F.2k.2c)

	return( scale*mean(crps) )
}



crps.reg <- function(par, obs, enspop, ensmean, ensmeandiff, par.climo, ensmean.avg)  {
	miss <- is.na(obs)

	mu.cl    <- par.climo[1]
	sigma.cl <- par.climo[2]
	shift.cl <- par.climo[3]

	log.arg <- par[2] + par[3]*enspop[!miss] + par[4]*ensmean[!miss]/ensmean.avg
	mu      <- mu.cl*log1p(expm1(par[1])*log.arg)/par[1]
	sigma   <- par[5]*sigma.cl*sqrt(mu/mu.cl) + par[6]*ensmeandiff[!miss]
	shift   <- shift.cl

	scale <- sigma^2/mu
	shape <- (mu/sigma)^2

	y <- obs[!miss]
	c.std <- -shift/scale
	y.std <- (y-shift)/scale

	F.k.c <- pgamma(c.std, shape=shape)
	F.kp1.c <- pgamma(c.std, shape=shape+1)
	F.2k.2c <- pgamma(2*c.std, shape=2*shape)
	B.05.kp05 <- beta(0.5,shape+0.5)

	F.k.y <- pgamma(y.std, shape=shape)
	F.kp1.y <- pgamma(y.std, shape=shape+1)

	crps <- y.std*(2*F.k.y-1) - shape*(2*F.kp1.y-1+F.k.c^2-2*F.kp1.c*F.k.c) - c.std*F.k.c^2 - (shape/pi)*B.05.kp05*(1-F.2k.2c)

	return( mean(scale*crps) )
}




crps.normal <- function(par, obs, ensmeanano, ensvar, obs.climo)  {
	miss <- is.na(obs)

	mu <- obs.climo[!miss] + par[1]*ensmeanano[!miss]
	sigma <- sqrt( par[2] + par[3]*ensvar[!miss] )

   	obs.stdz <- (obs[!miss]-mu)/sigma
	crps <- sigma * ( obs.stdz*(2*pnorm(obs.stdz)-1) + 2*dnorm(obs.stdz) - 1/sqrt(pi) ) 

	return( mean(crps) )
}




plot.reliability <- function(n.day, x.day, y.day, log.freq=FALSE, N.boot=1000, n.max=NULL, threshold=NULL, cleadb=NULL, cleade=NULL)
{
   hist.freq <- apply(n.day,2,sum)
   if (is.null(n.max)) {
      n.max <- max(hist.freq)
   }
   rg.hist <- c(0,n.max)
   l <- nrow(n.day)

   n <- apply(n.day,2,sum)
   r <- apply(x.day,2,sum)/n
   f <- apply(y.day,2,sum)/n

   boot.perm <- matrix(sample(1:l, N.boot*l, replace=TRUE), N.boot, l)
   n.boot <- apply(boot.perm, 1, FUN=function(ind) apply(n.day[ind,],2,sum))
   x.boot <- apply(boot.perm, 1, FUN=function(ind) apply(x.day[ind,],2,sum))
   y.boot <- apply(boot.perm, 1, FUN=function(ind) apply(y.day[ind,],2,sum))
   r.boot <- x.boot / n.boot
   f.boot <- y.boot / n.boot

   f.q90 <- matrix(r,2,length(breaks)-1,byrow=TRUE) + apply(f.boot-r.boot, 1, quantile, prob=c(0.05,0.95), na.rm=TRUE)

#   var.ok <- apply(f.q90,2,diff) < 0.8
   var.ok <- hist.freq > 20
   plot(r[var.ok], f[var.ok], type="b", ylim=c(0,1), xlim=c(0,1), cex=0.8, xlab="forecast probabilities",
      ylab="observed frequencies", pch=19, col="blue", cex.lab=1.2, lwd=2)
   abline(a=0, b=1, col=1, lwd=1)
   for (k in (1:(length(breaks)-1))[var.ok])  {
      lines(rep(r[k],2), f.q90[,k], col="blue", lwd=1)
      lines(r[k]+c(-0.01,0.01), rep(f.q90[1,k],2), col="blue", lwd=2)
      lines(r[k]+c(-0.01,0.01), rep(f.q90[2,k],2), col="blue", lwd=2)
   }
#   title(paste("Reliability Diagram for  'precipitation >", threshold, "'"), cex=0.9)
   title(paste("Lead time: ",cleadb,"-",cleade,"h,  threshold:", threshold, sep=""), cex=0.9)

   obs.sum.day <- apply(y.day,1,sum)
   n.sum.day <- apply(n.day,1,sum)

#   N.boot <- apply(n.boot,2,sum)
#   obs.mean.boot <- apply(boot.perm, 1, FUN=function(ind) mean(obs.sum.day[ind])) / apply(boot.perm, 1, FUN=function(ind) mean(n.sum.day[ind]))
#   REL.boot <- apply(n.boot*(f.boot-r.boot)^2,2,sum,na.rm=TRUE) / N.boot
#   RES.boot <- apply(n.boot*(f.boot-obs.mean.boot)^2,2,sum,na.rm=TRUE) / N.boot
#   UNC.boot <- obs.mean.boot*(1-obs.mean.boot)
#   BS.boot <- REL.boot - RES.boot + UNC.boot
#   n.boot[n.boot==0] <- NA
#   corr.term.boot <- apply(n.boot/(n.boot-1)*f.boot*(1-f.boot),2,sum,na.rm=TRUE) / N.boot
#   REL.boot.corr <- REL.boot - corr.term.boot
#   RES.boot.corr <- RES.boot - corr.term.boot + UNC.boot/(N.boot-1)
#   UNC.boot.corr <- (UNC.boot*N.boot)/(N.boot-1)

#   REL.q90 <- quantile(REL.boot.corr,prob=c(0.05,0.95),na.rm=TRUE)
#   RES.q90 <- quantile(RES.boot.corr,prob=c(0.05,0.95),na.rm=TRUE)
#   UNC.q90 <- quantile(UNC.boot.corr,prob=c(0.05,0.95),na.rm=TRUE)
#   BS.q90 <- quantile(BS.boot,prob=c(0.05,0.95),na.rm=TRUE)

   if(log.freq)  {
      hist.freq <- log10(hist.freq)
      rg.hist <- c(0,log10(n.max))
   }

   plot.hist <- function()  {
      at <- 0:floor(rg.hist)[2]
      barplot(hist.freq, ylim=rg.hist*1.05, xaxt="n", yaxt="n", xlab='', ylab='', main='', col=ifelse(var.ok,'blue','lightblue'))
      axis(4,las=2,at=at,labels=format(10^at,scientific=TRUE), cex.axis=0.8, hadj=0.35, tck=-0.03)
      box()
   }
   subplot( plot.hist(), cnvrt.coords(0.02,0.98,'plt')$usr, size=c(0.8,0.8), vadj=1, hadj=0 )

#   text(0.65, 0.21, "5%", pos=4, cex=1.1)
#   text(0.87, 0.21, "95%", pos=4, cex=1.1)
#   text(0.45, 0.14, "REL:", pos=4, cex=1.1)
#   text(0.60, 0.14, formatC(signif(REL.q90[1],3), digits=5, format="f", flag="-", drop0=TRUE), pos=4, cex=1.1)
#   text(0.82, 0.14, formatC(signif(REL.q90[2],3), digits=5, format="f", flag="-", drop0=TRUE), pos=4, cex=1.1)
#   text(0.45, 0.07, "RES:", pos=4, cex=1.1)
#   text(0.60, 0.07, formatC(signif(RES.q90[1],3), digits=5, format="f", flag="-", drop0=TRUE), pos=4, cex=1.1)
#   text(0.82, 0.07, formatC(signif(RES.q90[2],3), digits=5, format="f", flag="-", drop0=TRUE) ,pos=4, cex=1.1)
#   text(0.45, 0.00, "UNC:", pos=4, cex=1.1)
#   text(0.60, 0.00, formatC(signif(UNC.q90[1],3), digits=5, format="f", flag="-", drop0=TRUE), pos=4, cex=1.1)
#   text(0.82, 0.00, formatC(signif(UNC.q90[2],3), digits=5, format="f", flag="-", drop0=TRUE), pos=4, cex=1.1)
}





avg.rank <- function(obs,fcst)  {                     ## Average ranks
	x.ranks <- apply(cbind(obs,fcst),1,rank)
	x.preranks <- apply(x.ranks,1,mean)
	return(rank(x.preranks,ties="random")[1])
}


bd.rank <- function(obs,fcst)  {                      ## Band depth ranks
	d <- nrow(fcst)
	m <- ncol(fcst)+1
	x.ranks <- apply(cbind(obs,fcst),1,rank)
	x.preranks <- apply((m-x.ranks)*(x.ranks-1),1,mean) + m - 1
	return(rank(x.preranks,ties="random")[1])
}


bd.rank.ties <- function(obs,fcst)  {                                  ## Band depth ranks (variant permitting ties)
	d <- nrow(fcst)
	m <- ncol(fcst)+1
	x.ranks <- apply(cbind(obs,fcst),1,rank,ties.method='max')
	n.ties <- apply(cbind(obs,fcst),1,function(x) apply(outer(x,x,'=='),1,sum))
	x.preranks <- apply(x.ranks*(m-x.ranks)+(x.ranks-1)*n.ties,1,mean)
	return(rank(x.preranks,ties="random")[1])
}


variogram.score <- function(obs,fcst,wgt)  {
	obs.diff <- sqrt(as.vector(dist(obs,"manhattan")))
	if (nrow(fcst)<3)  {
		fcst.diff <- mean(apply(fcst,2,function(x) sqrt(as.vector(dist(x,"manhattan")))))
	} else {
		fcst.diff <- apply(apply(fcst,2,function(x) sqrt(as.vector(dist(x,"manhattan")))), 1, mean)
	}
	return( sum(wgt*(obs.diff-fcst.diff)^2) )
}


energy.score <- function(obs,fcst)  {
	T1 <- mean(sqrt(apply(sweep(fcst,1,obs,"-")^2,2,sum)))
	T2 <- mean(sqrt(apply(apply(fcst,1,function(x) as.vector(dist(x,"manhattan")))^2,1,sum)))*(1-1/ncol(fcst))
	return(T1-0.5*T2)
}


