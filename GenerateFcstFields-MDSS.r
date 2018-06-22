
###
#   Code for generating the MDSS-RO and MDSS-SDA ensemble forecast fields for the entire verification period
#
#      written by Michael Scheuerer, Feb 2018
#



library(ncdf4)
library(fields)
library(maps)

source("~/Desktop/QPF-T2M-MV/AuxiliaryFunctions.r")


month <- 1

years <- 2002:2013

lead.seq <- 8:11

nyears <- length(years)
nlt <- length(lead.seq)

filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/precip_06h_CCPA_2p5km_RRB.nc'
prec.nc <- nc_open(filename)
yyyymmddhh_begin <- ncvar_get(prec.nc, varid="yyyymmddhh_begin")
lon <- ncvar_get(prec.nc, varid="lons")-360
lat <- ncvar_get(prec.nc, varid="lats")
apcp.anal <- ncvar_get(prec.nc, varid="apcp_anal")
nc_close(prec.nc)

nda <- length(yyyymmddhh_begin)
nxa <- nrow(lon)
nya <- ncol(lon)

mask.conus <- matrix(is.na(map.where("state",lon-0.03,lat)), nxa, nya)
mask.domain <- mask.conus | (lon < -123.9845) | (lon > -122.1095) | (lat < 38.38457) | (lat > 39.78888)


filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/data/refcstv2_precip_ccpav3_066_to_072.nc'
prec.nc <- nc_open(filename)
lons.fcst <- ncvar_get(prec.nc, varid="lons_fcst")
lats.fcst <- ncvar_get(prec.nc, varid="lats_fcst")
nc_close(prec.nc)

coords.grid <- as.matrix(expand.grid(lons.fcst[5:8,1],lats.fcst[1,30:32]))[-1,]
nloc <- nrow(coords.grid)


###
#   Associate CCPA grid points with 1/2 degree grid and scale up

apcp.anal.upsc <- matrix(NA,nloc,nda)
ccpa.gpt.ind <- vector(nloc,mode='list')

idxa <- row(lon)
idya <- col(lat)
dim(apcp.anal) <- c(nxa*nya,nda)


dst2 <- outer(as.vector(lon),coords.grid[,1],'-')^2 + outer(as.vector(lat),coords.grid[,2],'-')^2
dst2[as.vector(mask.domain),] <- NA

I <- apply(dst2,1,function(x) ifelse(all(is.na(x)),NA,which.min(x)))
dim(I) <- c(nxa,nya)

for (iloc in 1:nloc)  {
	id.gpt <- !is.na(I) & (I==iloc)
	ccpa.gpt.ind[[iloc]] <- cbind(as.vector(idxa)[id.gpt],as.vector(idya)[id.gpt])
	apcp.anal.upsc[iloc,] <- apply(apcp.anal[id.gpt,],2,mean,na.rm=TRUE)
}
dim(apcp.anal) <- c(nxa,nya,nda)



###
#   Fit univariate forecast distributions


crps.reg <- function(par, obs, enspop, ensmean, ensmeandiff, par.climo)  {
	miss <- is.na(obs)

	mu.cl    <- par.climo[1]
	sigma.cl <- par.climo[2]
	shift.cl <- par.climo[3]

	log.arg <- par[2] + par[3]*enspop[!miss] + par[4]*ensmean[!miss]
	mu      <- mu.cl*log1p(expm1(par[1])*log.arg)/par[1]
	sigma   <- par[5]*sigma.cl*sqrt(mu/mu.cl) + par[6]*sigma.cl*ensmeandiff[!miss]
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


par.climo <- array(dim=c(nloc,nyears,nlt,3))
par.reg <- array(dim=c(nloc,nyears,nlt,6))

mu.fcst <- sigma.fcst <- shift.fcst <- array(dim=c(nloc,nyears,31,nlt))
obs.verif <- array(dim=c(nxa,nya,nyears,31,nlt))

train.ind.anal.bckp <- integer(nyears*nlt*91*(nyears-1))
dim(train.ind.anal.bckp) <- c(nyears,nlt,91*(nyears-1))


for (ilt in 1:nlt)  {
	cleadb <- formatC(lead.seq[ilt]*6,width=3,flag="0")    	  	# beginning of the accumulation period
	cleade <- formatC((lead.seq[ilt]+1)*6,width=3,flag="0")		# end of the accumulation period
	rho2 <- 1 + lead.seq[ilt]/4

	filename <- paste('~/Desktop/Russian-River-CaseStudy/data/refcstv2_precip_ccpav3_',cleadb,'_to_',cleade,'.nc',sep='')
	cat(paste('Loading ',filename,'\n'))
	prec.nc <- nc_open(filename)
	dates <- ncvar_get(prec.nc, varid="yyyymmddhh_init")
	dates.fcstb <- ncvar_get(prec.nc, varid="yyyymmddhh_fcstb")
	apcp.fcst <- ncvar_get(prec.nc, varid="apcp_fcst_ens")
	nc_close(prec.nc)

	nxf <- dim(apcp.fcst)[1]
	nyf <- dim(apcp.fcst)[2]
	nmb <- dim(apcp.fcst)[3]

	for (iyear in 1:nyears)  {
		print(c(iyear,nyears))

		mid.ind <- which( dates%/%1000000 != years[iyear] & (dates%/%10000)%%100 %in% month & (dates%/%100)%%100 == 15 )
		if(month==1 & iyear>1) mid.ind <- c(mid.ind,length(dates)+16)
		if(month==12 & iyear<nyears) mid.ind <- c(-16,mid.ind)
		train.ind <- as.vector(outer(seq(-45,45,1), mid.ind, '+'))
		train.ind <- train.ind[train.ind >= 1  &  train.ind <= length(dates)]
		train.ind.anal <- match(dates.fcstb[train.ind],yyyymmddhh_begin)
		n.train <- length(train.ind)

		verif.ind <- which( ((dates%/%10000)%%100) == month & (dates%/%1000000) == years[iyear] )
		verif.ind.anal <- match(dates.fcstb[verif.ind],yyyymmddhh_begin)
		n.verif <- length(verif.ind)

		train.ind.anal.bckp[iyear,ilt,1:n.train] <- train.ind.anal
		obs.verif[,,iyear,1:n.verif,ilt] <- apcp.anal[,,verif.ind.anal]

		obs.train <- apcp.anal.upsc[,train.ind.anal]
		fcst.train <- apcp.fcst[,,,train.ind]
		dim(fcst.train) <- c(nxf*nyf,nmb,n.train)
		fcst.verif <- apcp.fcst[,,,verif.ind]
		dim(fcst.verif) <- c(nxf*nyf,nmb,n.verif)

		cl.avg.anal <- rep(NA,nloc)
		cl.avg.fcst <- rep(NA,nxf*nyf)


		# --- estimate analysis and forecast climatological average

		for (iloc in 1:nloc)  {
			if (sum(!is.na(obs.train[iloc,]))<200) next
			cl.avg.anal[iloc] <- mean(obs.train[iloc,],na.rm=TRUE)
		}
		for (igrid in 1:(nxf*nyf))  {
			cl.avg.fcst[igrid] <- mean(fcst.train[igrid,,],na.rm=TRUE)
		}


		# --- Fit CSGD model to local climatology, adjusted forecasts, calculate ensemble statistics, and fit CSGD regression model

		for (iloc in 1:nloc)  {
	   	      # Fit observation climatology
			obs.mean <- mean(obs.train[iloc,obs.train[iloc,]>0],na.rm=TRUE)
			obs.pop <- mean(obs.train[iloc,]>0,na.rm=TRUE)
			sigma <- obs.mean

			if (obs.pop<0.005)  {
				par.climo[iloc,iyear,ilt,] <- c(0.0005,0.0182,-0.00049)		# Extremely dry location
				mu.fcst[iloc,iyear,1:n.verif,ilt] <- 0.0005			#  -> climatological forecast
				sigma.fcst[iloc,iyear,1:n.verif,ilt] <- 0.0182
				shift.fcst[iloc,iyear,1:n.verif,ilt] <- -0.00049
				next
			}
			for (mu in (40:1)*(sigma/40))  {
				shape <- (mu/sigma)^2
				scale <- mu/shape
				shift <- -qgamma(obs.pop,shape=shape,scale=scale,lower.tail=FALSE)
				if (shift > -mu/2) break
			}
			par0 <- c(mu,sigma,shift)

			if (obs.pop<0.02)  {							# Still dry
				par.climo[iloc,iyear,ilt,] <- par0				#  -> climatological forecast
				mu.fcst[iloc,iyear,1:n.verif,ilt] <- mu				#    with slighly better climatology
				sigma.fcst[iloc,iyear,1:n.verif,ilt] <- sigma
				shift.fcst[iloc,iyear,1:n.verif,ilt] <- shift
				next
			}
			par.climo[iloc,iyear,ilt,] <- optim(par0, crps.climo, obs=obs.train[iloc,], method="L-BFGS-B", lower=par0*c(0.5,0.5,2), upper=par0*c(2,2,0.1))$par

		      # Normalize (multiplicatively) the forecasts at all forecast grid points within a neighborhood of this location
			dst2 <- outer( (lons.fcst[,1]-coords.grid[iloc,1])^2, (lats.fcst[1,]-coords.grid[iloc,2])^2, '+')
			nbh.ind <- which(as.vector(dst2)<=rho2)
			nnbh <- length(nbh.ind)

			fcst.bc.train <- sweep(fcst.train[nbh.ind,,],1,cl.avg.fcst[nbh.ind],'/')
			fcst.bc.verif <- sweep(fcst.verif[nbh.ind,,],1,cl.avg.fcst[nbh.ind],'/')

		      # Calculate ensemble statistics with prediction error based weights
			ensmean.bc.train <- apply(fcst.bc.train,c(1,3),mean)
			obs.bc.train <- obs.train[iloc,] / cl.avg.anal[iloc]
			rmse <- apply(ensmean.bc.train, 1, function(x) sqrt(mean((obs.bc.train-x)^2,na.rm=TRUE)) )
			weights <- 1-((rmse-min(rmse))/diff(range(rmse)))
			weights <- weights/sum(weights)
			w <- array(weights/nmb, dim=c(nnbh,nmb))

			ensmean.train <- apply(sweep(fcst.bc.train, c(1,2), w, "*"), 3, sum, na.rm=TRUE)
			ensmeandiff.train <- apply(fcst.bc.train, 3, wgt.md, w=w, na.rm=TRUE)
			enspop.train <- apply(sweep(1*(fcst.bc.train>0), c(1,2), w, "*"), 3, sum, na.rm=TRUE)
			ensmean.verif <- apply(sweep(fcst.bc.verif, c(1,2), w, "*"), 3, sum, na.rm=TRUE)
			ensmeandiff.verif <- apply(fcst.bc.verif, 3, wgt.md, w=w, na.rm=TRUE)
			enspop.verif <- apply(sweep(1*(fcst.bc.verif>0), c(1,2), w, "*"), 3, sum, na.rm=TRUE)

			upper <- ifelse(rep(obs.pop<0.05,6), c(0.10,1.00,1.0,1.0,1.0,0.10), c(1.00,1.00,1.5,1.5,1.0,1.5))
			par0  <- ifelse(rep(obs.pop<0.05,6), c(0.05,0.05,1.0,0.8,0.5,0.05), c(0.05,0.05,1.0,0.8,0.5,0.1))
			if (iyear>1 & !all(is.na(par.reg[iloc,1:(iyear-1),ilt,1])))  {
				par0 <- apply(par.reg[iloc,1:(iyear-1),ilt,,drop=FALSE],4,mean,na.rm=TRUE)
			}

			par.opt <- optim(par0, crps.reg,
				obs = obs.train[iloc,],
				enspop = enspop.train,
				ensmean = ensmean.train,
				ensmeandiff = ensmeandiff.train,
				par.climo = par.climo[iloc,iyear,ilt,],
				method = "L-BFGS-B",
				control = list(maxit=ifelse(iyear==1,20,10)),
				lower = c(0.001,0.05,0.0,0.0,0.1,0.0),
				upper = upper)$par

			par.reg[iloc,iyear,ilt,] <- par.opt

			mu.cl    <- par.climo[iloc,iyear,ilt,1]
			sigma.cl <- par.climo[iloc,iyear,ilt,2]
			shift.cl <- par.climo[iloc,iyear,ilt,3]

			log.arg <- par.opt[2] + par.opt[3]*enspop.verif + par.opt[4]*ensmean.verif
			mu.fcst[iloc,iyear,1:n.verif,ilt] <- mu.cl*log1p(expm1(par.opt[1])*log.arg)/par.opt[1]
			sigma.fcst[iloc,iyear,1:n.verif,ilt] <- par.opt[5]*sigma.cl*sqrt(mu.fcst[iloc,iyear,1:n.verif,ilt]/mu.cl) + par.opt[6]*sigma.cl*ensmeandiff.verif
			shift.fcst[iloc,iyear,1:n.verif,ilt] <- shift.cl
		}
	}
}

rm(apcp.fcst)



###
#   Generate MDSS-SDA trajectories


d <- nloc*4
nmb <- 11
#nmb <- 55

k.seq <- c(rev(10+cumsum(rep(1:30,each=2))))
#k.seq <- c(rev(52+cumsum(c(rep(3:29,each=2),30))))

delta.lon <- 0.46875
delta.lat <- 0.46810

r.lon <- 1.5*delta.lon
r.lat <- 1.5*delta.lat


# Define basis functions for the multiplicative adjustment factor

dst.lon <- abs(outer(as.vector(lon),coords.grid[,1],'-'))
tricube.lon <- ifelse(dst.lon>r.lon,0,(1-(dst.lon/r.lon)^3)^3)
dst.lat <- abs(outer(as.vector(lat),coords.grid[,2],'-'))
tricube.lat <- ifelse(dst.lat>r.lat,0,(1-(dst.lat/r.lat)^3)^3)

basis.fcts <- tricube.lon*tricube.lat
basis.fcts[as.vector(mask.domain),] <- NA
dim(basis.fcts) <- c(nxa,nya,11)
basis.fcts <- sweep(basis.fcts,c(1,2),apply(basis.fcts,c(1,2),sum),'/')

dst.gefs <- as.matrix(dist(coords.grid))
n.nbh <- apply(dst.gefs>0 & dst.gefs<0.5, 1, sum)
L <- 1*(dst.gefs==0) - sweep(1*(dst.gefs>0 & dst.gefs<0.5),1,n.nbh,'/')


# Target function for fitting the coefficients of the multiplicative adjustment factor

adjm.fct.target <- function(alpha,M,xi)  {
	sum((M%*%alpha-xi)^2) + 0.25*sum((L%*%alpha)^2)
}

mdss.dates.fcstb <- integer(nyears*31*nmb)
dim(mdss.dates.fcstb) <- c(nyears,31,nmb)

ranks.mdss.ro <- integer(nxa*nya*nyears*31*nlt*11)
dim(ranks.mdss.ro) <- c(nxa,nya,nyears,31,nlt,11)

zeros.mdss.ro <- integer(nxa*nya*nyears*31*nlt)
dim(zeros.mdss.ro) <- c(nxa,nya,nyears,31,nlt)

apcp.mdss.traj <- array(dim=c(nxa,nya,nyears,31,nlt,nmb))


for (iyear in 1:nyears)  {
	print(c(iyear,nyears))

	use <- apply(!is.na(train.ind.anal.bckp[iyear,,]) & train.ind.anal.bckp[iyear,,]>0,2,all)
	n.train <- sum(use)

	apcp.obs.traj <- array(dim=c(nxa,nya,4,n.train))
	apcp.upsc.traj <- array(dim=c(nloc,4,n.train))

	for (ilt in 1:4)  {
		train.ind.anal <- train.ind.anal.bckp[iyear,ilt,use]
		apcp.obs.traj[,,ilt,] <- apcp.anal[,,train.ind.anal]
		apcp.upsc.traj[,ilt,] <- apcp.anal.upsc[,train.ind.anal]
	}

	miss <- apply(is.na(apcp.upsc.traj),3,any)
	if (any(miss))  {
		n.train <- sum(!miss)
		apcp.obs.traj <- apcp.obs.traj[,,,!miss] 
		apcp.upsc.traj <- apcp.upsc.traj[,,!miss] 
	}


	for (iday in 1:31)  {
		if (any(is.na(mu.fcst[,iyear,iday,]))) next
		e.fcst <- array(dim=c(d,nmb))

		shape <- (as.vector(mu.fcst[,iyear,iday,])/as.vector(sigma.fcst[,iyear,iday,]))^2
		scale <- as.vector(mu.fcst[,iyear,iday,])/shape
		shift <- as.vector(shift.fcst[,iyear,iday,])

		for (k in 1:nmb)  {
			e.fcst[,k] <- pmax(shift+qgamma((k-0.5)/nmb,scale=scale,shape=shape),0)
		}

		dim(apcp.upsc.traj) <- c(d,n.train)


	      # Select MDSS observation trajectories (using forecasts at 1/2 degree resolution)

		mdss.ind <- 1:n.train
		nss <- n.train

		for (k in k.seq)  {
			e.obs <- apcp.upsc.traj[,mdss.ind]
			o.obs <- t(apply(e.obs,1,order))

			eps.m <- matrix(NA,d,nss)
			eps.p <- matrix(NA,d,nss)

			for (j in 1:nss)  {
				diff.obs.fcst <- sweep(e.fcst,1,e.obs[,j],'-')
				eps.m[,j] <- apply(pmax(diff.obs.fcst,0),1,mean)
				eps.p[,j] <- apply(pmax(-diff.obs.fcst,0),1,mean)
			}
			eps.d <- eps.m - eps.p

			alpha <- ((1:(nss-1))-0.5)/(nss-1)
			T1 <- apply(eps.p,2,sum)
			dvg.loo <- rep(NA,nss)

			for (j in 1:nss)  {
				T2 <- sapply(1:d, function(i) mean(alpha*eps.d[i,o.obs[i,o.obs[i,]!=j]]))
				dvg.loo[j] <- 2*(mean(T1[-j])+sum(T2))
			}
			mdss.ind <- mdss.ind[tail(order(dvg.loo),k)]
			nss <- length(mdss.ind)
		}

		mdss.dates.fcstb[iyear,iday,] <- yyyymmddhh_begin[train.ind.anal.bckp[iyear,1,mdss.ind]]

		# matplot((1:d)-0.12,apcp.upsc.traj[,mdss.ind], type='p', pch=18, col=4, ylim=c(0,50))
		# matplot((1:d)+0.12, e.fcst, type='p', pch=18, col=2, add=TRUE)


	      # Generate 6-h MDSS forecast trajectories (post-processing at 1/2 degree resolution, spatial disaggregation via MDSS trajectories)

		e.obs <- apcp.upsc.traj[,mdss.ind]
		dim(e.obs) <- dim(e.fcst) <- c(nloc,4,nmb)
		e.obs.cpy <- e.obs

		apcp.hist.traj <- array(dim=c(nxa,nya,4,nmb))
		apcp.mdss.upsc <- array(dim=c(nloc,4,nmb))

		for (ilt in 1:4)  {

		       # Step 1: Use Schaake shuffle procedure to re-order forecasted quantiles at coarse scale

			for (iloc in 1:nloc)  {
				dst <- as.matrix(dist(coords.grid))[iloc,]
				wgt <- pmax(1-(dst/2)^2,0)
				wgt <- wgt/sum(wgt)

				apcp.mdss.ranks <- rep(NA,nmb)	
				zeros <- (e.obs[iloc,ilt,]<0.1)
				if (any(!zeros))  {
					apcp.mdss.ranks[!zeros] <- rank(e.obs[iloc,ilt,],ties.method="random")[!zeros]
				}
				for (i in which(zeros))  {
				      # replace zero observation with weighted average over a larger entire domain...
					e.obs[iloc,ilt,i] <- sum(e.obs.cpy[,ilt,i]*wgt)
				}
				if (any(zeros))  {
				      # ... and try ranking those upscaled values
					apcp.mdss.ranks[zeros] <- rank(e.obs[iloc,ilt,zeros],ties.method="random")
				}
				apcp.mdss.upsc[iloc,ilt,] <- e.fcst[iloc,ilt,apcp.mdss.ranks]
			}

		       # Step 2: Disaggregate these forecasts to fine scale using the historical templates and a spatially smooth adjustment function

			for (imb in 1:nmb)  {
				M <- matrix(NA,nloc,nloc)
				for (iloc in 1:nloc)  {
					gpt.ind <- ccpa.gpt.ind[[iloc]]
					apcp.hist.traj[cbind(gpt.ind,ilt,imb)] <- apcp.obs.traj[cbind(gpt.ind,ilt,mdss.ind[imb])]
					for (jloc in 1:nloc)  {
						M[iloc,jloc] <- mean(basis.fcts[cbind(gpt.ind,jloc)]*apcp.hist.traj[cbind(gpt.ind,ilt,imb)])
					}
				}
				xi <- apcp.mdss.upsc[,ilt,imb]
				alpha <- optim(rep(1,nloc), adjm.fct.target, M=M, xi=xi, method="L-BFGS-B", lower=rep(0,nloc))$par
				apcp.adjm.fctr <- apply(sweep(basis.fcts,3,alpha,'*'),c(1,2),sum)
				apcp.mdss.traj[,,iyear,iday,ilt,imb] <- apcp.adjm.fctr * apcp.hist.traj[,,ilt,imb]
			}

		       # Save out ranks of the historic ensemble for the reordering implementation of the MDSS method
			for (ix in 1:nxa)  {
				for (jy in 1:nya)  {
					if (mask.domain[ix,jy]) next
					zeros.mdss.ro[ix,jy,iyear,iday,ilt] <- sum(apcp.hist.traj[ix,jy,ilt,]==0)
					ranks.mdss.ro[ix,jy,iyear,iday,ilt,] <- rank(apcp.hist.traj[ix,jy,ilt,],ties.method="random")
				}
			}

		}
	}
}

save(lons.fcst, lats.fcst, mask.domain, apcp.mdss.traj, obs.verif, file=paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata",sep=""))
save(ranks.mdss.ro, zeros.mdss.ro, file=paste("~/Desktop/Russian-River-CaseStudy/data/HistoricRanks-MDSS-SDA-",month,".Rdata",sep=""))



# save(lons.fcst, lats.fcst, mask.domain, obs.verif, file=paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata",sep=""))




