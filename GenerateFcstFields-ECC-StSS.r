
###
#   Code for generating the StSS, ECC-Q, and ECC-mQ-SNP ensemble forecast fields for the entire verification period
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


par.climo <- array(dim=c(nxa,nya,nyears,nlt,3))
par.reg <- array(dim=c(nxa,nya,nyears,nlt,6))

mu.fcst <- sigma.fcst <- shift.fcst <- array(dim=c(nxa,nya,nyears,31,nlt))

ranks.ecc <- integer(nxa*nya*nyears*31*nlt*11)
dim(ranks.ecc) <- c(nxa,nya,nyears,31,nlt,11)
ranks.stss <- integer(nxa*nya*nyears*31*nlt*11)
dim(ranks.stss) <- c(nxa,nya,nyears,31,nlt,11)


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

		obs.train <- apcp.anal[,,train.ind.anal]
		fcst.train <- apcp.fcst[,,,train.ind]
		dim(fcst.train) <- c(nxf*nyf,nmb,n.train)
		fcst.verif <- apcp.fcst[,,,verif.ind]
		dim(fcst.verif) <- c(nxf*nyf,nmb,n.verif)

		cl.avg.anal <- matrix(NA,nxa,nya)
		cl.avg.fcst <- rep(NA,nxf*nyf)


	        # --- estimate analysis and forecast climatological average

		for (ix in 1:nxa)  {
			for (jy in 1:nya)  {
				if (mask.domain[ix,jy]) next
				if (sum(!is.na(obs.train[ix,jy,]))<200) next
				cl.avg.anal[ix,jy] <- mean(obs.train[ix,jy,],na.rm=TRUE)
			}
		}
		for (igrid in 1:(nxf*nyf))  {
			cl.avg.fcst[igrid] <- mean(fcst.train[igrid,,],na.rm=TRUE)
		}


		# --- Interpolate forecasts to CCPA grid, pick hostorical dates, and calculate ECC and StSS ranks

		loc.ip <- cbind(lon[!mask.domain],lat[!mask.domain])
		for (iday in 1:n.verif)  {
			jday <- ifelse(month==2 & iday==29, 28, iday)
			stss.ind <- which( ((dates%/%100)%%100)==jday & ((dates%/%10000)%%100)==month & (dates%/%1000000)!=years[iyear] )
			stss.ind.anal <- match(dates.fcstb[stss.ind],yyyymmddhh_begin)
			gefs.ip.verif <- apcp.hist.obs <- array(dim=c(nxa,nya,nmb))
			for (imb in 1:nmb)  {
				obj <- list(x=lons.fcst[,1], y=lats.fcst[1,], z=apcp.fcst[,,imb,verif.ind[iday]])
				gefs.ip.verif[,,imb][!mask.domain] <- interp.surface(obj, loc.ip)
				apcp.hist.obs[,,imb][!mask.domain] <- apcp.anal[,,stss.ind.anal[imb]][!mask.domain]
			}
			gefs.ip.verif.ranks <- apply(gefs.ip.verif,c(1,2),rank,ties.method="random")
			apcp.hist.obs.ranks <- apply(apcp.hist.obs,c(1,2),rank,ties.method="random")
			for (imb in 1:nmb)  {
				ranks.ecc[,,iyear,iday,ilt,imb] <- gefs.ip.verif.ranks[imb,,]
				ranks.stss[,,iyear,iday,ilt,imb] <- apcp.hist.obs.ranks[imb,,]
			}
		}


		# --- Fit CSGD model to local climatology, adjusted forecasts, calculate ensemble statistics, and fit CSGD regression model

		for (ix in 1:nxa)  {
			for (jy in 1:nya)  {
				if (mask.domain[ix,jy]) next

 			      # Fit observation climatology
				obs.mean <- mean(obs.train[ix,jy,obs.train[ix,jy,]>0],na.rm=TRUE)
				obs.pop <- mean(obs.train[ix,jy,]>0,na.rm=TRUE)
				sigma <- obs.mean

				if (obs.pop<0.005)  {
					par.climo[ix,jy,iyear,ilt,] <- c(0.0005,0.0182,-0.00049)	# Extremely dry location
					mu.fcst[ix,jy,iyear,1:n.verif,ilt] <- 0.0005			#  -> climatological forecast
					sigma.fcst[ix,jy,iyear,1:n.verif,ilt] <- 0.0182
					shift.fcst[ix,jy,iyear,1:n.verif,ilt] <- -0.00049
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
					par.climo[ix,jy,iyear,ilt,] <- par0				#  -> climatological forecast
					mu.fcst[ix,jy,iyear,1:n.verif,ilt] <- mu			#    with slighly better climatology
					sigma.fcst[ix,jy,iyear,1:n.verif,ilt] <- sigma
					shift.fcst[ix,jy,iyear,1:n.verif,ilt] <- shift
					next
				}
				par.climo[ix,jy,iyear,ilt,] <- optim(par0, crps.climo, obs=obs.train[ix,jy,], method="L-BFGS-B", lower=par0*c(0.5,0.5,2), upper=par0*c(2,2,0.1))$par


			      # Normalize (multiplicatively) the forecasts at all forecast grid points within a neighborhood of this location
				dst2 <- outer( (lons.fcst[,1]-lon[ix,jy])^2, (lats.fcst[1,]-lat[ix,jy])^2, '+')
				nbh.ind <- which(as.vector(dst2)<=rho2)
				nnbh <- length(nbh.ind)

				fcst.bc.train <- sweep(fcst.train[nbh.ind,,],1,cl.avg.fcst[nbh.ind],'/')
				fcst.bc.verif <- sweep(fcst.verif[nbh.ind,,],1,cl.avg.fcst[nbh.ind],'/')

			      # Calculate ensemble statistics with prediction error based weights
				ensmean.bc.train <- apply(fcst.bc.train,c(1,3),mean)
				obs.bc.train <- obs.train[ix,jy,] / cl.avg.anal[ix,jy]
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
				if (iyear>1 & !all(is.na(par.reg[ix,jy,1:(iyear-1),ilt,1])))  {
					par0 <- apply(par.reg[ix,jy,1:(iyear-1),ilt,,drop=FALSE],5,mean,na.rm=TRUE)
				}

				par.opt <- optim(par0, crps.reg,
					obs = obs.train[ix,jy,],
					enspop = enspop.train,
					ensmean = ensmean.train,
					ensmeandiff = ensmeandiff.train,
					par.climo = par.climo[ix,jy,iyear,ilt,],
					method = "L-BFGS-B",
					control = list(maxit=ifelse(iyear==1,20,10)),
					lower = c(0.001,0.05,0.0,0.0,0.1,0.0),
					upper = upper)$par

				par.reg[ix,jy,iyear,ilt,] <- par.opt

				mu.cl    <- par.climo[ix,jy,iyear,ilt,1]
				sigma.cl <- par.climo[ix,jy,iyear,ilt,2]
				shift.cl <- par.climo[ix,jy,iyear,ilt,3]

				log.arg <- par.opt[2] + par.opt[3]*enspop.verif + par.opt[4]*ensmean.verif
				mu.fcst[ix,jy,iyear,1:n.verif,ilt] <- mu.cl*log1p(expm1(par.opt[1])*log.arg)/par.opt[1]
				sigma.fcst[ix,jy,iyear,1:n.verif,ilt] <- par.opt[5]*sigma.cl*sqrt(mu.fcst[ix,jy,iyear,1:n.verif,ilt]/mu.cl) + par.opt[6]*sigma.cl*ensmeandiff.verif
				shift.fcst[ix,jy,iyear,1:n.verif,ilt] <- shift.cl
			}
		}
	}
	save(mask.domain, par.climo, par.reg, mu.fcst, sigma.fcst, shift.fcst, ranks.ecc, ranks.stss, file=paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata",sep=""))
}










###
#   Generate ECC-mQ-SNP trajectories


library(ncdf4)
library(fields)

month <- 1

years <- 2002:2013

lead.seq <- 8:11

nyears <- length(years)
nlt <- length(lead.seq)

filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/precip_06h_CCPA_2p5km_RRB.nc'
prec.nc <- nc_open(filename)
lon <- ncvar_get(prec.nc, varid="lons")-360
lat <- ncvar_get(prec.nc, varid="lats")
nc_close(prec.nc)

nxa <- nrow(lon)
nya <- ncol(lon)
nmb <- 11

mask.conus <- matrix(is.na(map.where("state",lon-0.03,lat)), nxa, nya)
mask.domain <- mask.conus | (lon < -123.9845) | (lon > -122.1095) | (lat < 38.38457) | (lat > 39.78888)

loc.ip <- cbind(lon[!mask.domain],lat[!mask.domain])
n.domain <- nrow(loc.ip)

filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/data/refcstv2_precip_ccpav3_066_to_072.nc'
prec.nc <- nc_open(filename)
lons.fcst <- ncvar_get(prec.nc, varid="lons_fcst")
lats.fcst <- ncvar_get(prec.nc, varid="lats_fcst")
nc_close(prec.nc)


lons.fcst.ss <- lons.fcst[4:9,1]
lats.fcst.ss <- lats.fcst[1,29:33]
nx.ss <- length(lons.fcst.ss)
ny.ss <- length(lats.fcst.ss)
grid.fcst.ss <- as.matrix(expand.grid(lons.fcst.ss,lats.fcst.ss))

delta.lon <- 0.46875
delta.lat <- 0.46810

r.lon <- 3*delta.lon
r.lat <- 3*delta.lat


# Define basis functions for simulating 'negative precipitation'

dst.lon <- as.matrix(dist(grid.fcst.ss[,1]))
tricube.lon <- ifelse(dst.lon>r.lon,0,(1-(dst.lon/r.lon)^3)^3)
dst.lat <- as.matrix(dist(grid.fcst.ss[,2]))
tricube.lat <- ifelse(dst.lat>r.lat,0,(1-(dst.lat/r.lat)^3)^3)
basis <- tricube.lon*tricube.lat
basis <- sweep(basis,1,apply(basis,1,sum),'/')


lambda <- 0.5

load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata",sep=""))

apcp.eccqm.snp.traj <- array(dim=c(nxa,nya,nyears,31,nlt,nmb))

for (ilt in 1:nlt)  {
	cleadb <- formatC(lead.seq[ilt]*6,width=3,flag="0")    	  	# beginning of the accumulation period
	cleade <- formatC((lead.seq[ilt]+1)*6,width=3,flag="0")		# end of the accumulation period

	filename <- paste('~/Desktop/Russian-River-CaseStudy/data/refcstv2_precip_ccpav3_',cleadb,'_to_',cleade,'.nc',sep='')
	cat(paste('Loading ',filename,'\n'))
	prec.nc <- nc_open(filename)
	dates <- ncvar_get(prec.nc, varid="yyyymmddhh_init")
	dates.fcstb <- ncvar_get(prec.nc, varid="yyyymmddhh_fcstb")
	apcp.fcst <- ncvar_get(prec.nc, varid="apcp_fcst_ens")
	nc_close(prec.nc)

	for (iyear in 1:nyears)  {
		print(c(iyear,nyears))
		verif.ind <- which( ((dates%/%10000)%%100) == month & (dates%/%1000000) == years[iyear] )
		n.verif <- length(verif.ind)

		for (iday in 1:n.verif)  {
			gefs.ss <- apcp.fcst[4:9,29:33,,verif.ind[iday]]
			shape.cal <- (mu.fcst[,,iyear,iday,ilt][!mask.domain]/sigma.fcst[,,iyear,iday,ilt][!mask.domain])^2
			scale.cal <- mu.fcst[,,iyear,iday,ilt][!mask.domain]/shape.cal
			shift.cal <- shift.fcst[,,iyear,iday,ilt][!mask.domain]
			qt.cal <- sweep(sweep(outer(shape.cal,((1:nmb)-0.5)/nmb,FUN=function(X,Y) qgamma(Y,shape=X)),1,scale.cal,'*'),1,shift.cal,'+')

		       # Simulate 'negative precipitation' where the raw ensemble forecasts are zero
			for (imb in 1:11)  {
				precip.sim.coef <- rep(0,nx.ss*ny.ss)
				zero.ind <- (as.vector(gefs.ss[,,imb])==0)
				n.zero <- sum(zero.ind)
				if (n.zero==0) next
				precip.sim.coef[zero.ind] <- runif(n.zero,-1,0)
				gefs.ss[,,imb][zero.ind] <- apply(sweep(basis[zero.ind,,drop=FALSE],2,precip.sim.coef,'*'),1,sum)
			}

		       # Interpolate augmented raw ensemble fields to analysis grid
			gefs.ip.augm <- matrix(NA,n.domain,nmb)
			for (imb in 1:nmb)  {
				obj <- list(x=lons.fcst.ss, y=lats.fcst.ss, z=gefs.ss[,,imb])
				gefs.ip.augm[,imb] <- interp.surface(obj, loc.ip)
			}

		       # Construct piecewise linear quantile-mapping function
			ecct.field <- matrix(0,n.domain,nmb)
			for (igp in 1:n.domain)  {
				l <- sum(qt.cal[igp,]<=0)
				if (l==nmb)  next					# all forecast quantiles are zero -> map to zero
				if (l==(nmb-1))  {					# only one non-zero forecast quantile -> map largest GEFS member to this quantile
					ind.max <- which.max(gefs.ip.augm[igp,])
					ecct.field[igp,ind.max] <- qt.cal[igp,nmb]
					next
				}							# otherwise: fit (regularized) linear regression spline
				ind.pos <- order(gefs.ip.augm[igp,])[l:nmb]
				n.pos <- length(ind.pos)
				x <- gefs.ip.augm[igp,ind.pos]
				if (length(unique(x))<n.pos)  {
					x <- x + runif(n.pos,0,0.01)	# resolve ties if necessary
				}
				y <- qt.cal[igp,l:nmb]
				A <- cbind(1,x,pmax(outer(x,x[c(-1,-n.pos)],'-'),0))
				coef <- solve(crossprod(A,A)+lambda*diag(c(0,0,pmax(1,y[c(-1,-n.pos)]))),crossprod(A,y))
				ecct.field[igp,ind.pos] <- pmax(A%*%coef,0)
			}

		       # Map the interpolated, augmented raw ensemble forecasts to the modified predictive CSGD sample
			for (imb in 1:nmb)  {
				apcp.eccqm.snp.traj[,,iyear,iday,ilt,imb][!mask.domain] <- ecct.field[,imb]
			}
		}
	}
}


save(apcp.eccqm.snp.traj, file=paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-mQ-SNP-",month,".Rdata",sep=""))


