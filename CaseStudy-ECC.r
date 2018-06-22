
###
#   Code for generating the ECC-Q and ECC-mQ-SNP ensemble forecast fields for the example shown in the paper and plot them
#
#      written by Michael Scheuerer, Feb 2018
#


library(ncdf4)
library(fields)
library(maps)

source("~/Desktop/QPF-T2M-MV/AuxiliaryFunctions.r")


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


filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/data/refcstv2_precip_ccpav3_066_to_072.nc'
prec.nc <- nc_open(filename)
lons.fcst <- ncvar_get(prec.nc, varid="lons_fcst")
lats.fcst <- ncvar_get(prec.nc, varid="lats_fcst")
nc_close(prec.nc)

coords.grid <- as.matrix(expand.grid(lons.fcst[5:8,1],lats.fcst[1,30:32]))[-1,]
nloc <- nrow(coords.grid)

mask.conus <- matrix(is.na(map.where("state",lon-0.03,lat)), nxa, nya)
mask.domain <- mask.conus | (lon < -123.9845) | (lon > -122.1095) | (lat < 38.38457) | (lat > 39.78888)


#verif.ind.anal <- 13095:13098     # original example
#verif.ind.anal <- 15942:15945     # really heavy event
verif.ind.anal <- 11758:11761

#lead.seq <- 9:12
lead.seq <- 8:11

nmb <- 11

month <- (yyyymmddhh_begin[verif.ind.anal[1]]%/%10000)%%100
year <- yyyymmddhh_begin[verif.ind.anal[1]]%/%1000000

ensmean.train <- ensmeandiff.train <- enspop.train <- obs.train <- array(dim=c(nxa,nya,1002,4))
ensmean.verif <- ensmeandiff.verif <- enspop.verif <- array(dim=c(nxa,nya,4))
gefs.verif <- array(dim=c(nloc,4,11))
gefs.ip.verif <- array(dim=c(nxa,nya,4,11))

par.climo <- array(dim=c(nxa,nya,4,3))


for (ilt in 1:4)  {
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

	mid.ind <- which( dates%/%1000000 != year & (dates%/%10000)%%100 %in% month & (dates%/%100)%%100 == 15 )
	if(month==1) mid.ind <- c(mid.ind,length(dates)+16)
	if(month==12) mid.ind <- c(-16,mid.ind)
	train.ind <- as.vector(outer(seq(-45,45,1), mid.ind, '+'))
	train.ind <- train.ind[train.ind >= 1  &  train.ind <= length(dates)]
	train.ind.anal <- match(dates.fcstb[train.ind],yyyymmddhh_begin)
	n.train <- length(train.ind)
	verif.ind <- match(yyyymmddhh_begin[verif.ind.anal[ilt]],dates.fcstb)
	
	obs.train[,,1:n.train,ilt] <- apcp.anal[,,train.ind.anal]
	fcst.train <- apcp.fcst[,,,train.ind]
	dim(fcst.train) <- c(nxf*nyf,nmb,n.train)
	fcst.verif <- apcp.fcst[,,,verif.ind]
	dim(fcst.verif) <- c(nxf*nyf,nmb)

	cl.avg.anal <- matrix(NA,nxa,nya)
	cl.avg.fcst <- rep(NA,nxf*nyf)


        # --- estimate analysis and forecast climatological average

	for (ix in 1:nxa)  {
		for (jy in 1:nya)  {
			if (mask.domain[ix,jy]) next
			if (sum(!is.na(obs.train[ix,jy,,ilt]))<200) next
			cl.avg.anal[ix,jy] <- mean(obs.train[ix,jy,,ilt],na.rm=TRUE)
		}
	}
	for (igrid in 1:(nxf*nyf))  {
		cl.avg.fcst[igrid] <- mean(fcst.train[igrid,,],na.rm=TRUE)
	}


	# --- Interpolate forecasts to CCPA grid

	loc.ip <- cbind(lon[!mask.domain],lat[!mask.domain])
	for (imb in 1:nmb)  {
		gefs.verif[,ilt,imb] <- apcp.fcst[5:8,30:32,imb,verif.ind][-1]
		obj <- list(x=lons.fcst[,1], y=lats.fcst[1,], z=apcp.fcst[,,imb,verif.ind])
		gefs.ip.verif[,,ilt,imb][!mask.domain] <- interp.surface(obj, loc.ip)
	}


	# --- Fit CSGD model to local climatology, adjusted forecasts, and calculate ensemble statistics

	for (ix in 1:nxa)  {
		print(c(ix,nxa))
		for (jy in 1:nya)  {
			if (mask.domain[ix,jy]) next
			obs.mean <- mean(obs.train[ix,jy,obs.train[ix,jy,,ilt]>0,ilt],na.rm=TRUE)
			obs.pop <- mean(obs.train[ix,jy,,ilt]>0,na.rm=TRUE)
			sigma <- obs.mean

			for (mu in (40:1)*(sigma/40))  {
				shape <- (mu/sigma)^2
				scale <- mu/shape
				shift <- -qgamma(obs.pop,shape=shape,scale=scale,lower.tail=FALSE)
				if (shift > -mu/2) break
			}
			par0 <- c(mu,sigma,shift)
			par.climo[ix,jy,ilt,] <- optim(par0, crps.climo, obs=obs.train[ix,jy,,ilt], method="L-BFGS-B", lower=par0*c(0.5,0.5,2), upper=par0*c(2,2,0.1))$par

		      # Normalize (multiplicatively) the forecasts at all forecast grid points within a neighborhood of this location
			dst2 <- outer( (lons.fcst[,1]-lon[ix,jy])^2, (lats.fcst[1,]-lat[ix,jy])^2, '+')
			nbh.ind <- which(as.vector(dst2)<=rho2)
			nnbh <- length(nbh.ind)

			fcst.bc.train <- sweep(fcst.train[nbh.ind,,],1,cl.avg.fcst[nbh.ind],'/')
			fcst.bc.verif <- sweep(fcst.verif[nbh.ind,],1,cl.avg.fcst[nbh.ind],'/')

		      # Calculate ensemble statistics with prediction error based weights
			ensmean.bc.train <- apply(fcst.bc.train,c(1,3),mean)
			obs.bc.train <- obs.train[ix,jy,1:n.train,ilt] / cl.avg.anal[ix,jy]

			rmse <- apply(ensmean.bc.train, 1, function(x) sqrt(mean((obs.bc.train-x)^2,na.rm=TRUE)) )
			weights <- 1-((rmse-min(rmse))/diff(range(rmse)))
			weights <- weights/sum(weights)
			w <- array(weights/nmb, dim=c(nnbh,nmb))

			ensmean.train[ix,jy,1:n.train,ilt] <- apply(sweep(fcst.bc.train, c(1,2), w, "*"), 3, sum, na.rm=TRUE)
			ensmeandiff.train[ix,jy,1:n.train,ilt] <- apply(fcst.bc.train, 3, wgt.md, w=w, na.rm=TRUE)
			enspop.train[ix,jy,1:n.train,ilt] <- apply(sweep(1*(fcst.bc.train>0), c(1,2), w, "*"), 3, sum, na.rm=TRUE)
			ensmean.verif[ix,jy,ilt] <- sum(fcst.bc.verif*w, na.rm=TRUE)
			ensmeandiff.verif[ix,jy,ilt] <- wgt.md(fcst.bc.verif, w=w, na.rm=TRUE)
			enspop.verif[ix,jy,ilt] <- sum(1*(fcst.bc.verif>0)*w, na.rm=TRUE)
		}
	}
}

save(mask.domain, par.climo, ensmean.train, ensmeandiff.train, enspop.train, ensmean.verif, ensmeandiff.verif, enspop.verif, obs.train, file="~/Desktop/Russian-River-CaseStudy/HighRes-Statistics.Rdata")





## Fit CSGD regression model at CCPA grid resolution


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


par.reg <- array(dim=c(nxa,nya,4,6))
mu.fcst <- sigma.fcst <- shift.fcst <- array(dim=c(nxa,nya,4))

for (ilt in 1:4)  {
	cat(paste('Lead time',ilt,'\n'))
	for (ix in 1:nxa)  {
		for (jy in 1:nya)  {
			if (mask.domain[ix,jy]) next

			par0 <- c(0.1,0.1,1.0,1.0,0.5,0.5)

			opt.res <- optim(par0, crps.reg,
				obs = obs.train[ix,jy,,ilt],
				enspop = enspop.train[ix,jy,,ilt],
				ensmean = ensmean.train[ix,jy,,ilt],
				ensmeandiff = ensmeandiff.train[ix,jy,,ilt],
				par.climo = par.climo[ix,jy,ilt,],
				method = "L-BFGS-B",
				lower = c(0.001, 0.05, 0.0, 0.0, 0.1, 0.0),
				upper = c(1.0,    1.0, 1.5, 1.5, 1.0, 1.5))

			par.reg[ix,jy,ilt,] <- opt.res$par

			mu.cl    <- par.climo[ix,jy,ilt,1]
			sigma.cl <- par.climo[ix,jy,ilt,2]
			shift.cl <- par.climo[ix,jy,ilt,3]

			log.arg <- opt.res$par[2] + opt.res$par[3]*enspop.verif[ix,jy,ilt] + opt.res$par[4]*ensmean.verif[ix,jy,ilt]
			mu.fcst[ix,jy,ilt] <- mu.cl*log1p(expm1(opt.res$par[1])*log.arg)/opt.res$par[1]
			sigma.fcst[ix,jy,ilt] <- opt.res$par[5]*sigma.cl*sqrt(mu.fcst[ix,jy,ilt]/mu.cl) + opt.res$par[6]*sigma.cl*ensmeandiff.verif[ix,jy,ilt]
			shift.fcst[ix,jy,ilt] <- shift.cl
		}
	}
}

shape.fcst <- (mu.fcst/sigma.fcst)^2
scale.fcst <- (sigma.fcst^2)/mu.fcst
csgd.fcst.mean <- shift.fcst*pgamma(-shift.fcst,shape.fcst,scale=scale.fcst,lower.tail=F) + mu.fcst*pgamma(-shift.fcst,shape.fcst+1,scale=scale.fcst,lower.tail=F)




##  Generate 6-h ECC-Q forecast trajectories

eccq.ranks <- integer(nxa*nya*4*nmb)
dim(eccq.ranks) <- c(nxa,nya,4,nmb)

apcp.eccq.traj <- array(dim=c(nxa,nya,4,nmb))
apcp.fcst.qt <- array(dim=c(nxa,nya,4,nmb))

for (ilt in 1:4)  {
	cat(paste('Lead time',ilt,'\n'))		# Generate ECC-Q ensemble
	for (ix in 1:nxa)  {
		for (jy in 1:nya)  {
			if (mask.domain[ix,jy]) next
			shape <- (mu.fcst[ix,jy,ilt]/sigma.fcst[ix,jy,ilt])^2
			scale <- mu.fcst[ix,jy,ilt]/shape
			shift <- shift.fcst[ix,jy,ilt]
			fcst.qt <- pmax(shift+qgamma(((1:nmb)-0.5)/nmb,scale=scale,shape=shape),0)
			apcp.fcst.qt[ix,jy,ilt,] <- fcst.qt
			eccq.ranks[ix,jy,ilt,] <- rank(gefs.ip.verif[ix,jy,ilt,],ties.method="random")
			apcp.eccq.traj[ix,jy,ilt,] <- fcst.qt[eccq.ranks[ix,jy,ilt,]]
		}
	}
}




# Plot interpolated GEFS fields and ECC-Q ensemble

library(fields)
library(maps)

load("~/Desktop/Russian-River-CaseStudy/FcstFields-ECC.Rdata")

nmb <- 11
lead.seq <- 8:11
zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
bwr <- colorRampPalette(colors=c("blue","white","red"))
breaks <- seq(0,zmax,,length(colors)+1)
colors.rank <- c('yellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','violet','violetred2','red3','red4')

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
bwr <- colorRampPalette(colors=c("blue","white","red"))
breaks <- seq(0,zmax,,length(colors)+1)
colors.rank <- c('yellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','violet','violetred2','red3','red4')

for (k in 1:nmb)  {
	pdf(paste("~/Desktop/Russian-River-CaseStudy/ECC/CaseStudy-FcstField",k,"-ECC-Q.pdf",sep=""), width=18, height=14)
	split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
	split.screen(c(4,4), screen=1) -> ind
	for (ilt in 1:4)  {
		screen(ind[ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		gefs.member <- matrix(NA,4,3)
		gefs.member[-1] <- gefs.verif[,ilt,k]
		ylab <- ifelse(ilt==1,paste('GEFS member',k),'')
		image(lons.fcst[5:8,1], lats.fcst[1,30:32], gefs.member, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		points(lons.fcst[5,1], lats.fcst[1,30], pch='X', cex=4)
		title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[4+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Interpolated GEFS member',k),'')
		poly.image(lon, lat, gefs.ip.verif[,,ilt,k], breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[8+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Rank of int. GEFS member',k),'')
		poly.image(lon, lat, eccq.ranks[,,ilt,k], breaks=seq(0.5,11.5,1), col=colors.rank, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[12+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('ECC-Q member',k),'')
		poly.image(lon, lat, pmin(apcp.eccq.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		map("state", add=TRUE)
		box()
	}
	split.screen(c(2,1), screen=2) -> ind
	screen(ind[1])
	image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.1,.8), col=colors, axis.args=axis.args)
	screen(ind[2])
	image.plot(zlim=c(1,11), legend.only=TRUE, smallplot=c(.27,.42,.2,.9), col=colors.rank, axis.args=list(at=1:11,labels=1:11), legend.lab='rank', legend.cex=1.2)
	close.screen(all=TRUE)
	dev.off()
}

order(apply(apcp.eccq.traj,4,mean,na.rm=TRUE))








##  Generate 6-h ECC-mQ-SNP forecast trajectories

loc.ip <- cbind(lon[!mask.domain],lat[!mask.domain])
n.domain <- sum(!mask.domain)

lambda <- 0.5

delta.lon <- 0.46875
delta.lat <- 0.46810

r.lon <- 3*delta.lon
r.lat <- 3*delta.lat


gefs.verif.snp <- array(dim=c(11,4,nmb))
gefs.ip.verif.snp <- apcp.eccmq.snp.traj <- array(dim=c(nxa,nya,4,nmb))

for (ilt in 1:4)  {						# Generate ECC-mQ-SNP ensemble
	cleadb <- formatC(lead.seq[ilt]*6,width=3,flag="0")
	cleade <- formatC((lead.seq[ilt]+1)*6,width=3,flag="0")
	filename <- paste('~/Desktop/Russian-River-CaseStudy/data/refcstv2_precip_ccpav3_',cleadb,'_to_',cleade,'.nc',sep='')
	cat(paste('Loading ',filename,'\n'))
	prec.nc <- nc_open(filename)
	dates.fcstb <- ncvar_get(prec.nc, varid="yyyymmddhh_fcstb")
	lons.fcst <- ncvar_get(prec.nc, varid="lons_fcst")
	lats.fcst <- ncvar_get(prec.nc, varid="lats_fcst")
	apcp.fcst <- ncvar_get(prec.nc, varid="apcp_fcst_ens")
	nc_close(prec.nc)

	shape.cal <- (mu.fcst[,,ilt][!mask.domain]/sigma.fcst[,,ilt][!mask.domain])^2
	scale.cal <- mu.fcst[,,ilt][!mask.domain]/shape.cal
	shift.cal <- shift.fcst[,,ilt][!mask.domain]
	qt.cal <- sweep(sweep(outer(shape.cal,((1:nmb)-0.5)/nmb,FUN=function(X,Y) qgamma(Y,shape=X)),1,scale.cal,'*'),1,shift.cal,'+')

	verif.ind <- match(yyyymmddhh_begin[verif.ind.anal[ilt]],dates.fcstb)
	lons.fcst.ss <- lons.fcst[4:9,1]
	lats.fcst.ss <- lats.fcst[1,29:33]
	nx.ss <- length(lons.fcst.ss)
	ny.ss <- length(lats.fcst.ss)
	grid.fcst.ss <- as.matrix(expand.grid(lons.fcst.ss,lats.fcst.ss))


       # Define basis functions for simulating 'negative precipitation'
	dst.lon <- as.matrix(dist(grid.fcst.ss[,1]))
	tricube.lon <- ifelse(dst.lon>r.lon,0,(1-(dst.lon/r.lon)^3)^3)
	dst.lat <- as.matrix(dist(grid.fcst.ss[,2]))
	tricube.lat <- ifelse(dst.lat>r.lat,0,(1-(dst.lat/r.lat)^3)^3)
	basis <- tricube.lon*tricube.lat
	basis <- sweep(basis,1,apply(basis,1,sum),'/')
	gefs.ss <- apcp.fcst[4:9,29:33,,verif.ind]


       # Simulate 'negative precipitation' where the raw ensemble forecasts are zero
	for (imb in 1:11)  {
		precip.sim.coef <- rep(0,nx.ss*ny.ss)
		zero.ind <- (as.vector(gefs.ss[,,imb])==0)
		n.zero <- sum(zero.ind)
		if (n.zero==0) {
			gefs.verif.snp[,ilt,imb] <- gefs.ss[2:5,2:4,imb][-1]
			next
		}
		precip.sim.coef[zero.ind] <- runif(n.zero,-1,0)
		gefs.ss[,,imb][zero.ind] <- apply(sweep(basis[zero.ind,,drop=FALSE],2,precip.sim.coef,'*'),1,sum)
		gefs.verif.snp[,ilt,imb] <- gefs.ss[2:5,2:4,imb][-1]
	}


       # Interpolate augmented raw ensemble fields to analysis grid
	gefs.ip.augm <- matrix(NA,n.domain,nmb)
	for (imb in 1:nmb)  {
		obj <- list(x=lons.fcst.ss, y=lats.fcst.ss, z=gefs.ss[,,imb])
		gefs.ip.augm[,imb] <- interp.surface(obj, loc.ip)
		gefs.ip.verif.snp[,,ilt,imb][!mask.domain] <- gefs.ip.augm[,imb]
	}

       # Construct piecewise linear quantile-mapping function
	eccmq.snp.field <- matrix(0,n.domain,nmb)
	for (igp in 1:n.domain)  {
		l <- sum(qt.cal[igp,]<=0)
		if (l==nmb)  next					# all forecast quantiles are zero -> map to zero
		if (l==(nmb-1))  {					# only one non-zero forecast quantile -> map largest GEFS member to this quantile
			ind.max <- which.max(gefs.ip.augm[igp,])
			eccmq.snp.field[igp,ind.max] <- qt.cal[igp,imb]
			next
		}							# otherwise: fit (regularized) linear regression spline
		ind.pos <- order(gefs.ip.augm[igp,])[l:nmb]
		n.pos <- length(ind.pos)
		x <- gefs.ip.augm[igp,ind.pos]
		y <- qt.cal[igp,l:nmb]
		A <- cbind(1,x,pmax(outer(x,x[c(-1,-n.pos)],'-'),0))
		coef <- solve(crossprod(A,A)+lambda*diag(c(0,0,pmax(1,y[c(-1,-n.pos)]))),crossprod(A,y))
		eccmq.snp.field[igp,ind.pos] <- pmax(A%*%coef,0)
	}
	for (imb in 1:nmb)  {
		apcp.eccmq.snp.traj[,,ilt,imb][!mask.domain] <- eccmq.snp.field[,imb]
	}
}

save(lons.fcst, lats.fcst, mask.domain, csgd.fcst.mean, gefs.verif, gefs.ip.verif, gefs.verif.snp, gefs.ip.verif.snp, par.reg, mu.fcst, sigma.fcst, shift.fcst, apcp.fcst.qt, eccq.ranks, apcp.eccq.traj, apcp.eccmq.snp.traj, file="~/Desktop/Russian-River-CaseStudy/FcstFields-ECC.Rdata")







# Plot interpolated GEFS fields and ECC-mQ-SNP ensemble

library(fields)
library(maps)

load("~/Desktop/Russian-River-CaseStudy/FcstFields-ECC.Rdata")

nmb <- 11
lead.seq <- 8:11
zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- c("#BEBEBE","#E9E9E9",colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51))
bwr <- colorRampPalette(colors=c("blue","white","red"))
breaks <- c(-2.784314,-1.392157,seq(0,zmax,,length(colors)-1))
colors.rank <- c('yellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','violet','violetred2','red3','red4')

for (k in 1:nmb)  {
	pdf(paste("~/Desktop/Russian-River-CaseStudy/ECC/CaseStudy-FcstField",k,"-ECC-mQ-SNP.pdf",sep=""), width=18, height=14)
	split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
	split.screen(c(4,4), screen=1) -> ind
	for (ilt in 1:4)  {
		screen(ind[ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		gefs.member <- matrix(NA,4,3)
		gefs.member[-1] <- gefs.verif.snp[,ilt,k]
		ylab <- ifelse(ilt==1,paste('GEFS member',k,'(& SNP)'),'')
		image(lons.fcst[5:8,1], lats.fcst[1,30:32], gefs.member, breaks=breaks, zlim=c(-2.784314,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		points(lons.fcst[5,1], lats.fcst[1,30], pch='X', cex=4)
		title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[4+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Interpolated GEFS member',k),'')
		poly.image(lon, lat, gefs.ip.verif.snp[,,ilt,k], breaks=breaks, zlim=c(-2.784314,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[8+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Rank of int. GEFS member',k),'')
		poly.image(lon, lat, eccq.ranks[,,ilt,k], breaks=seq(0.5,11.5,1), col=colors.rank, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[12+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('ECC-mQ-SNP member',k),'')
		poly.image(lon, lat, pmin(apcp.eccmq.snp.traj[,,ilt,k],zmax), breaks=breaks[-c(1,2)], zlim=c(0,zmax), col=colors[-c(1,2)], xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		map("state", add=TRUE)
		box()
	}
	split.screen(c(2,1), screen=2) -> ind
	screen(ind[1])
	image.plot(zlim=c(-2.784314,zmax), legend.only=TRUE, smallplot=c(.27,.42,.1,.8), col=colors, axis.args=axis.args)
	screen(ind[2])
	image.plot(zlim=c(1,11), legend.only=TRUE, smallplot=c(.27,.42,.2,.9), col=colors.rank, axis.args=list(at=1:11,labels=1:11), legend.lab='rank', legend.cex=1.2)
	close.screen(all=TRUE)
	dev.off()
}

order(apply(apcp.eccmq.snp.traj,4,mean,na.rm=TRUE))





















###
#   Plot GEFS ensemble mean, CSGD quantiles, and verifying CCPA analysis

library(fields)
library(maps)

load("~/Desktop/Russian-River-CaseStudy/FcstFields-ECC.Rdata")

#verif.ind.anal <- 13095:13098     # original example
#verif.ind.anal <- 15942:15945     # really heavy event
verif.ind.anal <- 11758:11761 

nmb <- 11
lead.seq <- 8:11
zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
breaks <- seq(0,zmax,,length(colors)+1)

date.string <- c("Jan 19, 2010, 0Z - 6Z", "Jan 19, 2010, 6Z - 12Z", "Jan 19, 2010, 12Z - 18Z", "Jan 19, 2010, 18Z - 0Z")

gefsmean.verif <- apply(gefs.verif, c(1,2), mean)


pdf("~/Desktop/Russian-River-CaseStudy/CaseStudy-GEFS-CCPA.pdf", width=18, height=10.5)
split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
split.screen(c(3,4), screen=1) -> ind
for (ilt in 1:4)  {
	screen(ind[ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	gefsmean <- matrix(NA,4,3)
	gefsmean[-1] <- gefsmean.verif[,ilt]
	ylab <- ifelse(ilt==1,'GEFS ensemble mean','')
	image(lons.fcst[5:8,1], lats.fcst[1,30:32], gefsmean, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	points(lons.fcst[5,1], lats.fcst[1,30], pch='X', cex=4)
	title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[4+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'CSGD mean','')
	poly.image(lon, lat, csgd.fcst.mean[,,ilt], breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[8+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'Analyzed field','')
	obs <- apcp.anal[,,verif.ind.anal[ilt]]
	obs[mask.domain] <- NA
	poly.image(lon, lat, obs, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	title(date.string[ilt], cex.main=1.5)
	map("state", add=TRUE)
	box()
}
screen(2)
image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.15,.85), col=colors, axis.args=axis.args)
close.screen(all=TRUE)
dev.off()









# Separate plots for ECC-Q and ECC-mQ-SNP (presentation-friendly format)

library(fields)
library(maps)

load("~/Desktop/Russian-River-CaseStudy/FcstFields-ECC.Rdata")

nmb <- 11
lead.seq <- 8:11
zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
bwr <- colorRampPalette(colors=c("blue","white","red"))
breaks <- seq(0,zmax,,length(colors)+1)
colors.rank <- c('yellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','violet','violetred2','red3','red4')

k <- 2


pdf("~/Desktop/Russian-River-CaseStudy/ECC/CaseStudy-FcstField2-ECC-Q.pdf", width=18, height=10.5)
split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
split.screen(c(3,4), screen=1) -> ind
for (ilt in 1:4)  {
	screen(ind[ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,paste('Interpolated GEFS member',k),'')
	poly.image(lon, lat, gefs.ip.verif[,,ilt,k], breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[4+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,paste('Rank of int. GEFS member',k),'')
	poly.image(lon, lat, eccq.ranks[,,ilt,k], breaks=seq(0.5,11.5,1), col=colors.rank, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[8+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,paste('ECC-Q member',k),'')
	poly.image(lon, lat, pmin(apcp.eccq.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
split.screen(c(2,1), screen=2) -> ind
screen(ind[1])
image.plot(zlim=c(1,11), legend.only=TRUE, smallplot=c(.27,.42,.1,.8), col=colors.rank, axis.args=list(at=1:11,labels=1:11), legend.lab='rank', legend.cex=1.2)
screen(ind[2])
image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.2,.9), col=colors, axis.args=axis.args)
close.screen(all=TRUE)
dev.off()



pdf("~/Desktop/Russian-River-CaseStudy/ECC/CaseStudy-FcstField2-ECC-mQ-SNP.pdf", width=18, height=10.5)
split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
split.screen(c(3,4), screen=1) -> ind
for (ilt in 1:4)  {
	screen(ind[ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	gefs.member <- matrix(NA,4,3)
	gefs.member[-1] <- gefs.verif[,ilt,k]
	ylab <- ifelse(ilt==1,paste('GEFS member',k),'')
	image(lons.fcst[5:8,1], lats.fcst[1,30:32], gefs.member, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	points(lons.fcst[5,1], lats.fcst[1,30], pch='X', cex=4)
	title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[4+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,paste('ECC-Q member',k),'')
	poly.image(lon, lat, pmin(apcp.eccq.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[8+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,paste('ECC-T member',k),'')
	poly.image(lon, lat, pmin(apcp.ecct.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
screen(2)
image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.15,.85), col=colors, axis.args=axis.args)
close.screen(all=TRUE)
dev.off()


