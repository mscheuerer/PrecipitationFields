
###
#   Code for generating the MDSS-SDA ensemble forecast fields for the example shown in the paper and plot them
#
#      written by Michael Scheuerer, Feb 2018
#


library(ncdf4)
library(fields)
library(maps)

source("~/Desktop/QPF-T2M-MV/AuxiliaryFunctions.r")


# Load analyzed data

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


# Load coordinates of forecast grid

filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/data/refcstv2_precip_ccpav3_066_to_072.nc'
prec.nc <- nc_open(filename)
lons.fcst <- ncvar_get(prec.nc, varid="lons_fcst")
lats.fcst <- ncvar_get(prec.nc, varid="lats_fcst")
nc_close(prec.nc)

coords.grid <- as.matrix(expand.grid(lons.fcst[5:8,1],lats.fcst[1,30:32]))[-1,]
nloc <- nrow(coords.grid)


# Mask out all analysis grid points over the ocean or outside the study domain

mask.conus <- matrix(is.na(map.where("state",lon-0.03,lat)), nxa, nya)
mask.domain <- mask.conus | (lon < -123.9845) | (lon > -122.1095) | (lat < 38.38457) | (lat > 39.78888)


# Associate CCPA grid points with 1/2 degree grid and scale up

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



#verif.ind.anal <- 13095:13098     # original example
#verif.ind.anal <- 15942:15945     # really heavy event
verif.ind.anal <- 11758:11761 

#lead.seq <- 9:12
lead.seq <- 8:11

month <- (yyyymmddhh_begin[verif.ind.anal[1]]%/%10000)%%100
year <- yyyymmddhh_begin[verif.ind.anal[1]]%/%1000000



## Calculate ensemble statistics for the CSGD method


ensmean.train <- ensmeandiff.train <- enspop.train <- obs.train <- array(dim=c(nloc,1002,4))
ensmean.verif <- ensmeandiff.verif <- enspop.verif <- gefsmean.verif <- array(dim=c(nloc,4))

par.climo <- array(dim=c(nloc,4,3))


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
	
	obs.train[,1:n.train,ilt] <- apcp.anal.upsc[,train.ind.anal]
	fcst.train <- apcp.fcst[,,,train.ind]
	dim(fcst.train) <- c(nxf*nyf,nmb,n.train)
	fcst.verif <- apcp.fcst[,,,verif.ind]
	dim(fcst.verif) <- c(nxf*nyf,nmb)

	cl.avg.anal <- rep(NA,nloc)
	cl.avg.fcst <- rep(NA,nxf*nyf)


        # --- estimate analysis and forecast climatological average

	for (iloc in 1:nloc)  {
		if (sum(!is.na(obs.train[iloc,,ilt]))<200) next
		cl.avg.anal[iloc] <- mean(obs.train[iloc,,ilt],na.rm=TRUE)
	}
	for (igrid in 1:(nxf*nyf))  {
		cl.avg.fcst[igrid] <- mean(fcst.train[igrid,,],na.rm=TRUE)
	}


	# --- Loop through all locations, adjusted forecasts, calculate ensemble statistics, and fit CSGD model to local climatology

	for (iloc in 1:nloc)  {
   	      # Fit observation climatology
		obs.mean <- mean(obs.train[iloc,obs.train[iloc,,ilt]>0,ilt],na.rm=TRUE)
		obs.pop <- mean(obs.train[iloc,,ilt]>0,na.rm=TRUE)
		sigma <- obs.mean

		for (mu in (40:1)*(sigma/40))  {
			shape <- (mu/sigma)^2
			scale <- mu/shape
			shift <- -qgamma(obs.pop,shape=shape,scale=scale,lower.tail=FALSE)
			if (shift > -mu/2) break
		}
		par0 <- c(mu,sigma,shift)
		par.climo[iloc,ilt,] <- optim(par0, crps.climo, obs=obs.train[iloc,,ilt], method="L-BFGS-B", lower=par0*c(0.5,0.5,2), upper=par0*c(2,2,0.1))$par

	      # Normalize (multiplicatively) the forecasts at all forecast grid points within a neighborhood of this location
		dst2 <- outer( (lons.fcst[,1]-coords.grid[iloc,1])^2, (lats.fcst[1,]-coords.grid[iloc,2])^2, '+')
		nbh.ind <- which(as.vector(dst2)<=rho2)
		nnbh <- length(nbh.ind)

		fcst.bc.train <- sweep(fcst.train[nbh.ind,,],1,cl.avg.fcst[nbh.ind],'/')
		fcst.bc.verif <- sweep(fcst.verif[nbh.ind,],1,cl.avg.fcst[nbh.ind],'/')

	      # Calculate ensemble statistics with prediction error based weights
		ensmean.bc.train <- apply(fcst.bc.train,c(1,3),mean)
		obs.bc.train <- obs.train[iloc,1:n.train,ilt] / cl.avg.anal[iloc]
		rmse <- apply(ensmean.bc.train, 1, function(x) sqrt(mean((obs.bc.train-x)^2,na.rm=TRUE)) )
		weights <- 1-((rmse-min(rmse))/diff(range(rmse)))
		weights <- weights/sum(weights)
		w <- array(weights/nmb, dim=c(nnbh,nmb))

#		plot(lons.fcst[nbh.ind], lats.fcst[nbh.ind], pch=19, cex=100*weights)

		ensmean.train[iloc,1:n.train,ilt] <- apply(sweep(fcst.bc.train, c(1,2), w, "*"), 3, sum, na.rm=TRUE)
		ensmeandiff.train[iloc,1:n.train,ilt] <- apply(fcst.bc.train, 3, wgt.md, w=w, na.rm=TRUE)
		enspop.train[iloc,1:n.train,ilt] <- apply(sweep(1*(fcst.bc.train>0), c(1,2), w, "*"), 3, sum, na.rm=TRUE)
		ensmean.verif[iloc,ilt] <- sum(fcst.bc.verif*w, na.rm=TRUE)
		ensmeandiff.verif[iloc,ilt] <- wgt.md(fcst.bc.verif, w=w, na.rm=TRUE)
		enspop.verif[iloc,ilt] <- sum(1*(fcst.bc.verif>0)*w, na.rm=TRUE)
#		readline("Press return to continue")
	}
	gefsmean.verif[,ilt] <- apply(apcp.fcst[5:8,30:32,,verif.ind],c(1,2),mean)[-1]
}





## Fit CSGD regression model at 1/2 degree resolution


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


par.reg <- array(dim=c(nloc,4,6))
mu.fcst <- sigma.fcst <- shift.fcst <- array(dim=c(nloc,4))

for (ilt in 1:4)  {
	for (iloc in 1:nloc)  {

		par0 <- c(0.1,0.1,1.0,1.0,0.5,0.5)

		opt.res <- optim(par0, crps.reg,
			obs = obs.train[iloc,,ilt],
			enspop = enspop.train[iloc,,ilt],
			ensmean = ensmean.train[iloc,,ilt],
			ensmeandiff = ensmeandiff.train[iloc,,ilt],
			par.climo = par.climo[iloc,ilt,],
			method = "L-BFGS-B",
			lower = c(0.001, 0.05, 0.0, 0.0, 0.1, 0.0),
			upper = c(1.0,    1.0, 1.5, 1.5, 1.0, 1.5))

		par.reg[iloc,ilt,] <- opt.res$par
#		print(round(opt.res$par,3))

		mu.cl    <- par.climo[iloc,ilt,1]
		sigma.cl <- par.climo[iloc,ilt,2]
		shift.cl <- par.climo[iloc,ilt,3]

		log.arg <- opt.res$par[2] + opt.res$par[3]*enspop.verif[iloc,ilt] + opt.res$par[4]*ensmean.verif[iloc,ilt]
		mu.fcst[iloc,ilt] <- mu.cl*log1p(expm1(opt.res$par[1])*log.arg)/opt.res$par[1]
		sigma.fcst[iloc,ilt] <- opt.res$par[5]*sigma.cl*sqrt(mu.fcst[iloc,ilt]/mu.cl) + opt.res$par[6]*sigma.cl*ensmeandiff.verif[iloc,ilt]
		shift.fcst[iloc,ilt] <- rep(shift.cl,length(ensmean.verif[iloc,ilt]))
	}
}






##  Create space-time dependence template


apcp.cl <- array(dim=c(nxa,nya,4))
apcp.obs.traj <- array(dim=c(nxa,nya,4,n.train))
apcp.upsc.traj <- array(dim=c(nloc,4,n.train))

for (ilt in 1:4)  {
	train.ind.anal <- match(dates[train.ind],yyyymmddhh_begin)+lead.seq[ilt]
	apcp.obs.traj[,,ilt,] <- apcp.anal[,,train.ind.anal]
	apcp.cl[,,ilt] <- apply(apcp.obs.traj[,,ilt,],c(1,2),mean,na.rm=TRUE)
	apcp.upsc.traj[,ilt,] <- apcp.anal.upsc[,train.ind.anal]
}


d <- nloc*4
nmb <- 11

#k.seq <- c(rev(20+5*cumsum(0:18)),15,11)
k.seq <- c(rev(10+cumsum(rep(1:30,each=2))))

e.fcst <- array(dim=c(d,nmb))

shape <- (as.vector(mu.fcst)/as.vector(sigma.fcst))^2
scale <- as.vector(mu.fcst)/shape
shift <- as.vector(shift.fcst)

for (k in 1:nmb)  {
	e.fcst[,k] <- pmax(shift+qgamma((k-0.5)/nmb,scale=scale,shape=shape),0)
}

dim(apcp.upsc.traj) <- c(d,n.train)


##  1. Select MDSS observation trajectories (using forecasts at 1/2 degree resolution)

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

# matplot((1:d)-0.12,apcp.upsc.traj[,mdss.ind], type='p', pch=18, col=4, ylim=c(0,50))
# matplot((1:d)+0.12, e.fcst, type='p', pch=18, col=2, add=TRUE)




##   2. Generate 6-h MDSS forecast trajectories (post-processing at 1/2 degree resolution, spatial disaggregation via MDSS trajectories)

e.obs <- apcp.upsc.traj[,mdss.ind]
dim(e.obs) <- dim(e.fcst) <- c(nloc,4,nmb)
apcp.hist.upsc <- e.obs

apcp.hist.traj <- array(dim=c(nxa,nya,4,nmb))
apcp.mdss.traj <- array(dim=c(nxa,nya,4,nmb))
apcp.adjm.fctr <- array(dim=c(nxa,nya,4,nmb))
apcp.mdss.upsc <- array(dim=c(nloc,4,nmb))


delta.lon <- 0.46875
delta.lat <- 0.46810

r.lon <- 1.5*delta.lon
r.lat <- 1.5*delta.lat

dst.lon <- abs(outer(as.vector(lon),coords.grid[,1],'-'))
tricube.lon <- ifelse(dst.lon>r.lon,0,(1-(dst.lon/r.lon)^3)^3)
dst.lat <- abs(outer(as.vector(lat),coords.grid[,2],'-'))
tricube.lat <- ifelse(dst.lat>r.lat,0,(1-(dst.lat/r.lat)^3)^3)

basis.fcts <- tricube.lon*tricube.lat
basis.fcts[as.vector(mask.domain),] <- NA
dim(basis.fcts) <- c(nxa,nya,nloc)
basis.fcts <- sweep(basis.fcts,c(1,2),apply(basis.fcts,c(1,2),sum),'/')

#par(mfrow=c(1,3))
#image.plot(basis.fcts[,,1])
#image.plot(basis.fcts[,,7])
#image.plot(basis.fcts[,,11])


dst.gefs <- as.matrix(dist(coords.grid))
n.nbh <- apply(dst.gefs>0 & dst.gefs<0.5, 1, sum)
L <- 1*(dst.gefs==0) - sweep(1*(dst.gefs>0 & dst.gefs<0.5),1,n.nbh,'/')


adjm.fct.target <- function(alpha,M,xi)  {
	sum((M%*%alpha-xi)^2) + 0.25*sum((L%*%alpha)^2)
}


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
			e.obs[iloc,ilt,i] <- sum(apcp.hist.upsc[,ilt,i]*wgt)
		}
		if (any(zeros))  {
		      # ... and try ranking those upscaled values
			apcp.mdss.ranks[zeros] <- rank(e.obs[iloc,ilt,zeros],ties.method="random")
		}
		apcp.mdss.upsc[iloc,ilt,] <- e.fcst[iloc,ilt,apcp.mdss.ranks]
	}


       # Step 2: Disaggregate these forecasts to fine scale using the historic templates and a spatially smooth adjustment function

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

		alpha <- optim(rep(1,nloc), adjm.fct.target, M=M, xi=xi, method="L-BFGS-B", lower=rep(0,nloc), upper=rep(10,nloc))$par
		print(round(alpha,3))

		apcp.adjm.fctr[,,ilt,imb] <- apply(sweep(basis.fcts,3,alpha,'*'),c(1,2),sum)
		apcp.mdss.traj[,,ilt,imb] <- apcp.adjm.fctr[,,ilt,imb] * apcp.hist.traj[,,ilt,imb]
	}
}

mdss.dates <- yyyymmddhh_begin[match(dates[train.ind[mdss.ind]],yyyymmddhh_begin)]

save(lons.fcst, lats.fcst, gefsmean.verif, mask.domain, apcp.mdss.upsc, apcp.hist.upsc, apcp.adjm.fctr, apcp.hist.traj, apcp.mdss.traj, file="~/Desktop/Russian-River-CaseStudy/FcstFields-MDSS-SDA.Rdata")




# Plot historic fields, adjustment factor, and modified MDSS field

load("~/Desktop/Russian-River-CaseStudy/FcstFields-MDSS-SDA.Rdata")

nmb <- 11
lead.seq <- 8:11
zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
bwr <- colorRampPalette(colors=c("blue","white","red"))
axis.args.adjf <- list(at=c(0,0.1,1,4,10)^0.3,labels=c(0,0.1,1,4,10),cex.axis=1.3)
breaks <- seq(0,zmax,,length(colors)+1)

for (k in 1:nmb)  {
	pdf(paste("~/Desktop/Russian-River-CaseStudy/MDSS/CaseStudy-FcstField",k,"-MDSS-SDA.pdf",sep=""), width=18, height=17.5)
	split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
	split.screen(c(5,4), screen=1) -> ind
	for (ilt in 1:4)  {
		screen(ind[ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Historic field',k,'(upscaled)'),'')
		fcst.field <- matrix(NA,4,3)
		fcst.field[-1] <- apcp.hist.upsc[,ilt,k]
		image(lons.fcst[5:8,1], lats.fcst[1,30:32], fcst.field, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		points(lons.fcst[5,1], lats.fcst[1,30], pch='X', cex=4)
		title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[4+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('MDSS member',k,'(coarse)'),'')
		fcst.field <- matrix(NA,4,3)
		fcst.field[-1] <- apcp.mdss.upsc[,ilt,k]
		image(lons.fcst[5:8,1], lats.fcst[1,30:32], fcst.field, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		points(lons.fcst[5,1], lats.fcst[1,30], pch='X', cex=4)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[8+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,'Adjustment function','')
		poly.image(lon, lat, apcp.adjm.fctr[,,ilt,k]^0.3, zlim=c(0,2), col=bwr(50), xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[12+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Historic field',k),'')
		poly.image(lon, lat, pmin(apcp.hist.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[16+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('MDSS-SDA member',k),'')
		poly.image(lon, lat, pmin(apcp.mdss.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		map("state", add=TRUE)
		box()
	}
	split.screen(c(3,1), screen=2) -> ind
	screen(ind[1])
	image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.05,.85), col=colors, axis.args=axis.args)
	screen(ind[2])
	image.plot(zlim=c(0,2), legend.only=TRUE, smallplot=c(.27,.42,.1,.9), col=bwr(50), axis.args=axis.args.adjf)
	screen(ind[3])
	image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.15,.95), col=colors, axis.args=axis.args)
	close.screen(all=TRUE)
	dev.off()
}

order(apply(apcp.mdss.traj,4,mean,na.rm=TRUE))












# Plot historic fields, adjustment factor, and modified MDSS field  (two sepearate plots)

load("~/Desktop/Russian-River-CaseStudy/FcstFields-MDSS-SDA.Rdata")

nmb <- 11
lead.seq <- 8:11
zmax <- 71

k < 11

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
bwr <- colorRampPalette(colors=c("blue","white","red"))
axis.args.adjf <- list(at=c(0,0.1,1,4,10)^0.3,labels=c(0,0.1,1,4,10),cex.axis=1.3)
breaks <- seq(0,zmax,,length(colors)+1)


pdf("~/Desktop/Russian-River-CaseStudy/MDSS/CaseStudy-FcstField11-MDSS-coarse.pdf", width=18, height=10.5)
split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
split.screen(c(3,4), screen=1) -> ind
for (ilt in 1:4)  {
	screen(ind[ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,paste('Historic field',k,'(upscaled)'),'')
	fcst.field <- matrix(NA,4,3)
	fcst.field[-1] <- apcp.hist.upsc[,ilt,k]
	image(lons.fcst[5:8,1], lats.fcst[1,30:32], fcst.field, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	points(lons.fcst[5,1], lats.fcst[1,30], pch='X', cex=4)
	title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[4+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'Coarse-scale MDSS member','')
	fcst.field <- matrix(NA,4,3)
	fcst.field[-1] <- apcp.mdss.upsc[,ilt,k]
	image(lons.fcst[5:8,1], lats.fcst[1,30:32], fcst.field, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	points(lons.fcst[5,1], lats.fcst[1,30], pch='X', cex=4)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[8+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'Adjustment function','')
	poly.image(lon, lat, apcp.adjm.fctr[,,ilt,k]^0.3, zlim=c(0,2), col=bwr(50), xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
split.screen(c(2,1), screen=2) -> ind
screen(ind[1])
image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.05,.85), col=colors, axis.args=axis.args)
screen(ind[2])
image.plot(zlim=c(0,2), legend.only=TRUE, smallplot=c(.27,.42,.1,.9), col=bwr(50), axis.args=axis.args.adjf, legend.lab='adjustment factor', legend.cex=1.2, legend.line=3)
close.screen(all=TRUE)
dev.off()




pdf("~/Desktop/Russian-River-CaseStudy/MDSS/CaseStudy-FcstField11-MDSS-fine.pdf", width=18, height=10.5)
split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
split.screen(c(3,4), screen=1) -> ind
for (ilt in 1:4)  {
	screen(ind[ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'Adjustment function','')
	poly.image(lon, lat, apcp.adjm.fctr[,,ilt,k]^0.3, zlim=c(0,2), col=bwr(50), xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[4+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,paste('Historic field',k),'')
	poly.image(lon, lat, pmin(apcp.hist.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[8+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'Fine-scale MDSS member','')
	poly.image(lon, lat, pmin(apcp.mdss.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
split.screen(c(2,1), screen=2) -> ind
screen(ind[1])
image.plot(zlim=c(0,2), legend.only=TRUE, smallplot=c(.27,.42,.1,.9), col=bwr(50), axis.args=axis.args.adjf, legend.lab='adjustment factor', legend.cex=1.2, legend.line=3)
screen(ind[2])
image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.15,.95), col=colors, axis.args=axis.args)
close.screen(all=TRUE)
dev.off()



