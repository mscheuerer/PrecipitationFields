
###
#   Code for generating all of the verification metrics shown in the paper and for ploting them
#
#      written by Michael Scheuerer, Feb 2018
#


library(ncdf4)
library(maps)

source("~/Desktop/QPF-T2M-MV/AuxiliaryFunctions.r")


## Load CCPA fields

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


## Load coordinates of forecast grid points

filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/data/refcstv2_precip_ccpav3_066_to_072.nc'
prec.nc <- nc_open(filename)
lons.fcst <- ncvar_get(prec.nc, varid="lons_fcst")
lats.fcst <- ncvar_get(prec.nc, varid="lats_fcst")
nc_close(prec.nc)

coords.grid <- as.matrix(expand.grid(lons.fcst[5:8,1],lats.fcst[1,30:32]))[-1,]
nloc <- nrow(coords.grid)


## Calculate index that associates each analysis grid point with the nearest forecast grid point 

dst2 <- outer(as.vector(lon),coords.grid[,1],'-')^2 + outer(as.vector(lat),coords.grid[,2],'-')^2
dst2[as.vector(mask.domain),] <- NA

I <- apply(dst2,1,function(x) ifelse(all(is.na(x)),NA,which.min(x)))
#dim(I) <- c(nxa,nya)



years <- 2002:2013
nyears <- length(years)

lead.seq <- 8:11
nlt <- length(lead.seq)

month.seq <- c(10,11,12,1,2,3,4,5)




###
#   Calculate Brier skill scores for marginal forecast distributions


obs.all <- array(dim=c(12,nxa,nya,nyears,31,nlt))

for (month in c(9,month.seq,6))  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata",sep=""))
	obs.all[month,,,,,] <- obs.verif
}



mdss.sda.crps <- stss.crps <- eccmq.snp.crps <- array(dim=c(12,12,31,4))
clim.crps <- array(dim=c(12,4))

for (month in month.seq)  {

	month.window <- c(12,1:12,1)[month+(0:2)]

	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-mQ-SNP-",month,".Rdata",sep=""))
	nmb <- dim(ranks.ecc)[6]

	for (ilt in 1:4)  {
		for (iyear in 1:12)  {
			for (iday in 1:31)  {
				use <- !is.na(mu.fcst[,,iyear,iday,ilt]) & !apply(is.na(apcp.mdss.traj[,,iyear,iday,ilt,]),c(1,2),any) & !apply(is.na(apcp.ecct.traj[,,iyear,iday,ilt,]),c(1,2),any)
				if(sum(use)==0) next
				shape <- (mu.fcst[,,iyear,iday,ilt]/sigma.fcst[,,iyear,iday,ilt])^2
				scale <- mu.fcst[,,iyear,iday,ilt]/shape
				shift <- shift.fcst[,,iyear,iday,ilt]

				mdss.sda.crps.day <- stss.crps.day <- eccmq.snp.crps.day <- matrix(NA,nxa,nya)

				for (ix in 1:nxa)  {
					for (jy in 1:nya)  {
						if (!use[ix,jy]) next
						fcst.qt <- pmax(shift[ix,jy]+qgamma(((1:nmb)-0.5)/nmb,scale=scale[ix,jy],shape=shape[ix,jy]),0)
						T1 <- mean(abs(fcst.qt-obs.verif[ix,jy,iyear,iday,ilt]))
						T2 <- - 0.5*gini.md(fcst.qt)
						stss.crps.day[ix,jy] <- T1 + T2
						T1 <- mean(abs(apcp.mdss.traj[ix,jy,iyear,iday,ilt,]-obs.verif[ix,jy,iyear,iday,ilt]))
						T2 <- - 0.5*gini.md(apcp.mdss.traj[ix,jy,iyear,iday,ilt,])
						mdss.sda.crps.day[ix,jy] <- T1 + T2
						T1 <- mean(abs(apcp.eccqm.snp.traj[ix,jy,iyear,iday,ilt,]-obs.verif[ix,jy,iyear,iday,ilt]))
						T2 <- - 0.5*gini.md(apcp.eccqm.snp.traj[ix,jy,iyear,iday,ilt,])
						eccmq.snp.crps.day[ix,jy] <- T1 + T2
					}
				}

				mdss.sda.crps[month,iyear,iday,ilt] <- mean(mdss.sda.crps.day, na.rm=TRUE)
				stss.crps[month,iyear,iday,ilt] <- mean(stss.crps.day, na.rm=TRUE)
				eccmq.snp.crps[month,iyear,iday,ilt] <- mean(eccmq.snp.crps.day, na.rm=TRUE)
			}
		}

		clim.crps.month <- matrix(NA,nxa,nya)

		for (ix in 1:nxa)  {
			for (jy in 1:nya)  {
				if (mask.domain[ix,jy]) next
				clim.crps.month[ix,jy] <- 0.5*gini.md(as.vector(obs.verif[ix,jy,,,ilt]),na.rm=TRUE)
			}
		}
		clim.crps[month,ilt] <- mean(clim.crps.month, na.rm=TRUE)
	}
}


round(1-apply(stss.crps,1,mean,na.rm=TRUE)/apply(clim.crps,1,mean,na.rm=TRUE),3)[month.seq]
round(1-apply(mdss.sda.crps,1,mean,na.rm=TRUE)/apply(clim.crps,1,mean,na.rm=TRUE),3)[month.seq]
round(1-apply(eccmq.snp.crps,1,mean,na.rm=TRUE)/apply(clim.crps,1,mean,na.rm=TRUE),3)[month.seq]


# StSS		0.342 0.249 0.295 0.259 0.299 0.276 0.268 0.300
# MDSS-SDA	0.336 0.238 0.283 0.247 0.286 0.263 0.257 0.287
# ECC-mQ-SNP	0.339 0.238 0.287 0.248 0.292 0.267 0.260 0.291




matplot((1-cbind(apply(stss.crps,1,mean,na.rm=TRUE),apply(mdss.sda.crps,1,mean,na.rm=TRUE),apply(eccmq.snp.crps,1,mean,na.rm=TRUE))/apply(clim.crps,1,mean,na.rm=TRUE))[month.seq,], lty=1, type='l', axes=FALSE, ylab='', ylim=c(0,0.35))
legend(4,0.15, col=1:3, legend=c('StSS / ECC-Q / MDSS-QRO','MDSS-SDA','ECC-mQ-SNP'), lwd=1, cex=1.5)
axis(1, at=1:8, label=month.string[month.seq], cex.axis=1.3)
axis(2, cex.axis=1.3)
box()





###
#   Calculate the subgrid scale observed and forecast maximum precipitation amounts


obs.verif.max <- array(dim=c(nloc,12,12,31))
mdss.ro.fcst.max <- mdss.sda.fcst.max <- stss.fcst.max <- array(dim=c(nloc,12,12,31,11))
eccq.fcst.max <- eccmq.snp.fcst.max <- array(dim=c(nloc,12,12,31,11))

for (month in c(9,month.seq,6))  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata",sep=""))
	# nmb <- dim(apcp.mdss.traj)[6]
	nmb <- 11
	dim(obs.verif) <- c(nxa*nya,12,31,4)
	dim(apcp.mdss.traj) <- c(nxa*nya,12,31,4,nmb)
	for (iloc in 1:nloc)  {
		id.gpt <- which(!is.na(I) & (I==iloc))
		obs.verif.max[iloc,month,,] <- apply(obs.verif[id.gpt,,,],c(2,3),max)
		mdss.sda.fcst.max[iloc,month,,,] <- apply(apcp.mdss.traj[id.gpt,,,,],c(2,3,5),max)
	}
}

for (month in month.seq)  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-mQ-SNP-",month,".Rdata",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/HistoricRanks-MDSS-SDA-",month,".Rdata",sep=""))
	nmb <- dim(ranks.ecc)[6]
	dim(ranks.ecc) <- c(nxa*nya,12,31,4,nmb)
	dim(ranks.stss) <- c(nxa*nya,12,31,4,nmb)
	dim(ranks.mdss.ro) <- c(nxa*nya,12,31,4,nmb)
	dim(apcp.eccqm.snp.traj) <- c(nxa*nya,12,31,4,nmb)
	dim(mu.fcst) <- c(nxa*nya,12,31,4)
	dim(sigma.fcst) <- c(nxa*nya,12,31,4)
	dim(shift.fcst) <- c(nxa*nya,12,31,4)	
	for (iloc in 1:nloc)  {
		id.gpt <- which(!is.na(I) & (I==iloc))
		ngpt <- length(id.gpt)
		for (iyear in 1:12)  {
			for (iday in 1:31)  {
				if (any(is.na(mu.fcst[id.gpt,iyear,iday,]))) next
				apcp.eccq.traj <- apcp.stss.traj <- apcp.mdss.qro.traj <- array(dim=c(ngpt,4,nmb))
				for (igpt in 1:ngpt)  {
					for (ilt in 1:4)  {
						shape <- (mu.fcst[id.gpt[igpt],iyear,iday,ilt]/sigma.fcst[id.gpt[igpt],iyear,iday,ilt])^2
						scale <- mu.fcst[id.gpt[igpt],iyear,iday,ilt]/shape
						shift <- shift.fcst[id.gpt[igpt],iyear,iday,ilt]
						fcst.qt <- pmax(shift+qgamma(((1:nmb)-0.5)/nmb,scale=scale,shape=shape),0)
						apcp.eccq.traj[igpt,ilt,] <- fcst.qt[ranks.ecc[id.gpt[igpt],iyear,iday,ilt,]]
						apcp.stss.traj[igpt,ilt,] <- fcst.qt[ranks.stss[id.gpt[igpt],iyear,iday,ilt,]]
						apcp.mdss.qro.traj[igpt,ilt,] <- fcst.qt[ranks.mdss.ro[id.gpt[igpt],iyear,iday,ilt,]]

					}
				}
				eccq.fcst.max[iloc,month,iyear,iday,] <- apply(apcp.eccq.traj,3,max)
				stss.fcst.max[iloc,month,iyear,iday,] <- apply(apcp.stss.traj,3,max)
				mdss.ro.fcst.max[iloc,month,iyear,iday,] <- apply(apcp.mdss.qro.traj,3,max)
			}
		}
		eccmq.snp.fcst.max[iloc,month,,,] <- apply(apcp.eccqm.snp.traj[id.gpt,,,,],c(2,3,5),max)
	}
}

save(obs.verif.max, mdss.sda.fcst.max, eccq.fcst.max, eccmq.snp.fcst.max, stss.fcst.max, mdss.ro.fcst.max, file="~/Desktop/Russian-River-CaseStudy/data/MaxFcstFields.Rdata")



###
#  Create reliability diagrams (for 11-member ensembles)


method.string <- "eccmq-snp"

if (method.string=="mdss-sda")  {
	fcst.max <- mdss.sda.fcst.max
}
if (method.string=="eccmq-snp")  {
	fcst.max <- eccmq.snp.fcst.max
}
if (method.string=="stss")  {
	fcst.max <- stss.fcst.max
}


month.string <- c("Jan","Feb","Mar","Apr","Mai","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
season.ind <- data.frame(DJ=c(12,1), FM=c(2,3), AM=c(4,5), ON=c(10,11))

breaks <- (0:11)/11


pdf(paste("~/Desktop/Russian-River-CaseStudy/reliability_",method.string,".pdf",sep=""), width=12, height=16)
par(mfrow=c(4,3))

for (irow in 1:4)  {

	n.pop  <- x.pop  <- y.pop  <- NULL
	n.pot.10 <- x.pot.10 <- y.pot.10 <- NULL
	n.pot.25 <- x.pot.25 <- y.pot.25 <- NULL

	# Calculate statistics for reliability diagrams

	for (month in season.ind[[irow]])  {

		for (iyear in 1:12)  {

			for (iday in 1:31)  {

				use <- !is.na(obs.verif.max[,month,iyear,iday]) & apply(!is.na(fcst.max[,month,iyear,iday,]),1,any)
				if (sum(use)==0) next

				exc.0p1 <- 1*(obs.verif.max[,month,iyear,iday]>0.1)
				exc.0p1[!use] <- NA
				pop <- apply(1*(fcst.max[,month,iyear,iday,]>0.1),1,mean)
				I.pop <- outer(pop, breaks[-length(breaks)], ">=") & outer(pop, breaks[-1], "<")
			    	n.pop <- rbind(n.pop, apply(I.pop, 2, sum))
		 		x.pop <- rbind(x.pop, apply(pop*I.pop, 2, sum, na.rm=TRUE))
      				y.pop <- rbind(y.pop, apply(exc.0p1*I.pop, 2, sum, na.rm=TRUE))
				
				exc.10 <- 1*(obs.verif.max[,month,iyear,iday]>10)
				exc.10[!use] <- NA
				pot.10 <- apply(1*(fcst.max[,month,iyear,iday,]>10),1,mean)
				I.pot.10 <- outer(pot.10, breaks[-length(breaks)], ">=") & outer(pot.10, breaks[-1], "<")
			    	n.pot.10 <- rbind(n.pot.10, apply(I.pot.10, 2, sum))
		 		x.pot.10 <- rbind(x.pot.10, apply(pot.10*I.pot.10, 2, sum, na.rm=TRUE))
      				y.pot.10 <- rbind(y.pot.10, apply(exc.10*I.pot.10, 2, sum, na.rm=TRUE))

				exc.25 <- 1*(obs.verif.max[,month,iyear,iday]>25)
				exc.25[!use] <- NA
				pot.25 <- apply(1*(fcst.max[,month,iyear,iday,]>25),1,mean)
				I.pot.25 <- outer(pot.25, breaks[-length(breaks)], ">=") & outer(pot.25, breaks[-1], "<")
			    	n.pot.25 <- rbind(n.pot.25, apply(I.pot.25, 2, sum))
		 		x.pot.25 <- rbind(x.pot.25, apply(pot.25*I.pot.25, 2, sum, na.rm=TRUE))
      				y.pot.25 <- rbind(y.pot.25, apply(exc.25*I.pot.25, 2, sum, na.rm=TRUE))
			}
		}
	}

	par(mfg=c(irow,1))
	plot.reliability (n.pop, x.pop, y.pop, log.freq=TRUE, N.boot=5000, threshold="0.1mm", cleadb=48, cleade=72)
	par(mfg=c(irow,2))
	plot.reliability (n.pot.10, x.pot.10, y.pot.10, log.freq=TRUE, N.boot=5000, threshold="10mm", cleadb=48, cleade=72)
	par(mfg=c(irow,3))
	plot.reliability (n.pot.25, x.pot.25, y.pot.25, log.freq=TRUE, N.boot=5000, threshold="25mm", cleadb=48, cleade=72)
}
dev.off()




###
#  Calculate and plot Brier skill scores for maximum precipitation

thresh <- c(0.1,10,25)

bs <- array(dim=c(12,6,3))

for (month in month.seq)  {

	month.window <- c(12,1:12,1)[month+(0:2)]
	use <- apply(!is.na(mdss.sda.fcst.max[,month,,,]),c(1,2,3),any) & apply(!is.na(eccq.fcst.max[,month,,,]),c(1,2,3),any) & apply(!is.na(stss.fcst.max[,month,,,]),c(1,2,3),any)
	obs.verif.max[,month,,][!use] <- NA

	for (ith in 1:3)  {
		exc <- 1*(obs.verif.max[,month,,]>thresh[ith])
		pot.clim <- apply(obs.verif.max[,month.window,,]>thresh[ith],1,mean,na.rm=TRUE)
		bs[month,6,ith] <- mean(sweep(exc,1,pot.clim,'-')^2,na.rm=TRUE)

		pot.mdss.sda <- apply(1*(mdss.sda.fcst.max[,month,,,]>thresh[ith]),c(1,2,3),mean)
		pot.mdss.ro <- apply(1*(mdss.ro.fcst.max[,month,,,]>thresh[ith]),c(1,2,3),mean)
		pot.eccq <- apply(1*(eccq.fcst.max[,month,,,]>thresh[ith]),c(1,2,3),mean)
		pot.eccmq.snp <- apply(1*(eccmq.snp.fcst.max[,month,,,]>thresh[ith]),c(1,2,3),mean)
		pot.stss <- apply(1*(stss.fcst.max[,month,,,]>thresh[ith]),c(1,2,3),mean)

		bs[month,1,ith] <- mean((pot.eccq-exc)^2,na.rm=TRUE)
		bs[month,2,ith] <- mean((pot.eccmq.snp-exc)^2,na.rm=TRUE)
		bs[month,3,ith] <- mean((pot.stss-exc)^2,na.rm=TRUE)
		bs[month,4,ith] <- mean((pot.mdss.sda-exc)^2,na.rm=TRUE)
		bs[month,5,ith] <- mean((pot.mdss.ro-exc)^2,na.rm=TRUE)

	}
}

bss <- 1-sweep(bs[,1:5,],c(1,3),bs[,6,],'/')


month.string <- c("Jan","Feb","Mar","Apr","Mai","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
colors <- c('darkgreen','lightseagreen','navy','red','violet')

pdf("~/Desktop/Russian-River-CaseStudy/BrierSkillScores.pdf", width=10, height=4)
par(mfrow=c(1,3), mar=c(3,3,3,1))
matplot(bss[month.seq,1:5,1], type='b', lty=1, pch=c(2,5,1,4,6), main='Brier skill scores, threshold: 0.1 mm', xlab='', ylab='Brier skill score', axes=FALSE, ylim=c(-0.65,0.65), col=colors, cex.main=1.4)
abline(h=0, lty=2)
axis(1, at=1:8, label=month.string[month.seq], cex.axis=1.3)
axis(2, cex.axis=1.3)
box()
matplot(bss[month.seq,1:5,2], type='b', lty=1, pch=c(2,5,1,4,6), main='Brier skill scores, threshold: 10 mm', xlab='', ylab='Brier skill score', axes=FALSE, ylim=c(-0.65,0.65), col=colors, cex.main=1.4)
abline(h=0, lty=2)
axis(1, at=1:8, label=month.string[month.seq], cex.axis=1.3)
axis(2, cex.axis=1.3)
box()
legend(2.3, -0.14, lwd=1, c("StSS","MDSS-RO","MDSS-SDA","ECC-Q","ECC-mQ-SNP"), col=colors[c(3,5,4,1,2)], pch=c(1,6,4,2,5), lty=1, cex=1.3)
matplot(bss[month.seq,1:5,3], type='b', lty=1, pch=c(2,5,1,4,6), main='Brier skill scores, threshold: 25 mm', xlab='', ylab='Brier skill score', axes=FALSE, ylim=c(-0.65,0.65), col=colors, cex.main=1.4)
abline(h=0, lty=2)
axis(1, at=1:8, label=month.string[month.seq], cex.axis=1.3)
axis(2, cex.axis=1.3)
box()
dev.off()






###
#   Calculate fraction exceedance ranks (for entire domain)


obs.verif.all <- array(dim=c(nxa,nya,12,12,31,4))

for (month in c(9,month.seq,6))  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata",sep=""))
	obs.verif.all[,,month,,,] <- obs.verif
}

thresh <- array(dim=c(nxa,nya,4,5))

for (month in month.seq)  {
	month.window <- c(12,1:12,1)[month+(0:2)]
	thresh[,,,1] <- 0.1
	thresh[,,,2] <- 10
	thresh[,,,3] <- 25
	thresh[,,,4] <- apply(obs.verif.all[,,month.window,,,],c(1,2,6),quantile,prob=0.995,na.rm=TRUE)
	thresh[,,,5] <- apply(obs.verif.all[,,month.window,,,],c(1,2,6),quantile,prob=0.999,na.rm=TRUE)
}

dim(thresh) <-c(nxa*nya,4,5)



obs.verif.frac <- array(dim=c(12,12,31,4,5))
mdss.sda.fcst.frac <- mdss.ro.fcst.frac <- stss.fcst.frac <- eccq.fcst.frac <- eccmq.snp.fcst.frac <- array(dim=c(12,12,31,4,11,5))

for (month in month.seq)  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-SDA-",month,".Rdata",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-mQ-SNP-",month,".Rdata",sep=""))
	nxa <- dim(apcp.mdss.traj)[1]
	nya <- dim(apcp.mdss.traj)[2]
	nmb <- dim(apcp.mdss.traj)[6]
	dim(obs.verif) <- c(nxa*nya,12,31,4)
	dim(apcp.mdss.traj) <- c(nxa*nya,12,31,4,nmb)
	dim(apcp.eccqm.snp.traj) <- c(nxa*nya,12,31,4,nmb)
	domain.ind <- which(!mask.domain)
	for (ith in 1:5)  {
		obs.verif.frac[month,,,,ith] <- apply(sweep(obs.verif[domain.ind,,,],c(1,4),thresh[domain.ind,,ith],'>'),c(2,3,4),mean)
		mdss.sda.fcst.frac[month,,,,,ith] <- apply(sweep(apcp.mdss.traj[domain.ind,,,,],c(1,4),thresh[domain.ind,,ith],'>'),c(2,3,4,5),mean)
		eccmq.snp.fcst.frac[month,,,,,ith] <- apply(sweep(apcp.eccqm.snp.traj[domain.ind,,,,],c(1,4),thresh[domain.ind,,ith],'>'),c(2,3,4,5),mean)
	}
}

for (month in month.seq)  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/HistoricRanks-MDSS-",month,".Rdata",sep=""))
	nxa <- dim(ranks.ecc)[1]
	nya <- dim(ranks.ecc)[2]
	nmb <- dim(ranks.ecc)[6]
	domain.ind <- which(!mask.domain)
	dim(ranks.ecc) <- c(nxa*nya,12,31,4,nmb)
	dim(ranks.mdss.ro) <- c(nxa*nya,12,31,4,nmb)
	dim(ranks.stss) <- c(nxa*nya,12,31,4,nmb)
	dim(mu.fcst) <- c(nxa*nya,12,31,4)
	dim(sigma.fcst) <- c(nxa*nya,12,31,4)
	dim(shift.fcst) <- c(nxa*nya,12,31,4)
	for (iyear in 1:12)  {
		for (iday in 1:31)  {
			for (ilt in 1:4)  {
				if (all(is.na(mu.fcst[domain.ind,iyear,iday,ilt]))) next
				shape <- (mu.fcst[domain.ind,iyear,iday,ilt]/sigma.fcst[domain.ind,iyear,iday,ilt])^2
				scale <- mu.fcst[domain.ind,iyear,iday,ilt]/shape
				shift <- shift.fcst[domain.ind,iyear,iday,ilt]
				for (imb in 1:nmb)  {
					apcp.eccq.traj <- pmax(shift+qgamma((ranks.ecc[domain.ind,iyear,iday,ilt,imb]-0.5)/nmb,scale=scale,shape=shape),0)
					eccq.fcst.frac[month,iyear,iday,ilt,imb,] <- apply(sweep(thresh[domain.ind,ilt,],1,apcp.eccq.traj,'<'),2,mean)
					apcp.stss.traj <- pmax(shift+qgamma((ranks.stss[domain.ind,iyear,iday,ilt,imb]-0.5)/nmb,scale=scale,shape=shape),0)
					stss.fcst.frac[month,iyear,iday,ilt,imb,] <- apply(sweep(thresh[domain.ind,ilt,],1,apcp.stss.traj,'<'),2,mean)
					apcp.mdss.ro.traj <- pmax(shift+qgamma((ranks.mdss.ro[domain.ind,iyear,iday,ilt,imb]-0.5)/nmb,scale=scale,shape=shape),0)
					mdss.ro.fcst.frac[month,iyear,iday,ilt,imb,] <- apply(sweep(thresh[domain.ind,ilt,],1,apcp.mdss.ro.traj,'<'),2,mean)
				}
			}
		}
	}
}

save(thresh, obs.verif.frac, mdss.ro.fcst.frac, mdss.sda.fcst.frac, stss.fcst.frac, eccq.fcst.frac, eccmq.snp.fcst.frac, file="~/Desktop/Russian-River-CaseStudy/data/FTE-domain.Rdata")



mdss.sda.exc.rank <- mdss.ro.exc.rank <- stss.exc.rank <- eccq.exc.rank <- eccmq.snp.exc.rank <- array(dim=c(12,12,31,4,5))
mdss.sda.exc.crps <- mdss.ro.exc.crps <- stss.exc.crps <- eccq.exc.crps <- eccmq.snp.exc.crps <- array(dim=c(12,12,31,4,5))
clim.exc.crps <- array(dim=c(12,4,5))

for (month in month.seq)  {
	for (ilt in 1:4)  {
		for (ith in 1:5)  {
			for (iyear in 1:12)  {
				for (iday in 1:31)  {
					obs.mdss.sda.frac <- c(obs.verif.frac[month,iyear,iday,ilt,ith],mdss.sda.fcst.frac[month,iyear,iday,ilt,,ith])
					mdss.sda.exc.rank[month,iyear,iday,ilt,ith] <- ifelse(all(obs.mdss.sda.frac==0),NA,rank(obs.mdss.sda.frac,ties='random')[1])
					T1 <- mean(abs(mdss.sda.fcst.frac[month,iyear,iday,ilt,,ith]-obs.verif.frac[month,iyear,iday,ilt,ith]))
					T2 <- - 0.5*gini.md(mdss.sda.fcst.frac[month,iyear,iday,ilt,,ith])
					mdss.sda.exc.crps[month,iyear,iday,ilt,ith] <- T1 + T2

					obs.mdss.ro.frac <- c(obs.verif.frac[month,iyear,iday,ilt,ith],mdss.ro.fcst.frac[month,iyear,iday,ilt,,ith])
					mdss.ro.exc.rank[month,iyear,iday,ilt,ith] <- ifelse(all(obs.mdss.ro.frac==0),NA,rank(obs.mdss.ro.frac,ties='random')[1])
					T1 <- mean(abs(mdss.ro.fcst.frac[month,iyear,iday,ilt,,ith]-obs.verif.frac[month,iyear,iday,ilt,ith]))
					T2 <- - 0.5*gini.md(mdss.ro.fcst.frac[month,iyear,iday,ilt,,ith])
					mdss.ro.exc.crps[month,iyear,iday,ilt,ith] <- T1 + T2

					obs.eccq.frac <- c(obs.verif.frac[month,iyear,iday,ilt,ith],eccq.fcst.frac[month,iyear,iday,ilt,,ith])
					eccq.exc.rank[month,iyear,iday,ilt,ith] <- ifelse(all(obs.eccq.frac==0),NA,rank(obs.eccq.frac,ties='random')[1])
					T1 <- mean(abs(eccq.fcst.frac[month,iyear,iday,ilt,,ith]-obs.verif.frac[month,iyear,iday,ilt,ith]))
					T2 <- - 0.5*gini.md(eccq.fcst.frac[month,iyear,iday,ilt,,ith])
					eccq.exc.crps[month,iyear,iday,ilt,ith] <- T1 + T2

					obs.eccmq.snp.frac <- c(obs.verif.frac[month,iyear,iday,ilt,ith],eccmq.snp.fcst.frac[month,iyear,iday,ilt,,ith])
					eccmq.snp.exc.rank[month,iyear,iday,ilt,ith] <- ifelse(all(obs.eccmq.snp.frac==0),NA,rank(obs.eccmq.snp.frac,ties='random')[1])
					T1 <- mean(abs(eccmq.snp.fcst.frac[month,iyear,iday,ilt,,ith]-obs.verif.frac[month,iyear,iday,ilt,ith]))
					T2 <- - 0.5*gini.md(eccmq.snp.fcst.frac[month,iyear,iday,ilt,,ith])
					eccmq.snp.exc.crps[month,iyear,iday,ilt,ith] <- T1 + T2

					obs.stss.frac <- c(obs.verif.frac[month,iyear,iday,ilt,ith],stss.fcst.frac[month,iyear,iday,ilt,,ith])
					stss.exc.rank[month,iyear,iday,ilt,ith] <- ifelse(all(obs.stss.frac==0),NA,rank(obs.stss.frac,ties='random')[1])
					T1 <- mean(abs(stss.fcst.frac[month,iyear,iday,ilt,,ith]-obs.verif.frac[month,iyear,iday,ilt,ith]))
					T2 <- - 0.5*gini.md(stss.fcst.frac[month,iyear,iday,ilt,,ith])
					stss.exc.crps[month,iyear,iday,ilt,ith] <- T1 + T2
				}
			}
			clim.exc.crps[month,ilt,ith] <- 0.5*gini.md(as.vector(obs.verif.frac[month,,,ilt,ith]),na.rm=TRUE)
		}
	}
}


thresh.name = c('0.1 mm','10 mm','25 mm','99.5% clim.','99.9% clim.')

#pdf("~/Desktop/Russian-River-CaseStudy/FTE-Histograms.pdf", width=12, height=5.6)
pdf("~/Desktop/Russian-River-CaseStudy/FTE-Histograms.pdf", width=15, height=5.6)
#par(mfrow=c(5,4), mar=c(0.5,5,3,0))
par(mfrow=c(5,5), mar=c(0.5,5,3,0))
for(ith in 1:5)  {
	hist(as.vector(stss.exc.rank[,,,,ith]), col='lightblue', breaks=seq(0.5,12.5,1), main=ifelse(ith==1,'StSS',''), axes=FALSE, ylab=thresh.name[ith], cex.lab=1.2, cex.main=1.5)
	abline(h=sum(!is.na(stss.exc.rank[,,,,ith]))/12, lty=2)
	hist(as.vector(mdss.ro.exc.rank[,,,,ith]), col='lightblue', breaks=seq(0.5,12.5,1), main=ifelse(ith==1,'MDSS-RO',''), axes=FALSE, ylab='', cex.main=1.5)
	abline(h=sum(!is.na(mdss.ro.exc.rank[,,,,ith]))/12, lty=2)
	hist(as.vector(mdss.sda.exc.rank[,,,,ith]), col='lightblue', breaks=seq(0.5,12.5,1), main=ifelse(ith==1,'MDSS-SDA',''), axes=FALSE, ylab='', cex.main=1.5)
	abline(h=sum(!is.na(mdss.sda.exc.rank[,,,,ith]))/12, lty=2)
	hist(as.vector(eccq.exc.rank[,,,,ith]), col='lightblue', breaks=seq(0.5,12.5,1), main=ifelse(ith==1,'ECC-Q',''), axes=FALSE, ylab='', cex.main=1.5)
	abline(h=sum(!is.na(eccq.exc.rank[,,,,ith]))/12, lty=2)
	hist(as.vector(eccmq.snp.exc.rank[,,,,ith]), col='lightblue', breaks=seq(0.5,12.5,1), main=ifelse(ith==1,'ECC-mQ-SNP',''), axes=FALSE, ylab='', cex.main=1.5)
	abline(h=sum(!is.na(eccmq.snp.exc.rank[,,,,ith]))/12, lty=2)
}
dev.off()


crpss <- matrix(NA,5,5)
colnames(crpss) <- c("StSS","MDSS-RO","MDSS-SDA","ECC-Q","ECC-mQ-SNP")
rownames(crpss) <- thresh.name

crpss[,1] <- 1 - apply(stss.exc.crps,5,mean,na.rm=TRUE) / apply(clim.exc.crps,3,mean,na.rm=TRUE)
crpss[,2] <- 1 - apply(mdss.ro.exc.crps,5,mean,na.rm=TRUE) / apply(clim.exc.crps,3,mean,na.rm=TRUE)
crpss[,3] <- 1 - apply(mdss.sda.exc.crps,5,mean,na.rm=TRUE) / apply(clim.exc.crps,3,mean,na.rm=TRUE)
crpss[,4] <- 1 - apply(eccq.exc.crps,5,mean,na.rm=TRUE) / apply(clim.exc.crps,3,mean,na.rm=TRUE)
crpss[,5] <- 1 - apply(eccmq.snp.exc.crps,5,mean,na.rm=TRUE) / apply(clim.exc.crps,3,mean,na.rm=TRUE)


round(crpss,3)


	     StSS   MDSS-RO  MDSS-SDA   ECC-Q    ECC-mQ-SNP
0.1 mm       0.338   0.439    0.478     0.462      0.473
10 mm        0.149   0.225    0.225     0.218      0.211
25 mm        0.065   0.126    0.126     0.116      0.122
99.5% clim.  0.120   0.206    0.207     0.196      0.191
99.9% clim.  0.088   0.161    0.160     0.145      0.142








###
#  Calculate frequencies of days with random reordering

nyears <- length(years)
nlt <- 4

zeros.qt <- integer(nxa*nya*nyears*31*nlt)
dim(zeros.qt) <- c(nxa,nya,nyears,31,nlt)

pop.avg <- array(dim=c(12,nyears,31,nlt))
frac.gt1.stss <- frac.gt1.ecc <- frac.gt1.mdss <- array(dim=c(12,nyears,31,nlt))

for (month in month.seq)  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/Zeros-ECC-StSS-",month,".Rdata",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/HistoricRanks-MDSS-",month,".Rdata",sep=""))

	for (iyear in 1:12)  {
		for (iday in 1:31)  {
			for (ilt in 1:4)  {
				pop.fcst <- array(dim=c(nxa,nya))
				for (ix in 1:nxa)  {
					for (jy in 1:nya)  {
						if (mask.domain[ix,jy]) {
							zeros.qt[ix,jy,iyear,iday,ilt] <- NA
							next
						}
						shape <- (mu.fcst[ix,jy,iyear,iday,ilt]/sigma.fcst[ix,jy,iyear,iday,ilt])^2
						scale <- mu.fcst[ix,jy,iyear,iday,ilt]/shape
						shift <- shift.fcst[ix,jy,iyear,iday,ilt]
						fcst.qt <- pmax(shift+qgamma(((1:nmb)-0.5)/nmb,scale=scale,shape=shape),0)
						zeros.qt[ix,jy,iyear,iday,ilt] <- sum(fcst.qt==0)
						pop.fcst[ix,jy] <- pgamma(-shift,scale=scale,shape=shape,lower.tail=FALSE)
					}
				}
				pop.avg[month,iyear,iday,nlt] <- mean(pop.fcst,na.rm=TRUE)
				frac.gt1.stss[month,iyear,iday,nlt] <- mean(zeros.stss[,,iyear,iday,ilt]-zeros.qt[,,iyear,iday,ilt]>1,na.rm=TRUE)
				frac.gt1.ecc[month,iyear,iday,nlt] <- mean(zeros.ecc[,,iyear,iday,ilt]-zeros.qt[,,iyear,iday,ilt]>1,na.rm=TRUE)
				frac.gt1.mdss[month,iyear,iday,nlt] <- mean(zeros.mdss.ro[,,iyear,iday,ilt]-zeros.qt[,,iyear,iday,ilt]>1,na.rm=TRUE)
			}
		}
	}
}




par(mfrow=c(1,5), mar=c(1,1,1,3))
image.plot(lon, lat, pop.fcst, xlim=c(-123.9845,-122.1095), ylim=c(38.38457,39.78888), zlim=c(0,1), axes=FALSE)
image.plot(lon, lat, zeros.qt[,,iyear,iday,ilt], xlim=c(-123.9845,-122.1095), ylim=c(38.38457,39.78888), zlim=c(0,11), axes=FALSE)
image.plot(lon, lat, zeros.stss[,,iyear,iday,ilt]-zeros.qt[,,iyear,iday,ilt], xlim=c(-123.9845,-122.1095), ylim=c(38.38457,39.78888), zlim=c(-11,11), axes=FALSE, col=cm.colors(51))
image.plot(lon, lat, zeros.ecc[,,iyear,iday,ilt]-zeros.qt[,,iyear,iday,ilt], xlim=c(-123.9845,-122.1095), ylim=c(38.38457,39.78888), zlim=c(-11,11), axes=FALSE, col=cm.colors(51))
image.plot(lon, lat, zeros.mdss.ro[,,iyear,iday,ilt]-zeros.qt[,,iyear,iday,ilt], xlim=c(-123.9845,-122.1095), ylim=c(38.38457,39.78888), zlim=c(-11,11), axes=FALSE, col=cm.colors(51))


month <- 11

par(mfrow=c(1,3), mar=c(3,3,1,1))
plot(as.vector(pop.avg[month,,,]), as.vector(frac.gt1.stss[month,,,]), pch=19, main='StSS', col='navy')
plot(as.vector(pop.avg[month,,,]), as.vector(frac.gt1.ecc[month,,,]), pch=19, main='ECC', col='navy')
plot(as.vector(pop.avg[month,,,]), as.vector(frac.gt1.mdss[month,,,]), pch=19, main='MDSS-QRO', col='navy')








###
#  Calculate band depth ranks (not useful here because it doesn't permit stratification according to precipitation level)


mdss.bd.rank <- stss.bd.rank <- ecc.bd.rank <- array(dim=c(12,12,31,4))

for (month in month.seq)  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-",month,".Rdata",sep=""))
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata",sep=""))
	nxa <- dim(apcp.mdss.traj)[1]
	nya <- dim(apcp.mdss.traj)[2]
	nmb <- dim(apcp.mdss.traj)[6]
	domain.ind <- which(!mask.domain)
	dim(obs.verif) <- c(nxa*nya,12,31,4)
	dim(apcp.mdss.traj) <- c(nxa*nya,12,31,4,nmb)
	dim(ranks.ecc) <- c(nxa*nya,12,31,4,nmb)
	dim(ranks.stss) <- c(nxa*nya,12,31,4,nmb)
	dim(mu.fcst) <- c(nxa*nya,12,31,4)
	dim(sigma.fcst) <- c(nxa*nya,12,31,4)
	dim(shift.fcst) <- c(nxa*nya,12,31,4)
	for (iyear in 1:12)  {
		print(iyear)
		for (iday in 1:31)  {
			for (ilt in 1:4)  {
				if (any(is.na(mu.fcst[domain.ind,iyear,iday,ilt])) | any(is.na(obs.verif[domain.ind,iyear,iday,ilt]))) next
				shape <- (mu.fcst[domain.ind,iyear,iday,ilt]/sigma.fcst[domain.ind,iyear,iday,ilt])^2
				scale <- mu.fcst[domain.ind,iyear,iday,ilt]/shape
				shift <- shift.fcst[domain.ind,iyear,iday,ilt]
				apcp.ecc.traj <- apcp.stss.traj <- matrix(NA,length(domain.ind),nmb)
				for (imb in 1:nmb)  {
					apcp.ecc.traj[,imb] <- pmax(shift+qgamma((ranks.ecc[domain.ind,iyear,iday,ilt,imb]-0.5)/nmb,scale=scale,shape=shape),0)
					apcp.stss.traj[,imb] <- pmax(shift+qgamma((ranks.stss[domain.ind,iyear,iday,ilt,imb]-0.5)/nmb,scale=scale,shape=shape),0)
				}
				if (any(obs.verif[domain.ind,iyear,iday,ilt]>0) | any(apcp.ecc.traj>0))  {
					ecc.bd.rank[month,iyear,iday,ilt] <- bd.rank.ties(obs.verif[domain.ind,iyear,iday,ilt],apcp.ecc.traj)
				}
				if (any(obs.verif[domain.ind,iyear,iday,ilt]>0) | any(apcp.stss.traj>0))  {
					stss.bd.rank[month,iyear,iday,ilt] <- bd.rank.ties(obs.verif[domain.ind,iyear,iday,ilt],apcp.stss.traj)
				}
				if (any(obs.verif[domain.ind,iyear,iday,ilt]>0) | any(apcp.mdss.traj[domain.ind,iyear,iday,ilt,]>0))  {
					mdss.bd.rank[month,iyear,iday,ilt] <- bd.rank.ties(obs.verif[domain.ind,iyear,iday,ilt],apcp.mdss.traj[domain.ind,iyear,iday,ilt,])
				}
			}
		}
	}
}

save(mdss.bd.rank, stss.bd.rank, ecc.bd.rank, file="~/Desktop/Russian-River-CaseStudy/data/BandDepthRank.Rdata")



pdf("~/Desktop/Russian-River-CaseStudy/BD-Histograms.pdf", width=12, height=1.5)
par(mfrow=c(1,3), mar=c(0,5,3,0))
hist(as.vector(ecc.bd.rank), col='lightblue', breaks=seq(0.5,12.5,1), main='Ensemble Copula Coupling', axes=FALSE, ylab='')
abline(h=sum(!is.na(ecc.bd.rank))/12, lty=2)
hist(as.vector(stss.bd.rank), col='lightblue', breaks=seq(0.5,12.5,1), main='Schaake Shuffle', axes=FALSE, ylab='')
abline(h=sum(!is.na(stss.bd.rank))/12, lty=2)
hist(as.vector(mdss.bd.rank), col='lightblue', breaks=seq(0.5,12.5,1), main='Minimum Divergence Schaake Shuffle', axes=FALSE, ylab='')
abline(h=sum(!is.na(mdss.bd.rank))/12, lty=2)
dev.off()







# Same for 55-member ensembles
# note: results were not very interesting, not included in the paper

mdss55.fcst.mean <- mdss55.fcst.max <- stss55.fcst.mean <- stss55.fcst.max <- array(dim=c(nloc,12,12,31,55))
ecc55.fcst.mean <- ecc55.fcst.max <- array(dim=c(nloc,12,12,31,55))

for (month in month.seq)  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-55-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-55-",month,".Rdata",sep=""))
	nmb <- dim(apcp.mdss.traj)[6]
	dim(apcp.mdss.traj) <- c(nxa*nya,12,31,4,nmb)
	for (iloc in 1:nloc)  {
		id.gpt <- which(!is.na(I) & (I==iloc))
		mdss55.fcst.mean[iloc,month,,,] <- apply(apcp.mdss.traj[id.gpt,,,,],c(2,3,5),mean)
		mdss55.fcst.max[iloc,month,,,] <- apply(apcp.mdss.traj[id.gpt,,,,],c(2,3,5),max)
	}
}

for (month in month.seq)  {
	cat(paste("Loading ~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata\n",sep=""))
	load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata",sep=""))
	dim(ranks.ecc) <- c(nxa*nya,12,31,4,11)
	dim(mu.fcst) <- c(nxa*nya,12,31,4)
	dim(sigma.fcst) <- c(nxa*nya,12,31,4)
	dim(shift.fcst) <- c(nxa*nya,12,31,4)
	ranks.stss <- integer(nxa*nya*12*31*4*55)
	dim(ranks.stss) <- c(nxa*nya,12,31,4,55)

	for (ilt in 1:4)  {							# Calculate 55-member Schaake shuffle ranks
		print(c(ilt,4))
		cleadb <- formatC(lead.seq[ilt]*6,width=3,flag="0")
		cleade <- formatC((lead.seq[ilt]+1)*6,width=3,flag="0")
		filename <- paste('~/Desktop/Russian-River-CaseStudy/data/refcstv2_precip_ccpav3_',cleadb,'_to_',cleade,'.nc',sep='')
		prec.nc <- nc_open(filename)
		dates <- ncvar_get(prec.nc, varid="yyyymmddhh_init")
		dates.fcstb <- ncvar_get(prec.nc, varid="yyyymmddhh_fcstb")
		nc_close(prec.nc)
		for (iyear in 1:12)  {
			for (iday in 1:31)  {
				jday <- ifelse(month==2 & iday==29, 28, iday)
				stss.ind.yrs <- which( ((dates%/%100)%%100)==jday & ((dates%/%10000)%%100)==month & (dates%/%1000000)!=years[iyear] )
				stss.ind.anal <- as.vector(outer(-2:2,match(dates.fcstb[stss.ind.yrs],yyyymmddhh_begin),'+'))
				apcp.hist.obs <- array(dim=c(nxa,nya,55))
				for (imb in 1:55)  {
					apcp.hist.obs[,,imb][!mask.domain] <- apcp.anal[,,stss.ind.anal[imb]][!mask.domain]
				}
				apcp.hist.obs.ranks <- apply(apcp.hist.obs,c(1,2),rank,ties.method="random")
				for (imb in 1:55)  {
					ranks.stss[,iyear,iday,ilt,imb] <- apcp.hist.obs.ranks[imb,,]
				}
			}
		}
	}
	for (iloc in 1:nloc)  {
		print(c(iloc,nloc))
		id.gpt <- which(!is.na(I) & (I==iloc))
		ngpt <- length(id.gpt)
		for (iyear in 1:12)  {
			for (iday in 1:31)  {
				if (any(is.na(mu.fcst[id.gpt,iyear,iday,]))) next
				apcp.ecc.traj <- apcp.stss.traj <- array(dim=c(ngpt,4,55))
				for (igpt in 1:ngpt)  {
					for (ilt in 1:4)  {
						shape <- (mu.fcst[id.gpt[igpt],iyear,iday,ilt]/sigma.fcst[id.gpt[igpt],iyear,iday,ilt])^2
						scale <- mu.fcst[id.gpt[igpt],iyear,iday,ilt]/shape
						shift <- shift.fcst[id.gpt[igpt],iyear,iday,ilt]
						for (rpt in 1:5) {
							fcst.smpl <- sort(pmax(shift+qgamma(runif(11),scale=scale,shape=shape),0))
							apcp.ecc.traj[igpt,ilt,11*rpt+(-10:0)] <- fcst.smpl[ranks.ecc[id.gpt[igpt],iyear,iday,ilt,]]
						}
						fcst.qt <- pmax(shift+qgamma(((1:55)-0.5)/55,scale=scale,shape=shape),0)
						apcp.stss.traj[igpt,ilt,] <- fcst.qt[ranks.stss[id.gpt[igpt],iyear,iday,ilt,]]
					}
				}
				ecc55.fcst.mean[iloc,month,iyear,iday,] <- apply(apcp.ecc.traj,3,mean)
				ecc55.fcst.max[iloc,month,iyear,iday,] <- apply(apcp.ecc.traj,3,max)
				stss55.fcst.mean[iloc,month,iyear,iday,] <- apply(apcp.stss.traj,3,mean)
				stss55.fcst.max[iloc,month,iyear,iday,] <- apply(apcp.stss.traj,3,max)
			}
		}
	}
}

save(mdss55.fcst.mean, ecc55.fcst.mean, stss55.fcst.mean, file="~/Desktop/Russian-River-CaseStudy/data/MeanFcstFields-55mb.Rdata")
save(mdss55.fcst.max, ecc55.fcst.max, stss55.fcst.max, file="~/Desktop/Russian-River-CaseStudy/data/MaxFcstFields-55mb.Rdata")




