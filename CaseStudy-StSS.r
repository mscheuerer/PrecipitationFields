
###
#   Code for generating the StSS and MDSS-RO ensemble forecast fields for the example shown in the paper and plot them
#
#      written by Michael Scheuerer, Feb 2018
#



library(ncdf4)
library(fields)
library(maps)


##  Load analyzed fields

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




##  Create space-time dependence template

load("~/Desktop/Russian-River-CaseStudy/FcstFields-ECC.Rdata")    # Univariate post-processing is the same as for ECC

gefsmean.verif <- apply(gefs.verif,c(1,2),mean)

stss.ind <- which(yyyymmddhh_begin%%1000000==11900 & yyyymmddhh_begin%/%1000000%in%c(2002:2009,2011:2013))    # Same mmddhh as in case study, different yyyy though
nmb <- length(stss.ind)

stss.ranks <- integer(nxa*nya*4*nmb)
dim(stss.ranks) <- c(nxa,nya,4,nmb)

apcp.stss.traj <- array(dim=c(nxa,nya,4,nmb))
apcp.hist.traj <- array(dim=c(nxa,nya,4,nmb))

for (ilt in 1:4)  {
	cat(paste('Lead time',ilt,'\n'))
	for (ix in 1:nxa)  {
		for (jy in 1:nya)  {
			if (mask.domain[ix,jy] | all(is.na(apcp.eccq.traj[ix,jy,ilt,]))) next
			apcp.hist.traj[ix,jy,ilt,] <- apcp.anal[ix,jy,stss.ind+ilt-1]
			fcst.qt <- sort(apcp.eccq.traj[ix,jy,ilt,])
			stss.ranks[ix,jy,ilt,] <- rank(apcp.hist.traj[ix,jy,ilt,],ties.method="random")
			apcp.stss.traj[ix,jy,ilt,] <- fcst.qt[stss.ranks[ix,jy,ilt,]]

		}
	}
}

save(lons.fcst, lats.fcst, gefsmean.verif, csgd.fcst.mean, mask.domain, stss.ranks, apcp.hist.traj, apcp.stss.traj, file="~/Desktop/Russian-River-CaseStudy/FcstFields-StSS.Rdata")






# Plot historic fields and StSS field

library(fields)
library(maps)

load("~/Desktop/Russian-River-CaseStudy/FcstFields-StSS.Rdata")

nmb <- 11
lead.seq <- 8:11
zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
breaks <- seq(0,zmax,,length(colors)+1)
colors.rank <- c('yellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','violet','violetred2','red3','red4')


for (k in 1:nmb)  {
	pdf(paste("~/Desktop/Russian-River-CaseStudy/StSS/CaseStudy-FcstField",k,"-StSS.pdf",sep=""), width=18, height=10.5)
	split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
	split.screen(c(3,4), screen=1) -> ind
	for (ilt in 1:4)  {
		screen(ind[ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Historic field',k),'')
		poly.image(lon, lat, pmin(apcp.hist.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[4+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Rank of historic field',k),'')
		poly.image(lon, lat, stss.ranks[,,ilt,k], breaks=seq(0.5,11.5,1), col=colors.rank, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[8+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('StSS member',k),'')
		poly.image(lon, lat, pmin(apcp.stss.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
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
}

order(apply(apcp.stss.traj,4,mean,na.rm=TRUE))











##  Create space-time dependence template for MDSS-RO implementation (fine scale post-processing, coarse-scale trajectory selection)

load("~/Desktop/Russian-River-CaseStudy/FcstFields-ECC.Rdata")    # Univariate post-processing is the same as for ECC
load("~/Desktop/Russian-River-CaseStudy/FcstFields-MDSS.Rdata")   # Load historic fields selected by MDSS algorithm

gefsmean.verif <- apply(gefs.verif,c(1,2),mean)

nmb <- 11

mdss.ranks <- integer(nxa*nya*4*nmb)
dim(mdss.ranks) <- c(nxa,nya,4,nmb)

for (ilt in 1:4)  {
	cat(paste('Lead time',ilt,'\n'))
	for (ix in 1:nxa)  {
		for (jy in 1:nya)  {
			if (mask.domain[ix,jy] | all(is.na(apcp.eccq.traj[ix,jy,ilt,]))) next
			fcst.qt <- sort(apcp.eccq.traj[ix,jy,ilt,])
			mdss.ranks[ix,jy,ilt,] <- rank(apcp.hist.traj[ix,jy,ilt,],ties.method="random")
			apcp.mdss.traj[ix,jy,ilt,] <- fcst.qt[mdss.ranks[ix,jy,ilt,]]

		}
	}
}




# Plot historic fields and MDSS-RO fields

library(fields)
library(maps)

lead.seq <- 8:11
zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
breaks <- seq(0,zmax,,length(colors)+1)
colors.rank <- c('yellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','violet','violetred2','red3','red4')


for (k in 1:nmb)  {
	pdf(paste("~/Desktop/Russian-River-CaseStudy/MDSS/CaseStudy-FcstField",k,"-MDSS-RO.pdf",sep=""), width=18, height=10.5)
	split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
	split.screen(c(3,4), screen=1) -> ind
	for (ilt in 1:4)  {
		screen(ind[ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Historic field',k),'')
		poly.image(lon, lat, pmin(apcp.hist.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		title(paste("lead time ", formatC(lead.seq[ilt]*6), " - ", formatC((lead.seq[ilt]+1)*6),"h", sep=''), cex.main=1.5)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[4+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('Rank of historic field',k),'')
		poly.image(lon, lat, mdss.ranks[,,ilt,k], breaks=seq(0.5,11.5,1), col=colors.rank, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		map("state", add=TRUE)
		box()
	}
	for (ilt in 1:4)  {
		screen(ind[8+ilt])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		ylab <- ifelse(ilt==1,paste('MDSS-RO member',k),'')
		poly.image(lon, lat, pmin(apcp.mdss.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
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
}










