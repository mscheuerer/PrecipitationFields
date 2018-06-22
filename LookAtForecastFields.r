
###
#   Code for displaying the ensemble forecast fields generated with 'GenerateFcstFields-ECC-StSS' and 'GenerateFcstFields-MDSS'
#
#      written by Michael Scheuerer, Mar 2018
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
nmb <- 11

mask.conus <- matrix(is.na(map.where("state",lon-0.03,lat)), nxa, nya)
mask.domain <- mask.conus | (lon < -123.9845) | (lon > -122.1095) | (lat < 38.38457) | (lat > 39.78888)


month <- 11
iyear <- 9
ilt <- 2



zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(ticks,"mm",sep=''),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
breaks <- seq(0,zmax,,length(colors)+1)



load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-mQ-SNP-",month,".Rdata",sep=""))

for (iday in 1:30)  {
	if (max(apcp.eccmq.snp.traj[,,iyear,iday,ilt,],na.rm=TRUE)<5) next
	split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
	split.screen(c(3,4), screen=1) -> ind
	for (k in 1:nmb)  {
		screen(ind[k])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		poly.image(lon, lat, pmin(apcp.eccmq.snp.traj[,,iyear,iday,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		title(paste('ECC-mQ-SNP member',k,',  lead = ',42+6*ilt,'h'))
		map("state", add=TRUE)
		box()
	}
	date.fcstb <- (2002:2013)[iyear]*1e6+month*1e4+iday*1e2
	verif.ind.anal <- match(date.fcstb,yyyymmddhh_begin)+7+ilt
	obs <- apcp.anal[,,verif.ind.anal]
	obs[mask.domain] <- NA
	screen(ind[12])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	poly.image(lon, lat, obs, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	title(paste('Analyzed field for ',yyyymmddhh_begin[verif.ind.anal]))
	map("state", add=TRUE)
	box()
	screen(2)
	image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.15,.85), col=colors, axis.args=axis.args)
	close.screen(all=TRUE)
	readline()
}





load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-MDSS-",month,".Rdata",sep=""))

for (iday in 1:30)  {
	if (max(apcp.mdss.traj[,,iyear,iday,ilt,k],na.rm=TRUE)<5) next
	split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
	split.screen(c(3,4), screen=1) -> ind
	for (k in 1:nmb)  {
		screen(ind[k])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		poly.image(lon, lat, pmin(apcp.mdss.traj[,,iyear,iday,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		title(paste('MDSS-SDA member',k,',  lead = ',42+6*ilt,'h'))
		map("state", add=TRUE)
		box()
	}
	date.fcstb <- (2002:2013)[iyear]*1e6+month*1e4+iday*1e2
	verif.ind.anal <- match(date.fcstb,yyyymmddhh_begin)+7+ilt
	obs <- apcp.anal[,,verif.ind.anal]
	obs[mask.domain] <- NA
	screen(ind[12])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'Analyzed field','')
	poly.image(lon, lat, obs, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	title(paste('Analyzed field for ',yyyymmddhh_begin[verif.ind.anal]))
	map("state", add=TRUE)
	box()
	screen(2)
	image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.15,.85), col=colors, axis.args=axis.args)
	close.screen(all=TRUE)
	readline()
}







load(paste("~/Desktop/Russian-River-CaseStudy/data/FcstFields-ECC-StSS-",month,".Rdata",sep=""))
load(paste("~/Desktop/Russian-River-CaseStudy/data/HistoricRanks-MDSS-",month,".Rdata",sep=""))

field <- matrix(NA,nxa,nya)

for (iday in 1:30)  {
	if (max(mu.fcst[,,iyear,iday,ilt],na.rm=TRUE)<1) next
	print(iday)
	split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
	split.screen(c(3,4), screen=1) -> ind
	shape <- (mu.fcst[,,iyear,iday,ilt]/sigma.fcst[,,iyear,iday,ilt])^2
	scale <- mu.fcst[,,iyear,iday,ilt]/shape
	shift <- shift.fcst[,,iyear,iday,ilt]
	for (k in 1:nmb)  {
		field[!mask.domain] <- pmax(shift[!mask.domain]+qgamma((ranks.mdss.ro[,,iyear,iday,ilt,k][!mask.domain]-0.5)/nmb,scale=scale[!mask.domain],shape=shape[!mask.domain]),0)
		screen(ind[k])
		par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
		poly.image(lon, lat, pmin(field,zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
		title(paste('MDSS-RO member',k,',  lead = ',42+6*ilt,'h'))
		map("state", add=TRUE)
		box()
	}
	date.fcstb <- (2002:2013)[iyear]*1e6+month*1e4+iday*1e2
	verif.ind.anal <- match(date.fcstb,yyyymmddhh_begin)+7+ilt
	obs <- apcp.anal[,,verif.ind.anal]
	obs[mask.domain] <- NA
	screen(ind[12])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'Analyzed field','')
	poly.image(lon, lat, obs, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	title(paste('Analyzed field for ',yyyymmddhh_begin[verif.ind.anal]))
	map("state", add=TRUE)
	box()
	screen(2)
	image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.15,.85), col=colors, axis.args=axis.args)
	close.screen(all=TRUE)
	readline()
}






