
###
#   Code for generating plots that illustrate the methodology in the paper
#
#      written by Michael Scheuerer, Apr 2018
#



library(maps)
library(fields)



load(paste('~/Desktop/QPF-T2M-MV/data/ensemble_forecasts_012_018.Rdata',sep=''))


apcp.fcst, t2m.fcst, lons.fcst, lats.fcst, dates




gridpoints.all <- as.matrix(expand.grid(lons.fcst[12:22,1],lats.fcst[1,30:38]))
gridpoints.rrb <- as.matrix(expand.grid(lons.fcst[15:18,1],lats.fcst[1,33:35]))[-9,]


pdf("~/Desktop/Russian-River-CaseStudy/StudyArea.pdf", width=5, height=5)
par(mar=c(0,0,0,0))
map("state", xlim=c(-124,-121.5), ylim=c(37.5,40), myborder=c(0.01,0.01))
points(-122.4, 37.77, pch=19, cex=1.8)
text(-122.95, 37.7, "San Francisco", cex=1.2)
for (i in 1:nrow(gridpoints.all))  {
	rect(gridpoints.all[i,1]-0.234375,gridpoints.all[i,2]-0.2340527,gridpoints.all[i,1]+0.234375,gridpoints.all[i,2]+0.2340527, lwd=0.5)
}
for (i in 1:nrow(gridpoints.rrb))  {
	rect(gridpoints.rrb[i,1]-0.234375,gridpoints.rrb[i,2]-0.2340527,gridpoints.rrb[i,1]+0.234375,gridpoints.rrb[i,2]+0.2340527, border=4, lwd=2)
}
dev.off()



library(ncdf4)

filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/precip_06h_CCPA_2p5km_RRB.nc'
prec.nc <- nc_open(filename)
#names(prec.nc$var)
yyyymmddhh_begin <- ncvar_get(prec.nc, varid="yyyymmddhh_begin")
lon <- ncvar_get(prec.nc, varid="lons")
lat <- ncvar_get(prec.nc, varid="lats")
apcp.anal <- ncvar_get(prec.nc, varid="apcp_anal")
nc_close(prec.nc)

nd <- length(yyyymmddhh_begin)
nx <- nrow(lon)
ny <- ncol(lon)

mask <- matrix(!(map.where("state",lon,lat)=="california"),nx,ny)


for (i in 13090:nd)  {
	if(yyyymmddhh_begin[i]%/%10000!=201012) next
	if(all(is.na(apcp.anal[,,i]))) next
	if(max(apcp.anal[,,i],na.rm=TRUE)<10) next
	image.plot(lon-360, lat, ifelse(mask,NA,apcp.anal[,,i]), col=tim.colors(50))
	map("state", add=TRUE)
	print(i)
	print(yyyymmddhh_begin[i])
	readline("")
}





###
#   Illustration of the basis functions


library(ncdf4)
library(fields)
library(maps)

filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/precip_06h_CCPA_2p5km_RRB.nc'
prec.nc <- nc_open(filename)
lon <- ncvar_get(prec.nc, varid="lons")-360
lat <- ncvar_get(prec.nc, varid="lats")
nc_close(prec.nc)

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


nmb <- 11

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
dim(basis.fcts) <- c(nxa,nya,nmb)
basis.fcts <- sweep(basis.fcts,c(1,2),apply(basis.fcts,c(1,2),sum),'/')

x <- seq(-1,1,,100)
y <- ifelse(abs(x)>r.lon,0,(1-(abs(x)/r.lon)^3)^3)

colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
breaks <- seq(0,1,,length(colors)+1)
screen.seq <- c(10,11,12,5,6,7,8,1,2,3,4)


pdf("~/Desktop/Russian-River-CaseStudy/MDSS-BasisFunctions.pdf", width=12, height=7)
split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
split.screen(c(3,4), screen=1) -> ind
for (imb in 1:11)  {
	screen(ind[screen.seq[imb]])
	par(mar=c(1,2,1,0), mgp=c(0.8,1,0))
	poly.image(lon, lat, basis.fcts[,,imb], breaks=breaks, zlim=c(0,1), col=colors, xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
#	box()
	for (i in 1:nrow(gridpoints.rrb))  {
		rect(coords.grid[i,1]-0.234375,coords.grid[i,2]-0.2340527,coords.grid[i,1]+0.234375,coords.grid[i,2]+0.2340527, border=1, lwd=1, lty=2)
	}
	rect(coords.grid[imb,1]-0.234375,coords.grid[imb,2]-0.2340527,coords.grid[imb,1]+0.234375,coords.grid[imb,2]+0.2340527, border=1, lwd=2)
}
screen(ind[9])
par(mar=c(3,4,1,0))
plot(x, y, type='l', col='navy', xlab='', ylab='', lwd=2, cex.axis=0.75)
screen(2)
image.plot(zlim=c(0,1), legend.only=TRUE, smallplot=c(.27,.42,.15,.85), col=colors)
close.screen(all=TRUE)
dev.off()








###
#   Illustration of the ECC-mQ-SNP quantile mapping


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

ilt <- 2


pdf("~/Desktop/Russian-River-CaseStudy/ECC-mapping.pdf", width=12, height=3.2)
par(mar=c(5,4,1,3), mfrow=c(1,3))
poly.image(lon, lat, pmin(apcp.eccq.traj[,,ilt,k],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='ECC-Q member 2, lead time 54 - 60h', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.5)
map("state", add=TRUE)
box()
points(lon[15:16,72], lat[15:16,72], pch=5, col=2, cex=1, lwd=2)
plot(gefs.ip.verif[15,72,ilt,], apcp.eccq.traj[15,72,ilt,], pch=18, xlim=c(-1.5,16), ylim=c(-2.5,35), xlab='interpolated GEFS forecasts', ylab='calibrated forecast sample', axes=FALSE, cex.lab=1.2)
axis(1, at=seq(0,15,5), labels=paste(seq(0,15,5),'mm'), cex.axis=1.2)
axis(2, at=seq(0,30,10), labels=paste(seq(0,30,10),'mm'), cex.axis=1.2)
box()
text(gefs.ip.verif[15,72,ilt,]+c(0.3,-0.1,-0.2,-0.3,0.6,0,0.5,0,0.5,0.3,0.3), apcp.eccq.traj[15,72,ilt,]+ifelse((1:11)%in%c(5,7,9),-1.2,1.55), labels=1:11, cex=1.2)
points(gefs.ip.verif[15,72,ilt,], apcp.eccmq.snp.traj[15,72,ilt,], pch=18, col=2)
lines(sort(gefs.ip.verif[15,72,ilt,]), sort(apcp.eccmq.snp.traj[15,72,ilt,]), col=2, lty=3)
points(rep(-1.5,11), apcp.eccq.traj[15,72,ilt,], pch=18)
points(rep(-1.0,11), apcp.eccmq.snp.traj[15,72,ilt,], pch=18, col=2)
points(gefs.ip.verif[15,72,ilt,], rep(-2.7,11), pch=18, col=8)
lines(c(0,16), c(0,0), lty=2)
lines(c(0,0), c(0,35), lty=2)
legend(1, 34, c('ECC-Q ensemble','ECC-mQ-SNP ensemble'), pch=18, col=c(1,2), cex=1.2)
plot(gefs.ip.verif[16,72,ilt,], apcp.eccq.traj[16,72,ilt,], pch=18, xlim=c(-1.5,16), ylim=c(-2.5,35), xlab='interpolated GEFS forecasts', ylab='calibrated forecast sample', axes=FALSE, cex.lab=1.2)
axis(1, at=seq(0,15,5), labels=paste(seq(0,15,5),'mm'), cex.axis=1.2)
axis(2, at=seq(0,30,10), labels=paste(seq(0,30,10),'mm'), cex.axis=1.2)
box()
text(gefs.ip.verif[16,72,ilt,]+c(0.3,0.1,-0.2,-0.3,0.6,-0.1,0.5,0,0.5,0.3,0.3), apcp.eccq.traj[16,72,ilt,]+ifelse((1:11)%in%c(5,7,9),-1.2,1.55), labels=1:11, cex=1.2)
points(gefs.ip.verif[16,72,ilt,], apcp.eccmq.snp.traj[16,72,ilt,], pch=18, col=2)
lines(sort(gefs.ip.verif[16,72,ilt,]), sort(apcp.eccmq.snp.traj[16,72,ilt,]), col=2, lty=3)
points(rep(-1.5,11), apcp.eccq.traj[16,72,ilt,], pch=18)
points(rep(-1.0,11), apcp.eccmq.snp.traj[16,72,ilt,], pch=18, col=2)
points(gefs.ip.verif[16,72,ilt,], rep(-2.7,11), pch=18, col=8)
lines(c(0,16), c(0,0), lty=2)
lines(c(0,0), c(0,35), lty=2)
legend(1, 34, c('ECC-Q ensemble','ECC-mQ-SNP ensemble'), pch=18, col=c(1,2), cex=1.2)
dev.off()





###
#   Illustration of the FTE metric


library(ncdf4)
library(fields)
library(maps)

load("~/Desktop/Russian-River-CaseStudy/FcstFields-ECC.Rdata")

ilt <- 3

verif.ind.anal <- 11758:11761 

date.string <- c("Jan 19, 2010, 0Z - 6Z", "Jan 19, 2010, 6Z - 12Z", "Jan 19, 2010, 12Z - 18Z", "Jan 19, 2010, 18Z - 0Z")


filename <- '/Users/mscheuerer/Desktop/Russian-River-CaseStudy/precip_06h_CCPA_2p5km_RRB.nc'
prec.nc <- nc_open(filename)
yyyymmddhh_begin <- ncvar_get(prec.nc, varid="yyyymmddhh_begin")
lon <- ncvar_get(prec.nc, varid="lons")-360
lat <- ncvar_get(prec.nc, varid="lats")
apcp.anal <- ncvar_get(prec.nc, varid="apcp_anal")
nc_close(prec.nc)



nmb <- 11
zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
breaks <- seq(0,zmax,,length(colors)+1)


pdf("~/Desktop/Russian-River-CaseStudy/FTE-example.pdf", width=12, height=4.5)
par(mfrow=c(1,2), mar=c(1,4,4,4))
obs <- apcp.anal[,,verif.ind.anal[ilt]]
obs[mask.domain] <- NA
image.plot(lon, lat, obs, breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, legend.cex=0.3, axis.args=axis.args)
title(paste('Analyzed field,',date.string[ilt]), cex.main=1.5)
map("state", add=TRUE)
box()
poly.image(lon, lat, 1*(obs>10), axes=FALSE, xlim=c(-124,-122.15), ylim=c(38.43,39.75), col=c('white','black'), xlab='', ylab='')
title('Exceedance of 25 mm precipitation', cex.main=1.5)
map("state", add=TRUE)
box()
dev.off()




pdf("~/Desktop/Russian-River-CaseStudy/FTE-example-ECC-T.pdf", width=8, height=5)
par(mfrow=c(3,4), mar=c(1,1,3,1))
for (k in 1:nmb)  {
	poly.image(lon, lat, 1*(apcp.ecct.traj[,,ilt,k]>10), zlim=c(0,1), col=c('white','black'), xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	title(paste('ECC-T member',k))
	map("state", add=TRUE)
	box()
}
obs <- apcp.anal[,,verif.ind.anal[ilt]]
obs[mask.domain] <- NA
poly.image(lon, lat, 1*(obs>10), col=c('white','black'), xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE)
title('Analyzed field')
map("state", add=TRUE)
box()
dev.off()




load("~/Desktop/Russian-River-CaseStudy/FcstFields-StSS.Rdata")

pdf("~/Desktop/Russian-River-CaseStudy/FTE-example-StSS.pdf", width=8, height=5)
par(mfrow=c(3,4), mar=c(1,1,3,1))
for (k in 1:nmb)  {
	poly.image(lon, lat, 1*(apcp.stss.traj[,,ilt,k]>10), zlim=c(0,1), col=c('white','black'), xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	title(paste('StSS member',k))
	map("state", add=TRUE)
	box()
}
obs <- apcp.anal[,,verif.ind.anal[ilt]]
obs[mask.domain] <- NA
poly.image(lon, lat, 1*(obs>10), col=c('white','black'), xlab='', ylab='', xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE)
title('Analyzed field')
map("state", add=TRUE)
box()
dev.off()









###
#  Plot of StSS and ECC-mQ-SNP field for proposal


library(ncdf4)
library(fields)
library(maps)


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


load("~/Desktop/Russian-River-CaseStudy/FcstFields-StSS.Rdata")
load("~/Desktop/Russian-River-CaseStudy/FcstFields-MDSS.Rdata")

nmb <- 11
lead.seq <- 8:11
zmax <- 71

ticks <- seq(0,70,10)
axis.args <- list(at=ticks,labels=paste(formatC(ticks,width=2),"mm"),cex.axis=1.3)
colors <- colorRampPalette(c('white','wheat','tan1','yellow','greenyellow','green1','springgreen3','springgreen4','skyblue1','royalblue','mediumblue','slateblue',
'violet','violetred2','red3','red4'))(51)
breaks <- seq(0,zmax,,length(colors)+1)



pdf("~/Desktop/CaseStudy-FcstField-StSS-MDSS.pdf", width=18, height=7)
split.screen( rbind(c(0,.9,0,1), c(.9,1,0,1)))
split.screen(c(2,4), screen=1) -> ind
for (ilt in 1:4)  {
	screen(ind[ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'StSS member 5','')
	poly.image(lon, lat, pmin(apcp.stss.traj[,,ilt,5],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
for (ilt in 1:4)  {
	screen(ind[4+ilt])
	par(mar=c(1,2,3,0), mgp=c(0.8,1,0))
	ylab <- ifelse(ilt==1,'MDSS member 11','')
	poly.image(lon, lat, pmin(apcp.mdss.traj[,,ilt,11],zmax), breaks=breaks, zlim=c(0,zmax), col=colors, xlab='', ylab=ylab, xlim=c(-124,-122.15), ylim=c(38.43,39.75), axes=FALSE, cex.lab=1.3)
	map("state", add=TRUE)
	box()
}
screen(2)
image.plot(zlim=c(0,zmax), legend.only=TRUE, smallplot=c(.27,.42,.15,.85), col=colors, axis.args=axis.args)
close.screen(all=TRUE)
dev.off()




