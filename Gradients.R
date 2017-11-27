# install.packages(c('RNetCDF', 'maps','mapdata','plotrix','akima','fields','OpenImageR','raster'), dependencies=T) # RUN THIS ONLY THE FIRST TIME
# note: Do you want to install from sources the package which needs compilation? YES

# DATA are now share usign DAT. Make sure to install DAT command line tools first. To install DAT, just open the OSX terminal and type:
# npm install -g dat
# if you have some weird error, type:
# sudo npm install -g dat


library(RNetCDF)
library(maps)
library(mapdata)
library(plotrix)
library(akima)
library(fields)
library(OpenImageR)
library(raster)

jet.colors <- colorRampPalette(c("blue4","royalblue4","deepskyblue3", "seagreen3", "yellow", "orangered2","darkred"))
jet.colors2 <- colorRampPalette(c("blue","white","red"))


#########################
### LOCATION OF FILES ###
#########################
path.to.data <- "~/Desktop/Gradients_Rcode/data/" # location of data for Ginger
path.to.data <- "~/Documents/DATA/SeaFlow/SF_GRADIENTS/Rcode-Ginger/data/" # location of data for Francois
setwd(path.to.data)

# Copy DAT data to the location
# system(paste("dat clone dat://bdbcddfe950911dd4abe607593681f030efb6821b93357adb4484bb7cf825735", path.to.data)) # need to do this only once, after that, data will sync automatically as long as DAT is open.

savepath <- "~/Desktop" # location of saved plots



#######################
### 1. SELECT DATA  ###
#######################

gradient <- 1 # value can be 1 or 2

out <- TRUE # FALSE # value can be TRUE (northward) or FALSE (southward)





####################
### 2. LOAD DATA ###
####################
if(gradient == 1){

  nc <- open.nc("wh300.nc")
  if(out){ start <- 1
           end <- 5633} ## OUT
  if(!out){ start <- 5633
            end <- 10785} ## BACK
 stat <- read.csv("stat1.csv")
 nut <- read.csv("Nitrogen_forODV.csv")
 if(out)  cur <- open.nc("oscar_vel8604.nc")
 if(!out) cur <- open.nc("oscar_vel8610.nc")
 if(out) sst <- open.nc("20160428090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc")
 if(!out) sst <- open.nc("20160501090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc")
 if(out) ctd <- read.csv("uCTD-OUT.csv")
 if(!out) ctd <- read.csv("uCTD-BACK.csv")
 sfl <- read.csv("time_lat_lon_temp_sal_Grad2.csv")
 }

if(gradient == 2){

  nc <- open.nc("os75bb.nc")
  if(out){ start <- 1
           end <- 3099} ## OUT
  if(!out){ start <- 3100
            end <- 4833} ## BACK
  stat <- read.csv("stat2.csv")
    stat[which(stat$lat < 20),'lat'] <- 33.25
  po4 <- read.delim('PO4__MGL1704.tsv')
  metals <- read.csv("High-resolution-NiCu.csv")
 if(out) cur <- open.nc("oscar_vel9001.nc")
 if(!out) cur <- open.nc("oscar_vel9006.nc")
 if(out) sst <- open.nc("20170531090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc")
 if(!out) sst <- open.nc("20170609090000-JPL-L4_GHRSST-SSTfnd-MUR-GLOB-v02.0-fv04.1.nc")
 sfl <- read.csv("time_lat_lon_temp_sal_Grad1.csv")
 }



#######################
### 3. PROCESS DATA ###
#######################

# ADCP
    dat <- read.nc(nc)
    v <- dat$v
    u <- dat$u
    lat <-dat$lat
    lon <-dat$lon
    depth <- dat$depth

    Lat <- rep(lat[start:end],each=dim(v)[1])
    Lon <- rep(lon[start:end],each=dim(v)[1])
    Depth <- rep(depth[,1], length(lat[start:end]))
    V <- unlist(list(v[,start:end]))
    U <- unlist(list(u[,start:end]))
    # V[is.na(V)] <- 0
    # U[is.na(U)] <- 0
    data <- data.frame(lon=Lon, lat =Lat,depth=Depth, v=V, u=U, s=sqrt(V^2 + U^2), theta=atan(V/U)*180/pi)
    data <- subset(data, depth > min(data$depth))
    data <- data[which(!is.na(data$lat)),]
    data[which(data$v > 0.6),'v'] <- 0.6

# OSCAR
    dat2 <- read.nc(cur)
    v2 <- dat2$v
      v2[is.na(v2)] <- 0
    u2 <- dat2$u
      u2[is.na(u2)] <- 0
    lat2 <- dat2$latitude
    lon2 <- dat2$longitude -360
    a <- expand.grid(lon2, lat2)


# SST
  ylat <- var.get.nc(sst, 'lat')
  xlon <- var.get.nc(sst, 'lon')
  z <- var.get.nc(sst, 'analysed_sst', unpack=TRUE)
      LonStartIdx <- min(which(xlon  > -175))
      LonEndIdx <- min(which(xlon > -140))
      LatStartIdx <- min(which(ylat  > 15))
      LatEndIdx <- min(which(ylat > 50))
    ylat <-  ylat[c(LatStartIdx:LatEndIdx)]
    xlon <-  xlon[c(LonStartIdx:LonEndIdx)]
    z <-  z[c(LonStartIdx:LonEndIdx),c(LatStartIdx:LatEndIdx)]
    z <- z - 273.15
    z[which(is.na(z))] <- min(z, na.rm=T)

# SeaFlow
if(out) stat2 <- stat[1:which(stat$lat == max(stat$lat, na.rm=T))[1],]
if(!out) stat2 <- stat[which(stat$lat == max(stat$lat, na.rm=T))[1]:nrow(stat),]

# CTD
data.sigma <- interp(ctd$lat, -ctd$pressure, ctd$sigma, duplicate="mean", nx=200)
data.sal <- interp(ctd$lat, -ctd$pressure, ctd$salinity, duplicate="mean", nx=200)
data.temp <- interp(ctd$lat, -ctd$pressure, ctd$temp, duplicate="mean", nx=200)

# PO4 (only Gradients 2.0)
stat$time <- as.POSIXct(stat$time, format = "%FT%T", tz = "GMT")
po4$time <- as.POSIXct(po4$time, format = "%FT%T", tz = "GMT")
id <- findInterval(po4$time, stat$time)
po4$lat <- stat[id,'lat']
po4$lon <- stat[id,'lon']
if(out) po4.2 <- po4[1:which(po4$lat == max(po4$lat, na.rm=T))[1],]
if(!out) po4.2 <- po4[which(po4$lat == max(po4$lat, na.rm=T))[1]:nrow(po4),]

# Fe/NO3 RATIO
metals <- metals[order(metals$Latitude),]
po4.2 <- po4.2[order(po4.2$lat),]
id <- findInterval(metals$Latitude, po4.2$lat)
metals$po4 <- po4.2[id, 'PO4']


# NO3 (only for Gradients 1.0)
nut <- nut[!is.na(nut$Nitrate..uM.),]
data.nut <- interp(nut$Latitude, -nut$Depth..m., nut$Nitrate..uM., duplicate="mean", nx=200)


# Underway SST
if(!out) sfl2 <- sfl[1:which(sfl$lat == max(sfl$lat, na.rm=TRUE))[1],]
if(out) sfl2 <- sfl[which(sfl$lat == max(sfl$lat, na.rm=TRUE))[1]:nrow(sfl),]

id <- findInterval(pro$time, sfl2$time)
pro$temp <- sfl2[id,'temp']

















# Note: EXECUTE only the plots of interest, do not run the entire code.


################
### 4. PLOTS ###
################

### CTD
if(out) png(paste0(savepath, "/Gradient-",gradient,"_uCTD-OUT.png"), width=114*2, height=114, res=600, units="mm", pointsize=8)
if(!out) png(paste0(savepath, "/Gradient-",gradient,"_uCTD-BACK.png"), width=114*2, height=114, res=600, units="mm", pointsize=8)

par(mfrow=c(3,1), mar=c(2,2,1,1))
image(data.sal, col=jet.colors(100), main="Salinity");contour(data.sal, add=T)
image(data.temp, col=jet.colors(100), main="Temperature");contour(data.temp, add=T)
image(data.sigma, col=jet.colors(100), main="Density");contour(data.sigma, add=T)

dev.off()





### ADCP
if(out) png(paste0(savepath, "/Gradient-",gradient,"_ADCP_OUT.png"), width=114*2, height=114/1.5, res=600, units="mm", pointsize=8)
if(!out) png(paste0(savepath, "/Gradient-",gradient,"_ADCP_BACK.png"), width=114*2, height=114/1.5, res=600, units="mm", pointsize=8)

par(mfrow=c(2,1), mar=c(2,2,1,6),oma=c(2,2,2,2),cex=1.2)
n <- 50
data2 <- subset(data, depth < 100)
param <- 'u'
plot(data2$lat, -data2$depth, pch=15, col=jet.colors2(n)[cut(data2[,param],n)],main=paste(param), cex=0.7, xlim=c(24,36))
ylim <- par('usr')[c(3,4)]
xlim <- par('usr')[c(1,2)]
color.legend(xlim[2]- 0.1 , ylim[1], xlim[2], ylim[2], legend=c("westward","eastward"), rect.col=jet.colors2(n), gradient='y',align='rb',cex=1.2)

param <- 'v'
plot(data2$lat, -data2$depth, pch=15, col=jet.colors2(n)[cut(data2[,param],n)], main=paste(param), cex=0.7, xlim=c(24,36))
ylim <- par('usr')[c(3,4)]
xlim <- par('usr')[c(1,2)]
color.legend(xlim[2]- 0.1 , ylim[1], xlim[2], ylim[2], legend=c("southward","northward"), rect.col=jet.colors2(n), gradient='y',align='rb',cex=1.2)

dev.off()



### MAPS of SST, CURRENTS + SEAFLOW

# location of the ZOOM region (rectangle)
y.rect <- c(-158.5,-157.5)
x.rect <- c(29,38)




if(out) png(paste0(savepath, "/Gradient-",gradient,"_SST-OSCAR-SeaFlowOUT.png"), width=114*2, height=114, res=600, units="mm", pointsize=8)
if(!out) png(paste0(savepath, "/Gradient-",gradient,"_SST-OSCAR-SeaFlowBACK.png"), width=114*2, height=114, res=600, units="mm", pointsize=8)

nf <- layout(matrix(c(1,1,1,2,3,4), 1,6, byrow=T))

par(cex=1, oma=c(1,1,1,3), mar=c(3,2,1,1))
if(out)  plot(1,1, pch=NA,asp=1, xlim=c(-165, -150), ylim=c(20,45), xlab="longitude", ylab='latitude', main=paste0('Gradient ', gradient, ".0 OUT"))#, ylim=c(20,30))
if(!out)  plot(1,1, pch=NA,asp=1, xlim=c(-165, -150), ylim=c(20,45), xlab="longitude", ylab='latitude', main=paste0('Gradient ', gradient, ".0 BACK"))#, ylim=c(20,30))
  image(x=xlon,y=ylat, z=z, xaxt='n', yaxt='n', xlab=NA, ylab=NA, col=jet.colors(100), add=TRUE)
  arrows(x0=a$Var1, y0=a$Var2, x1=a$Var1+as.vector(u2)*2, y1=a$Var2+as.vector(v2)*2, length=0.025, lwd=0.5)
  contour(x=xlon,y=ylat, z=z, levels=18, col='white',add=T)
  lines(c(-158,-158),range(stat$lat, na.rm=T), ,lwd=2)
  map('worldHires', border=NA, fill=TRUE, col='pink', add=TRUE)
  # ylim <- par('usr')[c(3,4)]
  # xlim <- par('usr')[c(1,2)]
  # color.legend(xlim[2]- 0.5 , ylim[1], xlim[2], ylim[2], legend=round(seq(0, max(z, na.rm=T), length.out=5)), rect.col=jet.colors(100), gradient='y',align='rb',cex=1.2)
  # mtext(substitute(paste("Temp (Deg C)")), side=4, line=1, cex=1.2)
  polygon(c(y.rect,rev(y.rect)),rep(x.rect, each=2),lwd=2, border='red3')


  df <- subset(stat2, pop == 'prochloro' &  lat < x.rect[2]+1 & lat > x.rect[1]-1 )
  plot(df$lon, df$lat, col=jet.colors(100)[cut(log10(df$abundance), 100)], pch=16,xlim=y.rect,ylim=x.rect,asp=1, xlab="longitude", ylab=NA, main="Prochloro")#, ylim=c(20,30))
  #contour(x=xlon,y=ylat, z=z, levels=18, lwd=2,col='darkgrey',,add=T)
  arrows(x0=a$Var1, y0=a$Var2, x1=a$Var1+as.vector(u2)*2, y1=a$Var2+as.vector(v2)*2,  length=0.025, lwd=0.7)
  # ylim <- par('usr')[c(3,4)]
  # xlim <- par('usr')[c(1,2)]
  # zlim <- range(log10(df$abundance),na.rm=T)
  # color.legend(xlim[2]- 0.1 , ylim[1], xlim[2], ylim[2], legend=round(10^seq(zlim[1],zlim[2], length.out=5)), rect.col=jet.colors(100), gradient='y',align='rb',cex=1.2)

  df <- subset(stat2, pop == 'synecho' &  lat < x.rect[2]+1 & lat > x.rect[1]-1)
  plot(df$lon, df$lat, col=jet.colors(100)[cut(log10(df$abundance), 100)], pch=16,xlim=y.rect,ylim=x.rect,asp=1, xlab="longitude", ylab=NA, main="Synecho")#, ylim=c(20,30))
  #contour(x=xlon,y=ylat, z=z, levels=18, lwd=2,col='darkgrey',add=T)
  arrows(x0=a$Var1, y0=a$Var2, x1=a$Var1+as.vector(u2)*2, y1=a$Var2+as.vector(v2)*2, length=0.025, lwd=0.7)
  # ylim <- par('usr')[c(3,4)]
  # xlim <- par('usr')[c(1,2)]
  # zlim <- range(log10(df$abundance),na.rm=T)
  # color.legend(xlim[2]- 0.1 , ylim[1], xlim[2], ylim[2], legend=round(10^seq(zlim[1],zlim[2], length.out=5)), rect.col=jet.colors(100), gradient='y',align='rb',cex=1.2)

  if(gradient ==2) df <- subset(stat2, pop == 'picoeuk' &  lat < x.rect[2]+1 & lat > x.rect[1]-1)
  if(gradient ==1) df <- subset(stat2, pop == 'picoeuks' &  lat < x.rect[2]+1 & lat > x.rect[1]-1)

  plot(df$lon, df$lat, col=jet.colors(100)[cut(log10(df$abundance), 100)], pch=16,xlim=y.rect,ylim=x.rect,asp=1, xlab="longitude", ylab=NA, main="Picoeuks")#, ylim=c(20,30))
  #contour(x=xlon,y=ylat, z=z, levels=18, lwd=2,col='darkgrey',add=T)
  arrows(x0=a$Var1, y0=a$Var2, x1=a$Var1+as.vector(u2)*2, y1=a$Var2+as.vector(v2)*2, length=0.025, lwd=0.7)
  # ylim <- par('usr')[c(3,4)]
  # xlim <- par('usr')[c(1,2)]
  # zlim <- range(log10(df$abundance),na.rm=T)
  # color.legend(xlim[2]- 0.1 , ylim[1], xlim[2], ylim[2], legend=NA, rect.col=jet.colors(100), gradient='y',align='rb',cex=1.2)
  mtext(substitute(paste("Abundance")), side=4, line=1, cex=1.2)

dev.off()




### NO3 profile
png(paste0(savepath, "/Gradient-",gradient,"_NO3.png"), width=114*2, height=114, res=600, units="mm", pointsize=8)

par(mfrow=c(2,1), mar=c(4,4,1,1))
# image(data.nut, col=jet.colors(100), main="NO3 (µmol L-1)");contour(data.nut, add=T)
# points(nut$Latitude, -nut$Depth..m., pch=16, col='grey')
plot(nut$Latitude, -nut$Depth..m., pch=16, col=jet.colors(100)[cut(nut$Nitrate..uM.,100)],cex=4, ylab='Depth (m)', xlab="Latitude")

dev.off()




plot(metals$Latitude, metals$Fe, xlim=xlim, yaxt='n', xaxt='n', xlab=NA, ylab=NA, type='p',pch=21,col=1, bg='grey')


### NUT CONC + SEAFLOW
xlim <- c(22,45)
syn <- subset(stat, pop == 'synecho')

if(out) png(paste0(savepath, "/Gradient-",gradient,"_NUT-SeaFlowOUT.png"), width=114*2, height=114, res=600, units="mm", pointsize=8)
if(!out) png(paste0(savepath, "/Gradient-",gradient,"_NUT-SeaFlowBACK.png"), width=114*2, height=114, res=600, units="mm", pointsize=8)

par(mar=c(5,5,2,5), oma=c(0,0,0,5))
plot(syn$lat, syn$abundance, xlim=xlim, ylim=c(0,150), xlab=NA, ylab=NA,type='p', pch=21, col=2, bg='orange')
mtext("Abundance (cells µL-1)", 2, line=3)
par(new=T)
plot(metals$Latitude, metals$Fe, xlim=xlim, yaxt='n', xaxt='n', xlab=NA, ylab=NA, type='p',pch=21,col=1, bg='grey')
axis(4)
mtext("Fe (nmol µL-1)", 4, line=3)
par(new=T)
plot(po4.2$lat, po4.2$PO4, xlim=xlim, col=3, yaxt='n', xaxt='n', xlab=NA, ylab=NA, type='l',lwd=2)
axis(4, line=5)
mtext("PO4 (nmol µL-1)", 4, line=7)
mtext('Latitude',side=1, line=3)
legend("topleft", c("Synechocococus abundance", "Fe concentration", "PO4 concentration"), pch=c(21,21,NA), pt.bg=c('orange','grey',NA), lwd=c(NA,NA,2), col=c(2,1,3), bty='n')

par(new=T)
plot(metals$Latitude, metals$po4/metals$Fe, xlim=xlim, yaxt='n', xaxt='n', xlab=NA, ylab=NA, type='p',pch=21,col=1, bg='purple')
axis(4, line=5)

dev.off()


### NUT RATIO + SEAFLOW
xlim <- c(22,45)
syn <- subset(stat, pop == 'synecho')
pico <- subset(stat, pop == 'picoeuk')

png(paste0(savepath, "/Gradient-",gradient,"_RATIO-SeaFlow.png"), width=114*2, height=114, res=600, units="mm", pointsize=8)

par(mar=c(5,5,2,5), oma=c(0,0,0,5))
plot(syn$lat, syn$abundance, xlim=xlim, ylim=c(0,150), xlab=NA, ylab=NA,type='p', pch=21, col=2, bg='orange')
par(new=T)
plot(pico$lat, pico$abundance, xlim=xlim, xlab=NA, ylab=NA,type='p', pch=21, col='darkgreen', bg='green', yaxt='n')
mtext("Abundance (cells µL-1)", 2, line=3)
par(new=T)
plot(metals$Latitude, metals$po4/metals$Fe, xlim=xlim, yaxt='n', xaxt='n', xlab=NA, ylab=NA, type='p',pch=21,col=1, bg='purple', ylim=c(0,2500))
axis(4)
mtext("PO4 / Fe ratio", 4, line=3)
mtext('Latitude',side=1, line=3)
legend("topleft", c("Synechocococus abundance", "Picoeukaryote abundance", "PO4 / Fe ratio"), pch=c(21,21,21),=c('orange','green','purple'), col=c(2,'darkgreen',3), bty='n')

dev.off()



### FE, CU and NI conc
par(mfrow=c(3,1))
plot(metals$Latitude, metals$Fe,pch=21,col=1, bg='grey')
plot(metals$Latitude, metals$Ni,pch=21,col=1, bg='red3')
plot(metals$Latitude, metals$Cu,pch=21,col=1, bg='seagreen3')

### FE / CU ratio
par(mfrow=c(2,1))
plot(metals$Latitude, metals$Cu/metals$Fe,pch=21,col=1, bg='gold3', main=c("Cu / Fe ratio"))
plot(metals$Latitude, metals$Fe/metals$Cu,pch=21,col=1, bg='orangered2', main=c("Fe / Cu ratio"))

### Underway SST + SAL
par(mfrow=c(2,1))
plot(sfl2$lat, sfl2$temp,col=1, bg='grey', type='l')
abline(h=18, col=2, lty=2)
plot(sfl2$lat, sfl2$sal,col=1, bg='grey', type='l') # SAL data for Gradient 1.0 not great
