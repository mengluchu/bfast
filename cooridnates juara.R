data.shape<-readShapePoly("C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/juara",IDvar="GEOCODIG_M",proj4string=CRS("+proj=longlat +ellps=aust_SA +no_defs "))
summary(data.shape)

#GDAL way:
shape=readOGR("C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/juara.shp", layer="juara") #will load the shapefile to your dataset.

extent(shape.utm)
coord.subdeter<-coordinates(subdeter)
spdeter<-SpatialPoints(coord.subdeter)
spmodist<-SpatialPoints(coord.changeinmodist)
spmodist51<-SpatialPoints(coordinates(changeinmot0.5.1))
spxy2d<-SpatialPoints(coordinates(xy2d))
#modis sinusoidal
#+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs
extent(deter1)
proj4string(xyc)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
proj4string(subset1)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'

proj4string(spmodist51)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
proj4string(around520)<-'+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
around520
#utm zone 21 south
tif1<-spTransform(tif1,CRS("+proj=utm +zone=21 +south"))
modis.t51<-spTransform(spmodist51,CRS("+proj=utm +zone=21 +south"))
plot(sb1,add=TRUE,color='blue')
shape.utm<-spTransform(shape,CRS("+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(shape.utm)
around520mo<-spTransform(around520,CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"))

dft1.utm<-spTransform(dft1,  CRS("+proj=utm +zone=21 +south"))

plot(shape.modis)
plot(dft1.utm,add=TRUE)

extent(dft1)
extent(shape.utm)
startm<-getxy(modisextent[1],modisextent[2],231.6564)
endm<-getxy(modisextent[3],modisextent[4],231.6564)
sexy<-as.data.frame(rbind(startm,endm))
sexy
            
ssexy<-SpatialPoints(coordinates(sexy))
coordinates(sx,sy)
?coordinates
plot(shape.utm)
plot(dc2,add=TRUE,col='green')
plot(mos,add=TRUE,col='pink')
plot(mot,add=TRUE,col='red')

plot(dft1.utm,add=TRUE,col='blue')
subset(StudentData, SchoolName=="Pine Tree Elementary" | Grade==3)
str(dc3)
extent(shape.utm)
xchange<-dc3@coords[,1]
ychange<-dc3@coords[,2]
xchangelat<-dc@coords[,1]
ychangelong<-dc@coords[,2]
xyll<-as.data.frame(cbind(xchangelat,ychangelong))
getcrMatrix <- function(colrowid.Matrix, pixelSize){
  #Returns the coords (MODIS synusoidal) of the center of the given pixel
  #SR-ORG:6974
  x <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  y <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  corner.ul.x <- -20015109.354
  corner.ul.y <- 10007554.677
  x<- floor( (colrowid.Matrix[,1] -corner.ul.x - pixelSize/2)/pixelSize) 
  y<- floor(-(colrowid.Matrix[,2] -corner.ul.y + pixelSize/2)/pixelSize) 
  
  cbind(x,y)
}
#y <- corner.ul.y - (pixelSize/2) - (colrowid.Matrix[,2] * pixelSize)
xy1<-cbind(xchange,ychange)
xy2<-getcrMatrix(xy1,231.6564)
xy2d<-as.data.frame(xy2)
xy2d[&xy2d$x<59079&xy2d$y>48210&xy2d$y<48360]
subdeter<-as.data.frame(xy2d[xy2d$x > 58930 & xy2d$x < 59079 & xy2d$y > 48210 & xy2d$y < 48360,])
summary(sub1)
#match<-subdeter$x%in%xychangeinmodist$xct 
coord.changeinmodist<-coordinates(xychangeinmodist)
coord.subdeter<-coordinates(subdeter)
?over
spmodist<-SpatialPoints(coord.changeinmodist)
spdeter<-SpatialPoints(coord.subdeter)

spmodist[!is.na(over(spmodist,spdeter))]
str(xychangeinmodist)
changeintrend<-which(test2nds[,,]!=0,arr.ind=TRUE)
subdeter$x
xct<-changeintrend[,1]+58929
yct<-changeintrend[,2]+48210
xychangeinmodist<-as.data.frame(cbind(xct,yct))
changeinmo<-getxyMatrix(xychangeinmodist,231.6564)
changeinseasonality<-which(test2ndt[,,]!=0,arr.ind=TRUE)
subdeter$x
xct<-changeinseasonality[,1]+58929
yct<-changeinseasonality[,2]+48210
xychangeinmodist<-as.data.frame(cbind(xct,yct))
changeinmos<-getxyMatrix(xychangeinmodist,231.6564)

saveGIF({
  
  ani.options(interval =1)
 

 
  jpeg('compare0.05t.jpg', height=20, width=10, res=1200,unit="in")

 
    plot(shape.utm)
 # title('change in deter (all)')
    points(dc2,col='green',cex=0.02)
 # plot(shape.utm,main='change in trend')
 # title('change in trend')
 # ani.pause()
    points(mos,col='gold',cex=0.02) #trend
 # plot(shape.utm,main='change in seasonality')
#  title('change in seasonality')
 # ani.pause(ani.options(interval =0.05))
    points(mot,col='red',cex=0.02) #seasonality
points(modis.t2,col='grey',cex=0.02) # pvalue 0.05
points(modis.s2,col='purple',cex=0.02) #pvalue 0.05
points(modis.t3,col='pink',cex=0.02) #p-value 0.01
points(modis.s3,col='cadetblue',cex=0.02) #pVALUE 0.01
plot(em1)
points(modis.t51,col='deepskyblue',cex=0.02)
points(modis.s51,col='darkgreen',cex=0.02)
 # plot(shape.utm,main='change in deter')
 # title('change in deter')
  points(sb1,col='black',cex=0.02) 
dev.off()

#ani.pause(ani.options(interval =2)) ## pause for a while ('interval')
 # ani.record()
}
, interval = 0.05, movie.name = "changec.gif", ani.width = 600, ani.height = 600)

extent(spmodiss)
extent(spmodist)
plot(shape.utm,add=TRUE)
plot(dc2,add=TRUE,col='green')
plot(mos,add=TRUE,col='pink')
plot(mo,add=TRUE,col='red')
plot(tif1)
plot(dft1.utm,add=TRUE,col='blue')
subset(StudentData, SchoolName=="Pine Tree Elementary" | Grade==3)