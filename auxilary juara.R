#### get a patch of points ##
x<-c(59138:59180)
y<-c(48712:48752)
##150by 150
x<-c(58930:59079)
y<-c(48210:48369)
#spatialPolygons
#42*40
x<-c(59139:59180)
y<-c(48713:48752)
cosp3<-cosp2[-1,]
cosp3[,1]+59138
c(59138:59180)
c(48712:48753)
#x<-c(round(ex42@xmin):round(ex42@xmax))
#y<-c(round(ex42@ymin):round(ex42@ymax))

x1<-rep(x,each=length(y))
y1<-rep(y,length(x))
x
xyd<-as.data.frame(cbind(x1,y1))
xyc<-getxyMatrix(xyd,231.6564)
xyc<-as.data.frame(xyc)
coordinates(xyc)<-~x+y
SpatialPoints(xyc)
as.character(round(unique(xyc@coords[,2])))
proj4string(xyc)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
xycs<-spTransform(xyc,CRS("+proj=utm +zone=21 +south  +ellps=WGS84 "))
plot(spatialPolygons)
deterpoints<-xycs[spatialPolygons,]
 
 
length(spatialPolygons)
deterpoints2<-c()
tl=0
for ( m in  1:155)
{
 
 # tl<-tl+(length(unique(xycs[spatialPolygons[m],][modis.mt52,]@coords))/2
  tl<-tl+(length(xycs[spatialPolygons[m],][is.na(over(xycs[spatialPolygons[m],],modis.mt52))]@coords))/2 
}
 

#in bfast 1175
length(over(deterpoints,modis.mt52)) #1898
length(which(!is.na(over(deterpoints,modis.mt52)))) #997
length(unique(deterpoints[is.na(over(deterpoints,modis.mt52))]@coords))/2 #901
length(unique(modis.mt52[is.na(over(modis.mt52,deterpoints))]@coords))/2 #5381 /4061
length(unique(modis.mt52[!is.na(over(modis.mt52,deterpoints))]@coords))/2 #1400/997
#number: 1898 tl:4250 the  way to account for the replication of the polygons would be to do the intercection one polygon by one polygon 
spatialPolygons[xycs,]
length(spatialPolygons[xycs,])
plot(xycs,add=TRUE)
points(deterpoints,add=TRUE,pch=2,cex=0.3)
plot(modis.mt52,col='red',pch=2,cex=0.2)
,add=TRUE)
plot(coordinates(xyd),col='green')
plot(xyc,pch=1)
cr520<-getcrMatrix(as.data.frame(xyc),231.6564)

match(cr520[,2],xyd[,2])
points(cr520,add=TRUE,col='red',pch='2',cex=0.5)
slot(xycs,'coords')[2,]
####
#time series 
w=1/46
#harmonic
co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2) #:fit with stl seasonality: p-value 0.9687
co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3) # fit with stl seasonality the best: p-value 0.9786
#co4 <- cos(2*pi*tl*w*4);si4 <- sin(2*pi*tl*w*4) # higher order will not improve:p-value 0.9779

season1<-lm(coredata(spt)~co+co2+co3+si+si2+si3) #almost perfect linear relationship with stl: season: R-squared:  0.9786
deseasonal<-coredata(spt)-fitted(season1)
plot(season1$residuals,typ='l')
plot(deseasonal,typ='l')
acf(coredata((lm(deseasonal~a1))$residuals),lag=100,main='8 days correlation 2 (deseasonalized and detrended 1,1)')
tsdiag(arima(coredata((lm(deseasonal~a1))$residuals),order=c(0,0,0)))

N = length(deseasonal)
I = abs(fft(deseasonal))/sqrt(N)^2
P = (4/N)*I # Scale periodogram
f = (0:floor(N/2))/N

plot(f, I[1:((N/2)+1)], type="o", xlab="frequency", ylab="")
