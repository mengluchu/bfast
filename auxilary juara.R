#### get a patch of points ##
x<-c(59138:59180)
y<-c(48712:48752)
##150by 150
x<-c(58930:59079)
y<-c(48210:48369)
#spatialPolygons
#42*40

 
#c(59138:59180)
#c(48712:48753)
#x<-c(round(ex42@xmin):round(ex42@xmax))
#y<-c(round(ex42@ymin):round(ex42@ymax))

## deter polygons to deter points
x<-c(59139:59180)
y<-c(48713:48752)
deterpolygontopoint<-function(x,y) {
x1<-rep(x,each=length(y))
y1<-rep(y,length(x))
 
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
plot(xycs,add=TRUE) 
return(deterpoints)
}
deterpoints<-deterpolygontopoint(x,y)
deterpoints
modis.mt52<-bfastchangepoint(test4nds)
load('test2ndt.Rdata')
load('test4nds.Rdata')
load('test5nds.Rdata')
load('test5ndt.Rdata')
 
#p.Vt <- sctest(efp(Vt ~ ti, h=h, type= "OLS-MOSUM"))
install_github("bfast-1","mengluchu")
bfast 
show1(test3nds)
 

show1<-function(arr ) {
filtered1<-neighborsum8(arr)
modis.mt52<-bfastchangepoint(filtered1)
 
deterpoints2<-c()
tl=0
for ( m in  1:155)
{
 
 # tl<-tl+(length(unique(xycs[spatialPolygons[m],][modis.mt52,]@coords))/2
  tl<-tl+(length(xycs[spatialPolygons[m],][is.na(over(xycs[spatialPolygons[m],],modis.mt52))]@coords))/2 
}
tl2=0
for ( m in  1:155)
{
  
  # tl<-tl+(length(unique(xycs[spatialPolygons[m],][modis.mt52,]@coords))/2
  tl2<-tl2+(length(xycs[spatialPolygons[m],][!is.na(over(xycs[spatialPolygons[m],],modis.mt52))]@coords))/2 
}
 
b9<-length(modis.mt52)
b10<-length(unique(modis.mt52@coords))/2
 
b1<-length( deterpoints) #1898
b2<-length(which(!is.na(over(deterpoints,modis.mt52)))) #997
b3<-length(unique(deterpoints[is.na(over(deterpoints,modis.mt52))]@coords))/2  

b4<-length(unique(deterpoints[!is.na(over(deterpoints,modis.mt52))]@coords))/2  
b5<-length(unique(modis.mt52[is.na(over(modis.mt52,deterpoints))]@coords))/2 #5381 /4061
b6<-length(unique(modis.mt52[!is.na(over(modis.mt52,deterpoints))]@coords))/2  
b7<-length(which(!is.na(over(modis.mt52,deterpoints)))) #1400/997
b8<-length(which(is.na(over(modis.mt52,deterpoints))))



TP = c(b7,b6, tl2)
FN=c(tl, b3)
FP=c( b8, b5)
P=TP[1:2]+FN
N=22500-P
TN=N-FP
ALL=22500
 
Precition =TP[1:2] /(TP[1:2]+FP)
Sensitivity = TP[1:2]/P 
Specifcity = TN/N
Accuracy= (TP[1:2]+TN)/ ALL
p1<-paste('total deter points:', b1, 'deterpoints in bfast:', tl2, 'deterpoints not in bfast', tl,'unique deter point not in bfast:',b3,
      'unique deter points in bfast',b4, 'unique bfast points not in deter', b5, 'unique bfast in deter',b6,'bfast points in deter:',b7,
      'bfast points not in deter',b8,'total bfast',b9, 'unique total bfast',b10,sep='        ')

p2<-paste( 'TP',TP[1],TP[2],TP[3],
       'FN',FN[1],FN[2],
       'FP', FP[1],FP[2],
       'P',P[1],P[2],
       'N',N[1],N[2],
       'TN',TN[1],TN[2],   
       'Precition', Precition[1],  Precition[2],
       'Sensitivity', Sensitivity[1],Sensitivity[2],
      ' Specificity', Specificity[1] ,Specificity[2],
      'accuracy',Accuracy[1], Accuracy[2],sep=',  ')
print(p1)
print(p2)

}
modis.mt52
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
