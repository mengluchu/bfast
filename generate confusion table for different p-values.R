
#load("pss4.Rdata")
load("pssww1.Rdata")
load("psswww2.Rdata")
load("pssww3.Rdata")
load("psswww4.Rdata")

load("waps1.Rdata")
load("waps2.Rdata")
load("waps3.Rdata")
load("waps4.Rdata")
####### functions implemented 
#bfastchangepoint()  bfast change to points
#deterpolygonpoint() deterpolygon to points
#show1() generate confusion matrix 
########


#install.packages("bfast", repos = "http://R-Forge.R-project.org")
# get confusion table

show1(ora150p005t)

#plot precition bfast detected point with different p value and aggregation
cusum150<-which(!is.na(t3darrbfamul),arr.ind=TRUE)
plot(cbind(cusum150[,1],cusum150[,2]))

which(!is.na(ora150p005t))
plot(bfastchangepoint(ora150p005t),col='blue',pch=1)
points(bfastchangepoint(a150p05t),col='red')
points(bfastchangepoint(a150p025t),col='brown4')
points(bfastchangepoint(a150p001t),col='skyblue2')
legend('bottomright',c('p-0.05 8-day','p-0.05 monthly','p-0.025 monthly','p-0.01 monthly'),col=c('blue','red','brown','skyblue'),pch=1)

# array index to MODIS coordinate 
#Returns the coords (MODIS synusoidal) of the center of the given pixel
#SR-ORG:6974

getxyMatrix <- function(colrowid.Matrix, pixelSize){
  
  x <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  y <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  corner.ul.x <- -20015109.354
  corner.ul.y <- 10007554.677
  x <- corner.ul.x + (pixelSize/2) + (colrowid.Matrix[,1] * pixelSize)
  y <- corner.ul.y - (pixelSize/2) - (colrowid.Matrix[,2] * pixelSize)
  cbind(x,y)
}
plot(bfastchangepoint(spacely40[,,5],x,y))
plot(bfastchangepoint(fspacely40[,,5],x,y),col='red')
extent(bcp2)
spfevi8<-bfastchangepoint(fevi8[,,1],x,y)

# plot changed points  
for(i in 1:20)
{
  
  if(length(which(!is.na(fspacelx42[,,i])))>0&&length(which(!is.na(f2spacely40[,,i])))>0) 
  {
    bcp<-bfastchangepoint(fspacelx42[,,i],x,y)
    bcpy<-bfastchangepoint(f2spacely40[,,i],x,y)
    #bcpa<-bfastchangepoint(spacelx150[,,i],x,y)
    #bcpya<-bfastchangepoint(spacely150[,,i],x,y)
    
    bcp2<-spTransform(bcp,CRS(usgspro))
    bcpy2<-spTransform(bcpy,CRS(usgspro))
    #bcp2a<-spTransform(bcpa,CRS(usgspro))
    #bcpy2a<-spTransform(bcpya,CRS(usgspro))
    
    plotRGB(crop(s6,epu4))
    #plotRGB(crop(s6,ex150usgs))
    points( bcp2,col=color2[i])
    points( bcpy2,col=color2[i],pch=2)
    #points( bcp2a,col=color2[i],pch=4)
    #points( bcpy2a,col=color2[i],pch=5)
  }
}
color2=rep(c('skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4'),1000)
x<-58929
y<-48219

#### get a patch of points ##
x<-c(59138:59180)
y<-c(48712:48752)
##150by 150
x<-c(58930:59079)
y<-c(48210:48359)

##10 by 10
x<-c(58930:58939)
y<-c(48210:48219)
##18 by 18
x<-c(58931:58948)
y<-c(48211:48228)
#1:28, 121:148
x<-c(58930:58958)
y<-c(48331:48358)
#121:148, 121:148
x<-c(59051:59078)
y<-c(48331:48358)
#spatialPolygons
#42*40
x<-c(58930:58971)
y<-c(48210:48249)

x<-c(58930:59079)
y<-c(48210:48211)
bfastchangepoint<- function(changearray,x,y)
{  
  change7<-which(!is.na(changearray ),arr.ind=TRUE) #0.05
  
  xct1<-change7[,1]+x[1]-1 # for the second 150 by 150 array
  xct2<-change7[,1] 
  #length(unique(alltct1))
  yct1<-change7[,2]+y[1]-1
  yct2<-change7[,2] 
  #alltct1<-change7[,3]
  
  dfallxyt<-as.data.frame(cbind(xct2,yct2))
  names(dfallxyt)<-c('x','y')
  
  
  coordinates(dfallxyt)<-~x+y #make the time value for searching
  
  #################
  xychangeinmodist<-as.data.frame(cbind(xct1,yct1))
  changeinmot0.5.1<-getxyMatrix(xychangeinmodist,231.6564)
  spmodist51<-SpatialPoints(coordinates(changeinmot0.5.1))
  proj4string(spmodist51)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
  
  modis.mt52<-spTransform(spmodist51,CRS("+proj=utm +zone=21 +south"))
  return(modis.mt52)
} 
bfastcusum<-bfastchangepoint(t3darrbfamul,x,y)

plot(bfastcusum)
pvaluepoint<- function(parray,x,y,xoff,yoff,pvalue)
{  
  change7<-which(parray<pvalue, arr.ind=TRUE) 
  
  xct1<-change7[,1]+x[1]-1+xoff # for the second 150 by 150 array
  xct2<-change7[,1] 
  #length(unique(alltct1))
  yct1<-change7[,2]+y[1]-1+yoff
  yct2<-change7[,2] 
  #alltct1<-change7[,3]
  
  dfallxyt<-as.data.frame(cbind(xct2,yct2))
  names(dfallxyt)<-c('x','y')
  
  
  coordinates(dfallxyt)<-~x+y #make the time value for searching
  
  #################
  xychangeinmodist<-as.data.frame(cbind(xct1,yct1))
  changeinmot0.5.1<-getxyMatrix(xychangeinmodist,231.6564)
  spmodist51<-SpatialPoints(coordinates(changeinmot0.5.1))
  proj4string(spmodist51)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
  
  modis.mt52<-spTransform(spmodist51,CRS("+proj=utm +zone=21 +south"))
  return(modis.mt52)
} 
wtssarcup0.05<-pvaluepoint(pssw1,x,y,1,1,0.05)#wavelet filter
wtssarcup0.005<-pvaluepoint(pssw1,x,y,1,1,0.005)
wtssarmop0.05<-pvaluepoint(pssw2,x,y,1,1,0.05)
wtssarmop0.025<-pvaluepoint(pssw2,x,y,1,1,0.025)
wtsCUSUM0.05<-pvaluepoint(pssw3,x,y,1,1,0.05) 
wtsMOSUM0.05<-pvaluepoint(pssw4,x,y,1,1,0.05)


h2tssarcup0.05<-pvaluepoint(pssw1,x,y,1,1,0.05) #h =0.2
h2tssarcup0.005<-pvaluepoint(pssw1,x,y,1,1,0.005)
h2tssarmop0.05<-pvaluepoint(pssw2,x,y,1,1,0.05)
h2tssarmop0.025<-pvaluepoint(pssw2,x,y,1,1,0.025)
h2tsCUSUM0.05<-pvaluepoint(pssw3,x,y,1,1,0.05) 
h2tsMOSUM0.05<-pvaluepoint(pssw4,x,y,1,1,0.05)
 
h5tssarcup0.05<-pvaluepoint(pssww1,x,y,1,1,0.05) #h=0.5
h5tssarcup0.005<-pvaluepoint(pssww1,x,y,1,1,0.005)
h5tssarmop0.05<-pvaluepoint(pssww2,x,y,1,1,0.05)
h5tssarmop0.025<-pvaluepoint(pssww2,x,y,1,1,0.025)
h5tsCUSUM0.05<-pvaluepoint(pssww3,x,y,1,1,0.05) 
h5tsMOSUM0.05<-pvaluepoint(pssww4,x,y,1,1,0.05)

h10tssarcup0.05<-pvaluepoint(pssww1,x,y,1,1,0.05) #h=1.0
h10tssarcup0.005<-pvaluepoint(pssww1,x,y,1,1,0.005)
h10tssarmop0.05<-pvaluepoint(psswww2,x,y,1,1,0.05)
h10tssarmop0.025<-pvaluepoint(psswww2,x,y,1,1,0.025)
h10tsCUSUM0.05<-pvaluepoint(pssww3,x,y,1,1,0.05) 
h10tsMOSUM0.05<-pvaluepoint(psswww4,x,y,1,1,0.05)

tssarcup0.05<-pvaluepoint(pss1,x,y,1,1,0.05) #trend and seasonality
tssarcup0.005<-pvaluepoint(pss1,x,y,1,1,0.005)
tssarmop0.05<-pvaluepoint(pss2,x,y,1,1,0.05)
tssarmop0.025<-pvaluepoint(pss2,x,y,1,1,0.025)
tsCUSUM0.05<-pvaluepoint(pss3,x,y,1,1,0.05) 
tsMOSUM0.05<-pvaluepoint(pss4,x,y,1,1,0.05)


sarp0.05<-pvaluepoint(p11,x,y,1,1,0.05)
sarp0.025<-pvaluepoint(p11,x,y,1,1,0.025)
sarpcusum0.05<-pvaluepoint(p1,x,y,1,1,0.05)
sarpcusum0.025<-pvaluepoint(p1,x,y,1,1,0.025)
sarpcusum0.01<-pvaluepoint(p1,x,y,1,1,0.01)
sarpcusum0.005<-pvaluepoint(p1,x,y,1,1,0.005)
MOSUM0.05<-pvaluepoint(p31,x,y,1,1,0.05)
CUSUM0.05<-pvaluepoint(p21,x,y,1,1,0.05)
plot(tssarp0.025,col='red',pch=12)
plot(rassp08671818)
library('rgdal')
deterpolygontopoint<-function(x,y,spatialPolygons) {
  
  x1<-rep(x,each=length(y))
  y1<-rep(y,length(x))
  
  xyd<-as.data.frame(cbind(x1,y1))
  xyc<-getxyMatrix(xyd,231.6564) # array coordinates to MODIS
  xyc<-as.data.frame(xyc)
  
  coordinates(xyc)<-~x+y
  SpatialPoints(xyc)
  as.character(round(unique(xyc@coords[,2])))
  
  proj4string(xyc)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
  xycs<-spTransform(xyc,CRS("+proj=utm +zone=21 +south  +ellps=WGS84 "))
  plot(spatialPolygons)
  proj4string(spatialPolygons)<-"+proj=utm +zone=21 +south  +ellps=WGS84 "
  deterpoints<-xycs[spatialPolygons,]
  points(xycs,col='yellow') 
  return(deterpoints)
}
plot(prodespoints00)
deterpoints<-deterpolygontopoint(x,y)
?spRbind
deterpointsc <- spChFIDs(deterpoints, as.character(paste('A',c(1:length(deterpoints)),sep='')))
prodesdeter0012<- spRbind(deterpoints,dif0012) #228 67 and 68
prodesdeter0012<-remove.duplicates(prodesdeter0012)
prodespoints00<-deterpolygontopoint(x,y,prodes2000)
prodespoints08<-deterpolygontopoint(x,y,sprodes2008)
prodespoints12<-deterpolygontopoint(x,y,sprodes2012)
prodespoints04<-deterpolygontopoint(x,y,sprodes2004)
save(prodesdeter0012,file='prodesdeter0012')

dif0012<-prodespoints00[which(is.na(over(prodespoints00,prodespoints12)))]
plot(dif0012)
load('dif0012.Rdata')
plot(prodespoints00[which(is.na(dif0008))])

plot(dif0008)
extent(prodes2000)
plot(deterpoints)
points(bfastcusum,col='red',pch=34)
load("spatialPolygons.R")
show1(a150p1s) # evi2 p-0.05
plot(xycs)
# extent(xycs)
#class       : Extent 
#xmin        : 367985.3 
#xmax        : 379129.5 
#ymin        : 8836861 
#ymax        : 8845887 
show1(sarp0.05)
arr<-t3darrbfamul
show1<-function(arr) {
  filtered1<-neighborsum8(arr)
  modis.mt52<-bfastchangepoint(filtered1,x,y)
  
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
  Specificity = TN/N
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
show3(sarp0.025)
show3(MOSUM0.05)
show3(CUSUM0.05)
show3<-function(changedpoints ) {
  
  modis.mt52<-changedpoints
  
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
  Specificity = TN/N
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
#prodes
plot(prodespoints00)

show2(t3darrbfamul,prodespoints00)
show2<-function(arr,spatialPolygons ) {
  filtered1<-neighborsum8(arr)
  modis.mt52<-bfastchangepoint(filtered1,x,y)
  
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
  Specificity = TN/N
  Accuracy= (TP[1:2]+TN)/ ALL
  p1<-paste('total prodes points:', b1, 'prodes points in bfast:', tl2, 'prodes points not in bfast', tl,'unique prodes point not in bfast:',b3,
            'unique prodes points in bfast',b4, 'unique bfast points not in prodes', b5, 'unique bfast in prodes',b6,'bfast points in prodes:',b7,
            'bfast points not in prodes',b8,'total bfast',b9, 'unique total bfast',b10,sep='        ')
  
  p2<-paste( 'TRUE POSITIVE',TP[1],TP[2],TP[3],
             'FALSE NEGATIVE',FN[1],FN[2],
             'FALSE POSITIVE', FP[1],FP[2],
             'POSITIVE',P[1],P[2],
             'NEGATIVE',N[1],N[2],
             'TRUE NEGATIVE',TN[1],TN[2],   
             'Precition', Precition[1],  Precition[2],
             'Sensitivity', Sensitivity[1],Sensitivity[2],
             'Specificity', Specificity[1] ,Specificity[2],
             'Accuracy',Accuracy[1], Accuracy[2],sep=', ')
  print(p1)
  print(p2)
  
}

tscusum<-show4(tsCUSUM0.05,dif0012)
show4(tsMOSUM0.05,dif0012)
show4(tssarcup0.05,dif0012)
show4(tssarcup0.025,dif0012)
show4(tssarmop0.05,dif0012)
show4(tssarmop0.025,dif0012)

show4(sarpcusum0.05,dif0012)

show4(sarpcusum0.01,dif0012)

#MOSUM0.05[which(!is.na(over(MOSUM0.05,prodespoints00)))]

#mask only forest area
MOSUM0.05<-MOSUM0.05[prodespoints00,]
CUSUM0.05<-CUSUM0.05[prodespoints00,]
sarp0.05<-sarp0.05[prodespoints00,]
sarp0.025<-sarp0.025[prodespoints00,]
sarpcusum0.05<-sarpcusum0.05[prodespoints00,]
sarpcusum0.005<-sarpcusum0.005[prodespoints00,]

tsMOSUM0.05<-tsMOSUM0.05[prodespoints00,]
tsCUSUM0.05<-tsCUSUM0.05[prodespoints00,]
tssarmop0.05<-tssarmop0.05[prodespoints00,]
tssarmop0.025<-tssarmop0.025[prodespoints00,]
tssarcup0.05<-tssarcup0.05[prodespoints00,]
tssarcup0.005<-tssarcup0.005[prodespoints00,]

h2tsMOSUM0.05<-h2tsMOSUM0.05[prodespoints00,]
h2tsCUSUM0.05<-h2tsCUSUM0.05[prodespoints00,]
h2tssarmop0.05<-h2tssarmop0.05[prodespoints00,]
h2tssarmop0.025<-h2tssarmop0.025[prodespoints00,]
h2tssarcup0.05<-h2tssarcup0.05[prodespoints00,]
h2tssarcup0.005<-h2tssarcup0.005[prodespoints00,]

h5tsMOSUM0.05<-h5tsMOSUM0.05[prodespoints00,]
h5tsCUSUM0.05<-h5tsCUSUM0.05[prodespoints00,]
h5tssarmop0.05<-h5tssarmop0.05[prodespoints00,]
h5tssarmop0.025<-h5tssarmop0.025[prodespoints00,]
h5tssarcup0.05<-h5tssarcup0.05[prodespoints00,]
h5tssarcup0.005<-h5tssarcup0.005[prodespoints00,]
par(mfrow=c(1,1))
prodespoints00[MOSUM0.05,]
plot(prodespoints00,col='orange')
points(tssarmop0.05,col='green',pch=4)
points(tssarmop0.025,col='skyblue4',pch=4)
points(sarp0.05 ,col='brown',pch=4)
points(tsMOSUM0.05,col='red2',pch=4)
points(MOSUM0.05,col='blue',pch=4)
points(prodesdeter0012,col='pink',pch=4)

svm<-show4(sarp0.05,MOSUM0.05)
plot(svm[,5:8])
points(MOSUM0.05[prodespoints00,],col='red')
groundtruth<-dif0012
groundtruth<-prodesdeter0012

m1<-show4(MOSUM0.05,groundtruth)
c1<-show4(CUSUM0.05,groundtruth)
sm11<-show4(sarp0.05,groundtruth)
sm21<-show4(sarp0.025,groundtruth)
cm11<-show4(sarpcusum0.05,groundtruth)
cm21<-show4(sarpcusum0.005,groundtruth)

m<-show4(tsMOSUM0.05,groundtruth)
c<-show4(tsCUSUM0.05,groundtruth)
sm1<-show4(tssarmop0.05,groundtruth)
sm2<-show4(tssarmop0.025,groundtruth)
cm1<-show4(tssarcup0.05,groundtruth)
cm2<-show4(tssarcup0.005,groundtruth)
##
hm<-show4(h2tsMOSUM0.05,groundtruth)
hc<-show4(h2tsCUSUM0.05,groundtruth)
hsm1<-show4(h2tssarmop0.05,groundtruth)
hsm2<-show4(h2tssarmop0.025,groundtruth)
hcm1<-show4(h2tssarcup0.05,groundtruth)
hcm2<-show4(h2tssarcup0.005,groundtruth)

hm5<-show4(h5tsMOSUM0.05,groundtruth)
hc5<-show4(h5tsCUSUM0.05,groundtruth)
hsm15<-show4(h5tssarmop0.05,groundtruth)
hsm25<-show4(h5tssarmop0.025,groundtruth)
hcm15<-show4(h5tssarcup0.05,groundtruth)
hcm25<-show4(h5tssarcup0.005,groundtruth)
allpd<-rbind(m1,c1,sm11,sm21,cm11,cm21) # trend  deseasonalized

all1<-rbind(m,c,sm1,sm2,cm1,cm2)       # trend and seasonality
all2<-rbind(hm,hc,hsm1,hsm2,hcm1,hcm2) #trend and seasonality h = 0.2
all5<-rbind(hm5,hc5,hsm15,hsm25,hcm15,hcm25)
all1
names1<-c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005")
#barplot(all1[,5:8],beside=TRUE, main="Precision, Sensitivity, Specificity and Accuracy",cex=0.8, cex.names = 0.7,    
#        legend.text=c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
#        args.legend = list(x = "topleft", bty = "n",cex=0.6))        

par(mfrow=c(1,2))
#rownames(all1)<-c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.025")
tex1='Trend + Seasonality'
tex2= 'Trend, Deseasonal'
tex1='h=0.15, ts'
tex2='h=0.5, ts'
tex2= 'Trend, Deseasonal'
barplot(all1[,1:4],beside=TRUE, main=tex1, cex.names = 0.7,    
        legend.text=c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
        args.legend = list(x = "topleft", bty = "n"))        

barplot(all5[,1:4],beside=TRUE, main=tex2, cex.names = 0.7,    
        legend.text=c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
        args.legend = list(x = "topleft", bty = "n"))        

barplot(all1[,5:8],beside=TRUE, main=tex1,cex=0.8, cex.names = 0.7,    
        legend.text=c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
        args.legend = list(x = "topleft", bty = "n",cex=0.6)) 

barplot(all5[,5:8],beside=TRUE, main=tex2,cex=0.8, cex.names = 0.7,    
        legend.text=c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
        args.legend = list(x = "topleft", bty = "n",cex=0.6))        


#rownames(all1)<-c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.025")

       
 

par(mfrow=c(2,2))

barplot(rbind(all1[,1],allpd[,1]),beside=TRUE, main="TRUE POSITIVE", cex.names = 0.7,    
        legend.text=c("PRODES ",'PRODES+DETER'),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        

axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.3)  

barplot(rbind(all1[,2],allpd[,2]),beside=TRUE, main="FALSE NEGATIVE", cex.names = 0.7,    
        legend.text=c("PRODES ",'PRODES+DETER'),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        

axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.3)  

barplot(rbind(all1[,3],allpd[,3]),beside=TRUE, main="FALSE POSITIVE", cex.names = 0.7,    
        legend.text=c("PRODES ",'PRODES+DETER'),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        

axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.3)  

barplot(rbind(all1[,4],allpd[,4]),beside=TRUE, main="TRUE NEGATIVE", cex.names = 1,    
        legend.text=c("PRODES ",'PRODES+DETER'),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        

axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  

barplot(rbind(all1[,5],allpd[,5]),beside=TRUE, main="PRECISION", cex.names = 0.7,    
        legend.text=c("PRODES ",'PRODES+DETER'),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        
axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  
barplot(rbind(all1[,6],allpd[,6]),beside=TRUE, main="SENSITIVITY", cex.names = 0.7,    
        legend.text=c("PRODES ",'PRODES+DETER'),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        
axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  
####

###
 
comp1<-all1
comp2<-all2
par(mfrow=c(2,2))

tex1='h=0.15 trend+seasonality'
tex2='h=0.5 trend+seasonlity'


tex1='deseasonlized, change in trend '
tex2='trend+seasonlity'
barplot(rbind(comp1[,1],comp2[,1]),beside=TRUE, main="TRUE POSITIVE", cex.names = 0.7,    
        legend.text=c(tex1 ,tex2),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        

axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  

barplot(rbind(comp1[,2],comp2[,2]),beside=TRUE, main="FALSE NEGATIVE", cex.names = 0.7,    
        legend.text=c(tex1 ,tex2),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        

axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  

barplot(rbind(comp1[,3],comp2[,3]),beside=TRUE, main="FALSE POSITIVE", cex.names = 0.7,    
        legend.text=c(tex1 ,tex2),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        

axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  

barplot(rbind(comp1[,4],comp2[,4]),beside=TRUE, main="TRUE NEGATIVE", cex.names = 0.7,    
        legend.text=c(tex1 ,tex2),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        

axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  

barplot(rbind(comp1[,5],comp2[,5]),beside=TRUE, main="PRECISION", cex.names = 0.7,    
        legend.text=c(tex1 ,tex2),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        
axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  

barplot(rbind(comp1[,6],comp2[,6]),beside=TRUE, main="SENSITIVITY", cex.names = 0.7,    
        legend.text=c(tex1 ,tex2),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        
axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  

barplot(rbind(comp1[,7],comp2[,7]),beside=TRUE, main="SPCIFICITY", cex.names = 0.7,    
        legend.text=c(tex1 ,tex2),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        
axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  

barplot(rbind(comp1[,8],comp2[,8]),beside=TRUE, main="ACCURACY", cex.names = 0.7,    
        legend.text=c(tex1 ,tex2),
        args.legend = list(x = "topright", bty = "n",cex=0.6),col=rainbow(2))        
axis(1,c(2,5,8,11,14,17),names1,lwd=1, cex.axis=0.5)  

prodesdeter0012<-remove.duplicates(prodesdeter0012)
plot(prodesdeter0012)
points(MOSUM0.05))

show4<-function(changedpoint,deterpoints ) {
  
  modis.mt52<-changedpoint
  
  result<-data.frame()
  deterpoints2<-c()
  
  
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
  
  TP = b7 
  FN=  b3
  FP=  b8
  P=TP+FN
  N=22500-P
  TN=N-FP
  ALL=22500
  
  Precision =TP/(TP+FP)
  Sensitivity = TP/P 
  Specificity = TN/N
  Accuracy= (TP+TN)/ ALL
  
  result<-cbind(TP,FN,FP,TN,Precision,Sensitivity,Specificity,Accuracy)
  names(result)<-c("TP","FN","FP","TN","Precision","Sensitivity","Specificity","Accuracy")
  
  p1<-paste('total prodes points:', b1,    'unique prodes point not in bfast:',b3,
            'unique prodes points in bfast',b4, 'unique bfast points not in prodes', b5, 'unique bfast in prodes',b6,
            'unique total bfast',b10, sep='        ')
  
  p2<-paste( 'TRUE POSITIVE',TP ,
             'FALSE NEGATIVE',FN ,
             'FALSE POSITIVE', FP ,
             'PRODES POSITIVE',P ,
             'PRODES NEGATIVE',N ,
             'TRUE NEGATIVE',TN   ,
             'Precision', Precision ,
             'Sensitivity', Sensitivity, 
             'Specificity', Specificity,
             'Accuracy',Accuracy, sep=', ')
  print(p1)
  print(p2)
  return(result)
  
}
################ auxilary functions ########################################
load('fevi8.Rdata')

a150p001t<-array(,c(150,150,167))
a150p001s<-array(,c(150,150,167))
a150p1t<-array(,c(150,150,167))
a150p1s<-array(,c(150,150,167))
date.trend2<-data.frame()
date.seasonal2<-data.frame()
############# filter changed points #########################
neighborsum8<-function(a)
{
  
  index1<-which(!is.na(a),arr.ind=TRUE)
  a[ which(is.na(a))]<-0
  
  for (i1 in (1:dim(index1)[1]))
  {
    
    i<-as.integer(index1[i1,1]) 
    
    j<-as.integer(index1[i1,2]) 
    t<-as.integer(index1[i1,3])
    
    ih<-i+1
    il<-i-1
    
    jh<-j+1
    jl<-j-1
    
    th<-t+1
    tl<-t-1
    
    if(ih>max(index1[,1]))
      ih=max(index1[,1])
    if(il<min(index1[,1]))
      il=min(index1[,1])
    
    if(jh>max(index1[,2]))
      jh=max(index1[,2])
    if(jl<min(index1[,2]))
      jl=min(index1[,2])
    
    if(th>max(index1[,3]))
      th=max(index1[,3])
    if(tl<min(index1[,3]))
      tl=min(index1[,3])
    
    
    asum<-sum(
      a[ih,j,t], a[il,j,t], a[i,jh,t], a[i,jl,t], a[ih,jh,t], a[il,jl,t], a[il,jh,t], a[ih,jl,t],
      a[ih,j,tl],a[il,j,tl],a[i,jh,tl],a[i,jl,tl],a[ih,jh,tl],a[il,jl,tl],a[il,jh,tl],a[ih,jl,tl],
      a[ih,j,th],a[il,j,th],a[i,jh,th],a[i,jl,th],a[ih,jh,th],a[il,jl,th],a[il,jh,th],a[ih,jl,th])
    if (asum==0)
    {
      a[i,j,t]<-NA
      
    }
  }
  a[ which(a==0)]<-NA
  return(a)   
}
################### run bfast in loop ##################################
for (i in 1:150)
{  
  for (j in 1:150)
  { 
    
    spt1<-fevi8[i,j,] 
    
    spt<-zoo( spt1,a1)
    monmean2 <- aggregate(spt, as.Date(as.yearmon(a1)), mean) #should aggregate in SciDB
    
    frequency(monmean2)<-12
    na.new <- function(x) ts(na.exclude(x), frequency = 12)
    
    stlmon<-stl(monmean2, na.action = na.new, s.window = "per")
    datamon <- ts(rowSums(stlmon$time.series)) 
    tsp(datamon) <- tsp(stlmon$time.series)
    
    fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=3,alpha=0.01) 
    
    print(i4)
    
    if(fitmon$nobp$Vt==FALSE)
    {
      
      date.trend2<-as.integer(fitmon$output[[1]]$Vt.bp)
      
      a150p001t[i,j, date.trend2]<-100
    }
    if(fitmon$nobp$Wt==FALSE)
    {
      date.seasonal2<-as.integer(fitmon$output[[1]]$Wt.bp)
      a150p001s[i,j, date.seasonal2]<-500
      
    }
    
    i4=i4+1
    
  }
}
#save( a150p05s,file='a150p05s.RData')
#save( a150p05t,file='a150p05t.RData')
#save( a150p025s,file='a150p025s.RData')
#save( a150p025t,file='a150p025t.RData')
#save( a150p001s,file='a150p001s.RData')
#save( a150p001t,file='a150p001t.RData')
getwd()
