library('rgdal')
library('raster')

load('prodespoints00.Rdata')
load('pdd.Rdata')
groundtruth<-pdd
load('tssarar1.Rdata') # corrected st model with original array cusum
load('tssarar2.Rdata')#mosum
load('tssarar3.Rdata')#sar cusum
load('tssarar4.Rdata')#sar mosum
load('tssarar5.Rdata')#ar cusum
load('tssarar6.Rdata')#ar mosum
 
color2=rep(c('skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4'),1000)

##150by 150
x<-c(58930:59079)
y<-c(48210:48359)

getxyMatrix <- function(colrowid.Matrix, pixelSize){
  
  x <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  y <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  corner.ul.x <- -20015109.354
  corner.ul.y <- 10007554.677
  x <- corner.ul.x + (pixelSize/2) + (colrowid.Matrix[,1] * pixelSize)
  y <- corner.ul.y - (pixelSize/2) - (colrowid.Matrix[,2] * pixelSize)
  cbind(x,y)
}
 
#spfevi8<-bfastchangepoint(fevi8[,,1],x,y)
 
bfastchangepoint<- function(changearray,x=c(58930:59079),y=c(48210:48359)) # map the changes from the array that stores the changes. multiple changes are mapped as one change point
{  
  change7<-which(!is.na(changearray ),arr.ind=TRUE) #0.05
  
  xct1<-change7[,1]+x[1]-1 # for the second 150 by 150 array
  xct2<-change7[,1] 

  yct1<-change7[,2]+y[1]-1
  yct2<-change7[,2] 

  dfallxyt<-as.data.frame(cbind(xct2,yct2))
  names(dfallxyt)<-c('x','y')
  
  
  coordinates(dfallxyt)<-~x+y #make the time value for searching

  xychangeinmodist<-as.data.frame(cbind(xct1,yct1))
  changeinmot0.5.1<-getxyMatrix(xychangeinmodist,231.6564)
  spmodist51<-SpatialPoints(coordinates(changeinmot0.5.1))
  proj4string(spmodist51)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
  
  modis.mt52<-spTransform(spmodist51,CRS("+proj=utm +zone=21 +south"))
  return(modis.mt52)
} 
 
 
pvaluepoint<- function(parray,x=c(58930:59079),y=c(48210:48359),xoff=1,yoff=1,pvalue=0.05) #map changes from an array of p-values
{  
  
  change7<-which(parray<pvalue, arr.ind=TRUE) 
  
  xct1<-change7[,1]+x[1]-1+xoff # for the second 150 by 150 array
  xct2<-change7[,1] 
 
  yct1<-change7[,2]+y[1]-1+yoff
  yct2<-change7[,2] 
 
  dfallxyt<-as.data.frame(cbind(xct2,yct2))
  names(dfallxyt)<-c('x','y')
  
  coordinates(dfallxyt)<-~x+y #make the time value for searching
 
  xychangeinmodist<-as.data.frame(cbind(xct1,yct1))
  changeinmot0.5.1<-getxyMatrix(xychangeinmodist,231.6564)
  spmodist51<-SpatialPoints(coordinates(changeinmot0.5.1))
  
  proj4string(spmodist51)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
 
  modis.mt52<-spTransform(spmodist51,CRS("+proj=utm +zone=21 +south"))
  # modis.mt52<-spTransform(modis.mt5,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  return(modis.mt52)
} 

show5<-function(changedpoint,deterpoints ) {
  
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
  N=19167-P
  TN=N-FP
  ALL=19167
  
  PRODUCER =TP/(TP+FP+FN)
  
  result<-cbind(TP,FN,FP,TN,PRODUCER)
  names(result)<-c("TP","FN","FP","TN","producer's accuracy")
  
  p1<-paste('total prodes points:', b1,    'unique prodes point not in bfast:',b3,
            'unique prodes points in bfast',b4, 'unique bfast points not in prodes', b5, 'unique bfast in prodes',b6,
            'unique total bfast',b10, sep='        ')
  
  p2<-paste( 'TRUE POSITIVE',TP ,
             'FALSE NEGATIVE',FN ,
             'FALSE POSITIVE', FP ,
             'PRODES POSITIVE',P ,
             'PRODES NEGATIVE',N ,
             'TRUE NEGATIVE',TN   ,
             'Precision', PRODUCER ,
             sep=', ')
  print(p1)
  print(p2)
  return(result)
  
}
generatecm<-function(groundtruth, pvaluemx,pv=0.05,x=c(58930:59079),y=c(48210:48359)) #  put functions together and generate confusion matrix (filtered with a mask)
{
  pvpoint<-pvaluepoint(pvaluemx,x,y,1,1,pv)
  fpvpoint<-pvpoint[prodespoints00,]
  showre<-show5(fpvpoint,groundtruth)
  return(showre)
}
generatecmbfast<-function(reference.sppoints=pdd, result.array,mask=prodespoints00,pv=0.05,x=c(58930:59079),y=c(48210:48359)) #  put functions together and generate confusion matrix (filtered with a mask)
{
  pvpoint<-bfastchangepoint(result.array,x,y )
  fpvpoint<-pvpoint[mask,]
  showre<-show5(fpvpoint,reference.sppoints)
  return(showre)
}

bfast1cm<-generatecmbfast(pdd,result.array=bfast1 )
edivi1cm<-generatecmbfast(pdd,result.array=edivisive1)

generatp<-function(groundtruth,pvaluemx,pv=0.05,x=c(58930:59079),y=c(48210:48359))  # filtered points
{
  
  pvpoint<-pvaluepoint(pvaluemx,x,y,1,1,pv)
  
  fpvpoint<-pvpoint[prodespoints00,]
 
  return(fpvpoint)
}
 
  
ptssarar1<-generatp(tssarar1,0.05)
ptssarar11<-generatp(tssarar1,0.005)
ptssarar2<-generatp(tssarar2,0.05)
ptssarar22<-generatp(tssarar2,0.025)
ptssarar3<-generatp(tssarar3,0.05)
ptssarar4<-generatp(tssarar4,0.05)
ptssarar5<-generatp(tssarar5,0.05)
ptssarar6<-generatp(tssarar6,0.05)

ttssarar1<-generatecm(tssarar1,0.05)
ttssarar11<-generatecm(tssarar1,0.2)
ttssarar2<-generatecm(tssarar2,0.05)
ttssarar22<-generatecm(tssarar2,0.1)
ttssarar3<-generatecm(tssarar3,0.05)
ttssarar4<-generatecm(tssarar4,0.05)
ttssarar5<-generatecm(tssarar5,0.05)
ttssarar6<-generatecm(tssarar6,0.05)

cts<-rbind(ttssarar4,ttssarar6,ttssarar2,ttssarar22，ttssarar3，ttssarar5，ttssarar1，ttssarar11)

cts2<-rbind(bfast1cm,edivi1cm)

barplot(cts[,1:4]/22500,beside=TRUE, main=tex1,    
        legend.text=c("bfast","e.divi"),
        
        cex.main=1,
        font.main=1,
        col=bpy.colors(2),
        #=rainbow(8,s = 1, v = 1, start = 0.05, end = 0.35),
        #legend.text=c("MOSUM","CUSUM","AR(1) MOSUM 0.05","AR(1) CUSUM 0.05","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
        args.legend = list(x = "topleft", bty = "n",cex=0.8))        


save(cts,file='cts.Rdata')
plot(1)
names1<-c("OLS-MOSUM p-value: 0.05","AR(1) OLS-MOSUM pvalue: 0.05","SAR OLS-MOSUM p-value: 0.05", "SAR OLS-MOSUM p-value: 0.1","OLS-CUSUM p-value: 0.05","AR(1) OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.2")
tex1='Confusion Matrix of Different Methods'
 
jpeg('Pontius Producer Accuracy ts.jpg', height=6, width=10, res=400,unit="in")

par(mfrow=c(1,2))
barplot(cts[,1:4]/22500,beside=TRUE, main=tex1,    
        legend.text=names1,
    
        cex.main=1,
        font.main=1,
        col=bpy.colors(8),
        #=rainbow(8,s = 1, v = 1, start = 0.05, end = 0.35),
        #legend.text=c("MOSUM","CUSUM","AR(1) MOSUM 0.05","AR(1) CUSUM 0.05","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
        args.legend = list(x = "topleft", bty = "n",cex=0.8))        

barplot(cts[,5],beside=TRUE, main='Pontius Producer\'s Accuracy',     
        col=bpy.colors(8),
        font.main=1,
        cex.main=1,
        #legend.text=c("MOSUM","CUSUM","SAR MOSUM 0.05", "SAR MOSUM 0.025","SAR CUSUM 0.05","SAR CUSUM 0.005"),
        #args.legend = list(x = "topleft", bty = "n",cex=0.6)) 
)
dev.off()
load("splamda.Rdata")

cor(as.vector(splamda),as.vector(tssplamda))
jpeg("s and t dependence2.jpg",width=6,height=3.5,res=1200,unit="in")
par(mfrow=c(1,2))
v.hist <- hist(splamda, plot=FALSE,breaks=20)
v.hist$counts <- v.hist$counts/sum(v.hist$counts)
plot(v.hist,xlab=expression(rho),ylab="probability",
 
     cex.main=1 ,
     xlim=c(0,1),ylim=c(0,1),main= expression(paste("Spatial Dependence  ", rho)) )
#,main= expression(paste(rho))
v.hist <- hist(tphi, plot=FALSE,breaks=20)
v.hist$counts <- v.hist$counts/sum(v.hist$counts)
plot(v.hist,xlab=expression(phi),ylab="probability",     cex.main=1 ,
     xlim=c(0,1),ylim=c(0,1),main= expression(paste("AR(1) Coefficient  ", phi)))

dev.off()



generatemap<-function (pva, groundtruth,titlename='map') 
{
  groundtruth<-spTransform(groundtruth,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  prodespoints00<-spTransform(prodespoints00,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  
  plot(prodespoints00,col='darkgreen',pch=16)#forest map
  pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  
  points(pva[which(is.na(over(pva,groundtruth)))],pch=16,col='orange') #false positive
  #faslse positive 
  points(groundtruth[which(is.na(over(groundtruth,pva)))],pch=16,col='skyblue')
  #false negative
  points(pva[groundtruth,],col='red',pch=16)# true positive
  #legend('bottomright',c('True  Positive', 'False Positive','False Negative','True  Negative'), col=c('red','orange','skyblue','darkgreen' ),pch=19,cex=1)
  title(titlename)
}


jpeg('tssarar.jpg', height=12, width=18, res=400,unit="in")
par(mfrow=c(2,2))

generatemap(ptssarar3,groundtruth, 'Change Detection with  OLS-CUSUM')
generatemap(ptssarar4,groundtruth, 'Change Detection with  OLS-MOSUM')
generatemap(ptssarar2,groundtruth, 'Change Detection with SAR OLS-MOSUM')
generatemap(ptssarar1,groundtruth, 'Change Detection with SAR OLS-CUSUM')

dev.off()

library(ggplot2)
library(gridExtra)

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}
 
#groundtruth<-spTransform(groundtruth,CRS("+proj=longlat +ellps=aust_SA +no_defs"))  
#prodespoints00<-spTransform(prodespoints00,CRS('+proj=longlat +ellps=aust_SA +no_defs'))
#pva<-spTransform(ptssarar1,CRS('+proj=longlat +ellps=aust_SA +no_defs'))
generatemapGGPLOT<-function (pva1,pva2,pva3,pva4, groundtruth,titlename='map') 
{
  pva<-pva1
  groundtruth<-spTransform(groundtruth,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  prodespoints00<-spTransform(prodespoints00,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
 
  FP<-pva[which(is.na(over(pva,groundtruth)))]
  FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
  TP<-pva[groundtruth,]
  PRODESNE<-prodespoints00[ which(is.na(over(prodespoints00,groundtruth))),]
  BFASTNE<-prodespoints00[ which(is.na(over(prodespoints00,pva))),]
  TN<-PRODESNE[BFASTNE,]

  FPx<-data.frame(FP)[,1]
  FPy<-data.frame(FP)[,2]
  FNx<-data.frame(FN)[,1]
  FNy<-data.frame(FN)[,2]
  TPx<-data.frame(TP)[,1]
  TPy<-data.frame(TP)[,2]
  ALLx<<-data.frame(prodespoints00)[,1]
  ALLy<-data.frame(prodespoints00)[,2]

gg1<-
  qplot()+ 
  geom_point( aes(x=ALLx, y=ALLy, color="TN"))+
  geom_point( aes(x=FPx, y=FPy, color="FP"))+
  geom_point( aes(x=FNx, y=FNy, color="FN" ))+
  geom_point( aes(x=TPx, y=TPy, color="TP" ))+
  xlab('X')+ylab('Y') + 
  ggtitle("OLS-CUSUM")+
  scale_colour_manual(values=c("TN"="darkgreen", "FP"="orange", "FN"="skyblue", "TP"="red"), 
                    name="CONFUSION MATRIX",      
                    breaks=c("TN","FP","FN","TP"),
                    labels=c("TRUE NEGATIVE", "FALSE POSITIVE", "FALSE NEGATIVE","TRUE POSITIVE"))+ 
theme(legend.position = "none")
############
pva<-pva2
groundtruth<-spTransform(groundtruth,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
prodespoints00<-spTransform(prodespoints00,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))

FP<-pva[which(is.na(over(pva,groundtruth)))]
FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
TP<-pva[groundtruth,]

FPx2<-data.frame(FP)[,1]
FPy2<-data.frame(FP)[,2]
FNx2<-data.frame(FN)[,1]
FNy2<-data.frame(FN)[,2]
TPx2<-data.frame(TP)[,1]
TPy2<-data.frame(TP)[,2]
ALLx2<<-data.frame(prodespoints00)[,1]
ALLy2<-data.frame(prodespoints00)[,2]

gg2<-
  qplot()+ 
  geom_point( aes(x=ALLx2, y=ALLy2, color="TN"))+
  geom_point( aes(x=FPx2, y=FPy2, color="FP"))+
  geom_point( aes(x=FNx2, y=FNy2, color="FN" ))+
  geom_point( aes(x=TPx2, y=TPy2, color="TP" ))+
  xlab('X')+ylab('Y') + 
  ggtitle("OLS-MOSUM")+
  scale_colour_manual(values=c("TN"="darkgreen", "FP"="orange", "FN"="skyblue", "TP"="red"), 
                      name="CONFUSION MATRIX",  
                      breaks=c("TN","FP","FN","TP"),
                      labels=c("TRUE NEGATIVE", "FALSE POSITIVE", "FALSE NEGATIVE","TRUE POSITIVE"))+ 
  theme(legend.position = "none")
############
pva<-pva3

pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))

FP<-pva[which(is.na(over(pva,groundtruth)))]
FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
TP<-pva[groundtruth,]

FPx3<-data.frame(FP)[,1]
FPy3<-data.frame(FP)[,2]
FNx3<-data.frame(FN)[,1]
FNy3<-data.frame(FN)[,2]
TPx3<-data.frame(TP)[,1]
TPy3<-data.frame(TP)[,2]
ALLx3<<-data.frame(prodespoints00)[,1]
ALLy3<-data.frame(prodespoints00)[,2]

gg3<-
  qplot()+ 
  geom_point( aes(x=ALLx3, y=ALLy3, color="TN"))+
  geom_point( aes(x=FPx3, y=FPy3, color="FP"))+
  geom_point( aes(x=FNx3, y=FNy3, color="FN" ))+
  geom_point( aes(x=TPx3, y=TPy3, color="TP" ))+
  xlab('X')+ylab('Y') + 
  ggtitle("SAR OLS-CUSUM")+
  scale_colour_manual(values=c("TN"="darkgreen", "FP"="orange", "FN"="skyblue", "TP"="red"), 
                      name="CONFUSION MATRIX",
                      breaks=c("TN","FP","FN","TP"),
                      labels=c("TRUE NEGATIVE", "FALSE POSITIVE", "TRUE NEGATIVE","TRUE POSITIVE"))+ 
  theme(legend.position = "none")
#############
pva<-pva4
 
pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))

FP<-pva[which(is.na(over(pva,groundtruth)))]
FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
TP<-pva[groundtruth,]

FPx4<-data.frame(FP)[,1]
FPy4<-data.frame(FP)[,2]
FNx4<-data.frame(FN)[,1]
FNy4<-data.frame(FN)[,2]
TPx4<-data.frame(TP)[,1]
TPy4<-data.frame(TP)[,2]
ALLx4<<-data.frame(prodespoints00)[,1]
ALLy4<-data.frame(prodespoints00)[,2]

gg4<-
  qplot()+ 
  geom_point( aes(x=ALLx4, y=ALLy4, color="TN"))+
  geom_point( aes(x=FPx4, y=FPy4, color="FP"))+
  geom_point( aes(x=FNx4, y=FNy4, color="FN" ))+
  geom_point( aes(x=TPx4, y=TPy4, color="TP" ))+
  xlab('X')+ylab('Y') + 
  ggtitle("SAR OLS-MOSUM")+
  scale_colour_manual(values=c("TN"="darkgreen", "FP"="orange", "FN"="skyblue", "TP"="red"), 
                      name="CONFUSION MATRIX",  
                      breaks=c("TN","FP","FN","TP"),
                      labels=c("TRUE NEGATIVE", "FALSE POSITIVE", "FALSE NEGATIVE","TRUE POSITIVE"))+ 
  theme(legend.position = "none")


grid_arrange_shared_legend(gg1,gg2,gg3,gg4)
}



generatemapGGRASTER<-function (pva1,pva2,pva3,pva4, groundtruth,titlename='map') 
{
  pva<-pva1
  groundtruth<-spTransform(groundtruth,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  prodespoints00<-spTransform(prodespoints00,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  
  FP<-pva[which(is.na(over(pva,groundtruth)))]
  FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
  TP<-pva[groundtruth,]
  PRODESNE<-prodespoints00[ which(is.na(over(prodespoints00,groundtruth))),]
  BFASTNE<-prodespoints00[ which(is.na(over(prodespoints00,pva))),]
  TN<-PRODESNE[BFASTNE,]
  
  TPV<-rep(1,length(data.frame(TP)[,1]))
  TNV<-rep(2,length(data.frame(TN)[,1]))
  FPV<-rep(3,length(data.frame(FP)[,1]))
  FNV<-rep(4,length(data.frame(FN)[,1]))
  #plot(raster(TPdf))
  TPdf<-SpatialPointsDataFrame( coordinates(data.frame(TP)),data=data.frame(TPV))
  TNdf<-SpatialPointsDataFrame( coordinates(data.frame(TN)),data=data.frame(TNV))
  FPdf<-SpatialPointsDataFrame( coordinates(data.frame(FP)),data=data.frame(FPV))
  FNdf<-SpatialPointsDataFrame( coordinates(data.frame(FN)),data=data.frame(FNV))
  a<-data.frame(TPdf)
  b<-data.frame(FPdf)
  c<-data.frame(FNdf)
  d<-data.frame(TNdf)
  
  names(c)<-names(a)
  names(b)<-names(a)
  names(d)<-names(a)
  abcd<-rbind(a,b,c,d)
  
  coordinates(abcd)<-~x+y
  gridded(abcd)<-TRUE
 
  r1<-raster(abcd)
 gg1<- gplot(r1)+ geom_tile(aes(fill = factor(value))) +
    scale_fill_manual(values=c("2"="seagreen4", "3"="palegoldenrod", "4"="skyblue", "1"="palevioletred1"), 
                      na.value="white",
                      name="CONFUSION MATRIX",  
                      breaks=c("2","3","4","1"),         
                      labels=c("TRUE NEGATIVE", "FALSE POSITIVE", "FALSE NEGATIVE","TRUE POSITIVE"))+
    coord_equal()+ 
    ggtitle("OLS-CUSUM")+    theme(legend.position = "none")
 
  ############
  pva<-pva2
  groundtruth<-spTransform(groundtruth,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  prodespoints00<-spTransform(prodespoints00,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  FP<-pva[which(is.na(over(pva,groundtruth)))]
  FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
  TP<-pva[groundtruth,]
  PRODESNE<-prodespoints00[ which(is.na(over(prodespoints00,groundtruth))),]
  BFASTNE<-prodespoints00[ which(is.na(over(prodespoints00,pva))),]
  TN<-PRODESNE[BFASTNE,]
  
  TPV<-rep(1,length(data.frame(TP)[,1]))
  TNV<-rep(2,length(data.frame(TN)[,1]))
  FPV<-rep(3,length(data.frame(FP)[,1]))
  FNV<-rep(4,length(data.frame(FN)[,1]))
  #plot(raster(TPdf))
  TPdf<-SpatialPointsDataFrame( coordinates(data.frame(TP)),data=data.frame(TPV))
  TNdf<-SpatialPointsDataFrame( coordinates(data.frame(TN)),data=data.frame(TNV))
  FPdf<-SpatialPointsDataFrame( coordinates(data.frame(FP)),data=data.frame(FPV))
  FNdf<-SpatialPointsDataFrame( coordinates(data.frame(FN)),data=data.frame(FNV))
  a<-data.frame(TPdf)
  b<-data.frame(FPdf)
  c<-data.frame(FNdf)
  d<-data.frame(TNdf)
  
  names(c)<-names(a)
  names(b)<-names(a)
  names(d)<-names(a)
  abcd<-rbind(a,b,c,d)
  
  coordinates(abcd)<-~x+y
  gridded(abcd)<-TRUE
  
  r2<-raster(abcd)
  gg2<-gplot(r2)+ geom_tile(aes(fill = factor(value))) +
    scale_fill_manual(values=c("2"="seagreen4", "3"="palegoldenrod", "4"="skyblue", "1"="palevioletred1"), 
                      na.value="white",
                      name="CONFUSION MATRIX",  
                      breaks=c("1","2","3","4"),         
                      labels=c("TRUE NEGATIVE", "FALSE POSITIVE", "FALSE NEGATIVE","TRUE POSITIVE"))+
    coord_equal()+ 
    ggtitle("OLS-MOSUM")+    theme(legend.position = "none")
  
   
  ############
  pva<-pva3
  
  pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
  
  FP<-pva[which(is.na(over(pva,groundtruth)))]
  FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
  TP<-pva[groundtruth,]
  PRODESNE<-prodespoints00[ which(is.na(over(prodespoints00,groundtruth))),]
  BFASTNE<-prodespoints00[ which(is.na(over(prodespoints00,pva))),]
  TN<-PRODESNE[BFASTNE,]
  
  TPV<-rep(1,length(data.frame(TP)[,1]))
  TNV<-rep(2,length(data.frame(TN)[,1]))
  FPV<-rep(3,length(data.frame(FP)[,1]))
  FNV<-rep(4,length(data.frame(FN)[,1]))
  #plot(raster(TPdf))
  TPdf<-SpatialPointsDataFrame( coordinates(data.frame(TP)),data=data.frame(TPV))
  TNdf<-SpatialPointsDataFrame( coordinates(data.frame(TN)),data=data.frame(TNV))
  FPdf<-SpatialPointsDataFrame( coordinates(data.frame(FP)),data=data.frame(FPV))
  FNdf<-SpatialPointsDataFrame( coordinates(data.frame(FN)),data=data.frame(FNV))
  a<-data.frame(TPdf)
  b<-data.frame(FPdf)
  c<-data.frame(FNdf)
  d<-data.frame(TNdf)
  
  names(c)<-names(a)
  names(b)<-names(a)
  names(d)<-names(a)
  abcd<-rbind(a,b,c,d)
  
  coordinates(abcd)<-~x+y
  gridded(abcd)<-TRUE
  
  r3<-raster(abcd)
 gg3<- gplot(r3)+ geom_tile(aes(fill = factor(value))) +
    scale_fill_manual(values=c("2"="seagreen4", "3"="palegoldenrod", "4"="skyblue", "1"="palevioletred1"), 
                      na.value="white",
                      name="CONFUSION MATRIX",  
                      breaks=c("1","2","3","4"),         
                      labels=c("TRUE NEGATIVE", "FALSE POSITIVE", "FALSE NEGATIVE","TRUE POSITIVE"))+
    coord_equal()+ 
    ggtitle("SAR OLS-CUSUM")+    theme(legend.position = "none")
  
     
  #############
  pva<-pva4
 
 pva<-spTransform(pva,CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'))
 
  FP<-pva[which(is.na(over(pva,groundtruth)))]
  FN<-groundtruth[which(is.na(over(groundtruth,pva)))]
  TP<-pva[groundtruth,]
  PRODESNE<-prodespoints00[ which(is.na(over(prodespoints00,groundtruth))),]
  BFASTNE<-prodespoints00[ which(is.na(over(prodespoints00,pva))),]
  TN<-PRODESNE[BFASTNE,]
  
  TPV<-rep(1,length(data.frame(TP)[,1]))
  TNV<-rep(2,length(data.frame(TN)[,1]))
  FPV<-rep(3,length(data.frame(FP)[,1]))
  FNV<-rep(4,length(data.frame(FN)[,1]))
  #plot(raster(TPdf))
  TPdf<-SpatialPointsDataFrame( coordinates(data.frame(TP)),data=data.frame(TPV))
  TNdf<-SpatialPointsDataFrame( coordinates(data.frame(TN)),data=data.frame(TNV))
  FPdf<-SpatialPointsDataFrame( coordinates(data.frame(FP)),data=data.frame(FPV))
  FNdf<-SpatialPointsDataFrame( coordinates(data.frame(FN)),data=data.frame(FNV))
  a<-data.frame(TPdf)
  b<-data.frame(FPdf)
  c<-data.frame(FNdf)
  d<-data.frame(TNdf)
  
  names(c)<-names(a)
  names(b)<-names(a)
  names(d)<-names(a)
  abcd<-rbind(a,b,c,d)
  
  coordinates(abcd)<-~x+y
  gridded(abcd)<-TRUE
  
  r4<-raster(abcd)
  gg4<-gplot(r4)+ geom_tile(aes(fill = factor(value))) +
    scale_fill_manual(values=c("2"="seagreen4", "3"="palegoldenrod", "4"="skyblue", "1"="palevioletred1"), 
                      na.value="white",
                      name="CONFUSION MATRIX",  
                      breaks=c("1","2","3","4"),         
                      labels=c("TRUE NEGATIVE", "FALSE POSITIVE", "FALSE NEGATIVE","TRUE POSITIVE"))+
    coord_equal()+ 
    ggtitle("SAR OLS-MOSUM")+    theme(legend.position = "none")
  
  grid_arrange_shared_legend(gg1,gg2,gg3,gg4)
}

 
 #c("red","darkgreen", "orange","skyblue")
#gplot(r)+ geom_tile(aes(fill=value)) +
#scale_fill_gradient2(low="red",high="green",
#                       limits=c(minValue(r),maxValue(r)), midpoint = 2.5)


jpeg('raster4-1.jpg', height=12, width=12, res=400,unit="in")

#generatemapGGPLOT (ptssarar3,ptssarar4,ptssarar1,ptssarar2,groundtruth,titlename='map') 
generatemapGGRASTER (ptssarar3,ptssarar4,ptssarar1,ptssarar2,groundtruth,titlename='map') 

dev.off()
save("ptssarar1",file="ptssarar1.Rdata")
save("ptssarar2",file="ptssarar2.Rdata")
save("ptssarar3",file="ptssarar3.Rdata")
save("ptssarar4",file="ptssarar4.Rdata")
pva3<-ptssarar3
pva4<-ptssarar4
pva1<-ptssarar1
pva2<-ptssarar2
,groundtruth,titlename='map') 
