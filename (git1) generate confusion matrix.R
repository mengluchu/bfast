library('rgdal')
library('raster')

load('prodespoints00.Rdata') # mask
load('pdd.Rdata')           # reference data
#groundtruth<-pdd            
load('tssarar1.Rdata') # corrected st model with original array cusum
load('tssarar2.Rdata')#mosum
load('tssarar3.Rdata')#sar cusum
load('tssarar4.Rdata')#sar mosum
load('tssarar5.Rdata')#ar cusum
load('tssarar6.Rdata')#ar mosum

#color2=rep(c('skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4'),1000)

##150by 150 # coordinate of MODIS array
#x<-c(58930:59079)
#y<-c(48210:48359)

##### MODIS array index to MODIS SINUSOIDAL coordinates ##############
getxyMatrix <- function(colrowid.Matrix, pixelSize,crs=CRS("+proj=utm +zone=21 +south")){
  
  x <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  y <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  corner.ul.x <- -20015109.354
  corner.ul.y <- 10007554.677
  x <- corner.ul.x + (pixelSize/2) + (colrowid.Matrix[,1] * pixelSize)
  y <- corner.ul.y - (pixelSize/2) - (colrowid.Matrix[,2] * pixelSize)
  cbind(x,y)
}


####### array of changepoints to spatial points  ############
changepoint2sp<- function(changearray,x=c(58930:59079),y=c(48210:48359)) # map the changes from the array that stores the changes. multiple changes are mapped as one change point
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
  if(!is.null(crs))
    modis.mt52<-spTransform(spmodist51,crs)
  return(modis.mt52)
} 
#example: spfevi8<-bfastchangepoint(fevi8[,,1],x,y)

######  array of p-value to spatial points   ####
pvaluepoint2sp<- function(parray,x=c(58930:59079),y=c(48210:48359),xoff=1,yoff=1,pvalue=0.05,crs=CRS("+proj=utm +zone=21 +south")) #map changes from an array of p-values
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
  
  if(!is.null(crs))
    modis.mt52<-spTransform(spmodist51, crs)
  
  return(modis.mt52)
} 
######### generate confusion matrix: we use pontias' producer's accuracy ###########################
gcm<-function(changepoint, reference.sppoint=pdd,totalp=19167){
  #changepoint<-pvpoint
  
  modis.mt52<-changepoint
  deterpoints<-reference.sppoint
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
  N=totalp-P
  TN=N-FP
  ALL=totalp
  
  PRODUCER =TP/(TP+FP+FN)
  
  result<-cbind(TP,FN,FP,TN,PRODUCER)
  names(result)<-c("TP","FN","FP","TN","producer's accuracy")
  
  p1<-paste('total reference points:', b1,    'unique reference points not in result changepoints:',b3,
            'unique referece points in result changepoints',b4, 'unique reference points not in result changepoints', b5, 'unique result changepoints in reference points',b6,
            'unique total result changepoints',b10, sep='        ')
  
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

### generate confusion matrix from result array, can use a mask to mask only interest area (sppoints) and set p-value#######

generatecmchange<-function( result.array,mask=prodespoints00,reference.sppoints=pdd,pv=0.05,x=c(58930:59079),y=c(48210:48359)) #  put functions together and generate confusion matrix (filtered with a mask)
{
  pvpoint<- changepoint2sp(result.array,x,y )
  if(!is.null(mask))
    pvpoint<-pvpoint[mask,]
  showre<-gcm(pvpoint,reference.sppoints)
  return(showre)
}

### generate confusion matrix from p-value array, can use a mask to mask only interest area (sppoints) and set p-value#######
generatecmpvalue<-function(result.array=pvaluemx ,reference.sppoints=pdd,mask=prodespoints00,pv=0.05,x=c(58930:59079),y=c(48210:48359)) #  put functions together and generate confusion matrix (filtered with a mask)
{
  
  pvpoint<-pvaluepoint2sp(result.array,x,y,1,1,pv)
  pvpoint<-pvpoint[mask,]
  showre<-gcm(pvpoint,reference.sppoints)
  return(showre)
}


generateppvalue<-function(result.array=pvaluemx ,reference.sppoints=pdd,mask=prodespoints00,pv=0.05,x=c(58930:59079),y=c(48210:48359)) 
{
  pvpoint<-pvaluepoint2sp(result.array,x,y,1,1,pv)
  pvpoint<-pvpoint[mask,]
  return(pvpoint)
}

### do things
ptssarar1<-generateppvalue(tssarar1,pv=0.05)
ptssarar11<-generateppvalue(tssarar1,pv=0.005)
ptssarar2<-generateppvalue(tssarar2,pv=0.05)
ptssarar22<-generateppvalue(tssarar2,pv=0.025)
ptssarar3<-generateppvalue(tssarar3,pv=0.05)
ptssarar4<-generateppvalue(tssarar4,pv=0.05)
ptssarar5<-generateppvalue(tssarar5,pv=0.05)
ptssarar6<-generateppvalue(tssarar6,pv=0.05)

ttssarar1<-generatecmpvalue(tssarar1,pdd,pv=0.05)
ttssarar11<-generatecmpvalue(result.array=tssarar1 ,pv=0.2)
ttssarar2<-generatecmpvalue(tssarar2,pdd,pv=0.05)
ttssarar22<-generatecmpvalue(tssarar2,pdd,pv=0.1)
ttssarar3<-generatecmpvalue(tssarar3,pdd,pv=0.05)
ttssarar4<-generatecmpvalue(tssarar4,pdd,pv=0.05)
ttssarar5<-generatecmpvalue(tssarar5,pdd,pv=0.05)
ttssarar6<-generatecmpvalue(tssarar6,pdd,pv=0.05)

cts<-rbind(ttssarar4,ttssarar6,ttssarar2,ttssarar22，ttssarar3，ttssarar5，ttssarar1，ttssarar11)

#jpeg('Pontius Producer Accuracy ts.jpg', height=6, width=10, res=400,unit="in")
names1<-c("OLS-MOSUM p-value: 0.05","AR(1) OLS-MOSUM pvalue: 0.05","SAR OLS-MOSUM p-value: 0.05", "SAR OLS-MOSUM p-value: 0.1","OLS-CUSUM p-value: 0.05","AR(1) OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.05","SAR OLS-CUSUM p-value: 0.2")
tex1='Confusion Matrix of Different Methods'

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
#dev.off()

