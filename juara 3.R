
 library(scidb)
library(xts)
library(zoo)
library(rgdal)
#libs <- c("rgdal", "maptools", "gridExtra")
#lapply(libs, require, character.only = TRUE)
 
############################ compare with deter data
#options(scidb.debug=TRUE)
deterchange<-read.csv("C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/deterchange.csv")
detime<-deterchange$DataHora
dc<-as.data.frame(deterchange)

coordinates(dc)<-c('Long','Lat')
proj4string(dc)<-"+proj=longlat +ellps=aust_SA +no_defs"

dc2<-spTransform(dc,CRS("+proj=utm +zone=21 +south"))
dc3<-spTransform(dc2,CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"))
 
xchangelat<-dc@coords[,1]
ychangelong<-dc@coords[,2]
xyll<-as.data.frame(cbind(xchangelat,ychangelong))

############### modis to array indices 
getcrMatrix <- function(colrowid.Matrix, pixelSize){
 
  x <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  y <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  corner.ul.x <- -20015109.354
  corner.ul.y <- 10007554.677
  x<- round( (colrowid.Matrix[,1] -corner.ul.x - pixelSize/2)/pixelSize) 
  y<- round(-(colrowid.Matrix[,2] -corner.ul.y + pixelSize/2)/pixelSize) 
  
  cbind(x,y)
}
 ###
 
cr520<-getcrMatrix(as.data.frame(around520mo),231.6564)
## sub deter: not a good way 
xy2d<-as.data.frame(xy2)
xy2d[&xy2d$x<59079&xy2d$y>48210&xy2d$y<48360]
subdeter<-as.data.frame(xy2d[xy2d$x > 58930 & xy2d$x < 59079 & xy2d$y > 48210 & xy2d$y < 48360,])

#match<-subdeter$x%in%xychangeinmodist$xct 
#dcsub<-subset(dc3,c(-6363482,-6328965,-1191758,-1160948))
e1 <- extent(dc3)
e2 <- extent(spmodist)
ext<-intersect(e1,e2)
subset1<-intersect(dc3, e2)# ??? not the e2

over(SpatialPoints(dc3),spmodist) # no over at all
####################################################### 
coord.changeinmodist<-coordinates(changeinmot)
coord.subdeter<-coordinates(subdeter)

spmodist<-SpatialPoints(coord.changeinmodist)
spmodiss<-SpatialPoints(coordinates(changeinmos))
spdeter<-SpatialPoints(coord.subdeter)
over(dc3,spmodiss)

spmodist[!is.na(over(spmodist,spdeter))]
 

load('test3nds.RData')
changeint<-which(test3nds[,,]!=0,arr.ind=TRUE)
#49:143
changeint  #3333+167 rows for ndvi; change in seasonality: 25
subdeter$x
xct1<-changeint[,1]+58929
yct1<-changeint[,2]+48210
xychangeinmodist<-as.data.frame(cbind(xct1,yct1))
changeinmot0.5.1<-getxyMatrix(xychangeinmodist,231.6564)
em1<-extent(modis.t51)
######################################################
changeint<-which(evi2t[,,]!=0,arr.ind=TRUE)
#49:143
changeint #1992 for evi 2 seasonality; 5843+3333 for trend
subdeter$x
xct1<-changeint[,1]+58929
yct1<-changeint[,2]+48210
xychangeinmodist<-as.data.frame(cbind(xct1,yct1))
changeinmot0.5.1<-getxyMatrix(xychangeinmodist,231.6564)
em1<-extent(modis.t52)
###################### all change in seasonality: points: can use project2() ########################################
changeinseasonality<-which(test3nds[,,49:143]!=0,arr.ind=TRUE)
 
xct<-changeinseasonality[,1]+58929
yct<-changeinseasonality[,2]+48210
xychangeinmodiss<-as.data.frame(cbind(xct,yct))
changeinmos0.5.1<-getxyMatrix(xychangeinmodiss,231.6564)
#58930,48210,6,59129,48409,643

save(xy2,file='xy2')
load('xy2')
ndviallr<-scidb('ndviinuseall')
bm1<-c()
timeused=0
bpt<-c()

#58828 59679 int64
#2  1  row_id 48103    948            502             5 48103 49050
xl<-length(xchange)
 
####################xy2: all the deter points coords ##
setwd('C:/Users/m_lu0002/Desktop/Climate/minnesota')
xy2d<-as.data.frame(xy2)
 load('xy2')
xy2
####################### run bfast on deter #######################
for (i1 in 1:791)
{  
 
    i<-xy2d$x[i1]
    j<-xy2d$y[i1]
    
  
 spt1<-eviallr[i,j,][]$evi2
    #summary(spt1)
    spt<-zoo(spt1,a1)
    monmean <- aggregate(spt, as.Date(as.yearmon(a1)), median) #should aggregate in SciDB
    summary(monmean)
    frequency(monmean)<-12
    na.new <- function(x) ts(na.exclude(x), frequency = 12)
    
    stlmon<-stl(monmean, na.action = na.new, s.window = "per")
    datamon <- ts(rowSums(stlmon$time.series)) 
    tsp(datamon) <- tsp(stlmon$time.series)
    
    fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=1) 

   
    jpeg(paste(i1,'evi2_monmedian.jpg'), height=4, width=7, res=400,unit="in")
    plot(fitmon,main='')
  title(main=paste("coordinates: ", "x:", xyll$xchangelat[i1],"y:",xyll$ychangelong[i1]),sub=paste('time in deter:',detime[i1]),cex=0.8)
    dev.off()
    print(i1)
}
 
spt1<-eviallr[i,j,][]$evi2
str(spt1)
spt<-zoo(spt1,a1)
monmean <- aggregate(spt, as.Date(as.yearmon(a1)), mean) #should aggregate in SciDB
monmean2<-monmean
#[1:46]
frequency(monmean2)<-12
na.new <- function(x) ts(na.exclude(x), frequency = 12)
time(monmean)
stlmon<-stl(monmean2, na.action = na.new, s.window = "per")
datamon <- ts(rowSums(stlmon$time.series)) 
tsp(datamon) <- tsp(stlmon$time.series)
bd<-breakpoints(datamon~1+ti) #45 99 124 
str(bd)
ti<-time(monmean)
l<-lm(datamon ~ breakfactor(bd)/ti)
plot(ts(fitted(l)))
bd
fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=1) # 45 104 129
fitmon$output
fitmon$output[[1]]$Tt
 
plot(fitmon,main='')
 
dev.off()
ndviallr
#simple way:
data.shape<-readShapePoly("C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/juara",IDvar="GEOCODIG_M",proj4string=CRS("+proj=longlat +ellps=aust_SA +no_defs "))
summary(data.shape)

#GDAL way:
shape=readOGR("C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/juara.shp", layer="juara") #will load the shapefile to your dataset.
summary(shape)
shape
#modis sinusoidal
#+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs
proj4string(dft1)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
proj4string(mo)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'

mo<-SpatialPoints(coordinates(changeinmo))
#utm zone 21 south
mot<-spTransform(mo,CRS("+proj=utm +zone=21 +south"))

shape.utm<-spTransform(shape,CRS("+proj=utm +zone=21 +south"))
shape.modis<-spTransform(shape.utm,CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"))
dft1.utm<-spTransform(dft1,  CRS("+proj=utm +zone=21 +south"))
plot(shape.modis)
plot(dft1.utm,add=TRUE)

extent(dft1)
extent(shape.utm)
#ogrListLayers("C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/juara.shp") 
#ogrInfo(dsn="C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/juara.shp", layer="juara")

#iquery('store(between(MOD09Q1_MENG_20140416,58828,48103,6,58828,48103,643),s1)')
#s1r<-scidb('s1')
#iquery('between(MOD09Q1_MENG_20140416,58828,48103,6,58828,48103,643)', return=TRUE)
#iquery("store(apply(s1,ndvi,1.0*(nir-red)/(nir+red)),ndvi2) ")
#iquery("apply(s1,ndvi,1.0*(nir-red)/(nir+red))",return=TRUE )
#1  0  col_id 58828    852            502             5 58828 59679 int64
#2  1  row_id 48103    948            502             5 48103 49050 int64
#3  2 time_id     0   9201              1             0     6   643 int64
modisextent<-c(58828, 48103,59679,49050)
#taking data out of scidb
#small testing array: iquery('store(between(MOD09Q1_MENG_20140416,58832,48110,6,58842,48120,12),inuse2)')
#inus2<-scidb('inuse2')
#inuse2n<-inus2[,,][]$red  # indext  second dimension first
#   i1 i2 i3 i4 (t1) 
#j1
#j2
#j3
iquery("store(apply(MOD09Q1_MENG_20140416,ndvi,1.0*(nir-red)/(nir+red)),ndviinuseall) ")
iquery('store(subarray(MOD09Q1_MENG_20140416,59229,48103,6,59628,48502,643),inuse5)')
iquery("store(apply(inuse5,ndvi,1.0*(nir-red)/(nir+red)),ndviinuse5) ")
inus4<-scidb('ndviinuse5')
#inuse4na<-array(,c(400,400,636))
ia<-c()
ia<-inus4[,,][]$ndvi

#inuse4na<-array(ia,c(400,400, 636))
 
#second 400*400
#get 100*100 array out of scidb
#iquery('store(between(MOD09Q1_MENG_20140416,58832,48110,6,58929,48209,643),inuse)')
############ 2nd array ##
iquery('store(subarray(MOD09Q1_MENG_20140416,58930,48210,6,59129,48409,643),inuse6)')
iquery("store(apply(inuse6,ndvi,1.0*(nir-red)/(nir+red)),ndviinuse6) ")
 
iquery("store(apply(inuse6,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),evi2inuse6) ")
iquery('dimensions(inuse6)',return=TRUE)

eviinus6<-scidb('evi2inuse6')
ib<-eviinus6[1:150,1:150,][]$evi2
evi6<-array(ib,c(150,150,636))
#inuse6na<-array(ia,c(200,200,636))
evi7<-aperm (evi6,c(2,1,3)) #rotate # new second array
save(evi7,file='evi2n.Rdata')
as.numeric(evi7[12,1,])-eviinus6[12,1,][]$evi2
 s
ib<-inus6[1:150,1:150,][]$ndvi
inuse6na<-array(ib,c(150,150,636))
#inuse6na<-array(ia,c(200,200,636))
inuse7nat<-aperm (inuse6na,c(2,1,3)) #rotate # new second array
save(inuse7nat,file='2n.Rdata')
as.numeric(inuse7nat[12,1,])-inus6[12,1,][]$ndvi
##
iquery('dimensions(evi2all1)',return=TRUE)
 
ib<-eviall[1:150,1:150,][]$ndvi
inuse6na<-array(ib,c(150,150,636))
#inuse6na<-array(ia,c(200,200,636))
inuse7nat<-aperm (inuse6na,c(2,1,3)) #rotate # new second array
save(inuse7nat,file='2n.Rdata')
as.numeric(inuse7nat[12,1,])-inus6[12,1,][]$ndvi
###############################################################################

#iquery("store(apply(inuse,ndvi,1.0*(nir-red)/(nir+red)),ndviinuse) ")

inus6<-scidb('ndviinuse6')

getxyMatrix <- function(colrowid.Matrix, pixelSize){
  #Returns the coords (MODIS synusoidal) of the center of the given pixel
  #SR-ORG:6974
  x <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  y <- vector(mode = "numeric", length = length(nrow(colrowid.Matrix)))
  corner.ul.x <- -20015109.354
  corner.ul.y <- 10007554.677
  x <- corner.ul.x + (pixelSize/2) + (colrowid.Matrix[,1] * pixelSize)
  y <- corner.ul.y - (pixelSize/2) - (colrowid.Matrix[,2] * pixelSize)
  
  cbind(x,y)
}
getxy(x)
################################################################################
testingt1<-inuse2na[,,1]
#rownames(testingt1)<-paste('x',c(1:100),sept="")
#colnames(testingt1)<-paste('y',c(1:100),sept="")
str(dft1)
dft1<-as.data.frame.table(inuse2na[,,1])
which(!is.na(testingt1))
x<-rep(58832:58931,98)
y<-rep(48110:48207,each=100)
xymatrix<-cbind(x,y)
xymatrix2<-getxyMatrix(changeinmodist,231.6564)
 
xt1<-xymatrix2[,1]
yt1<-xymatrix2[,2]
dft1['Var1']<-xt1
dft1['Var2']<-yt1
coordinates(dft1)<-c('Var1','Var2')
#############################################################################################
#object.size(array(1:20,c(2,10)))
# at last: this is still the fastest way: to get the (part of) scidb array into main memory and process in R
#inuse2na<-array(,c(100,98,638)) #first array
#inuse3na<-array(,c(200,200,636))
#58930,48210,6,59129,48409,643


# iquery("load_library('dense_linear_algebra')", release=1,resp=FALSE)
#a<-inus3[,,][]$ndvi
#save(a,file='ndvi200200.Rdata')
#58930,48210,6,59129,48409,643

inuse3na<-array(a,c(200,200, 636)) # remember 2 days data are missing. should be 636
load(inuse2na,file='inuse2na.Rdata')

inuse3nat[,1,1]-inus3[,48210,6][]$ndvi
save(inuse3na, file='RData')
inuse3nat<-aperm (inuse3na,c(2,1,3)) #rotate

 
# iquery('dimensions(inuse)',return=TRUE) 100*100*638
#No    name start length chunk_interval chunk_overlap   low  high  type
#1  0  col_id 58828    852            502             5 58832 58929 int64
#2  1  row_id 48103    948            502             5 48110 48209 int64
#3  2 time_id     0   9201              1             0     6   643 int64

# test for quality quality 1
#quality1<-c()
#quality1<-inus[,,][]$quality
#unique(quality1[which(quality1!=4096)])

#sc1<-iquery('apply(between(MOD09Q1_MENG_20140416,58828,48103,6,58828,48103,643),ndvi,1.0*(nir-red)/(nir+red))',return=TRUE)

#######get a1: all the time 

date<-time_id2date(s1r[]$time_id)
a<-c()
for (i in 1:635)
 a[i]<-as.Date(date[[i]])

a1<-as.Date(a) 
#save(a1, file="a1.saved") 
#save(ndviallr, file ='ndviallr.saved')

################################## bfast on each pixels   ######################################################
#ptm <- proc.time()

ndviallr<-scidb('ndviinuseall')
bm1<-c()
timeused=0
bpt<-c()
i2=1
i1=1
i3=1
bftimes<-c()
bftimet<-c()
bptVT<-c()
bptWT<-c()
date.trend<-data.frame(0,0,0,0,0)
date.trend2<-data.frame()
date.seasonal2<-data.frame()
date.seasonal<-data.frame(0,0,0,0,0,0)
date.seasonal4<-data.frame(0,0,0,0,0,0)
bptVT2<-c()
bptWT2<-c()
bptVT3<-c()
bptWT3<-c()
date.seasonal1<-0
date.seasonal3<-data.frame()
date.trend1<-0
nbs<-c()
nbt<-c()
i4=1
timeused<-0
 
 
ar520tt<-array(,c(43,41,167)) #havent run
ar520ss<-array(,c(43,41,167))
a2<-a1[39:635]
  
(59138:59180)
(48712:48753)
 
for (i in 1:42)
{  
  for (j in 1:40)
{ 
    for(i in 1:12){ 
xx<- cosp12[,1]
yy<- cosp12[,2]
xx2<-xx+59137  
yy2<-yy+48711

#xx<-artind2[,1]+1
#yy<-artind2[,2]+1
#xx2<-xx+59137  
#yy2<-yy+48711
xym1<-c()
a1<-c()

xym1<-cbind(xx2,yy2)
xym2<-getxyMatrix(xym1,231.6564)
xymc1<-coordinates(xym2) 
xyms1<-SpatialPoints(xymc1)

proj4string(xyms1)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
krigc1<-spTransform(xyms1,  CRS("+proj=utm +zone=21 +south"))
krigc2<-spTransform(xyms1,CRS("+proj=utm +zone=21  +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#krigc2<-spTransform(xyms1,CRS("+proj=utm +zone=21  +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plotRGB(crop(s6,epu4))

points(krigc2,col='gold')
##########################################################
  spt1<-as.numeric(itry2[21,39,])

  spt<-zoo(spt1,a1)
  monmean2 <- aggregate(spt, as.Date(as.yearmon(a1)), median) #should aggregate in SciDB

  frequency(monmean2)<-12
  na.new <- function(x) ts(na.exclude(x), frequency = 12)
  
  stlmon<-stl(monmean2, na.action = na.new, s.window = "per")
  datamon <- ts(rowSums(stlmon$time.series)) 
  tsp(datamon) <- tsp(stlmon$time.series)

  fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=1) 
 # plot(fitmon)
 
 print(i4)
#d.sea<-ifelse(!is.na(breakdates(fitmon$output[[1]]$bp.Wt)),as.character(breakdates(fitmon$output[[1]]$bp.Wt)),'none')
#d.tre<-ifelse(!is.na(breakdates(fitmon$output[[1]]$bp.Vt)),as.character(round(breakdates(fitmon$output[[1]]$bp.Vt),digit=2)),'none')
#d.trp<-paste(d.tre,collapse='; ')
 
#jpeg(paste(i4,'around520evi.jpg'), height=4, width=7, res=400,unit="in")
#plot(fitmon,main='')

#title(main=paste("coordinates: ", "x:", slot(xycs,'coords')[ij,1],"y:",slot(xycs,'coords')[ij,2]),sub=paste('time :','trend',d.trp,'seasonality:',d.sea,sep='  '),cex=0.8)
 #dev.off()
#ij=ij+1
  if(fitmon$nobp$Vt==FALSE)
  {
    bptVT[i2]<-i
   bptVT2[i2]<-j
   date.trend2<-as.integer(fitmon$output[[1]]$Vt.bp)
   date.trend<-rbind(date.trend, date.trend2)
   bptVT3[i2]<-i4

  ar520tt[bptVT[i2],bptVT2[i2],date.trend2]<-100
 
   nbt[i2]<-length(as.integer(fitmon$output[[1]]$Vt.bp))
   i2=i2+1
  }
  if(fitmon$nobp$Wt==FALSE)
  {
   bptWT[i3]<-i
   bptWT2[i3]<-j
  
   nbs[i3]<-length(as.integer(fitmon$output[[1]]$Wt.bp))
   date.seasonal2<-as.integer(fitmon$output[[1]]$Wt.bp)
   date.seasonal<-rbind(date.seasonal, date.seasonal2)
   bptWT3[i3]<-i4
ar520ss[bptWT[i3],bptWT2[i3],date.seasonal2]<-500
   i3=i3+1
   print(nbs)  
  }
  
  i4=i4+1

}
}
save(ar520tt,file='ar520tt.RData')
save(ar520ss,file='ar520ss.RData')

#####
i4=1
pv1t2<-array(,c(150,150))  
pv1s2<-array(,c(150,150))
for (i in 1:150)
{  
  for (j in 1:150)
  { 
  
    spt1<-evi7[i,j,] 
    
    spt<-zoo(as.numeric(spt1),a1)
    monmean2 <- aggregate(spt, as.Date(as.yearmon(a1)), mean) #should aggregate in SciDB
    
    frequency(monmean2)<-12
    na.new <- function(x) ts(na.exclude(x), frequency = 12)
    
    stlmon<-stl(monmean2, na.action = na.new, s.window = "per")
    datamon <- ts(rowSums(stlmon$time.series)) 
    tsp(datamon) <- tsp(stlmon$time.series)
    
    fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=3) 
    
    pv1t2[i,j ]<-fitmon$output[[1]]$pvaluet$p.value
    pv1s2[i,j ]<-fitmon$output[[1]]$pvalues$p.value
    print(i4)
    
  
    i4=i4+1
    
  }
}
save( pv1t2,file='pv1t2.RData')
save( pv1s2,file='pv1s2.RData')
 bfast
pv1t<-array(,c(150,150))  
pv1s<-array(,c(150,150))
 

load('evi2n.Rdata')
load("a1.Saved") 

 ##############

 
 
i4=1

ops<-array(,c(150,150))  
opt<-array(,c(150,150))
which(!is.na(opt), arr.ind=TRUE)
for (i in 1:1305)
{  

    spt1<-evi7[i1,j1,] 
    spt<-zoo(as.numeric(spt1),a1)
    monmean2 <- aggregate(spt, as.Date(as.yearmon(a1)), mean) #should aggregate in SciDB
    frequency(monmean2)<-12
    na.new <- function(x) ts(na.exclude(x), frequency = 12)
    stlmon<-stl(monmean2, na.action = na.new, s.window = "per")
    datamon <- ts(rowSums(stlmon$time.series)) 
    tsp(datamon) <- tsp(stlmon$time.series)
      
    fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=3) 
    #jpeg(paste(i4,'loc.p-value1.jpg '), height=4, width=7, res=400,unit="in")
    plot(fitmon,main='where p-value filtered') 
    #dev.off()
    print(i4)
    if(fitmon$nobp$Vt==FALSE)
    {
      opt[i1,j1 ]<-100
    }
    if(fitmon$nobp$Wt==FALSE)
    {
      ops[i1,j1 ]<-500       
    }
    i4=i4+1    
  }
save(opt,file='opt.RData')
save(ops,file='ops.RData')
 

                             #373 points changed in bfast why
############################################################ 
i4=1
#fops<-array(,c(150,150))  
#fopt<-array(,c(150,150))
pv1stest<-array(,c(150,150))
pv1ttest<-array(,c(150,150))
fops1<-array(,c(150,150))
fopt1<-array(,c(150,150))
for (i in 1:1305)
{  
  
 
  i1<-fil.plt[i,1]
  j1<-fil.plt[i,2]
  spt1<-evi7[i1,j1,] 
  
  spt<-zoo(as.numeric(spt1),a1)
  monmean2 <- aggregate(spt, as.Date(as.yearmon(a1)), mean) #should aggregate in SciDB
  
  frequency(monmean2)<-12
  na.new <- function(x) ts(na.exclude(x), frequency = 12)
  
  stlmon<-stl(monmean2, na.action = na.new, s.window = "per")
  datamon <- ts(rowSums(stlmon$time.series)) 
  tsp(datamon) <- tsp(stlmon$time.series)
  
  fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=3) 
  #fitmon <- bfast3(datamon,h=0.15, season="harmonic", max.iter=3,pvat=filtered.pv1t,pvas=filtered.pv1s,ijk=fil.pl2[i]) 
  #jpeg(paste(i4,'loc.fil.p-value1.jpg '), height=4, width=7, res=400,unit="in")
 # plot(fitmon,main='where p-value filtered') 
  #dev.off()
  print(i4)
  
  #if(fitmon$nobp$Vt==FALSE)
  #{
    
  #  fopt1[i1,j1 ]<-100
  #  pv1ttest[i1,j1 ]<-fitmon$output[[1]]$pvalues$p.value 
  
  #}
 
  if(fitmon$nobp$Wt==FALSE)
  {
    
    fops1[i1,j1 ]<-500
    pv1stest[i1,j1 ]<-fitmon$output[[1]]$pvaluet$p.value 
 
  }
  
  i4=i4+1
  
}
##########################################################
save(fopt1,file='fopt1.RData')
save(fops1,file='fops1.RData')

fpv1t<-pfilter1(pv1t2,0.05)
fpv1s<-pfilter1(pv1s2,0.05)

filtered.pv1t<-fpv1t[[1]]
filtered.pv1s<-fpv1s[[1]]
keeptrack.pv1t<-fpv1t[[2]]
keeptrack.pv1s<-fpv1s[[2]]
 
fil.plt<-which(keeptrack.pv1t==2.0,arr.ind=TRUE) #seasonality
fil.pls<-which(keeptrack.pv1s==2.0,arr.ind=TRUE) #trend
fil.plt1<-which(keeptrack.pv1t==2.0) #seasonality
fil.pls1<-which(keeptrack.pv1s==2.0)

match(which(  pv1s2[fil.pls1]<=0.05),letmeknow2)
length(which(  pv1s2[fil.pls1]<=0.05))
which(pv1s2[which(is.na(match(fil.pls1,letmeknow2)))]<=0.05)
length(which(  pv1s2[fil.pls1]<=0.05))
#filtered.pv1t[which(keeptrack.pv1t==2.0)]<0.05
#110              filtered p-value less than 0.05
length(which(pv1t2     <=0.05)) 
 #110 < 0.05

################################################################
pfilter1 <-function(a,val=0.05)
{
   b<-a
  row<-length(a[,1])-1
  col<-length(a[,2])-1
  for (i in 2:col)
    {
   for (j in 2:row)    
   {
   
    middle=a[i,j]
    ih<-i+1
    il<-i-1
    
    jh<-j+1
    jl<-j-1
 
    neimax<-max(c(
      a[ih,jh ], a[il,jh ], a[i,jh ], 
      a[ih,jl ], a[il,jl], a[i,jl], 
      a[ih,j],    a[il,j] ))
    
    neimin<-min(c(
      a[ih,jh ], a[il,jh ], a[i,jh ], 
      a[ih,jl ], a[il,jl], a[i,jl], 
      a[ih,j],    a[il,j]) )
    neimedian<-median(c(
      a[ih,jh ], a[il,jh ], a[i,jh ], 
      a[ih,jl ], a[il,jl], a[i,jl], 
      a[ih,j],    a[il,j]) )
    
    if (middle>=neimax+val | middle <= neimin-val )
    {
      b[i,j]<-2  
      a[i,j]<-neimedian
    }
  }
  }
  a[ which(a==0)]<-NA
  return(list(a,b))
  
}








setwd()
which(!is.na(ar520t),arr.ind=TRUE)
save(date.seasonal, file = "date.seasonal.RData")
#save
for(i6 in 1: 144)
#  testing1[bptVT[i6],bptVT2[i6],as.numeric(date.trend[i6+1,])]<-100
  testing6[bptWT[i6],bptWT2[i6],as.numeric(date.trend[i6+1,])]<-500
summary(testing6)
############################################################################################

for(i in 2:12)
  { 
  xx<- cosp12[i,1]
  
  yy<- cosp12[i,2]
  spt1<-as.numeric(arraystevid[xx,yy,])
 monmean1
  spt<-zoo(spt1,monmean1)
  #monmean2 <- aggregate(spt, as.Date(as.yearmon(a1)), median) #should aggregate in SciDB
  plot(spt)
 frequency(spt)<-12
  frequency(monmean2)<-12
 stlmon<-stl(spt, na.action = na.new, s.window = "per")
  na.new <- function(x) ts(na.exclude(x), frequency = 12)
  
 stlmon<-stl(monmean2, na.action = na.new, s.window = "per")
  datamon <- ts(rowSums(stlmon$time.series)) 
  tsp(datamon) <- tsp(stlmon$time.series)
 
  fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=1) 
   plot(fitmon)
  }
xx<- cosp3[,1]
yy<- cosp3[,2]
kric1<-as.data.frame(cbind(xx,yy))  
###########



  
save(date.seasonal, file = "date.seasonal.RData")
#save
#for(i6 in 1: 144)
  ##  testing1[bptVT[i6],bptVT2[i6],as.numeric(date.trend[i6+1,])]<-100
  #testing6[bptWT[i6],bptWT2[i6],as.numeric(date.trend[i6+1,])]<-500
#summary(testing6)
load('a1.saved')

list.files()
######################## original codes ######################################################################
 for (i in 58832:59400)
{  for (j in 48810: 48900)
 # i=58828
{ ptm <- proc.time()
   i=58828
j=48810
  spt1<-ndviallr[i,j,][]$ndvi
spt<-zoo(spt1,a1)
monmean <- aggregate(spt, as.Date(as.yearmon(a1)), mean)

#frequency(spt)<-46
frequency(monmean)<-12
na.new <- function(x) ts(na.exclude(x), frequency = 12)
#stl2<-stl(spt, na.action = na.new, s.window = "per")
stlmon<-stl(monmean, na.action = na.new, s.window = "per")
#datats <- ts(rowSums(stl2$time.series)) # sum of all the components (season,abrupt,remainder)
datamon <- ts(rowSums(stlmon$time.series)) 
tsp(datamon) <- tsp(stlmon$time.series)
fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=1) 
plot(fitmon)
fitmon$output
#tsp(datats) <- tsp(stl2$time.series) # assign correct time series attributes
#fit <- bfast(datats,h=0.15, season="harmonic", max.iter=1) 

proc.time() - ptm
print(fitmon$Time)
timeused1<-proc.time() - ptm
timeused<-timeused1+timeused
print(i4)
print(timeused)

plot(fitmon)

if(fitmon$nobp$Vt==FALSE)
{bptVT[i2]<-j
 bptVT2[i2]<-i4
# date.trend[i2]<-fitmon$output[[1]]$Vt.bp
 #date.trend2[i2]<-time(monmean[fitmon$output[[1]]$Vt.bp])
 nbt[i2]<-length(as.integer(fitmon$output[[1]]$Vt.bp))
 date.trend2<-as.integer(fitmon$output[[1]]$Vt.bp)
 date.trend<-rbind(date.trend, date.trend2)
bftimes[i2]<-fitmon$Time
i2=i2+1
}
if(fitmon$nobp$Wt==FALSE)
{bptWT[i3]<-j
 bptWT2[i3]<-i4
 #date.seasonal[i3,j]<-fitmon$output[[1]]$Wt.bp
 #date.seasonal2[i3,j]<-time(monmean[fitmon$output[[1]]$Wt.bp])
 date.seasonal3<-time(monmean[fitmon$output[[1]]$Wt.bp])
 nbs[i3]<-length(as.integer(fitmon$output[[1]]$Wt.bp))
 date.seasonal2<-as.integer(fitmon$output[[1]]$Wt.bp)
 date.seasonal<-rbind(date.seasonal, date.seasonal2)

 date.seasonal4<-rbind(date.seasonal4, date.seasonal3)


i3=i3+1
}

i4=i4+1

}
}
save(date.seasonal4, file="date.seasonal4.saved") 

for (q in 1:90)
{

 datase<- subset(as.Date(t(date.seasonal4))[(q*6+1) : (q*6+6)],!duplicated(as.Date(t(date.seasonal4))[(q*6+1) : (q*6+6)]))
 print(datase)

}

#as.Date(na.omit(date.seasonal))
<-fitmon$output[[1]]$Wt.bp
j=1
i3=1
nb1<-length(as.integer(fitmon$output[[1]]$Wt.bp))
date.seasonal<-as.integer(fitmon$output[[1]]$Wt.bp)
date.seasonal2<-c(1,1,1,1,1)
date.seasonal<-rbind(date.seasonal1, date.seasonal)
time(monmean[fitmon$output[[1]]$Wt.bp])
  #Named num [1:5] 26 51 81 110 142

###############################################################################################

################## missing data
b<-c()
for (i in 1:636)
  b[i]<-ifelse(as.integer((a1[i+1]-a1[i]))==8,NA,as.Date(a1[i]))
as.Date(na.omit(b))

#[1] "2000-12-26" "2001-06-10" "2    " "2002-12-27" "2003-12-27" "2004-12-26" "2005-12-27"
#[8] "2006-12-19" "2010-12-27" "2011-12-27" "2012-12-26"
#note the end of the year, the data in 06 18 2001 is missing, and the year 2007 2008 2009 is missing 

mdates1 <- seq(from=as.Date("2000-01-01"), to=as.Date("2000-02-11"), by=8)
mdates2<-c('2001-06-17','2006-12-27','2013-12-27')
mdates2<-as.Date(mdates2)
mdates3<-c(mdates1,mdates2)
#pad with na
empty <- zoo(,mdates3)
ndviwithna <- merge(ndviss1z, empty, all=TRUE)
########################################################################
       
       ########### time testing######
       i4=1
       timeused<-0
       proc.time() <-0
       for (i in 59450:59500)
       {  
         
         for (j in 48850: 48900)
           # i=58828
           for( i1 in 8:54)
           {
{
  ptm <-proc.time()
  timeused1<-proc.time() - ptm
  
  spt1<-ndviallr[i,,][]$ndvi
  
  timeused1<-proc.time() - ptm
  timeused<-timeused1+timeused
  
  print(i4)
  
  print(timeused)
  i4=i4+1
}
           }

######################## animation #############################################################
oopt = ani.options(interval = 0.2, nmax = 167)
## use a loop to create images one by one
for (i10 in 1:ani.options("nmax")) {
 testing12<-testing11[,,i10]
   as.matrix(testing12)
   testing12[1,1]<-2
   myImagePlot(testing12)
  ani.pause() ## pause for a while ('interval')
}

## restore the options
ani.options(oopt)
## see ?ani.record for an alternative way to set up an animation
### 2. Animations in HTML pages ###
saveHTML
saveHTML({
  ani.options(interval = 0.30, nmax = 167)
  par(mar = c(3, 3, 2, 0.5), mgp = c(2, 0.5, 0), tcl = -0.3, cex.axis = 0.8, cex.lab = 0.8,
      cex.main = 1)
  for (i10 in 1:ani.options("nmax")) {
    testing58<-testing5[,,i10]
    testing68<-testing6[,,i10]
    as.matrix(testing58)
    as.matrix(testing68)
    #testing57[1,1]<-2
    par(mfrow=c(1,2))
    image(testing58,main=paste(time(monmean)[i10],'Change in Trend'))
    image(testing68,main=paste(time(monmean)[i10],'Change in Seasonality'))
    ani.pause() ## pause for a while ('interval')
  }
}, img.name = "change in trend and seasonality", title = "change in trend", description = c("change in trend: for each pixel"))


#rollapply(Ts, width = 3, by = 2, FUN = mean, align = "left")

#####################################################

project1<-function(a)
{ 
  
    index1<-which(a!=0,arr.ind=TRUE)
    t<-as.integer(index1[,3])
    nr<-as.integer(table(t))
    i4=1
    i3=1
    tu<-unique(t)
 
     
    for (i2 in 1: (length(nr)-1))
    {
   
    i<-as.integer(index111[i3:(i3+nr[i2]),1]) 
    ii<-i+58831
    
    j<-as.integer(index111[i3:(i3+nr[i2]),2]) 
    jj<-j+48109
    
    
    
    i3=i3+nr[i2]

    a1<-c()
   
    a1<-(testing1[i,j,tu[i2]])
    a1<-as.data.frame(a1)
    xym1<-cbind(ii,jj)
    xym2<-getxyMatrix(xym1,231.6564)
    a1$x<-xym2[,1]
    a1$y<-xym2[,2]
    coordinates(a1)<-c('x','y')
 
    
    proj4string(a1)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
    #a1<-spTransform(a1,  CRS("+proj=utm +zone=21 +south"))
    #a2<-cbind(a2,a1)

    return(a1)
    }
        
}


project2<-function(a ,nr2 ,tu ,i2)
{  
 
 

  tu1<-tu
  #index1<-which(a!=0,arr.ind=TRUE)
 index1<-a
 
 
 i<-as.integer(index1[(1+nr2[i2]):( nr2[i2+1] ),1]) 
 
  ii<-i+59137  
    j<-as.integer(index1[(1+nr2[i2]):( nr2[i2+1] ),2]) 
    
  jj<-j+48711

  xym1<-c()
  a1<-c()
 
    xym1<-cbind(ii,jj)
    xym2<-getxyMatrix(xym1,231.6564)
    xymc1<-coordinates(xym2) 
    xyms1<-SpatialPoints(xymc1)
   
    proj4string(xyms1)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
    a1<-spTransform(xyms1,  CRS("+proj=utm +zone=21 +south"))
 a2<-spTransform(a1,  CRS("+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
    return(a2)
  }
 
   
 setwd("C:/Users/m_lu0002/Desktop/Climate/minnesota")
index2<-which(ar520t!=0,arr.ind=TRUE)
index2<-cosp4
t<-as.integer(index2[,3])

#extent(spatialPolygons)
nr<-as.integer(table(t))
tu<-unique(t)
nr1<-c()
nr1[2:(length(nr)+1)]<-nr
nr1[1]<-0 
nr2<-cumsum(nr1)

load('ar520t.Rdata')
for(i in 1:(length(tu)) )
{
   #testingip<-project2(ar520t,nr2,tu,i)
  testingip<-project2(cosp4,nr2,tu,i)
 print(testingip) 

plot(testingip,col='red')
}
#########################################################
plot(shape.utm)
plotRGB(crop(s5,ex4))
testxy<-cbind(401080,8725236)
testxyd<-as.data.frame(testxy)
testxyd
coordinates(testxyd)<-~V1+V2
#plot(testxyd)
points(testxyd,col='gold',add=TRUE,cex=3)
click()
cccc

ex4<-extent(394168.3 , 406183.6,8720344, 8731158)
plot(ex4)
setwd("C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/US")

rdate<-substr(band52 ,start=10,stop=16)
s6<-stack(band52[1],band42[1], band32[1])

rjul<-as.integer(substr(rdate,start=5,stop=nchar(rdate)))
ryear<-as.integer(substr(rdate,start=1,stop=4))
myear<-as.integer(substr(monmean3[tu],start=1, stop=4))
 
o1=1
o2=1
par(mar=c(1,1,1,4),oma=c(3,3,3,3))
load('monmean2.Rdata')
monmean3<-time(monmean2)
#save(monmean2,file='monmean.Rdata')
saveGIF({
  
  ani.options(interval =2,nmax = length(tu))

  for (i10 in 1:length(tu)) {
   # ani.options("nmax")
    testingip<-c()
   
    testingip<-project2(cosp4,nr2 ,tu,i10)
 
#    if(o1<length(band52))
#  {
#    if (myear[i10]>=ryear[o1])
 #     { 
  #      s6<-stack(band52[o1],band42[o1], band32[o1])
   #     o1=o1+1
  #      o2=o1-1
   #   }
   
  #  }
#  rdate2<-substr(band52[o2] ,start=10,stop=16)
 # rjul<-as.integer(substr(rdate2,start=5,stop=nchar(rdate)))
  #mon<- date.mdy(rjul, weekday = FALSE)$month
  #day<- date.mdy(rjul, weekday = FALSE)$day
  
  plotRGB(crop(s6,epu4))
     points(testingip,cex=1,pch=2,add=TRUE,col='skyblue')
    # text(dcu, Nr, cex=1,col='gold')
     title( paste('date:', monmean3[tu[i10]]),cex=0.7,outer=TRUE)
  #'  ', paste("raster date:", ryear[o2],mon,day,sep='-')),col='blue',
  } 
 
   
 

    ani.pause() ## pause for a while ('interval')
  }
 
, interval = 0.05, movie.name = "around520ima.gif", ani.width = 600, ani.height = 600)

time(monmean)[tu]
monmean
tu
str(testingip)

plot(1,1,main=1)

#######################################################################
a<-c()
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
 
testing5558<-array(,c(100,100,167))
testing5558<-neighborsum8(testing5)
summary(testing5558)

 aaaa<-neighborsum8(test4ndt)
  
length(aaaa[!is.na(aaaa)])


###
neighborsum2<-function(a)
  {
  
  index1<-which(a!=0,arr.ind=TRUE)
 

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
        a[ih,j,t],a[il,j,t],a[i,jh,t],a[i,jh,t],
        a[ih,j,tl],a[il,j,tl],a[i,jh,tl],a[i,jh,tl],
        a[ih,j,th],a[il,j,th],a[i,jh,th],a[i,jh,th])
  if (asum==0)
    {
    a[i,j,t]<-0.01
    print(c(i,j,t))

    print(asum)
    }
  }
  return(a)   
}
testing666<-array(,c(100,100,167))
testing666<-neighborsum2(testing66)
testing6.1<-testing66
summary(testing666)
summary(testing5)
length(which(testing6668==500)) # 1253 vs. 1027 #226 
                                #  113  vs. 76  #31     

#1253 Vs. 950  filtered 303  24.2%
which(testing666==500) #113  vs. 62   filtered 51   45.1%
31/113
dim(a)[1]
neighborsum(testing5)
index111<-which(!is.na(testing55),arr.ind=TRUE)
summary(index111[,3])

max(index111[,3])
a111<-array(1:125,dim<-c(5,5,5))
neighborsum

#Here's a simple way of getting a neighbour position avoiding worrying about the sides:

#int leftX = (x - 1 + width) % width;
#int rightX = (x + 1) % width;
#int aboveY = (y - 1 + height) % height;
#int belowY = (y + 1) % height;

testing56<-array(,c(100,100,167))
testing55<-testing5
testing66<-testing6

testing55.8<-floor(testing5558)
testing66.8<-floor(testing6668)
which(testing55.8==0.01)
summary(testing5)
testing5[which(is.na(testing5))]<-0
testing6[which(is.na(testing6))]<-0
testing56<-testing5+testing6

# same result as above
load('test4ndt.Rdata')
#testing56[1,1]<-2
summary(testing5)
testing55.8<-floor(testing5558)
test3ndt[!is.na(test3ndt[,,])]-
str(test4ndt[!is.na(test4ndt[,,])])

testing66.6<-test3nds
saveGIF({
  
  ani.options(interval =0.3,nmax = 160)
  for (i10 in 1:ani.options("nmax")) {
    
    testing588<-ar520t[,,i10]
    #testing598<-testing55.8[,,i10]
    
    as.matrix(testing588)
   # as.matrix(testing598)
    #testing57[1,1]<-2
    #par(mfrow=c(1,2))
  image(testing588,main=paste(time(monmean)[ i10 +11],'Change in trenf'))
  #image(testing598,main=paste(time(monmean)[i10],'filtered- 8 (24) neighbor'))
   ani.pause() ## pause for a while ('interval')
}
}
, interval = 0.05, movie.name = "520ndtrend.gif", ani.width = 600, ani.height = 600)

############################ animatino for residuals ###########################################################
str(test.dst3)
test.dst3<-dst(df)
saveGIF({
  
  ani.options(interval =0.3,nmax = 167)
  for (i10 in 1:ani.options("nmax")) {
    
    testing588<-test.dst3[,,i10]
    #testing598<-testing55.8[,,i10]
    
    as.matrix(testing588)
    # as.matrix(testing598)
    #testing57[1,1]<-2
 
    image(testing588,main=paste(time(monmean)[i10],'residuals'))
    #image(testing598,main=paste(time(monmean)[i10],'filtered- 8 (24) neighbor'))
    ani.pause() ## pause for a while ('interval')
  }
}
, interval = 0.05, movie.name = "anires.gif", ani.width = 600, ani.height = 600)

getwd()
############################# create images from array ########################################################

testing55<-testing5
testing66<-testing6
testing56<-testing5+testing6
testing56[1,1,1]<-2

myImagePlot(testing)}

############### plot matrix as image#################################
 
st = STIDF(geometry(cs), begin, as.data.frame(cs), end)
pt = SpatialPoints(cbind(7, 52), CRS(proj4string(cs)))
as.data.frame(st[pt,,1:5])

setwd("C:/Users/m_lu0002/Desktop/Climate/minnesota/deter")
deter040608<-read.shp('juara/Deter_20040608_pol.shp')
deter1=readOGR('juara/Deter_20040608_pol.shp', layer="Deter_20040608_pol") #will load the shapefile to your dataset.
 
proj4string(deter1)<-CRS('+init=epsg:4618')
#spplot(deter1)
detert1<-spTransform(deter1,CRS("+proj=utm +zone=21 +south"))
inter2<-crop(detert1,em1)
#em1
#monmean['2004-06-01'] #53
#time(monmean)

load('a1.saved')
#monmean <- aggregate(spt, as.Date(as.yearmon(a1)), mean) #should aggregate in SciDB
plot(shape.utm,add=TRUE)
plot(detert1,add=TRUE)
#get juara
inter<-intersect(detert1,shape.utm)
#get modis
#esexy<-extent(sssexy)
#inter2<-crop(detert1,esexy)
inter2<-crop(detert1,em1)
 
plot(inter3,add=TRUE,col='gold')
plot(inter2,add=TRUE,col='pink')
plot(sssexy,col='green') #extent of the whole array: esexy

plot(em1,add=TRUE)#extent of second array
#x<-c(coordinates(ssexy)[,1],coordinates(ssexy)[2,1],coordinates(ssexy)[1,1],coordinates(ssexy)[1,1])
#y<-c(rep(coordinates(ssexy)[,2],each=2),coordinates(ssexy)[1,2])
#exy<-cbind(x,y)

#2014 06

#save(monmean,file='monmean.Rdata')
setwd("C:/Users/m_lu0002/Desktop/Climate/minnesota")
load('test3ndt.Rdata')
load('monmean.Rdata')
monmean
str(test3ndt)
it<-grep('2003',time(monmean))

it[1]:it[12]
 
##############################deter polygoy and bfast ######################################
ex42<-extent(c( 366971.9, 403203.4, 8811367 ,8845732) ) # extent for the 150 by 150 array
plot(modis.mt52, add=TRUE,col='gold')
plot(ex42)
as(spatialPolygons, "SpatialGrid")
?SpatialGrids
load('mevi2t.Rdata')
load('monmean')
 
str(mevi2t) #150by150 media

#change6<-which(!is.na(mevi2t[,,45:167]),arr.ind=TRUE) # for a year from 2003 10

#for (vp in 45:167)
#{
#  load("mevi2t.Rdata")
bfastchangepoint<- function(changearray)
{  
change7<-which(!is.na(changearray[,, ]),arr.ind=TRUE) #0.05
 
xct1<-change7[,1]+58929
xct2<-change7[,1] 
#length(unique(alltct1))
yct1<-change7[,2]+48210
yct2<-change7[,2] 
alltct1<-change7[,3]
 
dfallxyt<-as.data.frame(cbind(xct2,yct2,alltct1))
names(dfallxyt)<-c('x','y','t')
 
coordinates(dfallxyt)<-~x+y #make the time value for searching
 
#################
xychangeinmodist<-as.data.frame(cbind(xct1,yct1))
changeinmot0.5.1<-getxyMatrix(xychangeinmodist,231.6564)
spmodist51<-SpatialPoints(coordinates(changeinmot0.5.1))
proj4string(spmodist51)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'

modis.mt52<-spTransform(spmodist51,CRS("+proj=utm +zone=21 +south"))
return(modis.mt52)
}
modis.mt52<-bfastchangepoint(test4ndt)
load('test4ndt.Rdata')

###################################################################
over(modis.mt52,deterpoints)


md1<-modis.mt52[deterpoints,]
md2<-length(deterpoints[modis.mt52,])
str(md1)
unique(slot(md1,"coords"))
deterpoints[modis.mt52,]
length(deterpoints)
################ get bfast time
 
m.in.d<- modis.mt52[spatialPolygons,] 
##
#if(!is.na(m.in.d@coords[1]))
#{
  ModisCRS<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
  modievi1<-spTransform(m.in.d,CRS(ModisCRS))

dfmoev1<-as.data.frame(slot(modievi1,'coords'))

comoev1<-getcrMatrix(dfmoev1,231.6564) 
comoev.x<-comoev1[,1] 
comoev.y<-comoev1[,2] 
comoev.x1<-comoev.x-58929
comoev.y1<-comoev.y-48210
#t<-which(!is.na(mevi2t[ comoev.x1,comoev.y1,]),arr.ind=TRUE)
dfcomo<-as.data.frame(cbind(comoev.x1,comoev.y1))
#
names(dfcomo)<-c('x','y')
coordinates(dfcomo)<-~x+y
tttt<-over(dfcomo,dfallxyt)
tttt1<-time(monmean)[tttt$t]
 load('monmean')
# get polygon time

o1<-over(modis.mt52,spatialPolygons)
#o2<-over(spatialPolygons,modis.mt52)
pind<-o1[which(!is.na(o1))]
length(o1)
 
tbfd3<-c()
pind<-o1[which(!is.na(o1))]
 
for (k in 1:length(pind))
  tbfd3[k]<-  slot(spatialPolygons[pind[k]] ,'polygons')[[1]]@ID

#str(modis.mt52[spatialPolygons,])
td<-substr(tbfd3,start=nchar(tbfd3)-9,stop=nchar(tbfd3))
 
dutd2<-as.Date(td, '%Y-%m-%d')
#timebf.in.d<-time(monmean)[dfallxyt[dfcomo,]$t]
#save(timebf.in.d,file='timebf.in.d.Rdata')

################## plot bfast and deter time 

#save(tttt1,file='tttt1')
load('tttt1.Rdata')
bf1<-as.integer(tttt1)
dt1<-as.integer(dutd2)
bdt<-as.data.frame(cbind(bf1,dt1))
names(bdt)<-c('x','y')
coordinates(bdt)<-~x+y

zbdt<-zerodist(bdt) 
tzbdt<-table(zbdt[,2])
tzbdt2<-c()
for(i in 1:length(tzbdt))
 tzbdt2[i]<-tzbdt[[i]]

str(zerodist(bdt) )
unre2<-as.numeric(dimnames(tzbdt)[[1]])
deter<-jitter(bdt[unre2]@coords[,2],amount=0)
bfast<-jitter(bdt[unre2]@coords[,1],amount=0)
deter<- bdt[unre2]@coords[,2] 
bfast<- bdt[unre2]@coords[,1] 
max(tzbdt2)
#jpeg('bvs.d.cex2.jpg', height=12, width=12, res=400,unit="in")
plot(bfast,deter,xlab='bfast time (seasonality)' ,pch=4,,cex=tzbdt2/5,ylab='deter time',ylim=c(11500,16000),xlim=c(11500,16000))
 hexbinplot(deter~bfast,aspect = 1,
shape=10,ylim=c(11500,16000),xlim=c(11500,16000))
plot(hex)
hexbinplot(hedb)
hedb<-hexbin(bfast,deter)
abline(a=0,b=1,col='red')
legvals<-c(1,5,10,15,20,25)
legend('bottomright',legend=legvals,pch=4,pt.cex=legvals/5,title='number of points',cex=0.8)
#dev.off()
#dbf2<-bf1[order(bf1)]
#dt2<-dt1[order(dt1)]
 
rin<-table(dt1)

rin2<-c()
for(i in 1:length(rin))
rin2[i]<-rin[[i]]
 
unre<-as.numeric(dimnames(rin)$dt1)
#rin3<-as.character(rin2)

bf1=rep(as.integer(time(monmean[vp])),length(rin))
 jpeg(paste(time(monmean)[vp],'det vs. bf-time.jpg '), height=4, width=7, res=400,unit="in")
plot(bf1,unre,xlab='bfast time (trend)',main=paste('bfast time:', time(monmean)[vp],sep=''),pch=1,cex=rin2/5,ylab='deter time (trend)',ylim=c(11500,16000),xlim=c(11500,16000))
#main=paste('bfast time:', time(monmean)[vp],sep='')
legvals<-c(1,5,10,15,20,25)
legend('bottomright',legend=legvals,pch=1,pt.cex=legvals/5,title='number of points',cex=0.8)
abline(a=0,b=1,col='red')

 dev.off()
}
}
extent(spatialPolygons)
######################################################
save(tbfd3,file='tbfd3.Rdata')
save(dutd2,file='dutd2.Rdata')
#timebf.in.d #bfast in deter
#dutd2 # deter polygons in bfast
 
################################################################
 
#inttimemevi2<-as.integer(timemevi2)
#indexst<-as.matrix(cbind(index(dfallxy),inttimemevi2))
#stmodis1<-STI (spmodist51,time(monmean),indexst,endTime = delta(time(monmean)))
#stmodis1<-STI (spmodist51,time(monmean),data.frame(time=inttimemevi2))
#               ,endTime = delta(time(monmean)))

#mevi2t[which(!is.na(over(dfallxy,dfcomo)))]
 )

#date.ddmmmyy(15279)
#testa1<-array(rep(0,6781),dim=c(150,150,117),dimnames=list(as.character(xct1),as.character(yct1),as.character(alltct1)))
 

 ################################ plot bfast points and deter polygon#######
 
color1=c('darkslategray3','darkslategray4','darkslategrey','darkturquoise','darkviolet','deeppink','blue3','blue4','blueviolet','brown','brown1','skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4')
color2=c('skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4')
 
#color=c('skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4')
setwd( "C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/deter")

deterall<-list.files(pattern = "\\.shp$")
 
deterg<-deterall
dl<-c()
for(i in 1:length(deterg))

{
  dl[i]<-ifelse(length(grep('_pol',deterg[i]))!=0,substr(deterg[i],start=7,stop=nchar(deterg[i])-8), substr(deterg[i],start=7,stop=nchar(deterg[i])-4))
 
  if( nchar(dl[i]) < 8 )
  dl[i]<-paste(dl[i],'01',sep='')   
}
dl[70]<-'20100430'
datedeter<-as.Date(dl,'%Y%m%d')
save(datedeter,file='datedeter.Rdata')
 
 
 
deterg<- deterall
dl<-substr(deterg,start=7,stop=nchar(deterg[i])-8) #time 
ldeter<-substr(deterg, start=1, stop=nchar(deterg)-4) #layer

#jpeg('around520evi2bfast061007vspoly2007all.jpg', height=12, width=20, res=800,unit="in")

plot(ex4,xlab='',ylab=''   )
plotRGB(crop(s5,ex4))
for (i in 1:length(deterg))
  {
 # j=i+0
  #deterall[j]
  
  deter1=readOGR(deterg[i], layer=ldeter[i]) #will load the shapefile to your dataset.
  proj4string(deter1)<-CRS('+init=epsg:4618')
 # detert1<-spTransform(deter1,CRS("+proj=utm +zone=21 +south"))
 detert1<-spTransform(deter1,CRS("+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
  inter5<-crop(detert1,epu4) #ex4
           
   if(!is.null(inter5))

    plot(inter5,col=color1[1],add=TRUE)
  }
plotRGB(crop(s6,epu4))
"+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "
################## put the time as ID #######################################
i=3
data.first <- readOGR(deterg[i], gsub(".shp","",deterg[i]))
proj4string(data.first)<-CRS('+init=epsg:4618')
data.first<-spTransform(data.first,CRS("+proj=utm +zone=21 +south"))
polygons<-crop(  data.first,ex42)
 
  for(j in 1:length(polygons)) 
  slot(slot(polygons,'polygons')[[j]],'ID')<-paste(j,datedeter[3])
polygons<-slot(polygons,'polygons') 
polygons    

for (i in 4:length(deterg)) {
     
     data.temp <- readOGR(deterg[i], gsub(".shp","",deterg[i]))
     proj4string(data.temp)<-CRS('+init=epsg:4618')
     data.temp<-spTransform(data.temp,CRS("+proj=utm +zone=21 +south"))
     cdata.temp<-crop(  data.temp,ex42)
     
     if (!is.null(cdata.temp))
     {
      
       for(j in 1:length(cdata.temp)) 
         {
         
           slot(slot(cdata.temp,'polygons')[[j]],'ID')<-paste(j,datedeter[i],sep='.')
           
         }
       polygons <- c(slot(cdata.temp, "polygons"),polygons)
       print(slot(polygons[[i]],'ID'))
     }
}
 
ss1<-c()
for (i in 1:155)
ss1[i]<- slot(polygons[[i]],'ID')
ss1
length(polygons)


proj4string(spatialPolygons)<-CRS("+proj=utm +zone=21 +south")
 
spatialPolygons <- SpatialPolygons(polygons)
setwd( "C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/deter") 
save(spatialPolygons,file='Spatial1.R')
load('Spatial1.R')
 

extent(spatialPolygons)
plot(spatialPolygons[!is.na(over(spatialPolygons, modis.mt52))])

tbf[k]<-slot(
 (modis.mt52[spatialPolygons,])[k] ,'polygons')[[1]]@ID
 
td<-substr(tdeter,start=nchar(tdeter)-9,stop=nchar(tdeter))
utd<-unique(td)
dutd2<-as.Date(td, '%Y-%m-%d')
 
#setwd("C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/deter")
tdeter<-c()
for (k in 1:length(spatialPolygons))
tdeter[k]<-slot(spatialPolygons[modis.mt52,][k],'polygons')[[1]]@ID
getwd()
save(tdeter,file='tdeter.R')
load('tdeter.R')
load('timedeter.R')
save(dutd,file='timedeter.R')
save(dutd2,file='timedeterrep.R')
plot(modis.mt52[spatialPolygons,]
spdf <- SpatialPolygonsDataFrame(spatialPolygons, data.frame(row.names=ss1[1:155]))
 plot(spdf)

spdf
writeOGR(spdf, dsn="combined4.shp", layer="combined4", driver="ESRI Shapefile")
 
com<-readOGR('combined.shp',layer='combined')
proj4string(com)<-CRS('+init=epsg:4618')
detert1<-spTransform(com,CRS("+proj=utm +zone=21 +south"))
################################################################
plot(com)
inter5<-crop(detert1,ex42)
#writeOGR(inter5, dsn="combinedcroped.shp", layer="combinedcroped", driver="ESRI Shapefile")
plot(inter5,color='gold')
#save(inter5, file='allpolygoncroped.Rdata')

#over(modis.mt52,inter5)
#plot(modis.mt52,col='red',pch=1,add=TRUE,cex=0.5)
for(i10 in 43:52)
{
  testingip<-project2(ar520t,nr2 ,tu,i10)  
  points(testingip,col='red' ,cex=1,add=TRUE)
}

dev.off()
#plot(modis.t52,col='gold',pch=2,add=TRUE,cex=0.6)
title('bfast monthly median trend Sep 03 to 04 vs. deter 04 to 05')
legend("topright", legend = dl, col = color1,ncol = 2, cex = 0.2, lwd = 4)
dev.off()
###################################################################
 
#deter1=readOGR('juara/Deter_20040729_pol.shp', layer="Deter_20040729_pol") #will load the shapefile to your dataset.
#deter1=readOGR('juara/Deter_20040825_pol.shp', layer="Deter_20040825_pol") #will load the shapefile to your dataset.
#use str for slot 
#slot(slot(slot(inter2,'polygons')[[2]],'Polygons')[[1]],'coords')

#yellow 0622 pink 0608 blue 0729 skyblue 0825
#orange May red June black July green August gray march to dec
color1=c('darkslategray3','darkslategray4','darkslategrey','darkturquoise','darkviolet','deeppink','blue3','blue4','blueviolet','brown','brown1','skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4')
color2=c('skyblue','red','gold','yellow','blue','antiquewhite2','darkolivegreen1','darkgreen','purple','pink','chartreuse1','deepskyblue4')
#monmean['2004-06-01'] #53
#time(monmean)
#monmean <- aggregate(spt, as.Date(as.yearmon(a1)), mean) #should aggregate in SciDB
#plot(shape.utm)

#inter<-intersect(detert1,shape.utm)

#inter2<-crop(detert1,em1)

plot(inter2,col='gold')
plot(inter2,add=TRUE,col='pink')
plot(inter4,add=TRUE,col='blue')
plot(inter5,add=TRUE,col='skyblue')
load('a1.saved')
aggregate(a1,'yearmon')
setwd( "C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/deter")
getwd()
deterall<-list.files(pattern = "\\.shp$")
deterall
#save(monmean,file='monmean')
#time(monmean)[48]
around520<-click(n=4)
around520<-as.data.frame(around520)
coordinates(around520)<-~x+y
SpatialPoints(around520)

#test2ndt  no-structral change
#test3ndt   p-value 0.05
#test4ndt           0.025
#test5ndt           0.001   
getwd()
setwd("C:\\Users\\m_lu0002\\Desktop\\Climate\\minnesota\\juara") 
rasterall<-list.files(pattern = "\\.tif$")
band5<-rasterall[grep('228_068_L2_BAND5',rasterall)]
band4<-rasterall[grep('228_068_L2_BAND4',rasterall)]
band3<-rasterall[grep('228_068_L2_BAND3',rasterall)]
raster(band5[1])
Nr<-slot(dc2,'data')$Nr
length(band3)

for(i in 1:38)
{
 i=20 
rdate<-substr(band5[i],start=14,stop=29)
s5<-stack(band5[i],band4[i], band3[i])

jpeg(paste('cropped3@',rdate,'.jpg'), height=10, width=10, res=1200,unit="in")
dev.off()
#plotRGB(s5,main='')
plotRGB(crop(s5,ex3))
#plot(slot(s5,'layers')[[1]])
#o1<-overlay(slot(s5,'layers')[[1]],slot(s5,'layers')[[2]], slot(s5,'layers')[[3]],fun=sum)
#pt1<-slot(dc2,'coords')[10:100,]
#points(pt1,col=color[3])
 

for(i in 1:(length(tu)) )
{
  testingip<-project2(ar520t,nr[i],tu[i])
  plot(testingip,col='red',add=TRUE)
  
}

points(dcu,cex=0.5)
text(dcu, Nr, cex=1,col='gold')
dcu
#pt1<-slot(dc2,'coords')[1:2,]
 
#points(pt1,col=color[1]) 
#text(pt1,'2',cex=0.5)
#legend('topright',col=c(color[3],color[1]),legend=c('Point 10','point 2'),lwd = 4,pch=1)
title(main = rdate,col='blue',cex=2)
dev.off()
}

setw
u=7
substr(rasterall)
s5<-stack(rasterall[24],rasterall[25], rasterall[26])
tiff('rgb20030905-11-2.tif', height=14, width=14, res=1200,unit="in")
plotRGB (s5, r=3, g=2, b=1,)
pt1<-slot(dc2,'coords')[10:11,]

i=3
points(pt1,col=color[i])
pt1<-slot(dc2,'coords')[1:2,]

i=1
points(pt1,col=color[i]) 
text(pt1,'2',cex=0.5)
dev.off()
ex1<-extent(385000,422860,8712000,8750000)
ex2<-extent(390000,420000,8714000,8748000)
ex3<-extent(394000,420000,8720000,8740000)
plotRGB(crop(s5,ex1))
list.files()
tif1<-raster('C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/LANDSAT_5_TM_20010627_228_068_L2_BAND3.tif')
tif2<-raster('C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/LANDSAT_5_TM_20010627_228_068_L2_BAND4.tif')
tif3<-raster('C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/LANDSAT_5_TM_20010627_228_068_L2_BAND5.tif')
s1<-stack(tif1,tif2,tif3)
tiff('rgb20110831-11-2.jpg', height=12, width=12, res=1200,unit="in")
plotRGB (s1, r=3, g=2, b=1,)
dev.off()
plot(tif3)
plot(shape.utm,add=TRUE)
plot(em1,add=TRUE)


pt1<-slot(dc2,'coords')[10:11,]
i=3
points(pt1,col=color[i])
pt1<-slot(dc2,'coords')[1:2,]
i=1
points(pt1,col=color[i]) 
dev.off()
 
getwd()
####get all USGS rasters #########
setwd( "C:/Users/m_lu0002/Desktop/Climate/minnesota/juara/US")
rasterall2<-list.files(pattern = "\\.TIF$")

band52<-rasterall2[grep('_B5',rasterall2)]
raster(band52[1])
band42<-rasterall2[grep('_B4',rasterall2)]
band32<-rasterall2[grep('_B3',rasterall2)]

################## try to reproject raster, so slow, not sceed

r521<-raster(band52[1])
newproj <- "+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0" 
proj4string(r521)<-'+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0' 
pr <- projectExtent(r521, newproj) 
projectRaster(s6,pr )
 
###### ## convert the extent to USGS CRS ## 
ep4<-as(ex4,'SpatialPolygons')
proj4string(ep4)<-'+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
epu4<-spTransform(ep4,CRS("+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
usgspro<-"+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "
 
### cropped data with point on ##
for (i in 1: length(band52))
{
i=1
  s6<-stack(band52[i],band42[i], band32[i])
  rdate<-substr(band52[i],start=10,stop=16)
  jpeg(paste(rdate,'USGS data around520.jpg'), height=12, width=12, res=800,unit="in")

plotRGB(crop(s6,epu4))
plot(krigc1)
s6
points(dcu,cex=0.5,pch=2)
text(dcu, Nr, cex=1,col='gold')
title( paste('date:',rdate),col='blue',cex=1,outer=TRUE)
dev.off()
}

## convert points to USGS projection
s7<-spTransform(band52[i],CRS("+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))
dcu<-spTransform(dc2,CRS("+proj=utm +zone=21 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
 
 
for(i10 in 41:51)
{
  testingip<-project2(ar520t,nr2 ,tu,i10)  
  points(testingip,col='red' ,cex=1,add=TRUE)
}
#epu4 the extent around point 520
load("evi52.Rdata")
load('testing5.Rdata')
load('inuse2na.Rdata')
load('xy2')
load('a1.saved')
load('ndviallr.saved')
load('ar520tt.RData')
load('monmean')
load('test3nds.Rdata')

setwd("C:/Users/m_lu0002/Desktop/Climate/minnesota")
