x<-1:100
x<-rep(x,98)
y<-1:98
y<-rep(y,each=100)
df<-cbind(x,y)


bfastresidual<-function(df) {
  
  test.dst<-array(,c(100,98,635))
  for (i in 1:length(df[,1]))
  {
    xx<-df[i,1]
    yy<-df[i,2]
    
    spt<-zoo(inuse2na[xx,yy,],a1)
    
    frequency(spt)<-46
    na.new <- function(x) ts(na.exclude(x), frequency = 46)
    
    stlmon<-stl(spt, na.action = na.new, s.window = "per")
    spt <- ts(rowSums(stlmon$time.series)) 
    tsp(spt) <- tsp(stlmon$time.series)
    
    fitspt <- bfast(spt,h=0.15, season="harmonic", max.iter=1) 
    output<-fitspt$output
    
    test.dst[xx,yy,]<-output[[1]]$Nt
    
  }
}

bfastresi.array<-bfastresidual(df)
tl=1:length(spt)
#tl<-time(monmean)
n.2f<-inuse2na[xi[1],yi[1],]
w=1/46
#harmonic
co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2) #:fit with stl seasonality: p-value 0.9687
co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3) # fit with stl seasonality the best: p-value 0.9786
#co4 <- cos(2*pi*tl*w*4);si4 <- sin(2*pi*tl*w*4) # higher order will not improve:p-value 0.9779
xi=seq(2,96,8);
yi=seq(2,96,8)
nxx2<-c(xi-1,xi,xi+1)
nyy2<-c(yi-1, yi,yi+1)

nx21<-c(xi-1,xi,xi+1, xi-1,xi,xi+1,xi-1,xi,xi+1)
ny21<-c(yi-1,yi-1,yi-1, yi,yi,yi,yi+1,yi+1,yi+1)

nx22<-c(xi-1,xi,xi+1, xi-1,xi+1,xi-1,xi,xi+1,xi,xi,xi,xi,xi,xi)
ny22<-c(yi-1,yi-1,yi-1, yi,yi,yi+1,yi+1,yi+1,yi,yi,yi,yi,yi,yi)

#########################################
setwd('C:/Users/m_lu0002/Desktop/Climate/minnesota')
load('inuse2na.Rdata')
trend.parameters<-function(df)
{
  test.dst<-array(,c(100,98,2))
  for (i in 1:length(df[,1]))
  {
    xx<-df[i,1]
    yy<-df[i,2]
    
    spt<-zoo(inuse2na[xx,yy,],a1)
    detrend<-as.numeric(coef(lm(coredata(spt)~a1)))
    
    test.dst[xx,yy,]<-detrend
    
  }
  return(test.dst)
}
tre.parameters<-trend.parameters(df)
str(tre.parameters)
summary(tre.parameters[,,1])
#########################################

seasonal.parameters<-function(df)
{
  test.dst<-array(,c(100,98,7))
  for (i in 1:length(df[,1]))
  {
    xx<-df[i,1]
    yy<-df[i,2]
    
    spt<-zoo(inuse2na[xx,yy,],a1)
    tl=1:length(spt)
    w=1/46
    #harmonic
    co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
    co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2) #:fit with stl seasonality: p-value 0.9687
    co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3) 
    season1<-lm(coredata(spt)~co+co2+co3+si+si2+si3) #almost perfect linear relationship with stl: season: R-squared:  0.9786
    parameter<-as.numeric(coef(season1))
    
    test.dst[xx,yy,]<-parameter
    
  }
  return(test.dst)
}
df
sea.parameters<-seasonal.parameters(df)
str(sea.parameters)
#function for detrend and deseasonalize
#fitmon$output[[1]]$Nt
dst<-function(df)
{
  test.dst<-array(,c(100,98,635))
  for (i in 1:length(df[,1]))
  {
    xx<-df[i,1]
    yy<-df[i,2]
    
    spt<-zoo(inuse2na[xx,yy,],a1)
    
    w=1/46
    #harmonic
    co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
    co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2) #:fit with stl seasonality: p-value 0.9687
    co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3) 
    season1<-lm(coredata(spt)~co+co2+co3+si+si2+si3) #almost perfect linear relationship with stl: season: R-squared:  0.9786
    deseasonal<-coredata(spt)-fitted(season1)
    detrend=lm(deseasonal~a1)$residuals #other way to detrend and deseasonalize? still compariable with the spatial neighbor?
    test.dst[xx,yy,]<-detrend
    
  }
  return(test.dst)
}

bfastdst<-function(df) # aggregated bfast residuals
{
  test.dst<-array(,c(42,40,167))
  for (i in 1:length(df[,1]))
  {
 
    xx<-df[i,1]
    yy<-df[i,2]
    
    spt<-zoo(itry2[xx,yy,],a1)
    
    monmean <- aggregate(spt, as.Date(as.yearmon(a1)), mean) #should aggregate in SciDB    
    frequency(monmean)<-12
    na.new <- function(x) ts(na.exclude(x), frequency = 12)
    
    stlmon<-stl(monmean, na.action = na.new, s.window = "per")
    datamon <- ts(rowSums(stlmon$time.series)) 
    tsp(datamon) <- tsp(stlmon$time.series)
    
    fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=1) 
    
    detrend<- fitmon$output[[1]]$Nt
    test.dst[xx,yy,]<-detrend
  
    
  }
  return(test.dst)
}

#test.dst2<-array(,c(100,98,635))
x<-1:42
load('a1.Saved')
load('itry')
x<-rep(x,40)
y<-1:40
y<-rep(y,each=42)
df<-cbind(x,y)
bfastrevi<-bfastdst(df)
save(bfastrevi,file='bfastrevi.Rdata')
 
str(bfastrevi)
test.dst3<-dst(df)
##### plot bfast residuals ############
#bfast.dst3<-array(,c(100,98,167))
#f.bfast.dst3<-filter.st.median(bfast.dst3)
t2sp<-f.bfast.dst3[,,20]

 load('bfast.dst3.Rdata')
length(which(abs(bfast.dst3)>1)) #20
#bfast.dst3<-bfastdst(df)
#save(bfast.dst3, file='bfast.dst3.RData')
t1sp<-as.data.frame(as.table(t2sp))
levels(t1sp$Var1)<-seq(1,100)
levels(t1sp$Var2)<-seq(1,98)
#t1sp$Freq<-log(t1sp$Freq) 
t1sp$Var1<-as.numeric(t1sp$Var1)
t1sp$Var2<-as.numeric(t1sp$Var2)
colnames(t1sp) <- c("x","y","ndvi")
coordinates(t1sp) <- ~x+y
proj4string(t1sp) <- CRS() 
jpeg(paste('bfast residuals map 2003 02','.jpg'), height=5, width=5, res=400,unit="in")
spplot(t1sp,main='bfast residual 2003 02')
dev.off()

gstat3d <- gstat(formula=ndvi~1, data=t1sp[!is.na(t1sp$ndvi),])
vgm3d <- variogram(gstat3d)
plot(vgm3d)

model3d <- fit.variogram(vgm3d,vgm(0.006, "Exp",5,0.002))
model3d
jpeg(paste('bfast residuals virogram 2002 02','.jpg'), height=5, width=5, res=400,unit="in")
plot(vgm3d, model3d,main='bfast residuals 2002 02')
dev.off()
#################################
#df<-cbind(nx2,ny2)
#test.dst2<-dst(df)
#test.dst2[df[2,1],df[2,2],]
#test.dst2[xi[1],yi[1],]
#ac<-acf(detrend,type='covariance',lag=1000)
sepa<-c('seasonality intercept','seasonality cos1','seasonality cos2','seasonality cos3','seasonality sin1','seasonality sin2','seasonality sin3')
for (i in 4:6)
{
  i=3
  t2sp<-sea.parameters[,,i]
  #str(t1sp)
  t1sp<-as.data.frame(as.table(t2sp))
  levels(t1sp$Var1)<-seq(1,100)
  levels(t1sp$Var2)<-seq(1,98)
  #t1sp$Freq<-log(t1sp$Freq) 
  t1sp$Var1<-as.numeric(t1sp$Var1)
  t1sp$Var2<-as.numeric(t1sp$Var2)
  colnames(t1sp) <- c("x","y","ndvi")
  coordinates(t1sp) <- ~x+y
  proj4string(t1sp) <- CRS()
  
  jpeg(paste('map seasonal parameter',sepa[i],'.jpg'), height=5, width=5, res=400,unit="in")
  spplot(t1sp, main= sepa[i])
  dev.off()
  dev.new()
  gstat3d <- gstat(formula=ndvi~1, data=t1sp[!is.na(t1sp$ndvi),])
  vgm3d <- variogram(gstat3d)
  
  plot(vgm3d,main= sepa[i])
  model3d <- fit.variogram(vgm3d,vgm(0.006, "Exp",5,0.002))
  model3d
  jpeg(paste('variogram seasonal parameter',sepa[i],'.jpg'), height=5, width=5, res=400,unit="in")
  plot(vgm3d, model3d,main=sepa[i])
  dev.off()
  
}

## 

tre.parameters
jpeg(paste('trendintercept.jpg'), height=4, width=7, res=400,unit="in")

spplot(t1sp, main='trend intercept')
dev.off()
t1sp<-tre.parameters[,,2]*1000
str(t1sp)
t1sp<-as.data.frame(as.table(t1sp))
levels(t1sp$Var1)<-seq(1,100)
levels(t1sp$Var2)<-seq(1,98)
#t1sp$Freq<-log(t1sp$Freq) 
t1sp$Var1<-as.numeric(t1sp$Var1)
t1sp$Var2<-as.numeric(t1sp$Var2)
colnames(t1sp) <- c("x","y","ndvi")
coordinates(t1sp) <- ~x+y
proj4string(t1sp) <- CRS() 

gstat3d <- gstat(formula=ndvi~1, data=t1sp[!is.na(t1sp$ndvi),])
vgm3d <- variogram(gstat3d)
plot(vgm3d)

model3d <- fit.variogram(vgm3d,vgm(6, "Exp",10,0.2))
model3d
plot(vgm3d, model3d)


dst<-function(df,dataarray,test.dst) # index, data array, output array
{

  for (i in 1:length(df[,1]))
  {
    xx<-df[i,1]
    yy<-df[i,2]
    
    spt<-zoo(dataarray[xx,yy,],a1)
    
    w=1/46
    tl=1:length(spt)
    #harmonic
    co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
    co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2) #:fit with stl seasonality: p-value 0.9687
    co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3) 
    season1<-lm(as.numeric(coredata(spt))~co+co2+co3+si+si2+si3) #almost perfect linear relationship with stl: season: R-squared:  0.9786
    deseasonal<- as.numeric(coredata(spt))-fitted(season1)
    detrend=lm(deseasonal~a1)$residuals #other way to detrend and deseasonalize? still compariable with the spatial neighbor?
    test.dst[xx,yy,]<-detrend
    
  }
  return(test.dst)
}

dstagg<-function(df,dataarray,test.dst) # index, data array, output array
{
  
  for (i in 1:length(df[,1]))
  {
   
    xx<-df[i,1]
    yy<-df[i,2]
    
    spt<-zoo(dataarray[xx,yy,],a1)
    monmean <- aggregate(spt, as.Date(as.yearmon(a1)), mean) 
    
    frequency(monmean)<-12
    w=1/12
    tl=1:length(monmean)
    #harmonic
    co <- cos(2*pi*tl*w); si <- sin(2*pi*tl*w)
    co2 <- cos(2*pi*tl*w*2);si2 <- sin(2*pi*tl*w*2) #:fit with stl seasonality: p-value 0.9687
    co3 <- cos(2*pi*tl*w*3);si3 <- sin(2*pi*tl*w*3) 
    season1<-lm(as.numeric(coredata(monmean))~co+co2+co3+si+si2+si3) #almost perfect linear relationship with stl: season: R-squared:  0.9786
    deseasonal<- as.numeric(coredata(monmean))-fitted(season1)
    detrend=lm(deseasonal~time(monmean))$residuals #other way to detrend and deseasonalize? still compariable with the spatial neighbor?
    test.dst[xx,yy,]<-detrend
  
  }
  return(test.dst)
}

monthmean<-function(df,dataarray, i) # time series
{
 
xx<-df[i,1]
yy<-df[i,2]

spt<-zoo(dataarray[xx,yy,],a1)
monmean <- aggregate(spt, as.Date(as.yearmon(a1)), mean) 
return(monmean)
}

monthmeanarr<-function(df,dataarray,monthmeanarr) # time series
{ 
  for (i in 1:length(df[,1]))
{
  xx<-df[i,1]
  yy<-df[i,2]
  
  spt<-zoo(dataarray[xx,yy,],a1)
  monmean <- aggregate(spt, as.Date(as.yearmon(a1)), mean) 
  monthmeanarr[xx,yy,]<-monmean
  }
  return(monthmeanarr)
}
a1
load('a1.Saved')
itry2mon<-monthmean(df, itry2)
itry2monarr<-monthmeanarr(df, itry2,test.dst2)
str(itry2monarr)
str(itry2mon)
x<-1:42
x<-rep(x,40)
y<-1:40
y<-rep(y,each=42)
df<-cbind(x,y)
test.dst2<-array(,c(42,40,167))
test.dst3<-array(,c(42,40,635))
a1
str(test.dst)
dst4240agg<-dstagg(df, itry2, test.dst2)
dst4240<-dst(df, itry2, test.dst3)
#save(dst4240agg, file='dst4240agg')
save(dst4240, file='dst4240')
str(dst4240)
bfastrevi<-bfastdst(df)