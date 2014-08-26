library(scidb)
library(xts)
library(zoo)
#options(scidb.debug=TRUE)

scidbconnect(host = "gis-ts.uni-muenster.de", port = 8080L)

iquery('store(between(MOD09Q1_MENG_20140416,58828,48103,6,58828,48103,643),s1)')
s1r<-scidb('s1')
#iquery('between(MOD09Q1_MENG_20140416,58828,48103,6,58828,48103,643)', return=TRUE)
#iquery("store(apply(s1,ndvi,1.0*(nir-red)/(nir+red)),ndvi2) ")
#iquery("apply(s1,ndvi,1.0*(nir-red)/(nir+red))",return=TRUE )
#1  0  col_id 58828    852            502             5 58828 59679 int64
#2  1  row_id 48103    948            502             5 48103 49050 int64
#3  2 time_id     0   9201              1             0     6   643 int64
#taking data out of scidb
 
#sc1<-iquery('apply(between(MOD09Q1_MENG_20140416,58828,48103,6,58828,48103,643),ndvi,1.0*(nir-red)/(nir+red))',return=TRUE)
date<-time_id2date(s1r[]$time_id)
a<-c()
for (i in 1:635)
  # a[i]<-date[[i]]
  a[i]<-as.Date(date[[i]])
a
a1<-as.Date(a)
str(a1)
#iquery("store(apply(MOD09Q1_MENG_20140416,ndvi,1.0*(nir-red)/(nir+red)),ndviall) ") #around 25 minutes
ndviallr<-scidb('ndviall')

save(a1, file="a1.saved") 
save(ndviallr, file ='ndviallr.saved')
load('a1.saved')
#load('ndviallr.saved')
a1
getwd()

################################## bfast on each pixels   ######################################################
#ptm <- proc.time()
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
testing5<-array(,c(100,100,167))
testing6<-array(,c(100,100,167))

for (i in 58832:58929)
{  
  for (j in 48110:48209)
{ 
  spt1<-ndviallr[i,j,][]$ndvi
  
  spt1<-ndviallr[58832,48170,][]$ndvi
  
  spt<-zoo(spt1,a1)
  monmean <- aggregate(spt, as.Date(as.yearmon(a1)), mean) #should aggregate in SciDB

  frequency(monmean)<-12
  na.new <- function(x) ts(na.exclude(x), frequency = 12)
  
  stlmon<-stl(monmean, na.action = na.new, s.window = "per")
  datamon <- ts(rowSums(stlmon$time.series)) 
  tsp(datamon) <- tsp(stlmon$time.series)

  fitmon <- bfast(datamon,h=0.15, season="harmonic", max.iter=1) 
  plot(fitmon)
  print(i4)

  
  if(fitmon$nobp$Vt==FALSE)
  {
    bptVT[i2]<-i-58829
   bptVT2[i2]<-j-48109
   date.trend2<-as.integer(fitmon$output[[1]]$Vt.bp)
   date.trend<-rbind(date.trend, date.trend2)
   bptVT3[i2]<-i4

   testing5[bptVT[i2],bptVT2[i2],date.trend2]<-100


  # date.trend[i2]<-fitmon$output[[1]]$Vt.bp
   #date.trend2[i2]<-time(monmean[fitmon$output[[1]]$Vt.bp])
   nbt[i2]<-length(as.integer(fitmon$output[[1]]$Vt.bp))
   i2=i2+1
  }
  if(fitmon$nobp$Wt==FALSE)
  {
   bptWT[i3]<-i-58829
   bptWT2[i3]<-j-48109
   
   nbs[i3]<-length(as.integer(fitmon$output[[1]]$Wt.bp))
   date.seasonal2<-as.integer(fitmon$output[[1]]$Wt.bp)
   date.seasonal<-rbind(date.seasonal, date.seasonal2)
   bptWT3[i3]<-i4
   testing6[bptWT[i3],bptWT2[i3],date.seasonal2]<-500
   i3=i3+1
   print(nbs)  
  }
  
  i4=i4+1

}
}
save(testing5,file='testing5.RData')
save(testing6,file='testing6.RData')
save(date.seasonal, file = "date.seasonal.RData")
#save
for(i6 in 1: 144)
#  testing1[bptVT[i6],bptVT2[i6],as.numeric(date.trend[i6+1,])]<-100
  testing6[bptWT[i6],bptWT2[i6],as.numeric(date.trend[i6+1,])]<-500
summary(testing6)
############################################################################################















######################## original codes ######################################################################
 for (i in 59300:59400)
{  for (j in 48820: 48900)
 # i=58828
{ ptm <- proc.time()

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


#Ts<-(c(1,1,1,2,3,3,2))
#rollapply(Ts, width = 3, by = 2, FUN = mean, align = "left")
#down vote
#?rollapply
# if (t==dim(a)[3])
#  asum<-sum(
#   a[i+1,j,t],a[i-1,j,t],a[i,j+1,t],a[i,j+1,t],
#  a[i+1,j,t-1],a[i-1,j,t-1],a[i,j+1,t-1],a[i,j+1,t-1]
# )

#if (i==dim(a)[2])
# asum<-sum(
#a[i-1,j,t],a[i,j-1,t],a[i,j+1,t],
#a[i-1,j,t-1],a[i,j-1,t-1],a[i,j+1,t-1],
#a[i-1,j,t+1],a[i,j-1,t+1],a[i,j+1,t+1])
#)

#######################################################################
a<-c()
neighborsum8<-function(a)
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
      a[ih,j,t],a[il,j,t],a[i,jh,t],a[i,jh,t],a[ih,jh,t],a[il,jl,t],a[il,jh,t],a[ih,jl,t],
      a[ih,j,tl],a[il,j,tl],a[i,jh,tl],a[i,jh,tl],a[ih,jh,tl],a[il,jl,tl],a[il,jh,tl],a[ih,jl,tl],
      a[ih,j,th],a[il,j,th],a[i,jh,th],a[i,jh,th],a[ih,jh,th],a[il,jl,th],a[il,jh,th],a[ih,jl,th])
    if (asum==0)
    {
      a[i,j,t]<-0.01
      print(c(i,j,t))
      
      print(asum)
    }
  }
  return(a)   
}
testing5558<-array(,c(100,100,167))
testing5558<-neighborsum8(testing5)
summary(testing5558)





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
index111[25,]
testing5[70,12,27]
testing5[83,2,143]
max(index111[,3])
a111<-array(1:125,dim<-c(5,5,5))
neighborsum
?apply
#Here's a simple way of getting a neighbour position avoiding worrying about the sides:

#int leftX = (x - 1 + width) % width;
#int rightX = (x + 1) % width;
#int aboveY = (y - 1 + height) % height;
#int belowY = (y + 1) % height;
### 3. GIF animations ###
#save(testing66,file='testing66.R')

testing56<-array(,c(100,100,167))
testing55<-testing5
testing66<-testing6
#testing666 filtered
#testing555 filtered
testing55.8<-floor(testing5558)
testing66.8<-floor(testing6668)
which(testing55.8==0.01)
summary(testing55.5)
testing5[which(is.na(testing5))]<-0
testing6[which(is.na(testing6))]<-0
testing56<-testing5+testing6

r <- raster(nrows=10, ncols=10)
adjacent(r, cells=c(1, 55), directions=8, pairs=TRUE)
a <- adjacent(r, cell = c(1,55,90), directions=4, sorted=TRUE)

a
r[c(1,55,90)] <- 1
r[a] <- 2
plot(r)
# same result as above
rook <- matrix(c(NA, 1, NA,
                 1, 0, 1,
                 NA, 1, NA), ncol=3, byrow=TRUE)
adjacent(r, cells = c(1,55,90), directions=rook, sorted=TRUE)
#testing56[is.na(testing56)]<-0
#testing56[1,1]<-2
summary(testing5)
testing55.8<-floor(testing5558)

testing66.6<-floor(testing666)
saveGIF({
  
  ani.options(interval =0.3,nmax = 167)
  for (i10 in 1:ani.options("nmax")) {
    testing588<-testing5[,,i10]
    testing598<-testing55.8[,,i10]
    
    as.matrix(testing588)
    as.matrix(testing598)
    #testing57[1,1]<-2
    par(mfrow=c(1,2))
  image(testing588,main=paste(time(monmean)[i10],'Change in trend'))
  image(testing598,main=paste(time(monmean)[i10],'filtered- 8 (24) neighbor'))
   ani.pause() ## pause for a while ('interval')
}
}
, interval = 0.05, movie.name = "bm_demo2.gif", ani.width = 600, ani.height = 600)




############################# create images from array ########################################################

testing55<-testing5
testing66<-testing6
testing56<-testing5+testing6
testing56[1,1,1]<-2

myImagePlot(testing)}
############### plot matrix as image#################################
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(mar = c(3,5,2.5,2))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7)
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}