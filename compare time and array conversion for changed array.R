bt<-time(monmean)
#save(bt,file='bt.Rdata')
 
load("Spatial.R")
load("tdall2.Rdata")
load("pp0.Rdata")
load("rpp.Rdata")
load("xyc.Rdata") # whole 22500 points

load('C:/Users/m_lu0002/Dropbox/mengluchu/bfast2/bt.Rdata')
load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast2/t3darrbfamul2.Rdata")
MODISCRS<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
UTM21S<-"+proj=utm +zone=21 +south"

SpatialPolygon<-spatialPolygons 
tdall<-tdall2
##
pp0<-SpatialPoints(matrix(c(0,0),c(1,2)))@coords
pp<-lapply(1:length(SpatialPolygon), function(i) polygontopoint(c(58930:59079),c(48210:48359),spatialPolygons[i]))
#number of points of each polygon
allpoints<-lapply(1:154, function(i)   pp[[i]]@coords )
###
rpp<-lapply(1:length(SpatialPolygon), function(i) length(pp[[i]]))

for(i in 1:155)
pp0<-rbind(pp0,xyc[pp[[i]],]@coords)
  
pp0<-data.frame(pp0)
names(pp0)<-c('x','y')
coordinates(pp0)<-~x+y  # deter points
#pp00<-data.frame(pp0)
pp0<-pp0[-1]
###

#save(pp0,file="pp0.Rdata") # deter points
#save(rpp,file="rpp.Rdata") #repeating
 
polygontopoint<-function(x=c(58930:59079),y=c(48210:48359),sppolygon,crs=MODISCRS) {
  #sppolygon<-prodes67808
  
  x1<-rep(x,each=length(y))
  y1<-rep(y,length(x))
  
  xyd<-as.data.frame(cbind(x1,y1))
  xyc<-getxyMatrix(xyd,231.6564)
  xyc<-as.data.frame(xyc)
  coordinates(xyc)<-~x+y
  SpatialPoints(xyc)
  #as.character(round(unique(xyc@coords[,2])))
  
  proj4string(xyc)<-crs
 
  proj4string(sppolygon)<-CRS("+proj=utm +zone=21 +south  +ellps=WGS84 ")
  sppolygon=spTransform(sppolygon,  CRS(crs) )
  points<-xyc[sppolygon,]
 
  
 return( points)
}
#p1<-polygontopoint(c(58930:59079),c(48210:48359),spatialPolygons)

deterpolytopointstoSTSDF<-function(tdall2=tdall,pp01=pp0,rpp2=rpp, crs=MODISCRS,attr.name="value",months=0.3)
{
  
  detertime<-tdall2
  dtime<-as.POSIXct(detertime )
  
  uniquetime<-unique(dtime)
  rp1<-table(dtime)
  tt2<-uniquetime[order(uniquetime)]
  
  lt<-length(uniquetime)
  
  tn<-unlist(lapply(lt:1, function(i) rep(i,rp1[i])))
  
  tn2<- unlist(lapply(1:length(tn), function(i) rep(tn[i],rpp2[[i]])))  #rpp2: the number of points
  
  sn2<-c(1:length(pp01))
  index1<-as.matrix(cbind( sn2 , tn2  ))
  
  data2<-data.frame(rep(1,length(pp01)))
  names( data2)<-attr.name
  proj4string(pp01)<-CRS('+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs')
  
  ltime<-  tt2+months*3600*24*30
  etime<- tt2-months*3600*24*30
  stsdf2<-STSDF(pp01,etime,index=index1,data2,ltime ) 
  
  #stpolygon<-spTransform(stsdf2, CRS(crs))
  
  return(stsdf2)
}

arraytoSTFDF2<-function(array,crs,attr.name="value", alltime=bt, x=c(58930:59079),y=c(48210:48359),months=0.3)
{
  
  t<-c(1:dim(array)[3])
  tt2<-unique(as.POSIXct(alltime[t]))  # array index to time
  x1<-rep(x,each=length(y))
  y1<-rep(y,length(x))
  xyd<-as.data.frame(cbind(x1,y1))
  xyc<-getxyMatrix(xyd,231.6564)
  xyc<-as.data.frame(xyc)
  coordinates(xyc)<-~x+y
  proj4string(xyc)<-crs
  
  data2<-as.vector(array) 
  names( data2)<-attr.name
  data2<-data.frame(data2)
  
  
  
  
  ltime<-  tt2+months*3600*24*30
  etime<-  tt2-months*3600*24*30
  stfdf1<-STFDF(xyc, etime, data2,ltime) 
  return(stfdf1)
}


changearraytoSTSDF<-function(array,crs,attr.name="value", alltime=bt, x=c(58930:59079),y=c(48210:48359),months=6)
{
  
  #itrydf<-as.data.frame.table(array) #y x t interate y ->x ->t: (x1,y1,t1)(x1,y2,t1)(x2,y1,t1)(2,2,t1)(1,1,t2)(1,2,t2)
  #aa2<-itrydf$Freq
  #array<-t3darrbfamul2
  #crs=MODISCRS
  #attr.name="value"
  # tt1<-dimnames(array)[3]
  # tt<-substr(tt1[[1]],start=2,stop=nchar(tt1[[1]]))
  # tt2<-as.Date(tt,format='%Y.%m.%d')
  ##
  #woud be easier to get the unordered index of time than get the right index of space
  change7<-which(!is.na(array),arr.ind=TRUE) #0.05
  change8<-which(!is.na(array))
  x11<-change7[,1]+x[1]-1 # for the second 150 by 150 array   
  y11<-change7[,2]+y[1]-1
  t<-change7[,3]
  
  xyd<-data.frame(cbind(x11,y11))
  xyc<-getxyMatrix(xyd,231.6564)
  xyc<-data.frame(xyc)
  names(xyc)<-c('x','y')
  xyc$ID<-c(1:length(x11))
  coordinates(xyc)<-~x+y 
  
  tt2<-unique(as.POSIXct(alltime[t]))  # array index to time
  
  zdxyc<-zerodist(xyc)     # non-unique spatial points   
  
  sel = lapply(2:length(zdxyc[,2]),  function(i) !identical( zdxyc[i-1,2], zdxyc[i,2]))
  sel<- unlist(sel)
  sel2<-c(TRUE,sel)
  sel2<-sel2[-length(sel2)]
  zuni = zdxyc[sel2,]
  length(zdxyc[,2])
 
  spl<-length(xyc)        
  
  uniquexyc<- remove.duplicates(xyc)
  dunixyc<-data.frame(uniquexyc)
  dunixyc$IDnew<-c(1:length(uniquexyc))
  
  noID<-data.frame(cbind(dunixyc$ID,dunixyc$IDnew))
  
  remo<- match(zuni[,1],dunixyc$ID)  # return the index of ID
  
  oldindex<-c(dunixyc$ID,zuni[,2])
  newindex<-c(dunixyc$IDnew,noID[remo,2])
  
  oldandnewindex<-data.frame(cbind(newindex,oldindex))
  names(oldandnewindex)<-c("new","old")
  oldnew<-oldandnewindex[order(oldindex),] 
  lt<-length(tt2)
  #spl<-length(newone2)
  #runi<-rank(unique())
  tn<-lapply(1:lt, function(i) rep(i,table(t)[i]))
  #sn<-lapply(1:spl, function(i) rep(runi[i],newone2[i]))
  index1<-as.matrix(cbind(oldnew$new,unlist(tn) ))
  
  data2<-na.omit(as.vector(array)) 
  names( data2)<-attr.name
  data2<-as.data.frame(data2)
  
  
  ltime<-  tt2+months*3600*24*30
  etime<-  tt2-months*3600*24*30
  spxyc<-as( uniquexyc,"SpatialPoints")
  proj4string( spxyc)<-crs
  stsdf1<-STSDF(spxyc, etime,index=index1,data2,ltime) 
  return(stsdf1)
} 
  
  #find the index of point with ID: uniquexyc)$ID = zdxyc[,1]
  # which(is.na( match( data.frame(xyc)$ID, data.frame(uniquexyc2)$ID)))- unique(zdxyc[,2])
  #  means the same way of removing points when the second pair is removed
  #removed points: unique(zdxyc[,2])
  #idofremove<-   xyc[zdxyc[,2],]@data$ID
  #idofcoresp<-   xyc[zdxyc[,1],]@data$ID
  #length(u1)
  #u1<-uID[ which(!is.na( match(data.frame(uniquexyc)$ID, zuni[,1])))]
  
  #matched<-which(!is.na( match( data.frame(uniquexyc)$ID,zuni[,1])))
  #matched1<-which(!is.na( match( zuni[,1],data.frame(uniquexyc)$ID)))
  #matched2<-which(is.na( match( zuni[,1],data.frame(uniquexyc)$ID)))
  #matched4<-which(is.na( match( zdxyc[,1][matched3],data.frame(uniquexyc)$ID))) #so we repeat this stupid thing
  
  #matched5<-which(!is.na( match(zdxyc[,2], zdxyc[,1][matched3][matched4])))
  #
  #length(matched2)
  #length(u2)
  #originalidtonewid
  #originalid<-zdxyc[,2]
  #correspondingid1<-zdxyc[,1]
  #
  #newid2<- zdxyc[,1][matched1] are u1
  #newid2<- zdxyc[,1][which(!is.na( match(zdxyc[,2], zdxyc[,1][matched2])))]#:duplicated points correspond to zdxyc[,1]
  #matched3<-which(!is.na( match(zdxyc[,2], zdxyc[,1][matched2]))) # these points in 1 can be found in 2
  #zdxyc[,2][matched] 
  #zdxyc[,2][matched3] # but if there are 3 duplicated points, not all matched 3 can be found
  
  #u2<-uID[ which(!is.na( match( data.frame(uniquexyc)$ID,zdxyc[,1][matched1])))]
  #u2<-uID[ which(!is.na( match( data.frame(uniquexyc)$ID,zdxyc[,1][matched3])))]
  
  #u3<-uID[ which(!is.na( match( data.frame(uniquexyc)$ID,zdxyc[,1][matched3][matched5])))]
  #c(zdxyc[,2][matched],
  #  u2<- which(!is.na( match( zdxyc[,1][matched3],data.frame(uniquexyc)$ID)))
  # table(zdxyc[,2])!=1
  #  zdxyc[,1][matched3]
  
  
  
  #  newid2<-
  #    zdxyc[,1][which(!is.na( match(zdxyc[,2], zdxyc[,1][matched2])))]
  
  #  uID[ which(!is.na( match( data.frame(uniquexyc)$ID,zdxyc[,1])))]
  
  #uID[ which(!is.na( match( data.frame(uniquexyc)$ID,rep1)))]
  
  # newone<-rep(0,spl)        # all spatial points
  #  allind<- 1:lxyc
  #  uni<- unique( zdxyc[,1] )
  #  table_r = rbind(label=uni , count=sapply(uni ,function(x)sum( zdxyc[,1]==x)))
  #  reppoints<-unique(c(zdxyc[,1], zdxyc[,2])) 
  #  uniqueind<- allind[which(is.na(match(allind,reppoints)))]
  #  allpoints<-c(unique(zdxyc[,1]),uniqueind)
  
  # table_r is the table() but the results are not sorted
  #  table_r[1,] euqal to uni          
  
  # data.frame(xyc)[1,]
  
  #sorted
  #uni<- unique(( zdxyc[,1][order(zdxyc[,1])])) # get unique index   
  #iind<-table( zdxyc[,1][order(zdxyc[,1])])  #get the repeatition
  #as.numeric(unlist(dimnames(table( zdxyc[,1])))) 
  
  #  newone[uni]<-table_r[2,]  # put in how many times the points are replicated
  # newone2<-newone[-zdxyc[,2]]  #create index and remove replicated corresponence
  #newone2
  #newone2<- newone2+1 # repeat 1 times is not to repeat ( n points repeat n+1 times)
  
 # too painful to build stsdf directly. 

length(which(!is.na(t3darrbfamul2)))
#not so computationally efficient, but easiest way
changestfdf<-arraytoSTFDF2(t3darrbfamul2, alltime=bt, x=c(58930:59079),y=c(48210:48359),MODISCRS,months=0.3)   #change array, points
stsdfchangepoints1<-as(changestfdf,"STIDF")
stsdfchangepoints2<-as(changestfdf,"STSDF")

deterpoinf<-deterpolytopointstoSTSDF(tdall ,pp01=pp0,rpp2=rpp,crs=MODISCRS,months=0.1)
changests<-changearraytoSTSDF(t3darrbfamul2, alltime=bt, x=c(58930:59079),y=c(48210:48359),MODISCRS,months=0.3)
plot(deterpoinf@sp)
over(deterpoinf,stsdfchangepoints1)

stsdfchangepoints2<-as(stsdfchangepoints1,"STSDF")
plot(stsdfchangepoints1)
plot(changests)

over(sts)
stsdfchangepoints1[1:20, ]
summary(changestfdf)
summary(stsdfchangepoints1)
as(stsdfchangepoints1,"Spatial")
#wholepatchsdf<-as(wholepatchdf,"STSDF")
#deterstfploy<-as(deterstfdf,"STF")
#deterstf<-as(stdeterpoints,"STF")
#changests<-as(changestsdf,"STF")
#wholepatch<-as(wholepatchdf,"STF")


from<-changestfdf
n = length(from@sp)
m = nrow(from@time)
index1 = cbind(rep(1:n, m), rep(1:m, each=n))
# copied from sp:

df1<-data.frame(cbind(index1,from@data))
todrop<-which(is.na(df1))

from[df1[,1],df1[,2]]
dffrom<-data.frame(from)
index1
sel = apply(sapply(from@data, is.na), 1, function(x) !all(x), arr.ind=TRUE)
index = index[sel,,drop=FALSE]
length(unlist(from@data))/22500
sapply(arr1[,],is.na)
from2=from[sel,drop=FALSE]
stsdfchangepoints1@index
from@data[sel,,drop=FALSE]
from[FALSE,]

sel = apply(sapply(from@data, is.na), 1, function(x) !all(x))
index = index[sel,,drop=FALSE]


stsdfchangepoints1@sp
stsdfchangepoints2<-as(stsdfchangepoints1,"STS")
sp1<-data.frame(stdfchangep@sp)
sp1$ID<-c(1:length(sp1[,1]))
stdfchangep@index[,1]
index[,1]
stdfchangep<-na.omit(stsdfchangepoints )
#deterstsdf<-deterpolygontoSTSDF(spatialPolygons,tdall,MODISCRS,months=24) #deter polygon
#wholepatchdf<-wholepointtoSTFDF(tdall,months=0)   # points

#deterpoinf[1,]
#changestfdf[1,]

#stpoints<-over( wholepatch,deterstf)
#aa<-wholepatchdf@sp[deterstsdf@sp,]  #polygon to points, get spatial points index
#stdeterpoints<-wholepatchdf[aa,] # polygon to points, pick up the overlapped spatial points
#aaa1<-which(!is.na(over(wholepatchsdf, geometry(deterstsdf))),arr.ind=TRUE)
  #x<- geometry(stsdfchangepoints2)

#compare the time for the spatially overlapped points
x<-geometry(changests)
y<-geometry(deterpoinf) #2125 deter points
deterpoinf@sp
 #str(x)
xspin<-na.omit(over(y@sp,x@sp)) #bfast in deter 1040
xspin1<-x@index[,1][xspin]

yspin<-na.omit(over(x@sp,y@sp)) # deter in bfast  (more replicated points?) 1370 # 873?
yspin1<-y@index[,1][yspin]
length(yspin)
plot(x@sp[y@sp,]) # mind they are not the same since there are duplicated spatial points # 1370
points(y@sp[x@sp,],col='red') # 1040

#check replicated points
#xspin2<-which(!is.na(over(x@sp,y@sp)))
#which(is.na(match(xspin2,xspin)))
  
#x[1995,]# this index is related to the index when creating STS
#when there are replicated spatial points with more time stamps, but i would expect more.
#should store the deter points in a right way (spatially replicated points should have same index)
#only consider the match of time when data are spatially correlated

#starting and end time for the interval alt 1 alt2  
#bfast change points
  alt1<-lapply(1:length(xspin1), function(i) time(x[xspin1[i],]@time))
  alt2<-lapply(1:length(xspin1), function(i) x[xspin1[i],]@endTime)
  int1<-cbind(unlist(alt1),  unlist(alt2))
  interval1<-Intervals(int1, closed = c(TRUE, FALSE))
  


#table(unlist(alt4)) #can check the times when the spatial points are met
#as.POSIXct( unlist(alt4), origin = "1960-01-01")

  # timeindex1 <- x@index[,2][xspin]
 
  alt3<-lapply(1:length(yspin1), function(i) start(y[yspin1[i],]))
  alt4<-lapply(1:length(yspin1), function(i) y[yspin1[i],]@endTime)
  int2<-cbind(as.numeric(alt3),as.numeric(alt4)) 
  interval2<-Intervals(int2, closed = c(TRUE, FALSE))
  # timeindex2<-y@index[,2][yspin]
   
  ret = interval_overlap(interval2,interval1)  # bfast in deter (deter 1370) /
  ret2 = interval_overlap(interval1,interval2) #deter in bfast
  #ret = interval_overlap(interval1[timeindex1], interval2[timeindex2])
  head( ret,n=100)
  length (unlist(ret))
length(alt4) #1370 deter change points
 table(interval2)
which(ret!=0)
alt2[unique( unlist(ret))]
x1<-x[,y@index[,2][unique( unlist(ret))]]  #times of bfast change points when there are change happening in deter
x1[spinx,]
length(x1)
x1@sp
x1@time
x1[xspin1,]
list(ret) length(ret)
alt1[unique( unlist(ret))]
spin
ret
bfastindeter<-unlist(ret)
plot( interval1@.Data[bfastindeter])
length(ret)

str(interval1) # 1049
aaa<-na.omit(over(wholepatchsdf, geometry(deterstsdf)))
stdeter1<-wholepatchsdf[aaa1]
 
aaa[1:length(aaa)]
length(aaa1 )
 
dif<-c()
for(i in 1:132)
 dif[i]<- rpp[[aaa2[i]]] - as.vector(table(aaa))[59]
which(dif!=0)
 
 
aaa2<-as.numeric(dimnames(table(aaa))$aaa)
 
 time (geometry(deterstsdf)@time)
 time (geometry(wholepatchdf)@time)
 

 
deterpolygontoSTSDF<-function(SpatialPolygon=spatialPolygons,tdall2=tdall,crs=MODISCRS,attr.name="value",months=6)
{
   detertime<-tdall2
   dtime<-as.POSIXct(detertime )
 
   uniquetime<-unique(dtime)
   rp1<-table(dtime)
   tt2<-uniquetime[order(uniquetime)]
   
   lt<-length(uniquetime)
 
   tn<-lapply(lt:1, function(i) rep(i,rp1[i]))
   sn<-c(1:155)
   index1<-as.matrix(cbind(unlist(sn),unlist(tn) ))
   data2<-data.frame(rep(1,length(SpatialPolygon)))
   names( data2)<-attr.name
   proj4string(SpatialPolygon)<-CRS("+proj=utm +zone=21 +south")
  
   ltime<-  tt2+months*3600*24*30
   etime<- tt2-months*3600*24*30
 stsdf2<-STSDF(SpatialPolygon,etime,index=index1,data2,ltime ) 
  
  stpolygon<-spTransform(stsdf2, CRS(crs))
  
  return(stpolygon)
}

wholepointtoSTFDF<-function( tdall2=tdall,x=c(58930:59079),y=c(48210:48359),crs=MODISCRS,attr.name="value",months=6)
{
  
  detertime<-tdall2
  dtime<-as.POSIXct(detertime )
  
  odtime<-order(dtime)
  tt2<-unique(dtime[odtime])
  
  x1<-rep(x,each=length(y))
  y1<-rep(y,length(x))
  
  xyd<-as.data.frame(cbind(x1,y1))
  xyc<-getxyMatrix(xyd,231.6564)
  xyc<-as.data.frame(xyc)
  coordinates(xyc)<-~x+y
  SpatialPoints(xyc)
  #as.character(round(unique(xyc@coords[,2])))
  ltime<-  tt2+months*3600*24*30
  etime<- tt2-months*3600*24*30
  proj4string(xyc)<-crs
  data2<-data.frame(rep(1,length(xyc)*length(tt2)))
  names( data2)<-attr.name
   stsdf2<-STFDF(xyc,etime,data2,ltime ) 
   
  return(stsdf2)
}

 
stsdf1[,time='2002-02-01']

x = as(stsdf1[,1:10], "xts")

changearraytoSTSDF<-function(array,crs,attr.name="value", alltime=bt, x=c(58930:59079),y=c(48210:48359),months=6)
{
  
  #itrydf<-as.data.frame.table(array) #y x t interate y ->x ->t: (x1,y1,t1)(x1,y2,t1)(x2,y1,t1)(2,2,t1)(1,1,t2)(1,2,t2)
  #aa2<-itrydf$Freq
  #array<-t3darrbfamul2
  #crs=MODISCRS
  #attr.name="value"
  # tt1<-dimnames(array)[3]
  # tt<-substr(tt1[[1]],start=2,stop=nchar(tt1[[1]]))
  # tt2<-as.Date(tt,format='%Y.%m.%d')
  ##
  #woud be easier to get the unordered index of time than get the right index of space
  change7<-which(!is.na(array),arr.ind=TRUE) #0.05
  change8<-which(!is.na(array))
  x11<-change7[,1]+x[1]-1 # for the second 150 by 150 array   
  y11<-change7[,2]+y[1]-1
  t<-change7[,3]
  
  xyd<-data.frame(cbind(x11,y11))
  xyc<-getxyMatrix(xyd,231.6564)
  xyc<-data.frame(xyc)
  names(xyc)<-c('x','y')
  xyc$ID<-c(1:length(x11))
  coordinates(xyc)<-~x+y 
  
  tt2<-unique(as.POSIXct(alltime[t]))  # array index to time
  
  zdxyc<-zerodist(xyc)     # non-unique spatial points   

  sel = lapply(2:length(zdxyc[,2]),  function(i) !identical( zdxyc[i-1,2], zdxyc[i,2]))
  sel<-unlist(sel)
  sel2<-  c(TRUE,sel)
  sel2<-sel2[-length(sel2)]
  zuni = zdxyc[sel2,]
 
  spl<-length(xyc)        
 
  uniquexyc<- remove.duplicates(xyc)
  dunixyc<-data.frame(uniquexyc)
  dunixyc$IDnew<-c(1:length(uniquexyc))
 
noID<-data.frame(cbind(dunixyc$ID,dunixyc$IDnew))
 
remo<- match(zuni[,1],dunixyc$ID)  # return the index of ID
 
oldindex<-c(dunixyc$ID,zuni[,2])
newindex<-c(dunixyc$IDnew,noID[remo,2])

oldandnewindex<-data.frame(cbind(newindex,oldindex))
names(oldandnewindex)<-c("new","old")
oldnew<-oldandnewindex[order(oldindex),] 
 
 
  #find the index of point with ID: uniquexyc)$ID = zdxyc[,1]
# which(is.na( match( data.frame(xyc)$ID, data.frame(uniquexyc2)$ID)))- unique(zdxyc[,2])
#  means the same way of removing points when the second pair is removed
#removed points: unique(zdxyc[,2])
#idofremove<-   xyc[zdxyc[,2],]@data$ID
#idofcoresp<-   xyc[zdxyc[,1],]@data$ID
#length(u1)
#u1<-uID[ which(!is.na( match(data.frame(uniquexyc)$ID, zuni[,1])))]

#matched<-which(!is.na( match( data.frame(uniquexyc)$ID,zuni[,1])))
#matched1<-which(!is.na( match( zuni[,1],data.frame(uniquexyc)$ID)))
#matched2<-which(is.na( match( zuni[,1],data.frame(uniquexyc)$ID)))
#matched4<-which(is.na( match( zdxyc[,1][matched3],data.frame(uniquexyc)$ID))) #so we repeat this stupid thing

#matched5<-which(!is.na( match(zdxyc[,2], zdxyc[,1][matched3][matched4])))
#
#length(matched2)
#length(u2)
#originalidtonewid
#originalid<-zdxyc[,2]
#correspondingid1<-zdxyc[,1]
#
#newid2<- zdxyc[,1][matched1] are u1
#newid2<- zdxyc[,1][which(!is.na( match(zdxyc[,2], zdxyc[,1][matched2])))]#:duplicated points correspond to zdxyc[,1]
#matched3<-which(!is.na( match(zdxyc[,2], zdxyc[,1][matched2]))) # these points in 1 can be found in 2
#zdxyc[,2][matched] 
#zdxyc[,2][matched3] # but if there are 3 duplicated points, not all matched 3 can be found

#u2<-uID[ which(!is.na( match( data.frame(uniquexyc)$ID,zdxyc[,1][matched1])))]
#u2<-uID[ which(!is.na( match( data.frame(uniquexyc)$ID,zdxyc[,1][matched3])))]

#u3<-uID[ which(!is.na( match( data.frame(uniquexyc)$ID,zdxyc[,1][matched3][matched5])))]
#c(zdxyc[,2][matched],
#  u2<- which(!is.na( match( zdxyc[,1][matched3],data.frame(uniquexyc)$ID)))
 # table(zdxyc[,2])!=1
#  zdxyc[,1][matched3]
  
  
  
#  newid2<-
#    zdxyc[,1][which(!is.na( match(zdxyc[,2], zdxyc[,1][matched2])))]

#  uID[ which(!is.na( match( data.frame(uniquexyc)$ID,zdxyc[,1])))]

  #uID[ which(!is.na( match( data.frame(uniquexyc)$ID,rep1)))]

 # newone<-rep(0,spl)        # all spatial points
#  allind<- 1:lxyc
#  uni<- unique( zdxyc[,1] )
#  table_r = rbind(label=uni , count=sapply(uni ,function(x)sum( zdxyc[,1]==x)))
#  reppoints<-unique(c(zdxyc[,1], zdxyc[,2])) 
#  uniqueind<- allind[which(is.na(match(allind,reppoints)))]
#  allpoints<-c(unique(zdxyc[,1]),uniqueind)
  
  # table_r is the table() but the results are not sorted
  #  table_r[1,] euqal to uni          
  
 # data.frame(xyc)[1,]
  
  #sorted
  #uni<- unique(( zdxyc[,1][order(zdxyc[,1])])) # get unique index   
  #iind<-table( zdxyc[,1][order(zdxyc[,1])])  #get the repeatition
  #as.numeric(unlist(dimnames(table( zdxyc[,1])))) 
  
#  newone[uni]<-table_r[2,]  # put in how many times the points are replicated
 # newone2<-newone[-zdxyc[,2]]  #create index and remove replicated corresponence
  #newone2
  #newone2<- newone2+1 # repeat 1 times is not to repeat ( n points repeat n+1 times)
  
  lt<-length(tt2)
  #spl<-length(newone2)
  #runi<-rank(unique())
  tn<-lapply(1:lt, function(i) rep(i,table(t)[i]))
  #sn<-lapply(1:spl, function(i) rep(runi[i],newone2[i]))
  index1<-as.matrix(cbind(oldnew$new,unlist(tn) ))
 
  data2<-na.omit(as.vector(array)) 
  names( data2)<-attr.name
  data2<-as.data.frame(data2)
  
  
  ltime<-  tt2+months*3600*24*30
  etime<-  tt2-months*3600*24*30
  spxyc<-as( uniquexyc,"SpatialPoints")
  proj4string( spxyc)<-crs
  stsdf1<-STSDF(spxyc, etime,index=index1,data2,ltime) 
  return(stsdf1)
} # too painful to build stsdf directly. 




 
arraytoSTFDF<-function(array,crs,attr.name, time)
{
  
  #itrydf<-as.data.frame.table(array) #y x t interate y ->x ->t: (x1,y1,t1)(x1,y2,t1)(x2,y1,t1)(2,2,t1)(1,1,t2)(1,2,t2)
  #aa2<-itrydf$Freq
 
  tt1<-dimnames(array)[3]
  tt<-substr(tt1[[1]],start=2,stop=nchar(tt1[[1]]))
  tt2<-as.Date(tt,format='%Y.%m.%d')
  ##
  x<-as.numeric(dimnames(array )[1][[1]])
  y<-as.numeric(dimnames(array )[2][[1]])
  x1<-as.numeric(c(x))
  x1<-rep(x,each=length(y))
  y1<-rep(y,length(x))
  xyd<-as.data.frame(cbind(x1,y1))
  xyc<-getxyMatrix(xyd,231.6564)
  xyc<-as.data.frame(xyc)
  
  data1<-as.data.frame.table(array) 
  data2<-as.data.frame(data1$Freq)
  names( data2)<-attr.name
  coordinates(xyc)<-~x+y
  
  SpatialPoints(xyc)
  
  proj4string(xyc)<-crs
  
  stfdf1<-STFDF(xyc,tt2,data2) 
  return(stfdf1)
}
