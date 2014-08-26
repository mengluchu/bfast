setwd( "C:/Users/m_lu0002/Desktop/Climate/minnesota")
load('evi2n.Rdata' ) #evi7, second array
load('a1.saved')
dst<-function(df)
{
  test.dst<-array(,c(100,98,635))
  for (i in 1:length(df[,1]))
  {
    xx<-df[i,1]
    yy<-df[i,2]
   
    spt<-zoo(evi7[xx,yy,],a1)
    
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

test.dst2<-array(,c(100,98,635))
x<-1:100
x<-rep(x,98)
y<-1:98
y<-rep(y,each=100)
df<-cbind(x,y)
test.dst2<-dst(df)

##get neighbors
getnei<-function(x,y,t,array)
{
  # x=41;y=25;t=31
  # array=arraystevid
  nxx2<-c(x-1,x,x+1)
  nyy2<-c(y-1, y,y+1)
 # t2<-c(rep(4,8),1,2,3,5,6,7)+(t-1)
  t2<-c(rep(0,8),-1,-2,-3)+(t)
 
 n.1<-c(array[nxx2,nyy2,t2[1]])
  
  n.2<-c(array[nxx2[2],nyy2[2],t2[9:11]])
  
  n.3<-c(n.1[-5],n.2)
  return(n.3)
}
# getnei(xi[1],yi[1],ti[1],arraystevid)
 
#artind<-which(!is.na(ar520t2),arr.ind=TRUE)

#xi<-artind[,1]
#yi<-artind[,2]
#ti<-artind[,3]


str(arraystevid)

#ti<-seq(1,10)

####

getcube<-function(x,y,t,array ) 
{ 
  air3d2<-as.data.frame(matrix(NA,1,4))
  names(air3d2)<-c('nei', 'nx22','ny22', 't2')
  
 # for (  j in 1: length(t)) 
  #{
   # ti<-t[j]
    
    for (i in 1: length(x)) 
    { ti<-t[i]
      xi<-x[i]
      yi<-y[i]
      
      nei<- getnei(xi,yi,ti,array)
      
      nx22<-c(xi-1,xi,xi+1, xi-1,xi+1,xi-1,xi,xi+1,xi,xi,xi)
      
      ny22<-c(yi-1,yi-1,yi-1, yi,yi,yi+1,yi+1,yi+1,yi,yi,yi)
#      t2<-c(rep(4,8),1,2,3,5,6,7)+(ti-1)
      t2<-c(rep(0,8),-1,-2,-3)+(ti)
      ## time as the third dimension, use 14 neighbors to predict
      air3d <- as.data.frame(cbind(nei,nx22,ny22,t2))
      air3d2<-rbind(air3d2,air3d)
    }
 # }
  air3d2<-na.omit(air3d2)
  colnames(air3d2) <- c("ndvi","x","y","t")
  
# coordinates(air3d2) <- ~x+y+t
 # proj4string(air3d2) = CRS() 
  return(air3d2)
}

#
getnei2<-function(x,y,t,array)
{

 
  nxx2<-c(x-1,x,x+1)
  nyy2<-c(y-1, y,y+1)
  # t2<-c(rep(4,8),1,2,3,5,6,7)+(t-1)
  t2<-c(-3:0) +(t)

  n.1<-c(array[nxx2,nyy2,t2 ])
 
  
 # n.3<- n.1[-5] 
  return(n.1)
}
 
getcube2<-function(x,y,t,array ) 
{ 
  air3d2<-as.data.frame(matrix(NA,1,4))
  names(air3d2)<-c('nei', 'nx22','ny22', 't2')
  
 
  for (i in 1: length(x)) 
  { ti<-t[i]
    xi<-x[i]
    yi<-y[i]
    
    nei<- getnei2(xi,yi,ti,array)
    
    nx22<-c(rep(c(xi-1,xi,xi+1, xi-1,xi, xi+1,xi-1,xi,xi+1),4))
 
 #nx22<-nx22[-5]
    ny22<-c(rep(c(yi-1,yi-1,yi-1, yi,yi,yi,yi+1,yi+1,yi+1),4))
 #ny22<-ny22[-5]
    t2<-c(rep(c(-3:0),each =9 ))+(ti)
# t2<-t2[-1]
    
 air3d <- as.data.frame(cbind(nei,nx22,ny22,t2))
    air3d2<-rbind(air3d2,air3d)
  }
 
  air3d2<-na.omit(air3d2)
  colnames(air3d2) <- c("ndvi","x","y","t")
  
  # coordinates(air3d2) <- ~x+y+t
  # proj4string(air3d2) = CRS() 
  return(air3d2)
}

#xi=seq(2,96,8)
#yi=seq(2,96,8)
#ti=seq(1,600,50)
 
c(59138:59180)
c(48712:48753)
### shrink the array

#artind<-which(!is.na(ar520t2),arr.ind=TRUE)
#artind1<-which(!is.na(ar520tt),arr.ind=TRUE)
 
#ar520t23<-ar520t2[-1,-1,-1:-3]
#ar520t24<-ar520t23[-41,-39,] # get rid of the border effect by shrinking the array
#
#str(ar520tt)

ar520tt1<-ar520tt[-43,-41,]
ar520tt23<-ar520tt1[-1,-1,-1:-3]
ar520tt24<-ar520tt23[-41,-39,]

artind2<-which(!is.na(ar520tt24),arr.ind=TRUE)
xi2<-artind2[,1]+1
yi2<-artind2[,2]+1
ti2<-artind2[,3]+3
xi2<-xi2[1:10]
yi2<-yi2[1:10]
ti2<-ti2[1:10]
#
ti3<-c(4:167)
ti2<-rep(ti3,each=length(xi2))
 length(xi2)
xi2<-rep(xi2, length(ti3))
yi2<-rep(yi2, length(ti3))
#unique(coordinates(cbind(xi2,yi2))) # 691 (1071)
#array1<-getcube(xi2,yi2,ti2,arraystevid)
 
#str(array1)
#save(arraystevid,file='arraystevid')
 
#coordinates(array1) <- ~x+y+t
 #proj4string(array1) = CRS() 

#gstat3d <- gstat(formula=ndvi~1, data=array1[!is.na(array1$ndvi),])
#vgm3d <- variogram(gstat3d,cutoff=2)
#vgm3d
#plot(vgm3d)
#model3d <- fit.variogram(vgm3d,vgm(0.014, "Gau",1.7,0.4))
#model3d 
 
#plot(vgm3d, model3d)
##predict
#rt<-rep(ti+3,each=length(xi))
#rt<-rep(ti2,each=length(xi2))
#rx<-rep(xi2, length(ti2))
#ry<-rep(yi2, length(ti2))
#length(ry)

i0=1:length(xi2) #1:144
ti<-c(6:15)
#orig<-lapply(i0, function(i0) test.dst2[rx[i0],ry[i0],rt[i0]])
 
orig<-lapply(i0, function(i0) arraystevid[xi2[i0],yi2[i0],ti2[i0]])
orig<-(as.numeric(orig))
#ngrid<-as.data.frame(cbind(rep(orig,each=14),rx,ry,rep(ti,each=length(xi))))
#colnames(ngrid) <- c("ndvi","x","y","t")
#ngrid<-as.data.frame(cbind((orig),xi2,yi2,rep(ti2,each=length(xi2))))
length(  orig)
xi4<-xi2+ 59137
yi4<-yi2 +48711 
xyi2<-as.data.frame(cbind(xi4,yi4))
xyi2m<-getxyMatrix(xyi2,231.6564)
xi3<- xyi2m[,1] 
yi3<- xyi2m[,2] 
#ti3<- monmean1[ti2-3]
ti3<-monmean1[-1:-3]
length(xi3)
#monmean1<-time(monmean)
ti4<- rep(monmean1[-1:-3],each=length(xi3))
ngrid<-as.data.frame(cbind(orig,xi3,yi3, ti4 ))
ngrid2<-as.data.frame(cbind(xi3,yi3 ))
 str(ngrid2)
colnames(ngrid) <- c("ndvi","x","y","t")
colnames(ngrid2) <- c("x","y")
coordinates(ngrid2)=~x+y

str(grid2)
gridded(ngrid2)<-TRUE
grid2<-STF(sp=as(ngrid,'SpatialPoints'),time=ti3)
str(grid2)
grid2[1,164]

#+t
save(ngrid,file='ngrid.Rdata')
save(grid2,file='grid2.Rdata')
proj4string(ngrid) =  CRS()

plot(ngrid)
################################################# 

cosp2<-cbind(0,0,0)
cosp3<-cbind(0,0,0)
cosp4<-cbind(0,0,0)
cosp5<-cbind(0,0,0)
cosp6<-cbind(0,0,0)
#for(i in 1:8)#{ # for (j in 1:8 )  #{
   for (k in 1:length(xi2))
   {
 
 
     array1<-getcube(xi2[k],yi2[k],ti2[k],arraystevid) 
  
coordinates(array1) <- ~x+y+t
 
proj4string(array1) = CRS()
v = function(x, y = x) { exp(-(spDists(coordinates(x),coordinates(y))/0.34)^2)}
k.1<- krige0(ndvi~1, data=array1[!is.na(array1$ndvi),],ngrid[g,], v, computeVar = TRUE)

crit<-qnorm(c(0.025,0.975),mean=k.1$pred, sd=sqrt(k.1$var))
#print(k.1$var)
#print(paste('kriging variance:', k.1$var))
zscore<-(ngrid[g,]$ndvi-k.1$pred )/sqrt(k.1$var)

#print(paste('zscore:',zscore))
#print(paste('critical value:',crit))
pval<- 2*pnorm( - abs(zscore))
 
  if(pval<0.1)
{
  print(paste('p-value:', pval,'coords',xi2[k],',',yi2[k],',',ti2[k]))
  cosp<-cbind(xi2[k],yi2[k],ti2[k])
  cosp3<-rbind(cosp3,cosp)
  }
if(pval<0.05)
{
 # print(paste('p-value:', pval,'coords',xi2[k],',',yi2[k],',',ti2[k]))
  cosp<-cbind(xi2[k],yi2[k],ti2[k])
  cosp4<-rbind(cosp4,cosp)
}
if(pval<0.2)
{
 # print(paste('p-value:', pval,'coords',xi2[k],',',yi2[k],',',ti2[k]))
  cosp<-cbind(xi2[k],yi2[k],ti2[k])
  cosp5<-rbind(cosp5,cosp)
}


#hx<-rnorm(30,mean=k.1$pred, sd=sqrt(k.1$var ))
#hist(seq(crit[1],crit[2],length=30),hx,breaks=16)
 
#xlim=c(min(crit[1],array1@data$ndvi)-0.05,max(crit[2],array1@data$ndvi+0.05))
#abline(v=crit[1],col='red')
#abline(v=crit[2],col='red')
g=g+1
   }
save(cosp5,file='cosp5')  #p-value is less than 0.2 2193 points 101
save(cosp4,file='cosp4')  #0.05 1374points  81
save(cosp3,file='cosp3')  # 0.1  1572 points     82    
#unique(coordinates(cbind(xi2,yi2))) # #bfast   691 1071



cosp9<-cbind(0,0,0)
cosp10<-cbind(0,0,0)
cosp11<-cbind(0,0,0)  


cosp12<-cbind(0,0,0)
cosp13<-cbind(0,0,0)
cosp14<-cbind(0,0,0)
error<-cbind(0,0,0)  
x<-cosp11[,1]
y<-cosp11[,2]
z<-cosp11[,3]
 
length(
  unique(
  coordinates(cbind(x,y,z))
  )
  ) 
/3
length(cosp9)
#1885
save(cosp6,file='cosp6')  #p-value is less than 0.1 64 points  53
save(cosp7,file='cosp7')  #0.05                     55 points  48
save(cosp8,file='cosp8')  #0,2                      72 poionts  59
# wrong 
# before filtering use 35 neighbors fitsummet
save(cosp9,file='cosp9')  #p-value is less than 0.1 6857 points  281
save(cosp10,file='cosp10')  #0.05                    6158 points  173
save(cosp11,file='cosp11')  #0,2                      8488 poionts  512

#after filtering
save(cosp12,file='cosp12')  #p-value is less than 0.1 1885 points  256
save(cosp13,file='cosp13')  # wrong                    1885

save(cosp14,file='cosp14')   #6068/504
unique(coordinates(cbind(cosp9[,1],cosp9[,2])))  # 101
cosp3<-cosp2[-1,]
insert.at <- function(a, pos, ...){
  dots <- list(...)
  stopifnot(length(dots)==length(pos))
  result <- vector("list",2*length(pos)+1)
  result[c(TRUE,FALSE)] <- split(a, cumsum(seq_along(a) %in% (pos+1)))
  result[c(FALSE,TRUE)] <- dots
  unlist(result)
}

 
g=1
for (k in 1:length(xi2))
{
 
  array1<-getcube2(xi2[k],yi2[k],ti2[k],farraystevi) 
  evi<-array1$ndvi
 eviorig<- evi[31]
  #evi<-insert.at(evi,4,NA)
 evi[31]<-NA
  ti2m<-monmean1[unique(array1$t)]
  tiend<-monmean1[unique(array1$t)+1]
  xi<-xi2[k]
  yi<-yi2[k]
  nx22<-c(xi-1,xi,xi+1, xi-1,xi, xi+1,xi-1,xi,xi+1)
 
  ny22<- c(yi-1,yi-1,yi-1, yi,yi,yi,yi+1,yi+1,yi+1) 
  xi4<-nx22+ 59137
  yi4<-ny22 +48711 
  xyi4<-as.data.frame(cbind(xi4,yi4))
  xyi2m<-getxyMatrix(xyi4,231.6564)
  xi5<-xyi2m[,1]
  yi5<-xyi2m[,2]
  ti3<- monmean1[ti2]

  array2<-as.data.frame(cbind(xi5,yi5))
#  ti2mo<-ti2m[order(ti2m)]
  tiend<-tiend[order(tiend)]
  #array1<-getcube(xi,yi,ti) 
  coordinates(array2) <- ~xi5+yi5
  gridded(array2)<-TRUE
 
  array3<-STFDF(sp=as(array2,'SpatialPoints'),time=ti2m,as.data.frame(evi))
 
 proj4string(array3)<-CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
proj4string(grid2)<-CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
  array4<-as( array3,'STSDF') 
 #endTime = tiend
tii<-(ti2[k]-3)

  kg2<-krigeST(formula=evi~1, array4,grid2[k,tii],fitSumMet,computeVar=TRUE)
if(kg2@data$var1.var>0.0001)
  {
#kg2@data$var1.pred
#eviorig
#array4@data$evi
 
zscore<-(eviorig-kg2@data$var1.pred )/sqrt(kg2@data$var1.var)

#print(paste('zscore:',zscore))
 
pval<- 2*pnorm( - abs(zscore))

if(pval<0.1)
{
  print(paste('p-value:', pval,'coords',xi2[k],',',yi2[k],',',ti2[k]))
  cosp<-cbind(xi2[k],yi2[k],ti2[k])
  cosp12<-rbind(cosp12,cosp)
}
if(pval<0.05)
{
  # print(paste('p-value:', pval,'coords',xi2[k],',',yi2[k],',',ti2[k]))
  cosp<-cbind(xi2[k],yi2[k],ti2[k])
  cosp13<-rbind(cosp13,cosp)
}
if(pval<0.2)
{
  # print(paste('p-value:', pval,'coords',xi2[k],',',yi2[k],',',ti2[k]))
  cosp<-cbind(xi2[k],yi2[k],ti2[k])
  cosp14<-rbind(cosp14,cosp)
}
}
else 
  {
    error1<-cbind(xi2[k],yi2[k],ti2[k])
    error<-rbind(error,error1)
print('error')
print(k)
print(error)
}
g=g+1
}

#  kg3<-krigeST(formula=evi4042~1, stfdf2,grid2,fitSumMet,nmax=20,computeVar=TRUE)
#how to select neighbors for estimation
#variance not computed
#grid2:   STF   all the points ST points to estimate ( where the change happened)
#stfdf2 : STDFT the whole array 40*42*167
#Error: cannot allocate vector of size 586.5 Gb 
#select a small cube and estimate one point:
#array3 : STDFT ndvi cube with t t-1 t-2 t-3 and the spatial neighbours 
#(36, the one to estimate has been set to NA) 
# if evi[5]<- NA  # how to do if only STFDF can be used
#Error in cbind(v0, X) : number of rows of matrices must match (see arg 2)
vstbr1 <-variogramST(evi ~1,array3)
#Error in `[.xts`(x@time, j) : subscript out of bounds
str(array3)
#proj4string(array3)<-CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

#crit<-qnorm(c(0.025,0.975),mean=k.1$pred, sd=sqrt(k.1$var))
#print(k.1$var)
#print(paste('kriging variance:', k.1$var))



#hx<-rnorm(30,mean=k.1$pred, sd=sqrt(k.1$var ))
#hist(seq(crit[1],crit[2],length=30),hx,breaks=16)

#xlim=c(min(crit[1],array1@data$ndvi)-0.05,max(crit[2],array1@data$ndvi+0.05))
#abline(v=crit[1],col='red')
#abline(v=crit[2],col='red')


 
 
#k.1$var
#k.1$pred

# k.2@data$var1.var
 

#k.2 <- krige(ndvi~1,  array1[!is.na(array1$ndvi),], ngrid,model=model3d)



k.2@data$var1.pred-k.1$pred # almost the same!
 
scale(k.1$pred)
#scatterplot3d(rx,  rt,k.1, main="Predicted values")
 
#ngrid <- as(ngrid, "SpatialPixelsDataFrame")

# variogram
gstat3d <- gstat(formula=ts~1, data=t1ts[!is.na(t1ts$ts),])
vgm3d <- variogram(gstat3d,cutoff=30)

## cross validation
#cv<-krige.cv(ndvi~1, array1[!is.na(array1$ndvi),],vgm(0.04, "Gau",22, 0.03 ),nfold=5,verbose=FALSE)
#bubble(cv)


save(k.1,file='k.1')
