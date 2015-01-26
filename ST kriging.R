itry2<-itry[-43,-41,-636:-637] # 3 dimensional x y t
save(itry2,file='itry2.Rdata')
 
load('monmean1.Rdata')
table(itry)
#itrya<-aperm (itry2,c(2,1,3))
bfastrevi<-dst4240agg
dst4240stfdf<-getstfdf(dst4240agg)
originalmonthly4042<-getstfdf(itry2monarr)
str(itry2monarr)
dst4240stfdf1<-getstfdf(dst4240)
length(which(abs(dst4240)<1))
getstfdf<-function(arr)
  {
#itrydf<-as.data.frame.table(bfastrevi) #y x t interate y ->x ->t: (x1,y1,t1)(x1,y2,t1)(x2,y1,t1)(2,2,t1)(1,1,t2)(1,2,t2)
  itrydf<-as.data.frame.table(arr)
  evi4042<-itrydf$Freq # bfast residuals 
evi4042<-as.data.frame(evi4042)
#evi4042<-as.data.frame(s2@data$evi4042)
names(evi4042)<-'evi4042'
#length(evi4042$evi4042)

#s3<-spTransform(s2,CRS("+proj=utm +zone=21 +south"))
#40*42
x<-c(59139:59180)
y<-c(48713:48752)
x1<-rep(x,each=length(y))
y1<-rep(y,length(x))
xyd<-as.data.frame(cbind(x1,y1))
xyc<-getxyMatrix(xyd,231.6564)
xyc<-as.data.frame(xyc)
coordinates(xyc)<-~x+y
SpatialPoints(xyc)
proj4string(xyc)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
#  xycs<-spTransform(xyc,CRS("+proj=utm +zone=21 +south  +ellps=WGS84 "))
stfdf<-STFDF( xyc,a1,evi4042) 
# stfdf<-STFDF( xyc,monmean1,evi4042) #bfast
return(stfdf)
}

#bfastre510stfdf<-STFDF( xyc,monmean1,evi4042) #bfast
 
arraystevid<-as.ndarray.STFDF(dst4240stfdf1,'evi4042')
 
farraystevi<-filter.st.median(arraystevid)
dst4240stfdf<-arraytoSTFDF(farraystevi)

length(which(dst4240stfdf1@data-dst4240stfdf@data!=0)) #5817

arraystevid<-farraystevi # saved filtered bfast residuals

 
arraystevid<-as.ndarray.STFDF(bfastre510stfdf,'evi4042') # array bfast residuals
farraystevi<-filter.st.median(arraystevid)
arraystevid<-farraystevi # filtered bfast residuals

#save(bfastre510stfdf,file='bfastre510stfdf.Rdata')

load('arraystevid')
arraystevid<-bfastre510stfdf # bfast residuals around point 520
save(bfastre510stfdf,file='bfastre510stfdf')
#40*42*167
## STFDF to array
as.ndarray.STFDF = function(x3,var1) {
 
  dfstevi<-as.data.frame(x3)
 
  x2<-spTransform(x3,CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")) 
  cr<-getcrMatrix(as.data.frame(x2@sp@coords),231.6564)
  
  x<-as.character(unique(cr[,1]))
  y<-as.character(unique(cr[,2]))
 
  a = array(NA, c(length(x),length(y),length(x3@time)))

  a[,,] =dfstevi[var1][,]
  
  dimnames(a) = list(x,y, make.names(index(x3@time)))
  return(a)
}

#arraystevi<-as.ndarray.STFDF(stfdf1,'evi4042')

###################################################### 
setMethod('brick', signature(x='array'), 
          function(x, xmn=0, xmx=1, ymn=0, ymx=1, crs=NA, transpose=FALSE) {
            dm <- dim(x)
            if (is.matrix(x)) {
              stop('cannot coerce a matrix to a RasterBrick')
            }
            if (length(dm) != 3) {
              stop('array has wrong number of dimensions (needs to be 3)')
            }
            b <- brick(xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, crs=crs, nl=dm[3])
            names(b) <- dimnames(x)[[3]]
            
            if (transpose) {
              dim(b) <- c(dm[2], dm[1], dm[3])
            } else {
              dim(b) <- dm
              # aperm etc suggested by Justin McGrath
              # https://r-forge.r-project.org/forum/message.php?msg_id=4312
              x = aperm(x, perm=c(2,1,3))
            }
            attributes(x) <- list()
            dim(x) <- c(dm[1] * dm[2], dm[3])
            setValues(b, x)
          }
)

arraytoSTFDF<-function(array  ) # array with modis index, time as attributes
{
  itrydf<-as.data.frame.table(array) #y x t interate y ->x ->t: (x1,y1,t1)(x1,y2,t1)(x2,y1,t1)(2,2,t1)(1,1,t2)(1,2,t2)
  aa2<-itrydf$Freq
  aa3<-as.data.frame(aa2)
  
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
 names( data2)<-'evi4042'
 coordinates(xyc)<-~x+y
  
  SpatialPoints(xyc)
  
  proj4string(xyc)<-'+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs'
   
  stfdf1<-STFDF(xyc,tt2,data2) 
  return(stfdf1)
}

filter.st.median<-function(a) #if evi is more than 1, replace with median
{
  
  index1<-which(abs(a)>1,arr.ind=TRUE)
  
  
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
    
    
    a[i,j,t]<-median(c(
      a[ih,j,t],a[il,j,t],a[i,jh,t],a[i,jh,t],
      a[ih,j,tl],a[il,j,tl],a[i,jh,tl],a[i,jh,tl],
      a[ih,j,th],a[il,j,th],a[i,jh,th],a[i,jh,th]))
    
  }
  return(a)   
}

 
load('bra520stfdf')
plot(stfdf2@sp,col='green')
setwd("C:/Users/m_lu0002/Desktop/Climate/minnesota")
arraystevid<-as.ndarray.STFDF(bfastre510stfdf,'evi4042') # array bfast residuals
arraystevid<-as.ndarray.STFDF(stfdf2,'evi4042')
arraystevid<-as.ndarray.STFDF(dst4240stfdf,'evi4042')
farraystevi<-filter.st.median(arraystevid)
fity2<-filter.st.median(itry2) #filtered 
 save(fity2,file='fity2.Rdata')
 
arraystevid<-farraystevi # filtered bfast residuals
#attributes(ar520tt1)<-attributes(arraystevid)
str(arraystevid)
dst4240stfdf1#filtered detrended stdft

stfdf3<-arraytoSTFDF(farraystevi)
#stfdfchange<-arraytoSTFDF(ar520tt1)
#geometry(stfdfchange)
summary(stfdf2@data)
kg2<-gstat(id=evi4042k,formula=evi4042~1, stfdf2,grid2,fitSumMet,nmax=9,maxdist=5)
,computeVar=TRUE)
# how to select neighbor
#variance not computed
str(grid2)
kg1@data$var1.pred

proj4string(grid2)<-CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
gridded(ngrid) <- TRUE
 
save(stfdf2,file='stfdf2')

load('stfdf2')
length(which(stfdf2@data-stfdf3@data!=0)) #monthly 1738
vstbr1 <-variogramST(evi4042~1,stfdf2,tlags=0:5,cutoff=2000) # not filtered bfast residuals
vstbr1.1 <-variogramST(evi4042~1,stfdf3,tlags=0:5,cutoff=2000) # filtered bfast residuals

vstbr2 <-variogramST(evi4042~1,dst4240stfdf,tlags=0:20,cutoff=2000) #filtered detrended and deseasonalized


arraystevid<-as.ndarray.STFDF(stfdf3,'evi4042')
farraystevi<-filter.st.median(arraystevid)


plot(vstbr2,wireframe=T)
save(vstbr2,file='vsbr2.Rdata')
load('vsbr2.Rdata')
#par(mfrow=c(1,1))
plot(vstbr1.1,wireframe=T,main='filtered 8-day')

sepa1<-vgmST('separable',space=vgm(0.008,'Exp',2000,0.002),time=vgm(.008,'Exp',20,0.002),sill=0.01)
fit1<-fit.StVariogram( vstbr1.1 , sepa1)
attr(fit1, "optim")$value
save(fit1,file='fit1.Rdata') # best fit model #0.0001660289
plot( vstbr1.1,fit1,wireframe=T, both=T)


ps1 <- vgmST("productSum",
             space=vgm(0.005,'Exp',2000 ,0),
             time=vgm(0.0002,'Exp',20,0),
             sill=15,
             nugget=10)
SumMet1 <- vgmST("sumMetric",
             space=vgm(0.008,'Exp',2000,0.003),
             time=vgm(0.008,'Sph',20,0.003),
             joint=vgm(0.002, "Exp", 10,0.003),
             stAni=0.2)

fitSumMetmf1 <- fit.StVariogram( vstbr1.1, SumMet1, lower=c(0,0.1,0,
                                                            0,0.1,0,
                                                            0,0.1,0,0.01))
attr(fitSumMetmf1, "optim")$value # Exp+Exp+Exp: 0.3623, Exp+Sph+Exp: 0.3602
plot(vstbr1.1,fitSumMetmf1,wireframe=T,all=T, scales=list(arrows=F),main='filtered monthly summet')
fitSumMetmf1

plot(vstbr1,SumMet1,wireframe=T)
plot( vstbr1.1,sepa1,wireframe=T, all=T)
fit1<-fit.StVariogram( vstbr1.1 , sepa1)
fitPs1<-fit.StVariogram( vstbr1.1 , ps1, lower=c(0.1,0.1,0.1,0.1,0.1,0.1))

# scale
#vstbr1.1$dist <- vstbr1.1$dist/100
#vstbr1.1$spacelag <- vstbr1.1$spacelag/100

# fit
fitSumMetmf1 <- fit.StVariogram( vstbr1.1, SumMet1, lower=c(0,0.1,0,
                                                        0,0.1,0,
                                                        0,0.1,0,0.01))
save(fitSumMet8day1,file='fitSumMet8day1.Rdata')
save(fitSumMetmf1,file='fitSumMetmf1.Rdata')
#,method = "L-BFGS-B",lower = c(10,0,0.01,0,1),upper = c(500,1,20,1,200))

attr(fit1, "optim")$value
attr(fitPs1, "optim")$value
attr(fitSumMetmf1, "optim")$value # Exp+Exp+Exp: 0.3623, Exp+Sph+Exp: 0.3602
 
plot(vstbr1.1,fit1,wireframe=T,all=T, scales=list(arrows=F),main='filtered 8-day separable')
plot(vstbr2,fitPs1,wireframe=T,all=T, scales=list(arrows=F))
plot(vstbr1.1,fitSumMetmf1,wireframe=T,all=T, scales=list(arrows=F),main='filtered 8-day summet')
fitSumMetmf1
# rescale
vstbr1$dist <- vstbr1$dist*100
vstbr1$spacelag <- vstbr1$spacelag*100

fitSumMet$space$range <- fitSumMet$space$range*100
fitSumMet$joint$range <- fitSumMet$joint$range*100
fitSumMet$stAni <- fitSumMet$stAni*100 

plot(vstbr1,fitSumMet,wireframe=T,both=T, scales=list(arrows=F))
plot(vstbr1,fitSumMet)

krige(log(zinc)~1, meuse, meuse.grid, model = m)
spplot(x["var1.pred"], main = "ordinary kriging predictions")

plot(vstbr1,wireframe=T)

str(stfdf2)
farraystevi[order(farraystevi)][1:10]

stfdf2@sp
plot(vst3, wireframe=T)
plot(vst, wireframe=T)

plot(vst2, wireframe=T)
getwd()
#save(farraystevi,file='farraystevi.Rdata')
#save(stfdf1,file='stfdf1.Rdata')
load('farraystevi.Rdata')
itry5<-arraystevi[-order(arraystevi)[1:52]]
str(itry5)
itry6<-itry5[-order(itry5, decreasing=TRUE)[1:30]]   
plot(arraystevi[1,3,])

dimnames(arraystevi)[3]
 
as.Date(tt1,format='%Y.%m.%d')



tt1[[1]] 
tt1
str(itry2)
dfstevi<-as.data.frame(stevi)
as(dfstevi, 'STFDF')
str(dfstevi)
order(stevi[,])
order(stevi[,]
      [])
stevi[cbind(1:3),c(1,2,1),]
#get rid of mal value
itry5<-evi4042[-order(evi4042)[1:52]]
itry6<-itry5[-order(itry5, decreasing=TRUE)[1:30]]   

min(itry6)
str(itry3) 
 # itry4<- ( rm.outlier(itry2,fill=TRUE,median=TRUE))
  max(itry4)
  str(itry4) 
  #itry3<-itry2[-order(itry2)[1:52]]
#itry4<-itry3[-order(itry3, decreasing=TRUE)[1:30]]   
str(itry3)
 
 
itry2 [9,7,6]
 
   evi_a53<-evi_a52[order(evi_a52)]

order(evi_a52)[1:52]
str(evi_a52)
20294/(43*41)
711537/(43*41)
evi_a52[20294]
evi_a52[41,21,12]
evi52re[8,7,5][]
#7683 data
#319409 271350 711537 893195   6754 899939 249052 249051 583789
evi_a52[41:50]
evi_a52[5,8,182]
evi_a52[15,24,404]
evi52re[15,24,404][]
str(evi_a52)
tail(evi_a53,n=100)  
evi_a52[1,2,1]
load('evi52.Rdata')
which.min(evi_a52)
evi_a52[973851] > -2
< (-1)
as.numeric(evi_a52)
length(evi_a52)
lapply(evi_a52, as.numeric)
str(evi_a52)

evt2<-as.data.frame(as.data.frame.table(evt))

head(evt2,n=100000)
evt2$Freq
           [21000:21059])
slot(xycs,'coords')
str(evt2)

evi_a52[1,2,1:5]
head(evi_a52[,,1])
order(evt2$Freq)[1]
summary(evt)
 stplot(stevi)
head(evt2$Freq[order(evt2$Freq)])

tail(evt2$Freq)
variogram(PM10~1, stevi)
         
str(stevi)
nrow(evt2)
a2<-a1[39:635]
a2
str(evi_a52)
length(a1)
length(xycs@coords)*length(monmean3)
evt2
as.data.frame.table(evi_a52[,,1])

str(evi_a52)
 str(evi_a52)

#slogo <- stack(system.file("external/rlogo.grd", package="raster")) 
#plot(slogo,3) #77 101

