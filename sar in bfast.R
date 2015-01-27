install_github("bfast2","mengluchu")
ti<-1:636
install_github("strucchange","mengluchu",build_vignettes = FALSE)
library(devtools)
library(strucchange)
library(bfast)


setwd("C:/Users/m_lu0002/Dropbox/mengluchu/bfast")
list.files()
load("a1.saved")
load("fevi8.Rdata")
save(fevi20b20,file="fevi20b20.Rdata")
fevi20b20<-fevi8[1:20,1:20,]
save(fevi20b20,file="fevi20b20.Rdata")
#fevi20b20<-fevi8[20:50,20:50,] # input example array
#fevi220b220<-fevi8[1:221,1:221,]
fevi9b9<-fevi8[1:9,1:9,]
fevi3b3<-fevi8[1:3,1:3,]
#fevi3b4<-fevi8[1:3,1:4,]

#output1<-bfmarray(fevi20b20,dates=a1,aggre='month',start=9) # bfast monitor
########## try spatial CAR bfast

#get STF

#array<-fevi3b3
#itrydf<-as.data.frame.table(array) #y x t interate y ->x ->t: (x1,y1,t1)(x1,y2,t1)(x2,y1,t1)(2,2,t1)(1,1,t2)(1,2,t2)
dim(fevi8)
# deseasonality
#for (i in 1:3)
#{
#  for ( j in 1:3)

#plot(bfast(ts(fevi3b3[i,j,], start=c(2000,1),frequency=46),h=0.15,max.iter=1,type = "OLS-CUSUM",) )
#}
eday <- as.Date("2000-01-30")           # date 
e8day <- seq(eday, length.out=636, by="8 days")

x<-c(1:3)
y<-c(1:3)
x1<-rep(x,length(y))
y1<-rep(y,each=length(x))
xyd<-as.data.frame(cbind(x1,y1))
coordinates(xyd)<-~x1+y1
 
 sarneighbor<-function(xi,yi) #xi by yi
{
eday <- as.Date("2000-01-30")           # date 
e8day <- seq(eday, length.out=636, by="8 days")
xyd<-expand.grid(x1=1:xi,y1=1:yi)
coordinates(xyd)<-~x1+y1
lecube<-xi*yi*636
aa3<-as.data.frame(c(1:lecube))
stfdf3b3<-STFDF(xyd,e8day,aa3) ## for creating neighbors only, aa3 could be any data?
cn<-cell2nb(xi,yi, type ="queen",torus =FALSE)
neigh1<-nbMult(cn, stfdf3b3, addT = FALSE, addST = FALSE) # only spatial neighbours are added for each time step
listcn636<-nb2listw(neigh1)
return(listcn636)
}

list2020<-sarneighbor(20,20)
str(list2020)
 
 
fevi20b20des<-apply(fevi20b20,c(1,2),function(x) (ts(x,start=c(2000,1),frequency=46)-stl(ts(x,start=c(2000,1),frequency=46),'per')$time.series[,"seasonal"]))
f2<-aperm(fevi20b20des,c(2,3,1))

aa2<-as.vector(f2)
tcm<-proc.time()
try2<-spautolm(aa2~c(rep(c(1:636),each=400)),family="SAR",method= "Matrix", listw=list2020)
proc.time()-tcm
#try3<-spautolm(aa2~c(rep(c(1:636),9)),family="SAR",method= "Matrix", listw=listcn636)  not exactly the same but less different than i thought..
#try2<-spautolm(aa2~c(rep(c(1:636),each=400)),family="CAR",method= "Matrix",listw=list2020) #a bit different than SAR
#so much faster than dense array! 1.3 s  22.06s #16.74    4.58   22.06
 
rn<-lapply(1:400,function(i) {residuals(try2)[seq(i,636*400-(400-i),400)]})
#get residuals for each time series
#rn
#ii<-5   # get the middle pixel (5 for 3*3 matrix)
ti<-1:636
p1<-array (,c(20 ,20 ))
p2<-array (,c(20 ,20 ))
p3<-array(,c(20,20 ))

for(i in 1:20)
  {
  for (j in 1:20)
  {
    ii<- (i-1)*20+j
fevi3b312t1<-ts(f2[i,j,],start=c(2000,1),frequency=46) # reconstruct the time series
p.Vt2 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM" )) 
p.Vt3 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-MOSUM" ))
p.Vt  <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]]))  )

p1[i,j]<- p.Vt$p.value # spautolm residuals
p2[i,j]<- p.Vt2$p.value # linear regression residuals
p3[i,j]<- p.Vt3$p.value
print(ii)
}
}
rn26<-rn[[22]]
cor(rn25,rn26)
f2[2,4,]<-fevi8[120,150,] 
f2[3,4,]<-fevi8[120,150,]
f2[1,4,]<-fevi8[120,150,]  
f2[4,1,]<-fevi8[120,150,]
f2[4,2,]<-fevi8[120,150,]
f2[4,3,]<-fevi8[120,150,] #is sar local? cor(rn24,rn22) #[1] 0.9998499 #is car local? cor(rn25,rn26) 1
 
#sar vs. car cor(rn25,rn22)[1] 0.9889769
save(rn22,file='rn22.Rdata')
length(which(p3<0.05))
stfdf3b3<-STFDF(xyd,e8day, aa3)  
cn<-cell2nb(3,3, type ="queen",torus =FALSE)
cn1<-nb2listw(cn)
neigh1<-nbMult(cn, stfdf3b3, addT = FALSE, addST = FALSE) # only spatial neighbours are added for eath time replicate
listcn636<-nb2listw(neigh1)
save(listcn636,file="listcn636.Rdata")
# get neighbor
p1<-array (,c(150 ,150 ))
p2<-array (,c(150 ,150 ))
p3<-array(,c(150,150 ))
for (i in 3:148 )
{
  for (j in 1:148 )
{

    fevi3b3<-fevi8[i:(i+2),j:(j+2),]
fevi3b312t<-apply(fevi3b3,c(1,2),function(x) (ts(x,start=c(2000,1),frequency=46)-stl(ts(x,start=c(2000,1),frequency=46),'per')$time.series[,"seasonal"]))
#dim(fevi3b312t)
f2<-aperm(fevi3b312t,c(2,3,1))
aa2<-as.data.frame.table(f2)$Freq
aa3<-as.data.frame(as.data.frame.table(f2)$Freq)
eday <- as.Date("2000-01-30")           # date 
e8day <- seq(eday, length.out=636, by="8 days")
 

try2<-spautolm(aa2~c(rep(c(1:636),each=9)),family="SAR",listw=listcn636)
#try1<-spautolm(aa2~c(rep(c(1:636),each=12)),family="CAR",listw=listcn636)
rn<-lapply(1:9,function(i) {residuals(try2)[seq(i,636*9-(9-i),9)]})
#get residuals for each time series


 
#p3<-array(,c(3,3))
#dim(f2) #deseasonalized 
    ii<-5    
    fevi3b312t1<-ts(f2[2,2,],start=c(2000,1),frequency=46)
    p.Vt2 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM" )) 
    p.Vt3 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-MOSUM" ))
    p.Vt <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]]))  )
    
    p1[i,j]<- p.Vt$p.value # spautolm residuals
    p2[i,j]<- p.Vt2$p.value # linear regression residuals
    p3[i,j]<- p.Vt3$p.value
    #p3[i,j]<-bfast(fevi3b312t,max.iter=1,level=0.05,h=0.15, type = "OLS-CUSUM")$output[[1]]$Vt.bp[1]
  }
}
p1
save(p1,file='100b100sar.Rdata')
save(p2,file='100b100cusum.Rdata')
save(p3,file='100b100mosum.Rdata')
length(which(p3<0.05))
summary(p1)
load('listcn636.Rdata')

summary(p2)
which(p1<0.05)
dev.new()
par(mfrow=c(1,1))
par()
for (j in 1:3) 
{
  for(i in 1:3)  
  { 
    ii<-i+(j-1)*3
    
    fevi3b312t<-ts(f2[i,j,],start=c(2000,1),frequency=46)
    p.Vt2 <- sctest(efp(fevi3b312t ~ ti, h = 0.15, type = "OLS-CUSUM" )) 
    p.Vt <- sctest(efp(fevi3b312t ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]]))  )
 
    p1[i,j]<- p.Vt$p.value # spautolm residuals
    p2[i,j]<-p.Vt2$p.value # linear regression residuals
    
    #p3[i,j]<-bfast(fevi3b312t,max.iter=1,level=0.05,h=0.15, type = "OLS-CUSUM")$output[[1]]$Vt.bp[1]
  }
}

#}
plot(f2[1,1,],typ='l')
p1
p2
#p1: two more breakpoints: i=1:2 j=1
i=2; j=2
ii<-i+(j-1)*2
fevi3b312<-ts(f2[i,j,],start=c(2000,1),frequency=46)
plot(stl(fevi3b312,"per"))
plot(efp(fevi3b312 ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]])))
plot(efp(fevi3b312 ~ ti, h = 0.15, type = "OLS-CUSUM"))
str(rn)

plot(rn[[5]],typ='l')
plot(residuals(lm(aa2[seq(i,636*9-(9-i),9)]~c(1:636))),typ='l')
i=5
plot(aa2[seq(i,636*9-(9-i),9)],typ='l')
cor(as.numeric(aa2), as.numeric(rn[[5]]))
seq(i,636*9-(9-i),9)]
# rec-cusum cusum mosum most breakpoints ->least breakpoints
#plot(bfast(fevi3b312,max.iter=1,level=0.05,h=0.15))
#p1 spatial
#p2 original
#tsa1<-ts(aa2[(seq(2,636*9-7,9))],start=c(2000,8),frequency=46,end=c(2013,45))
#remainder<- stl(tsa1,"period")$time.series[,"remainder"]
#plot(residuals(lm(remainder~c(rep(1:636,1)))),typ="l") # seasonality is still here ?
#plot(aa2[(seq(1,636*9-8,9))],typ="l")
#summary(residuals(tryll))
#plot(aa2[(seq(1,636*9-8,9))],typ="l")

#a<-seq(1, 636*9-8,9)
#b<-seq(2, 636*9-7,9)
#c<-seq(3, 636*9-6,9)
#d<-seq(4, 636*9-5,9)

#plot(fitted(try1)[a],typ="l")
#plot(fitted(try1)[b],typ="l")
#plot(fitted(try1)[c],typ="l")
#plot(fitted(try1)[d],typ="l")


#r1<-residuals(try1)[a]
#r2<-residuals(try1)[b]
#r3<-residuals(try1)[c]
#r4<-residuals(try1)[d]
#r5<-residuals(try1)[seq(5,636*9-4,9)]
#r6<-residuals(try1)[seq(6,636*9-3,9)]
#r7<-residuals(try1)[seq(7,636*9-2,9)]
#r8<-residuals(try1)[seq(8,636*9-1,9)]
#r9<-residuals(try1)[seq(9,636*9,9)]
#rr<-list(r1,r2,r3,r4,r5,r6,r7,r8,r9)



#fevi3b312<-ts(fevi3b3[1,2,],start=c(2000,1),frequency=46)

#p3<-bfast(fevi3b312,max.iter=1,level=0.05,h=0.15)$output[[1]]$bp.Vt

#edit(bfast)
#fevi3b3<-fevi8[40:42,40:42,]
#ti<-1:636




str(pvalue.efp(fevi3b312 ~ ti, h = 0.15,lim.process="Brownian motion",alt.boundary=FALSE)
    
    
    help(sctest)
    edit(efp)
    efp()
    ## convert time to time dimension
    
    mtime<-monthofyears(output1[1,,]) #timearr
    
    dimx<-dim(output1)[2]
    dimy<-dim(output1)[3]
    dimt<-max(mtime[!is.na(mtime)])
    
    newarray<-array(,c(dimx,dimy,dimt)) #newarr
    
    mag<- output1[2,,]   # variablearr (magnitude)
    
    #load("C:/Users/m_lu0002/Dropbox/mengluchu/bfast/app1.R")
    
    t3darr<-timeasdim(newarray,mtime,mag) # 3 dimensional array with magnitude 
    
    which(!is.na(t3darr),arr.ind=TRUE)
    image(t3darr[,,176]) #test: bfastmonitor magnitude
    
    #test bfast: save only time
    t1<-fevi20b20
    efp()
    dimt=23*12
    dimx<-dim(output3)[1]
    dimy<-dim(output3)[2]
    
    newarray<-array(,c(dimx,dimy,dimt)) #newarr
    
    mag<- array(100,c(dimx,dimy))
    
    t3dbfaarr<-timeasdimbfa(newarray,output3,mag)
    which(!is.na(t3dbfaarr),arr.ind=TRUE)
    image(t3dbfaarr[,,126])
    ########################################### bfast with save other values
    p<-proc.time()
    
    
    t1<-fevi20b20[1:21,1:21,]
    output2 <-bfaarray2(t1,dates=a1,aggre='month',season="harmonic",max.iter=1,level=0.05)
    
    str(output2)
    allbreakpoints<-output2[1:6,,] # breakpoint
    allmagnitudes<-output2[7:12,,] # magnitude
    #which(allbreakpoints!=0)
    mtime<-allbreakpoints #timearr
    
    dimx<-dim(output2)[2]
    dimy<-dim(output2)[3]
    dimt<-max(mtime[ mtime!=0])
    
    newarray<-array(,c(dimx,dimy,dimt)) #newarr
    
    mag<- allmagnitudes  # variablearr (magnitude)
    
    t3darrbfamul<-tasd.bfa.mul(newarray,mtime,mag) # 3 dimensional array with magnitude 
    proc.time()-p
    #t3darrbfamul[which(is.na(t3darrbfamul))]<-0
    
    which(!is.na(t3darrbfamul),arr.ind=TRUE)
    image(t3darrbfamul[,,126]) #test: bfastmonitor magnitude
    t3darrbfamul[1,17,140]
    #
    
    setwd('C:/Users/m_lu0002/Desktop/Climate/minnesota')
    xy2d<-as.data.frame(xy2)
    load('xy2')
    xy2
    
    ####################### run bfast on deter #######################
    xcm<-c()
    ycm<-c()
    z<-array(,c(150,150))
    i2=1
    for (i1 in 1:791)
    {  
      
      i<-xy2d$x[i1]-58929
      j<-xy2d$y[i1]-48209
      #
      
       
      
      if( i<150 &&  j<150 && i>0 && j>0)
      {
  
        fets<-ts(fevi8[i,j,],start=c(2000,1),frequency=46) 
        b1<- bfast(fets,h=0.15,type= "OLS-MOSUM", max.iter=1)
        if (b1$nobp$Vt == FALSE ) {      
        xcm[i2]<-i
        ycm[i2]<-j
        i2=i2+1  
        print(i)
      } else {
        z[i,j]<-0
  print (0)
      }}}
    
  save(ycm,file='ycm.Rdata')
    length(which(z==0)) # CUSUM 58 #MOSUM 37 
  length(xcm)
   
    i2=3
    fets<-ts(fevi8[xc[i2],yc[i2],],start=c(2000,1),frequency=46)
    b1<- bfast(fets,h=0.15,type= "OLS-CUSUM", max.iter=1)
    plot(b1)
    
    #28 cusum -- mosum 49 x y total: 86
  result3<-data.frame(0,0,0,0,0)
  names(result3)<-c('intercept',"trend","lambda",'LL','LL0')
    p11<-c()
    p21<-c()
  p1<-c()
  p2<-c()
  p3<-c()
  
    #for (i1 in 2:49 )
    #{
      xcms<-c()
      ycms<-c()
      z<-array(,c(150,150))
      i2=1
      for (i1 in 1:791)
      {  
        print(i1)
        i<-xy2d$x[i1]-58929
        j<-xy2d$y[i1]-48209
        #             
        if( i<150 &&  j<150 && i>1 && j>1)
        #i=xcm [i1];j=ycm[i1]
        #i=i1;j=i1
       {
        fevi3b3<-fevi8[(i-1):(i+1),(j-1):(j+1),]
        fevi3b312t<-apply(fevi3b3,c(1,2),function(x) (ts(x,start=c(2000,1),frequency=46)-stl(ts(x,start=c(2000,1),frequency=46),'per')$time.series[,"seasonal"]))
        #dim(fevi3b312t)
        f2<-aperm(fevi3b312t,c(2,3,1))
        aa2<-as.data.frame.table(f2)$Freq
        aa3<-as.data.frame(as.data.frame.table(f2)$Freq)
        eday <- as.Date("2000-01-30")           # date 
        e8day <- seq(eday, length.out=636, by="8 days")
        
        
        try2<-spautolm(aa2~c(rep(c(1:636),each=9)),family="SAR",listw=listcn636)
        
 
      #  result3<-rbind(result3,c(try2$fit$coefficients[1],try2$fit$coefficients[2],try2$lambda,try2$LL,try2$LL0)) 
     
  #write.csv(result2,file='resultfromsar2.csv')
        #try1<-spautolm(aa2~c(rep(c(1:636),each=12)),family="CAR",listw=listcn636)
        rn<-lapply(1:9,function(i) {residuals(try2)[seq(i,636*9-(9-i),9)]})
        #get residuals for each time series  
        #p3<-array(,c(3,3))
        #dim(f2) #deseasonalized 
        ii<-5    
        fevi3b312t1<-ts(f2[2,2,],start=c(2000,1),frequency=46)
        p.Vt2 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM" )) 
        p.Vt3 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-MOSUM" )) 
        p.Vt <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]]))  )
      # jpeg(paste(i1,'randomSAR.jpg'), height=4, width=7, res=400,unit="in")
       # plot(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]])),main='')
      #title(paste('SAR:','(',xcm[i1],',', ycm[i1],')'),sub=paste('i=',i1,'','coordinates:','(',xcm[i1]+58929,'',ycm[i1]+48209,')'),cex=0.05 )
 #title(paste('SAR:','(',i,',', j,')'),sub=paste('i=',i1,'','coordinates:','(',i+58929,'',j+48209,')'),cex=0.05 )
  #     dev.off()
   #    jpeg(paste(i1,'randomoriginal.jpg'), height=4, width=7, res=400,unit="in")
    #   plot(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM"),main='')
 #title(paste('original:','(',i,',', j,')'),sub=paste('i=',i1,'','coordinates:','(',i+58929,'',j+48209,')'),cex=0.05 )     
 #title(paste('original:','(',xcm[i1],',', ycm[i1],')'),sub=paste('i=',i1,'','coordinates:','(',xcm[i1]+58929,'',ycm[i1]+48209,')'),cex=0.05 )
  #     dev.off()  
     
        p1[i1]<- p.Vt$p.value # spautolm residuals
        p2[i1]<-p.Vt2$p.value # linear regression residuals
        p3[i1]<-p.Vt3$p.value # MOSUM
       
 #p3[i,j]<-bfast(fevi3b312t,max.iter=1,level=0.05,h=0.15, type = "OLS-CUSUM")$output[[1]]$Vt.bp[1]
      }}
 save(p3,file='p3.R')
 length(which(!is.na(p3)))   
 p3[which(!is.na(p3))] 
 p2[which(!is.na(p2))] 
 length(which(p1[which(!is.na(p1))]<0.05)) 
 length(which(p2[which(!is.na(p2))]<0.05)) 
 length(which(p3[which(!is.na(p3))]<0.05)) 
 ycm[48]
  ip1<- which(p1<0.05)
  ip2<-which(p2<0.05)
length(match(ip1,ip2)) #28 vs.44
#there are 49 points detected by OLS MOSUM, 28 points detected by OLS CUSUM, 44points detected by OLS CUSUM SAR  
i=2; j=2
ii<-i+(j-1)*2
fevi3b312<-ts(f2[i,j,],start=c(2000,1),frequency=46)
 p1 - p2
plot(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]])))
plot(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM"))
    save(p1,file='100b100sar.Rdata')
    save(p2,file='100b100p.Rdata')
    summary(p1)
fets<-ts(fevi8[i,j,],start=c(2000,1),frequency=46)
fets<-ts(result3$lambda[4:639],start=c(2000,1),frequency=46)
fets<-ts(result3[,3][638:1273],start=c(2000,1),frequency=46)
plot(fets)
b1<- bfast(fets,h=0.15,type= "OLS-MOSUM", max.iter=1)
plot(b1)
str(b1)
b1
######## sar on each time step to see if the lambda changed###
result3<-data.frame(0,0,0,0,0)
names(result3)<-c('intercept',"trend","lambda",'LL','LL0')
  
i2=41
     i=xcm [i2];j=ycm[i2]
    #i=i1;j=i1
    fevi3b3<-fevi8[(i-1):(i+1),(j-1):(j+1),]
  fevi3b312t<-apply(fevi3b3,c(1,2),function(x) (ts(x,start=c(2000,1),frequency=46)-stl(ts(x,start=c(2000,1),frequency=46),'per')$time.series[,"seasonal"]))
  #dim(fevi3b312t)
  f2<-aperm(fevi3b312t,c(2,3,1))
  aa2<-as.data.frame.table(f2)$Freq
  aa3<-as.data.frame(as.data.frame.table(f2)$Freq)
  eday <- as.Date("2000-01-30")           # date 
  e8day <- seq(eday, length.out=636, by="8 days")
  
for (i1 in 0:635)
{  
  try2<-spautolm(aa2[(1+i1*9):(9+i1*9)]~1,family="SAR",listw=cn1)  
  result3<-rbind(result3,c(try2$fit$coefficients[1],try2$fit$coefficients[2],try2$lambda,try2$LL,try2$LL0)) 
print(i1)
}
which.max( result3[,3][2:637] )
  #write.csv(result2,file='resultfromsar2.csv')