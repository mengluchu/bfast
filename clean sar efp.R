
load("fevi20b20.Rdata") #fevi20by20: 20 by 20 subarray of evi values

library(devtools)
library(strucchange)
library(spacetime) 
library(spdep)
install_github("strucchange","mengluchu",build_vignettes = FALSE)
# get neighbor
eday <- as.Date("2000-01-30")           # date 
e8day <- seq(eday, length.out=636, by="8 days")
xyd<-expand.grid(x1=1:3,y1=1:3)
coordinates(xyd)<-~x1+y1
lecube<-3*3*636
aa3<-as.data.frame(c(1:lecube))
stfdf3b3<-STFDF(xyd,e8day,aa3) ## for creating neighbors only, aa3 could be any data?
cn<-cell2nb(3,3, type ="queen",torus =FALSE)
neigh1<-nbMult(cn, stfdf3b3, addT = FALSE, addST = FALSE) # only spatial neighbours are added for each time step
listcn636<-nb2listw(neigh1)

#p1<-array(,c(18,18))
#p2<-array(,c(18,18))
#p3<-array(,c(18,18))
#for(i in 1:18)
# {
#   for (j in 1:18)
#   {
i=1;j=18  #i, j can be from 1 to 18    
fevi3b3<-fevi20b20[i:(i+2),j:(j+2),]
fevi3b312t<-apply(fevi3b3,c(1,2),function(x) (ts(x,start=c(2000,1),frequency=46)-stl(ts(x,start=c(2000,1),frequency=46),'per')$time.series[,"seasonal"]))
f2<-aperm(fevi3b312t,c(2,3,1))
aa2<-as.vector(f2)

try2<-spautolm(aa2~c(rep(c(1:636),each=9)),family="SAR",method= "Matrix", listw=listcn636)
#try3<-spautolm(aa2~c(rep(c(1:636),9)),family="SAR",method= "Matrix", listw=listcn636)  not exactly the same but less different than i thought..
#try1<-spautolm(aa2~c(rep(c(1:636),each=9)),family="CAR",method= "Matrix",listw=listcn636) a bit different than SAR
#so much faster than dense array!

rn<-lapply(1:9,function(i) {residuals(try2)[seq(i,636*9-(9-i),9)]})
#get residuals for each time series

ii<-5   # get the middle pixel (5 for 3*3 matrix)
ti<-1:636

fevi3b312t1<-ts(f2[2,2,],start=c(2000,1),frequency=46) # reconstruct the time series
p.Vt2 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM" )) 
p.Vt3 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-MOSUM" ))
p.Vt  <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = "OLS-CUSUM", spatial1=as.numeric(rn[[ii]]))  )

#p1[i,j]<-p.Vt$p.value # spautolm residuals
#p2 [i,j]<-p.Vt2$p.value # linear regression residuals
#p3 [i,j]<-p.Vt3$p.value # MOSUM
#p-values of efp    
#bfast(fevi3b312t,max.iter=3,level=0.05,h=0.15, type = "OLS-CUSUM")$output[[1]]$Vt.bp[1] #bfast result
#}
#}
#test correlation of two residuals
resorig<-residuals(lm(f2[2,2,]~ti))
plot(resorig,typ='l') #residuals of time series
lines(rn[[5]],col='red') # residuals of SAR
cor(resorig,rn[[5]])


