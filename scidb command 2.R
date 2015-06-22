library(scidb)
install.packages("bfast", repos="http://R-Forge.R-project.org")

iquery("save(MOD09Q1_JUARA, '/home/menglu/shared/JUARA1.bin', -2, '(int16, int16, uint16)')")


iquery("
store(
  project(
      apply(      
        repart(
            apply(
              subarray(
              
              MOD09Q1_JUARA,58828,48103,6,59679,49050,643),
              evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),
 <red:int16,nir:int16,quality:uint16,evi2:double>[col_id=0:851,3,2,row_id=0:947,3,2,time_id=0:637,638,0]),
              the_i,double(col_id),the_j,double(row_id),the_t,double(time_id)),
              the_i, the_j,the_t, evi2),
              allarray1)")
       
 
iquery("store(apply(try221b221,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),evi2t221221) ")

iquery("store(repart(evi2t2121, <red:int16,nir:int16,quality:uint16,evi2:double>
       [col_id=0:20,3,0,row_id=0:20,3,0,time_id=0:637,638,0]),evi2t2121nooverlap3)");

iquery("store(apply(evi2t2121nooverlap3,the_i,double(col_id),the_j,double(row_id),the_t,double(time_id)),tryijt2) ")


iquery("store(project(tryijt2, the_i, the_j,the_t, evi2), evi210b210ijt3)")

iquery("show(allarray1)",return=TRUE)

library("bfast", lib.loc="C:/Users/m_lu0002/Documents/R/win-library/3.1")
scidbconnect("gis-bigdata.uni-muenster.de", 49974, "meng", "bigdata2014")
library(scidb)
scidbconnect("gis-bigdata.uni-muenster.de", 49972, "menglu", "ms3A5Bfn9PCycMr")
scidblist()
iquery('show(all150150o)',return=TRUE)
"try221b221"
iquery("store(repart(evi21b21ijt3, <red:int16,nir:int16,quality:uint16,evi2:double>
       [col_id=0:20,3,2,row_id=0:20,3,2,time_id=0:637,638,0]),evi21b21o2)");
iquery("store(subarray(allarray1,0,0,0,848,944,635),all848944o)")
a<-scidb('all150150o')
which(is.na(a[,,][]))
length(
  a[1:3,1:3,633][]
 
a<-scidb("MOD09Q1_JUARA")
which(is.na(
  dim(
    length(a[58930,48210,6:643][]$red)
  )
  )
  )

#all848944o
#evi21b21o
#all150150
build(<v:int64>[k=0:9999,1000,0],random())
iquery("store(build(<val:double> [time=0:300,100,20],time),ttrc)")
iquery("show(ewmann)",return=TRUE)
a<-scidb("ewmann")
a[1,0][]

iquery("store(r_exec(apply(ttrc,t,double(time)),
             'output_attrs=1',
             'expr=
    
    chunk<-as.double(length(t))
    
  
    list(chunk)'),
      ewmann)
")


iquery("subarray(try221b221,58828,48103,46,58835,48110,56)", `return` = TRUE, afl = TRUE, iterative = FALSE, n = 10000)

time1<-proc.time()
iquery("store(r_exec(all150150o1,         
       'output_attrs=2',
       'expr= 
       require(bfast)
       require(strucchange)
       require(spdep)
       load(\"/home/menglu/listcn636.Rdata\")
       sarsarbfastt1<-array(, 636)
       sarsarbfasts1<-array(, 636)
       
       eviarray<-array(evi2,c(636,3,3))          # scidb array to r array                
       eviarray2<-aperm(eviarray,c(2,3,1))
       bfa <-try(bfastsar(eviarray2, season= \"harmonic\",max.iter=1 ))       
       
    if(class(bfa) == \"try-error\") 
    {
     sarsarbfastt1[1:636]<- 1
      sarsarbfasts1[1:636]<- 1
      
    }   else if (bfa$nobp$Vt == FALSE) 
    {
       bkpttr <- as.integer(bfa$output[[1]]$Vt.bp) #   of 16 days e.g. year=16*50/365
       sarsarbfastt1[ bkpttr]<-100
    }  else if (bfa$nobp$Wt == FALSE) 
    {
       bkptse <- as.integer(bfa$output[[1]]$Wt.bp)
       sarsarbfasts1[ bkptse]<-200
    } else  {
       sarsarbfastt1[1:636]<-2
       sarsarbfastt1[1:636]<-2
    }     
       sarsarbfastt1[which(is.na(sarsarbfastt1))]<-0 
       sarsarbfasts1[which(is.na(sarsarbfasts1))]<-0 
       sarsarbfastt1<-as.data.frame.table (sarsarbfastt1)$Freq  
       sarsarbfasts1<-as.data.frame.table (sarsarbfasts1)$Freq  
      list( sarsarbfastt1,sarsarbfasts1)'),rexe4)",return=TRUE) 

time3<-proc.time()-time1
time3
user    system   elapsed 
4.93      3.72 164359.90 
a<-scidb("rexe3")
a[1,1][]

iquery("store(r_exec(all2121o,         
       'output_attrs=2',
       'expr= 
        require(bfast)
        require(strucchange)
        require(spdep)
        load(\"/home/menglu/listcn636.Rdata\")
        sarsarbfastt1<-array(,c(3,3,636))
        sarsarbfasts1<-array(,c(3,3,636))
        
        eviarray<-array(evi2,c(636,3,3))          # scidb array to r array                
        eviarray2<-aperm(eviarray,c(2,3,1))
        bfa <-bfastsar(eviarray2, season= \"harmonic\",max.iter=1 )       

        if (bfa$nobp$Vt == FALSE) 
        {
        bkpttr <- as.integer(bfa$output[[1]]$Vt.bp) #   of 16 days e.g. year=16*50/365
        sarsarbfastt1[2,2,bkpttr]<-100
        } 
        
        if (bfa$nobp$Wt == FALSE) 
        {
         bkptse <- as.integer(bfa$output[[1]]$Wt.bp)
         sarsarbfasts1[2,2,bkptse]<-200
        } else  {
         bkpttr<-0
         bkptse<-0
                }     
 sarsarbfastt1[which(is.na(sarsarbfastt1))]<-0 
 sarsarbfasts1[which(is.na(sarsarbfasts1))]<-0 
 sarsarbfastt1<-as.data.frame.table (sarsarbfastt1)$Freq  
 sarsarbfasts1<-as.data.frame.table (sarsarbfasts1)$Freq  
      sarsarbfastt1,sarsarbfasts1)'),rexe3)",return=TRUE) 






load("/home/menglu/listcn636.Rdata")
        sarsarbfastt1<-array(,c(148,148,636))
        sarsarbfasts1<-array(,c(148,148,636))
        636*9
        eviarray<-array(rnorm(5724),c(636,3,3))          # scidb array to r array                
        eviarray2<-aperm(eviarray,c(2,3,1))
        bfas <-bfastsar(eviarray2, season= "harmonic",max.iter=1 )       

        if (bfa$nobp$Vt == FALSE) 
        {
        bkpttr <- as.integer(bfa$output[[1]]$Vt.bp) #   of 16 days e.g. year=16*50/365
        sarsarbfastt1[the_i,the_j,bkpttr]<-100
        } 
        
        if (bfa$nobp$Wt == FALSE) 
        {
         bkptse <- as.integer(bfa$output[[1]]$Wt.bp)
         sarsarbfasts1[the_i,the_j,bkptse]<-200
        } else  {
         bkpttr<-0
         bkptse<-0
                }     
 sarsarbfastt1[which(is.na(sarsarbfastt1))]<-0 
 sarsarbfasts1[which(is.na(sarsarbfasts1))]<-0 




iquery("store(r_exec(t7,         
       'output_attrs=1',
       'expr= 
       require(bfast)
       require(strucchange)
       require(spdep)
       load(\"/home/menglu/listcn636.Rdata\")
 
        output2 <-bfastsar(array(rnorm(5764),c(3,3,636)), season= \"harmonic\",max.iter=1 )       
      # output2<-bfast(ts(1:636,start=2001,frequency=46)  , season= \"harmonic\",max.iter=1  )
        a<-as.double(output2[1]$Yt)
       list(a  )'),rexe2)",return=TRUE)       

scidbconnect("gis-bigdata.uni-muenster.de", 49962, "menglu", "m7KLKHYKdrQ8fQM") #enterprise

x <- as.scidb(matrix(rnorm(5000*20),nrow=5000))
y <- as.scidb(rnorm(5000))
M <- glm.fit(x, y)
iquery("list('aggregates')",return=TRUE)

coef(M)[]
scidblist()
iquery("load('t1.txt')")

iquery("create array JUARA <red:int16,nir:int16,quality:uint16> [col_id=58828:59679,502,5,row_id=48103:49050,502,5,time_id=0:9200,1,0]")

iquery("store(build <exposure:double> [i=1:8,1,0,j=1:8,1,0]),trya)")
 
iquery("show(MOD09Q1_JUARA )",return=TRUE)

iquery("load(JUARA,'/home/scidb/shared/JUARA.bin','(int16, int16, uint16)')") 

wiquery("show(t1)",return=TRUE)

iquery("scan(t1)",return=TRUE)
iquery("window(TEST_ARRAY,1,1,1,1,median(num))",return=TRUE)
iquery("list('aggregates')",return=TRUE)
# Using glm (similar to glm)
# First, let's run a standard glm example from Dobson (1990) Page 93:
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
d.AD <- data.frame(treatment, outcome, counts)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),data=d.AD)
summary(glm.D93)


sl<-
  scidblist()
sl[-2]
 
x <- as.scidb(d.AD) # couldnt omit name 
scidbremove( sl[-2],force=TRUE)
?rm
iquery("scan(x)",return=TRUE)
glm.D93_sci = glm(counts ~ outcome + treatment, family = poisson(), data=a2)
summary(glm.D93_sci)

persist(glm.D93_sci)
save(glm.D93_sci, file="my_model.rdata")
#ms3A5Bfn9PCycMr
#VxSxrk5cnXed74N
#MOD09Q1_JUARA
scidblist()
scidbversion()


iquery("load_library('dense_linear_algebra')", release=1,resp=FALSE)
iquery("load_library('linear_algebra')", release=1,resp=FALSE)
#iquery("load_library('collate')", release=1,resp=FALSE)
iquery("build(<z:double>[i=1:10 ,1,0],i*0.1+0.9)",return=TRUE)


x = iquery(" r_exec(build(<z:double>[i=1:100,100,0],0.1*i+0.9),
'expr=x<-c(0.01:0.08)
require(raster)
r1 <- raster(nrows=108, ncols=21, xmn=0, xmx=10)
r1<-setValues(r1,c(1:ncell(r1)))
 
xxx<-\\\'raster1\\\'
 
 try<-values(r1)*0.1
list(try)')",return=TRUE);
##
 

#######
iquery(" r_exec(build(<z:double>[i=1:100,100 ,0],0.1*i+0.9),
'expr=
require(strucchange)
tsz<-ts(z,frequency=12)
 
btsz<-efp(tsz~1 )
wt<-coef(btsz)[1]
 
list( wt)')",return=TRUE)
##
 
iquery(" r_exec(build(<z:double>[i=1:100,100 ,0],0.1*i+0.9),
'expr= ts(seq())
require(zoo)
library(bfast)
NDVIb <- as.ts(zoo(som$NDVI.b, som$Time))
monb <- bfastmonitor(NDVIb, history=\\\'ROC\\\',start = c(2010, 13))
bre<-monb$model$residuals

list(as.double(bre))')
       ",return=TRUE)

plot(aa$expr_value_0,typ='l',main='bfastmonitor residuals')
############################################################
iquery('dimensions(MOD09Q1_MENG_20140416)',return=TRUE)
iquery('store(repart(smalltest1, <red:int16,nir:int16,quality:uint16> [col_id=0:3,1,0,row_id=0:3,1,0,time_id=0:636,637,0]),target2)')
iquery("store(apply(target2,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),smallevi2) ")
iquery("store(build(< z:int64> [i=1:8,2,0,j=1:8,2,0], i ),t4)")
iquery("store(apply(t4,z1,z+0.5),t5)")


iquery(" store(r_exec(t6,
'output_attrs=2',
'expr=
output1=z1-0.001;
location=c(rep(2.2',64));
list(output1,location)'),t8)
       ")

iquery(" store(r_exec(t6,
'output_attrs=2',
'expr=
output1=z1-0.001; #length(z1)=4
output2=c(1:4) ;
list(output1,output2*1.0)'),t9)
       ")
iquery("scan(t9)",return=TRUE)

iquery("store(project(t5,z1),t6)",return=TRUE)
iquery("store(build(< z:double> [i=1:800,2,0], i ),t7)")
#aggregate time series
#iquery("regrid(t7,100,avg(z)) ",return=TRUE)
iquery("dimensions(MOD09Q1_MENG_20140416)",return=TRUE)

require(bfast)
ts1= ts( evi2  ,frequency=46)


monb <- bfastmonitor(ts1, start = 3)
bre<-monb$model$residuals
as.double(bre
s1<-scidb("smallevi1")
ts1<-as.double(s1[][]$evi2[1:493])
ts2<-ts(ts1,frequency=46)
monb <- bfastmonitor(ts2, start = 8)
bre<-monb$model$residuals


 

# if you dont want to deal with chunk, set chunck size same as dimension
# test the instance and chunks
# it doesnt process sequentially, it sometimes put several chunks in one instance. 
# In the output you will have no idea which result conrrespond to which chunk.
# if theres overlap, then looks more complex

iquery(" r_exec(build(<z:double>[i=1:8,2,0,j=1:8,2,0],  i+j*0.1 ),
'expr=  
 x1<-z+0.0001
 list(as.double(x1))')
       ",return=TRUE)

iquery(" r_exec(build(<z:double>[i=1:8,2,0,j=1:8,2,0],  i+j*0.1 ),
'expr=  
 x1<-z+0.0001
 list(as.double(x1))')
       ")
iquery("show(z)",return=TRUE)
####################################################################################
iquery("store( redimension(apply(project(sort(apply(build(<v:int64>[k=0:9999,1000,0],random()),p,k)),p),m,n),
  <p:int64> [m=0:*,1000,0]),B)",return=TRUE)
       iquery("show(B)",return=TRUE)
btsz<-bfast(tsz,h=0.15, max.iter=1)
plot(btsz)
z <- zoo(cbind(foo = rnorm(5), bar = rnorm(5)))
z$foo
z$xyz <- zoo(rnorm(3), 2:4)
 
z4 <- zoo(11:15, complex(real = c(1, 3, 4, 5, 6), imag = c(0, 1, 0, 0, 1)))
merge(z4, lag(z4))
merge(z5, lag(z5))

# index values relative to 2001Q1

z5 <- zoo(11:15, letters[1:5])
merge(z5, lag(z5))
str(btsz)
plot(simts)
as.double(simts$time.series)
library(zoo)
lm(z~1)$residuals)
r1 <- raster(nrows=108, ncols=21, xmn=0, xmx=10)
r1<-setValues(r1,c(1:ncell(r1)))
ncell(r1)
plot(r1)
values(r1)
x<-c(0.01:0.08)
y<-runif(10)
z<-seq(1.0,19,0.1)
z
ts(z,frequency=3)
sum(x^2+y^2)/250
lm(z~1)$coefficient
x <- as.scidb(matrix(rnorm(5000*20),nrow=5000))
y <- as.scidb(rnorm(5000))
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
d.AD <- data.frame(treatment, outcome, counts)
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson(),data=d.AD)
summary(glm.D93)
ff <- log(Volume) ~ log(Height) + log(Girth)
utils::str(m <- model.frame(ff, trees))
mat <- model.matrix(ff, m)

dd <- data.frame(a = gl(3,4), b = gl(4,1,12)) # balanced 2-way
options("contrasts")
model.matrix(~ a+b , dd)
model.matrix(~ a + b, dd, contrasts = list(a = "contr.sum"))
model.matrix(~ a + b, dd, contrasts = list(a = "contr.sum", b = "contr.poly"))
m.orth <- model.matrix(~a+b, dd, contrasts = list(a = "contr.helmert"))
crossprod(m.orth) # m.orth is  ALMOST  orthogonal
# Compare with:
d.AD_sci = as.scidb(d.AD)
glm.D93_sci = glm_scidb(counts ~ outcome + treatment, family = poisson(), data=d.AD_sci)
summary(glm.D93_sci)
M <- glm.fit(x, y)
coef(M)[]
M<-glm_scidb(x,y)
M <- glm.fit(x, y,weights=NULL,family=gaussian())
library('devtools')
install_github("SciDBR","paradigm4")
iquery('dimensions(MOD09Q1_MENG_20140416)',return=TRUE)
scidblist()
#ndviallr<-scidb('ndviinuseall')
#eviallr<-scidb('evi2all1')
iquery('store(subarray(MOD09Q1_MENG_20140416,58930,48210,6,59129,48409,643),inuse6)')
iquery("store(apply(inuse6,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),evi2inuse6) ")
iquery("store(apply(inuse6,ndvi,1.0*(nir-red)/(nir+red)),ndviinuse6) ")
iquery('load_library("linear_algebra")')
iquery("store(apply(inuse6,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),evi2inuse6) ")
iquery('dimensions(MOD09Q1_JUARA)',return=TRUE)
#150by150
quality<-scidb('inuse6')
length(quality[,,][]$quality==4096)
quality150<-quality[1:150,1:150,][]$quality
save(quality150, file='quality150')

load('quality150')
length(which(quality150==4096))/length(quality150)
malquality<-quality150[quality150!=4096]
length(malquality)
table(malquality)
length(quality150)/22500
arrquality150<-array(quality150,c(150,150,637))
qin1<-which(arrquality150!=4096,arr.ind=TRUE)
x<-qin1[,1]+58930-1
y<-qin1[,2]+48210-1
z<-qin1[,3]
coor.quality<-coordinates(getxyMatrix(cbind(x,y),231.6564))
length(coordinates(cbind(x,y)))
plot(coor.quality)
plot(x,y)
length(y)
22772/2
36964/2
plot(modis.mt52,col='purple',add=TRUE)
spplot(coor.quality)
 
 
plot(1)
dev.new()
hexbinplot(y~x)
qualitysp<-coordinates(cbind(x,y))
##3rd array: 41 by 43 around 520 for validation 
iquery('store(subarray(MOD09Q1_MENG_20140416,59138,48712,6,59180,48752,643),around6)')
iquery("store(apply(around6,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),aevi5) ")
iquery('dimensions(MOD09Q1_MENG_20140416)',return=TRUE)
iquery('store(subarray(around6,1,1,1,4,4,638),smalltest1)')

# create evi around 52
evi52re<-scidb('aevi5')

 
id<-c()
id<-evi52re[1,1,][]$evi2
 
itry<-array(,c(43,41,637))
#i0<-1:41;i01<-rep(i0,43);j0<-1:43;j01<-rep(j0,each=41)
  #tr1<-lapply(c(i01,j01), function(j01,i01) evi52re[j01,i01,][]$evi2)

for(i in 0:42)
{
  for (j in 0:40)
  {
    itry[i,j,]<-evi52re[i,j,][]$evi2
  }
}
itry[9,7,6]-evi52re[9,7,5][]$evi2
 
 
#str(evi_a52)
length(ic)/43/41
#evi52re[1,1,2][]$evi2
eviar510<-array(ic,c(41,43,638))
ib[1:10]
head(eviar510)
evi52re[1:10,1,1][]$evi2
head(ic)
evi_a520<-aperm (evi510,c(2,1,3)) #rotate # new second array
save(itry,file='evi510.Rdata')
save(itry,file='itry')
str(itry2)
setwd("C:/Users/m_lu0002/Desktop/Climate/minnesota")
load('evi510.Rdata')
evi52re[,1,1][]$evi2-as.numeric(eviar510[,1,2])
#############evi all
iquery("store(apply(MOD09Q1_MENG_20140416,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),evi2all1) ") #around 25 minutes
#2.5*((NIR-Red)/(NIR+2.4*Red+1)
###
##the 1st 100 by 100 array
# iquery('dimensions(inuse)',return=TRUE) 100*100*638
#No    name start length chunk_interval chunk_overlap   low  high  type
#1  0  col_id 58828    852            502             5 58832 58929 int64
#2  1  row_id 48103    948            502             5 48110 48209 int64
#3  2 time_id     0   9201              1             0     6   643 int64
###41 43
c(59138:59180)
c(48712:48753)


if(!file.exists("EQY_US_NYSE_TRADE_20130404.zip"))
{
  download.file("ftp://ftp.nyxdata.com/Historical%20Data%20Samples/TAQ%20NYSE%20Trades/EQY_US_NYSE_TRADE_20130404.zip",
                "EQY_US_NYSE_TRADE_20130404.zip")
}

# Load just the variable names first
u  = unz("EQY_US_NYSE_TRADE_20130404.zip",filename="nysetrades20130404")
N = names(read.table(u,sep="|",stringsAsFactors=FALSE,header=TRUE,nrows=1))
# Now load a small part of the data (later in the day)...
u = unz("EQY_US_NYSE_TRADE_20130404.zip",filename="nysetrades20130404")
x = read.table(u,sep="|",stringsAsFactors=FALSE,header=FALSE,nrows=1e4,skip=1e6)
names(x) = N

# Let's limit our tiny example to three bank stocks
x = x[x$Symbol %in% c("JPM","MS","C"),]

# Assemble some of this into an xts object (requires the xts package)
library("xts")
options(digits.secs=5)
x$t = as.POSIXct(x$SourceTime/1000,origin="2013-4-4",tz="EDT")
x$price = x$PriceNumerator/10^x$PriceScaleCode
idx = x$Symbol=="C"
p = cbind(C=as.xts(x[idx,"price"],order.by=x[idx,"t"]))
idx = x$Symbol=="MS"
p = cbind(p,MS=as.xts(x[idx,"price"],order.by=x[idx,"t"]))
idx = x$Symbol=="JPM"
p = cbind(p,JPM=as.xts(x[idx,"price"],order.by=x[idx,"t"]))

# Let's cut this down to a really small example:
p = p[!duplicated(index(p)),]
p = p[1:20,]

# These data are both sparse (in time) and contain missing values.
print(p)

x$time = as.numeric(x$t)*1000
x$symbol = as.integer(as.integer(factor(x$Symbol,levels=c("C","MS","JPM"))) - 1)
# R doesn't have 64-bit integers, so we need to cast the big time values from
# doubles...
raw = project(bind(as.scidb(x[1:20,c("time","price","symbol")],types=c("double","double","int64")),"t","int64(time)"), c("t","symbol","price"))
# Redimension the raw data into a (sparse) time x symbol array:
prices = redimension(raw, dim=c("t","symbol"))

# Although the prices array is very sparse because its time coordinate axis is
# in milliseconds, it contains the same data as our R xts array p.
# Let's view the schema of that:
print(schema(prices))

 
load_library('r_exec');


 

 
############################ r exec############################
scidbremove("evi2t2121",force=TRUE)
iquery('store(subarray(MOD09Q1_MENG_20140416,58930,48210,6,58950,48230,643),try21b21)')
iquery('store(subarray(MOD09Q1_MENG_20140416,58930,48210,6,58950,48230,643),try21b21)')


iquery("store(apply(try21b21,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),evi2t2121) ")

iquery("store(repart(evi2t2121, <red:int16,nir:int16,quality:uint16,evi2:double>
       [col_id=0:20,3,0,row_id=0:20,3,0,time_id=0:637,638,0]),evi2t2121nooverlap3)");

iquery("store(apply(evi2t2121nooverlap3,the_i,double(col_id),the_j,double(row_id),the_t,double(time_id)),tryijt2) ")


iquery("store(project(tryijt2, the_i, the_j,the_t, evi2), evi21b21ijt3)")
iquery("store(subarray(evi21b21ijt,0,0,0,2,2,635),test3by3)" )
iquery("store(subarray(evi21b21ijt,0,0,0,1,1,635),test2by2)" )
###########################################################################
#evi21b21ijt  :  array with attribute i j t for further retrieving array###
###########################################################################

iquery("scan(test2by2)",return=TRUE)
iquery("scan(evi21b21ijt)",n= 800, return=TRUE)
iquery("r_exec(MENG_test, 
'output_attrs=3', 
'expr=list(the_i, the_j, value + 100)')");
 
iquery("show(evi21b21ijt)",return=TRUE)
#evi21b21ijt 21 by 21
#test3by3 3 by 3


 

ptm <- proc.time()
 

iquery("store(r_exec(evi21b21ijt3,         
'output_attrs=4',
'expr= eviarray<-array(evi2,c(638,3,3))
eviarray2<-aperm(eviarray,c(2,3,1))
list()
require(bfast)

eviarray<-array(evi2,c(638,3,3))
eviarray2<-aperm(eviarray,c(2,3,1))

eday <- as.Date(\"2000-01-30\")
e8day <- seq(eday, length.out=635, by=\"8 days\")
output2 <-bfaarray2(eviarray2,dates=e8day,aggre= \"month\",season= \"harmonic\",max.iter=1,level=0.05)

    
allbreakpoints<-output2[1:6,,]  
allmagnitudes<-output2[7:12,,]  

mtime<-allbreakpoints  

dimx<-dim(output2)[2]
dimy<-dim(output2)[3]
dimt<-591

newarray<-array(,c(dimx,dimy,dimt))  

mag<- allmagnitudes   

t3darrbfamul<-tasd.bfa.mul(newarray,mtime,mag)
t3darrbfamul[which(is.na(t3darrbfamul))]<-0

t3df<-as.data.frame.table (t3darrbfamul)
 
t3dff<-as.data.frame.table (t3darrbfamul)$Freq

list( the_j,the_i,the_t,t3dff )'),rexec21by21)",return=TRUE)       

proc.time() - ptm
#list( the_j,the_i,the_t)')",return=TRUE)    
     
#start time 11:39    
  # iquery("show(rexec21by21)",return=TRUE)  
#levels(t3df$Var1)<-c(1:length(levels(t3df$Var1)))
#levels(t3df$Var2) <-c(1:length(levels(t3df$Var2)))
#levels(t3df$Var3)<-c(1:length(levels(t3df$Var3)))
#as.double(t3df$Var1),as.double(t3df$Var2),as.double(t3df$Var3)
iquery("store(r_exec(MENGtest2, 'output_attrs=3', 'expr=list(thei, thej, value + 100)'), MENGtest3)")

iquery("store(attribute_rename(project(unpack(rexec21by21, tmpDim), 
        expr_value_0,expr_value_1,expr_value_2,expr_value_3),
        expr_value_0, thei, expr_value_1, thej, expr_value_2, thet, expr_value_3,value), re2121)");


#Restore the dimensions to int64
iquery("store(project(apply(re2121, i, int64(thei),
       j, int64(thej), t, int64(thet)), i, j, t, value), rexec2121)");

#Restore to 2D array
iquery("store(redimension(rexec2121, <value:double> [i=0:20,3,0,j=0:20,3,0,t=0:637,638,0]), rexec2121re)");

iquery("scan(rexec2121)",n=700,return=TRUE)
iquery("scan(try21b21)",n=700,return=TRUE)

test<-scidb("MOD09Q1_MENG_20140416")
test[58930,48210,570:589][]
########################################################################################
NDVIb <- as.ts(zoo(som$NDVI.b, som$Time))
monb <- bfastmonitor(NDVIb, history=\\\'ROC\\\',start = c(2010, 13))
bre<-monb$model$residuals

scidblist()
scidbremove("MENG_test_res",force=TRUE);

remove(MENG_test3);

remove(MENG_test4);

--Reduce to single dimension; 
get rid of inst and n, rename attributes
store(attribute_rename(project(unpack(MENG_test1, tmpDim), expr_value_0,expr_value_1,expr_value_2), expr_value_0, the_i, expr_value_1, the_j, expr_value_2, value), MENG_test2);



--Restore the dimensions to int64
store(project(apply(MENG_test2, i, int64(the_i), j, int64(the_j)), i, j, value), MENG_test3);


--Restore to 2D array
store(redimension(MENG_test3, <value:double> [i=0:3,4,0,j=0:3,4,0]), MENG_test4);


iquery(" r_exec(evi21b21ijt3, 'output_attrs=4', 
'expr=list(as.double(length(evi2)),as.double(length(the_i)),as.double(length(the_j)),as.double(length(the_t)))')", n=10, return=TRUE)
591*3*3

iquery(" r_exec(evi21b21ijt3, 'output_attrs=1', 
'expr=list(as.double(the_t))')", n=700, return=TRUE)
 
iquery("dimensions(MOD09Q1_MENG_20140416)",n=10,return=TRUE)
iquery("show(evi21b21ijt3)",n=10,return=TRUE)
length(c(1,NA,1))'

    iquery("store(build(<val:double NULL DEFAULT null>[i=0:9,9,0,j=0:9,9,0],random()%2), COMWAY)")

 
    iquery("insert(project(apply(apply(apply(join(COMWAY, window(COMWAY, 1, 1, 1, 1, sum(val))), sum, val_sum - val),factor,iif((val = 0 and sum != 3), 0,iif((val = 0 or sum = 3), 1,iif(sum < 2, -1,iif(sum > 3,-1,0))))    ),newval, double(val + factor)), newval), COMWAY)")

 
   
                     iquery("scan(COMWAY)",return=TRUE)

