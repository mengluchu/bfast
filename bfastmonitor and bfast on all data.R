#The changes can be detected by comparing the observations and the distribution of the Kriging estimations. install.packages("bfast", repos = "http://R-Forge.R-project.org")
load('fevi8.Rdata')
ora150p005t<-array(,c(150,150,637))
ora150p005s<-array(,c(150,150,637))
ora150p1t<-array(,c(150,150,637))
ora150p1s<-array(,c(150,150,637))
ora150p0025t<-array(,c(150,150,637))
ora150p0025s<-array(,c(150,150,637))
length(which(!is.na(ora150p005s)))

date.trend2<-data.frame()
date.seasonal2<-data.frame()
i4=1

or4240p005t<-array(,c(42,40,637))
or4240p005s<-array(,c(42,40,637))
or4240p0025t<-array(,c(42,40,637))
or4240p0025s<-array(,c(42,40,637))
or4240p1t<-array(,c(42,40,637))    
or4240p1s<-array(,c(42,40,637)) # takes far more time!
for (i in 1:42)
{  
  for (j in 1:40)
  { 
    
    #spt1<-fevi8[i,j,] 
    spt1<-itry2[i,j,] 
    
    spt<-ts( spt1,start=c(2000,7),end=c(2013,44),frequency=46)
    
    
    fitmon <- bfast(spt,h=0.15, season="harmonic", max.iter=3,level=1) 
    
    
    
    if(fitmon$nobp$Vt==FALSE)
    {
      
      date.trend2<-as.integer(fitmon$output[[1]]$Vt.bp)
      
      or4240p1t[i,j, date.trend2]<-100
    }
    if(fitmon$nobp$Wt==FALSE)
    {
      date.seasonal2<-as.integer(fitmon$output[[1]]$Wt.bp)
      or4240p1s[i,j, date.seasonal2]<-500
      
    }
    
    i4=i4+1
    print(i4)
  }
  
}

length(which(!is.na(ora150p005s)))
save(ora150p005s, file='ora150p005s.Rdata')
save(ora150p005t, file='ora150p005t.Rdata')

save(or4240p005t, file='or4240p005t.Rdata')
save(or4240p005s, file='or4240p005s.Rdata')

save(or4240p1t, file='or4240p1t.Rdata')
save(or4240p1s, file='or4240p1s.Rdata')
length(which(!is.na(or4240p005s)))
############################################# monitoring #########

monitor01t<-array(,c(150,150,637))
monitor01s<-array(,c(150,150,637))
monitor01t08<-array(,c(150,150,637))
i=1
j=1
for (i in 1:150)
{  
  for (j in 1:150)
  { 
    print(cbind(i,j))
    #134
    
    #spt1<-fevi8[i,j,] 
    spt1<-fevi8[i,j,] 
    
    # spt<-ts( spt1,start=c(2000,7),end=c(2013,44),frequency=46)
    spt<-ts( spt1, frequency=46)
    
    #bmonit<- bfastmonitor(spt,h=0.25,start=8,history='ROC') # start from 2007 index 322-368-414
    bmonit<- bfastmonitor(spt,h=0.25,start=9,history='ROC')    #  2008
    timeindex<-round((bmonit$breakpoint-1)*46)
    
    if(!is.na(timeindex))
      
    {
      
      moni<-window(spt, start=bmonit$history[1],end= bmonit$monitor[2]) 
      moni01<-bfast01(moni)
      cl<-bfast01classify(moni01, alpha = 0.05, pct_stable = 0.25)
      cl$flag_type
      monitor01t08[i,j, timeindex]<- cl$flag_type
       plot(bmonit)
    }
  }
}
save(monitor01t08,file='monitor1t08.Rdata')
summary(monitor01t08)
save(monitor1t,file='monitor1tt.Rdata')
save(monitor1t08,file='monitor1t08.Rdata')
length(which(!is.na(monitor1t08)))  #10396 #monitoring from 2007: 9898 # monitoring from 2008:  10219 
scatterplot3d(which(!is.na(monitor1t08),arr.ind=TRUE))
which(!is.na(monitor1t),arr.ind=TRUE)[,3]
str(monitor1t)
monitor1t
plot(fitmon)
if(fitmon$nobp$Vt==FALSE)
{
  
  date.trend2<-as.integer(fitmon$output[[1]]$Vt.bp)
  
  or4240p1t[i,j, date.trend2]<-100
}
if(fitmon$nobp$Wt==FALSE)
{
  date.seasonal2<-as.integer(fitmon$output[[1]]$Wt.bp)
  or4240p1s[i,j, date.seasonal2]<-500
  
}
o1<-bfast(ndvi,h=0.15,season="dummy", max.iter=1)
o1$output[[1]]$Vt.bp
z <- array(1:24, dim = 2:4)

zseq <- apply(z, 1:2, function(x) list(seq_len(max(x)),sum(x)))
zseq
str(app1)
str(zseq)
i4=i4+1
 seq_len(50)
print(i4)
}

}

length(which(!is.na(ora150p005s)))
save(ora150p005s, file='ora150p005s.Rdata')
save(ora150p005t, file='ora150p005t.Rdata')

bmonit<- bfastmonitor(datamon,h=0.25,start=2007,history='BP')

a<-structure(1:6, dim = 2:3)
b<-c(1,1)
dim(a)

a[,4,]<-c(1:10)
str(a)
library(bfast)
library(zoo)
NDVIb <- as.ts(zoo(som$NDVI.b, som$Time))
plot(NDVIb)
NDVIb

library('devtools')
install_github("strucchange","mengluchu")
monb <- bfastmonitor(NDVIb, start = c(2010, 13))
monb
bmonit<- bfastmonitor(NDVIb,h=0.25,start= c(2010, 13),history='ROC') #works
bmonit<- bfastmonitor(NDVIb,h=0.25,start=c(2010, 13),history='ROC',type = "OLS-CUSUM") #MOSUM has limitations, cusum doesnt
bmonit<- bfastmonitor(NDVIb,h=0.15,start=c(2010, 13),history='ROC',type = "OLS-MOSUM") # detect change later
plot(bmonit)

moni<-window(NDVIb, start=c(bmonit$history[1]), bmonit$monitor[2]) 
moni01<-bfast01(moni)
cl<-bfast01classify(moni01, alpha = 0.05, pct_stable = 0.25)
cl$flag_type
 
data("Nile")
plot(Nile)
plot(as.zoo(a))
## test the null hypothesis that the annual flow remains constant
## over the years
fs.nile <- Fstats(Nile ~ 1)
str(fs.nile)
sctest(fs.nile)
summary(monb$model)
plot(monb$mefp, functional = NULL)

monb<-bfastmonitor(NDVIb, start = c(2010, 13),
                   history = c(2002, 7), order = 1, plot = TRUE,formula = response ~ harmon)
#select 2nd order seasonality, set stable period
str(bmonit)
bmonit$breakpoint
bmonit$mefp$breakpoint

plot(bmonit)
str(fitmon$output)

install.packages("bfast", repos="http://R-Forge.R-project.org")
spacelx42<-array(,c(42,40,637))    
spacely40<-array(,c(42,40,637)) # takes far more time!
fspacelx42<-array(,c(42,40,637))    
f2spacely40<-array(,c(42,40,637)) 

fspacelx150<-array(,c(150,150,637))    
fspacely150<-array(,c(150,150,637)) 
for (i in 1:636)
{  
  for (j in 1:42)
  { 
    
    spt1<-fity2[j,,i] 
    ti<-time(spt1)
    p.val<-sctest(efp(spt1 ~ ti, h = h, type =  "OLS-CUSUM"))$p.value
    
    #plot(breakpoints(spt1 ~ ti))
    
    # fm1<-lm(spt1~breakfactor(breakpoints(spt1 ~ ti)))
    # Tt<- ts(fitted(fm1))
    #par(mfrow=c(2,1))
    #plot(spt1)
    #plot(Tt)
    
    if(p.val<0.05)
    {
      
      bp<-breakpoints(spt1 ~ ti)$breakpoints    
      f2spacely40[j,bp,i]<-100
    }
    
    i4=i4+1
    print(i4)
  }
  
}

for (i in 1:636)
{  
  for (j in 1:150)
  { 
    
    spt1<-fevi8[,j,i] 
    ti<-time(spt1)
    p.val<-sctest(efp(spt1 ~ ti, h = h, type =  "OLS-CUSUM"))$p.value
    
    #plot(breakpoints(spt1 ~ ti))
    
    # fm1<-lm(spt1~breakfactor(breakpoints(spt1 ~ ti)))
    # Tt<- ts(fitted(fm1))
    #par(mfrow=c(2,1))
    #plot(spt1)
    #plot(Tt)
    
    if(p.val<0.05)
    {
      
      bp<-breakpoints(spt1 ~ ti)$breakpoints    
      fspacelx150[bp,j,i]<-100
    }
    
    i4=i4+1
    print(i4)
  }
  
}

for (i in 1:636)
{  
  for (j in 1:150)
  { 
    
    spt1<-fevi8[j,,i] 
    ti<-time(spt1)
    p.val<-sctest(efp(spt1 ~ ti, h = h, type =  "OLS-CUSUM"))$p.value
    
    #plot(breakpoints(spt1 ~ ti))
    
    # fm1<-lm(spt1~breakfactor(breakpoints(spt1 ~ ti)))
    # Tt<- ts(fitted(fm1))
    #par(mfrow=c(2,1))
    #plot(spt1)
    #plot(Tt)
    
    if(p.val<0.05)
    {
      
      bp<-breakpoints(spt1 ~ ti)$breakpoints    
      fspacely150[j,bp,i]<-500
    }
    
    i4=i4+1
    print(i4)
  }
  
}
save(fspacelx150, file='fspacelx150.Rdata')
save(fspacely150, file='fspacely150.Rdata')

save(fspacelx42, file='fspacelx42.Rdata')
save(f2spacely40, file='fspacely40.Rdata')

match(aaa,fspacelx42)
fspacelx42-aaa
which(!is.na(f2spacely40))
save(ora150p005t, file='ora150p005t.Rdata')
length(which(!is.na(fspacely150)))

save(ora150p005s, file='ora150p005s.Rdata')
save(ora150p005t, file='ora150p005t.Rdata')

