 
#explored scidb overlay, scidb overlay overlay on the diagnal pixel as well.in the following example it will replicat each pixel 
#4 6 4
#6 9 6
#4 6 4
#times
iquery("store(build(<evi2:double> 
  [col_id=0:2,1,1,row_id=0:2,1,1,time_id=0:9,10,0] , col_id*100+row_id*10+time_id*0.1),tryoverlap9)")
iquery("store(apply(tryoverlap9,the_i,double(col_id),the_j,double(row_id),the_t,double(time_id)),tryoverlap10) ")
iquery("show(ev22ov)",return=TRUE)

#evi2t2121overlap4
#iquery('store(subarray(evi2t2121overlap4,0,0,0,2,2,10),ev22ov1)')

#x=c(58930:59079),y=c(48210:48359)
iquery("
store(
  project(
      apply(      
        repart(
            apply(
              subarray(
              
              MOD09Q1_JUARA,58930,48210,6,59079,48359,643),
              evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),
              <red:int16,nir:int16,quality:uint16,evi2:double>[col_id=0:149,1,1,row_id=0:149,1,1,time_id=0:637,638,0]),
              the_i,double(col_id),the_j,double(row_id),the_t,double(time_id)),
              the_i, the_j,the_t, evi2),
              repro1)")


iquery(" r_exec(repro1, 'output_attrs=8',

'expr= 
 
library(strucchange)
library(spdep)
library(nlme)

load(\"/home/menglu/fevi8.Rdata\") # coefficient matrix 
load(\"/home/menglu/listcn636.Rdata\") # neighbor


dim1<-length(unique(the_i)) 
dim2<-length(unique(the_j)) 
dim3<-length(unique(the_t))
tl=1:dim3
w=1/46       
co <- cos(2*pi*tl*w)
si <- sin(2*pi*tl*w)
co2 <- cos(2*pi*tl*w*2)
si2 <- sin(2*pi*tl*w*2)  
co3 <- cos(2*pi*tl*w*3)
si3 <- sin(2*pi*tl*w*3) 
newarray<-array(evi2,c(dim1,dim2,dim3))          # scidb array to r array                

 
      if((dim1+dim2)==6)
       {  
       
       fevi3b3<-newarray[2,2,]    
      
       fevi3b3[is.na(fevi3b3)] <- 0
       
       if(length(which(fevi3b3==0))<100)
          {
            aa2<-as.vector(fevi3b3)
            aa2[aa2==0]<-NA
            fevi3b3[fevi3b3==0]<-median(fevi3b3)
            fevi3b312t1<-ts(fevi3b3[2,2,],start=c(2000,1),frequency=46) # reconstruct the time series
 
       try2<-spautolm(aa2~. , data.frame(aa2,X),family=\"SAR\",method= \"Matrix\", listw=listcn636,na.action=na.exclude,zero.policy=TRUE)
 
       rn<-lapply(1:9,function(i) {residuals(try2)[seq(i,636*9-(9-i),9)]})
    
       
      ii<-5   # get the middle pixel (5 for 3*3 matrix)
   
      resar1<-coredata(residuals(gls(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,correlation=corAR1())))
       
       p.Vt1  <- sctest(efp(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,    h = 0.15, type = \"OLS-CUSUM\", spatial1=as.numeric(rn[[ii]]))  )
       p.Vt2  <- sctest(efp(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,    h = 0.15, type = \"OLS-MOSUM\", spatial1=as.numeric(rn[[ii]]))  )
       p.Vt3 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" )) 
      p.Vt4 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ))
      p.Vt5 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" ,spatial1=as.numeric(resar1)) ) 
      p.Vt6 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ,spatial1=as.numeric(resar1)) )
       
      spcusum1  <-p.Vt1$p.value # spautolm residuals CUSUM
      spmosum1  <-p.Vt2$p.value # spautolm residuals  MOSUM  
       cusum1    <-p.Vt3$p.value # CUSUM
       mosum1    <-p.Vt4$p.value # MOSUM
       cusumar1 <-p.Vt5$p.value # CUSUM ar 1
       mosumar1 <-p.Vt6$p.value # MOSUM ar1
        rcol<-(min(the_i)+1)
       rrow <-(min(the_j)+1)  
       }
}
list(spcusum1,spmosum1,cusum1,mosum1,cusumar1,mosumar1,rcol,rrow)
   
 ') ",
       return=TRUE
       )

 

 





iquery(" r_exec(repro1, 'output_attrs=1',
       
       'expr= 
       
       library(strucchange)
       library(spdep)
       library(nlme)
       
       load(\"/home/menglu/fevi8.Rdata\") # coefficient matrix 
       load(\"/home/menglu/listcn636.Rdata\") # neighbor
       
       
       dim1<-length(unique(the_i)) 
       dim2<-length(unique(the_j)) 
       dim3<-length(unique(the_t))
       tl=1:dim3
       w=1/46       
       co <- cos(2*pi*tl*w)
       si <- sin(2*pi*tl*w)
       co2 <- cos(2*pi*tl*w*2)
       si2 <- sin(2*pi*tl*w*2)  
       co3 <- cos(2*pi*tl*w*3)
       si3 <- sin(2*pi*tl*w*3) 
       newarray<-array(evi2,c(dim1,dim2,dim3))          # scidb array to r array                
       
       
       if(dim1<3 || dim2<3)
{
       
       fevi3b312t1<-ts(newarray[1,1,],start=c(2000,1),frequency=46) # reconstruct the time series with its own
       
       resar1<-coredata(residuals(gls(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,correlation=corAR1())))
       
       p.Vt3 <- sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" )) 
       p.Vt4 <- sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ))
       p.Vt5 <- sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" ,spatial1=as.numeric(resar1)) ) 
       p.Vt6 <- sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ,spatial1=as.numeric(resar1)) )
       
       spcusum1<-p.Vt3$p.value # spautolm residuals CUSUM
       spmosum1 <-p.Vt4$p.value # spautolm residuals  MOSUM  
       cusum1 <-p.Vt3$p.value # CUSUM
       mosum1<-p.Vt4$p.value # MOSUM
       cusumar1 <-p.Vt5$p.value # CUSUM ar 1
       mosumar1 <-p.Vt6$p.value # MOSUM ar1
       rcol<-min(the_i)*100
       rrow<-min(the_j)*100
} else {  
       
       fevi3b3<-newarray[2,2,]    
       fevi3b3[is.na(fevi3b3)] <- 0
       
       aa2<-as.vector(fevi3b3)
       aa2[aa2==0]<-NA
       
       fevi3b312t1<-ts(fevi3b3[2,2,],start=c(2000,1),frequency=46) # reconstruct the time series
       
       try2<-spautolm(aa2~. , data.frame(aa2,X),family=\"SAR\",method= \"Matrix\", listw=listcn636,na.action=na.exclude,zero.policy=TRUE)
       
       rn<-lapply(1:9,function(i) {residuals(try2)[seq(i,636*9-(9-i),9)]})
       #get residuals for each time series
       
       ii<-5   # get the middle pixel (5 for 3*3 matrix)
       
       resar1<-coredata(residuals(gls(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,correlation=corAR1())))
       
       p.Vt1  <- sctest(efp(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,    h = 0.15, type = \"OLS-CUSUM\", spatial1=as.numeric(rn[[ii]]))  )
       p.Vt2  <- sctest(efp(fevi3b312t1 ~ tl+co+co2+co3+si+si2+si3,    h = 0.15, type = \"OLS-MOSUM\", spatial1=as.numeric(rn[[ii]]))  )
       p.Vt3 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" )) 
       p.Vt4 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ))
       p.Vt5 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-CUSUM\" ,spatial1=as.numeric(resar1)) ) 
       p.Vt6 <-  sctest(efp(fevi3b312t1 ~  tl+co+co2+co3+si+si2+si3,   h = 0.15, type = \"OLS-MOSUM\" ,spatial1=as.numeric(resar1)) )
       
       spcusum1  <-p.Vt1$p.value # spautolm residuals CUSUM
       spmosum1  <-p.Vt2$p.value # spautolm residuals  MOSUM  
       cusum1    <-p.Vt3$p.value # CUSUM
       mosum1    <-p.Vt4$p.value # MOSUM
       cusumar1 <-p.Vt5$p.value # CUSUM ar 1
       mosumar1 <-p.Vt6$p.value # MOSUM ar1
        rcol<-(min(the_i)+1)
       rrow <-(min(the_j)+1)  
 
}
list(0.1)
   
 ') ",
       return=TRUE
       )
spcusum1,spmosum1,cusum1,mosum1,cusumar1,mosumar1,rcol,rrow


#tt2<-scidb("tt3")
#allx<-summary(tt2[][][])$x


min(1)








iquery('store(subarray(repro1,1,1,1,3,3,10),srepro2)')
#tryoverlay10
iquery("r_exec(srepro2, 'output_attrs=2',       
       'expr= 
        dim1<-length(evi2))*1.0 
        list(  dim1 )')",
        return=TRUE
       )



iquery("show(tryoverlap10)",return=TRUE)

min(c(1,2))
ddd<-scidb("ex2")

length(ddd[,,][1])
table(summary(ddd[,,][][])$x)

2548/4
3822/6
5733/9
#something with subarray is wrong?
summary(ddd[])
#problem1 if go with subarray then the chunk is not right