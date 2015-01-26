
iquery('store(subarray(MOD09Q1_JUARA,58828,48103,6,58857,48132,643),MODsmall1)') #30by30
#scidblist()
iquery("store(apply(MODsmall1,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),eviM2s1) ")

iquery(" apply(eviM2s,evif1, iff((evi2 <1) or( evi2 >-1), evi2, window(eviM2s,1,1,1,1,1,1, median(evi2 ))))")
iquery(" window(eviM2s1,1,1,1,1,1,1, median(evi2 ))",return=TRUE)
iquery("list('aggregates')",return=TRUE)
eviM2s1re<-scidb("eviM2s1")
mean(eviM2s1re)
iquery("store(repart(eviM2s, <red:int16,nir:int16,quality:uint16,evi2:double> 
       [col_id=0:848,3,2,row_id=0:944,3,2,time_id=0:637,638,0]),eviw5s)");
#(49048-48103) 
# change chunk size
iquery("store(apply(eviw5s,the_i,double(col_id),the_j,double(row_id),the_t,double(time_id)),eviw6s) ")


iquery("store(project(eviw6s, the_i, the_j,the_t, evi2), eviw7s)")
#iquery('store(subarray(eviw4,0,0,0,6,6,637),test1)')

iquery("dimensions(MODsmall1)",return=TRUE)
#iquery("dimensions(eviM1)",return=TRUE)
ptm <- proc.time()

iquery("store(r_exec(eviw5s,         
       'output_attrs=5',
       'expr= 
       #require(bfast)
       require(strucchange)
       require(spdep)
       load(listcn636.Rdata)
       eviarray<-array(evi2,c(637,3,3))          # scidb array to r array                
       eviarray2<-aperm(eviarray,c(2,3,1))       # rotate r array 
       eday <- as.Date(\"2000-01-30\")           # date 
       e8day <- seq(eday, length.out=637, by=\"8 days\")
      
  
       fevi3b312t<-apply(eviarray2,c(1,2),function(x) (ts(x,start=c(2000,1),frequency=46)-stl(ts(x,start=c(2000,1),frequency=46),'per')$time.series[,\"seasonal\"]))
      
       f2<-aperm(fevi3b312t,c(2,3,1))
       aa2<-as.data.frame.table(f2)$Freq
       aa3<-as.data.frame(as.data.frame.table(f2)$Freq)
   
       try2<-spautolm(aa2~c(rep(c(1:636),each=9)),family=\"SAR\",listw=listcn636)
 
       rn<-lapply(1:9,function(i) {residuals(try2)[seq(i,636*9-(9-i),9)]})
 
     
       fevi3b312t1<-ts(f2[2,2,],start=c(2000,1),frequency=46)
       p.Vt2 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = \"OLS-CUSUM\" )) 
       p.Vt3 <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = \"OLS-MOSUM\" )) 
       p.Vt <- sctest(efp(fevi3b312t1 ~ ti, h = 0.15, type = \"OLS-CUSUM\", spatial1=as.numeric(rn[[5]]))  )
       
       p1[i,j]<- p.Vt$p.value  
       p2[i,j]<-p.Vt2$p.value  
       p3[i,j]<-p.Vt3$p.value  
 
       
       
       list( p1,p2,p3,the_i, the_j)'), sarefp1)",return=TRUE)  

#the_j,the_i,the_t, 

proc.time() - ptm
install_github("strucchange","mengluchu",build_vignettes = FALSE)
