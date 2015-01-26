iquery('store(subarray(MOD09Q1_JUARA,58828,48103,6,59676,49047,643),MOD2)')
iquery("store(apply(MOD2,evi2,2.5*((nir*0.0001-red*0.0001)/(nir*0.0001+2.4*red*0.0001+1.0))),eviM2) ")

iquery("store(repart(eviM2, <red:int16,nir:int16,quality:uint16,evi2:double>
       [col_id=0:848,3,0,row_id=0:944,3,0,time_id=0:637,638,0]),eviw5)");
#(49048-48103) 

iquery("store(apply(eviw5,the_i,double(col_id),the_j,double(row_id),the_t,double(time_id)),eviw6) ")


iquery("store(project(eviw6, the_i, the_j,the_t, evi2), eviw7)")
#iquery('store(subarray(eviw4,0,0,0,6,6,637),test1)')

iquery("dimensions(eviw4)",return=TRUE)
#iquery("dimensions(eviM1)",return=TRUE)
ptm <- proc.time()

iquery("store(r_exec(eviw5,         
       'output_attrs=1',
       'expr= 
       require(bfast)
       eviarray<-array(evi2,c(637,3,3))          # scidb array to r array                
       eviarray2<-aperm(eviarray,c(2,3,1))       # rotate r array -> x y t
       eday <- as.Date(\"2000-01-30\")           # date 
       e8day <- seq(eday, length.out=637, by=\"8 days\")
       output2 <-bfaarray2(eviarray2,dates=e8day,aggre= \"month\",season= \"harmonic\",max.iter=1,level=0.05)       
       allbreakpoints<-output2[1:6,,]            #breakpoint time
       allmagnitudes<-output2[7:12,,]            #breakpoint magnitude
       
       mtime<-allbreakpoints  
       mag<- allmagnitudes
       
       dimx<-dim(output2)[2]
       dimy<-dim(output2)[3]
       dimt<-637
       
       newarray<-array(,c(dimx,dimy,dimt))         
       
       t3darrbfamul<-tasd.bfa.mul(newarray,mtime,mag) # time as a 3rd dimension
       t3darrbfamul[which(is.na(t3darrbfamul))]<-0    # na to 0
       
       t3df<-as.data.frame.table (t3darrbfamul)$Freq   # output double 
       
       
       list( t3df )'),eviwout)",return=TRUE)  
the_j,the_i,the_t, 

proc.time() - ptm
