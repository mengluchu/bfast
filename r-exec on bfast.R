## test r-exec 
## bfast input as double, restore to 3-d r array 
## bfast works on 3-d r array per-pixel
## output 3-d array with breakpoint magnitude (or other variables)
## output as double, restore to 3-d scidb array

iquery("show(evi21b21ijt3)",return=TRUE)
# evi21b21ijt3<the_i:double,the_j:double,the_t:double,evi2:double> [col_id=0:20,3,0,row_id=0:20,3,0,time_id=0:637,638,0]
# 21 by 21 by 591/(638)  array

ptm <- proc.time()

iquery("store(r_exec(evi21b21ijt3,         
       'output_attrs=4',
       'expr= 
        require(bfast)
        eviarray<-array(evi2,c(638,3,3))          # scidb array to r array                
        eviarray2<-aperm(eviarray,c(2,3,1))       # rotate r array 
        eday <- as.Date(\"2000-01-30\")           # date 
        e8day <- seq(eday, length.out=635, by=\"8 days\")
        output2 <-bfaarray2(eviarray2,dates=e8day,aggre= \"month\",season= \"harmonic\",max.iter=1,level=0.05)       
        allbreakpoints<-output2[1:6,,]            #breakpoint time
        allmagnitudes<-output2[7:12,,]            #breakpoint magnitude
       
        mtime<-allbreakpoints  
   
        dimx<-dim(output2)[2]
        dimy<-dim(output2)[3]
        dimt<-591
       
        newarray<-array(,c(dimx,dimy,dimt))  
       
        mag<- allmagnitudes          
        t3darrbfamul<-tasd.bfa.mul(newarray,mtime,mag) # time as a 3rd dimension
        t3darrbfamul[which(is.na(t3darrbfamul))]<-0    # na to 0
       
       t3df<-as.data.frame.table (t3darrbfamul)$Freq   # output double 
       
      
       
       list( the_j,the_i,the_t,t3df  )'),rexec21by212)",return=TRUE)       

proc.time() - ptm

# rebuild scidb array 
  
iquery("store(attribute_rename(project(unpack(rexec21by21, tmpDim), 
        expr_value_0,expr_value_1,expr_value_2,expr_value_3),
        expr_value_0, thei, expr_value_1, thej, expr_value_2, thet, expr_value_3,value), re2121)");


#Restore the dimensions (attribute double to int64, then redimension)
iquery("store(project(apply(re2121, i, int64(thei),
       j, int64(thej), t, int64(thet)), i, j, t, value), rexec2121)");

iquery("store(redimension(rexec2121, <value:double> [i=0:20,3,0,j=0:20,3,0,t=0:637,638,0]), rexec2121re)");
#rexec2121re is the 3d scidb array with bfast results 