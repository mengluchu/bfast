
#all array 1; get the whole array of evi2 with array dimensions also as attribute, array chunk overlap is 2

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

#no overlap array
iquery("store(repart(allarray1, <red:int16,nir:int16,quality:uint16,evi2:double>
       [col_id=0:851,3,0,row_id=0:947,3,0,time_id=0:637,638,0]),allarray1no)");

#iquery("store(subarray(allarray1no,0,0,0,851,947,635),all851947o)")
#get 150 by 150 array ( no overlap)
iquery("store(subarray(allarray1no,0,0,0,150,150,635),all150150o)")

iquery("store(r_exec(all150150o,         
       'output_attrs=6',
       'expr= 
       require(bfast)
       require(strucchange)
       
       eviarray<-array(evi2,c(635,3,3))          # scidb array to r array                
       eviarray2<-aperm(eviarray,c(2,3,1))       # rotate r array 
       eday <- as.Date(\"2000-01-30\")           # date 
       e8day <- seq(eday, length.out=635, by=\"8 days\")
       output2 <-bfaarray2(eviarray2,dates=e8day,aggre= \"month\",season= \"harmonic\",max.iter=1,level=0.05)       
       allbreakpoints<-output2[1:6,,]            #breakpoint time
       allmagnitudes<-output2[7:12,,]            #breakpoint magnitude
       breakpointclass<-output1[13:18,,] # class
       mtime<-allbreakpoints #timearr
       
       dimx<-dim(output2)[2]
       dimy<-dim(output2)[3]
       dimt<-635
       
       newarray<-array(,c(dimx,dimy,dimt)) #newarr
       
       t3darrbfamul1<-tasd.bfa.mul(newarray,mtime,allmagnitudes) # 3 dimensional array with magnitude     
       t3darrbfamul2<-tasd.bfa.mul(newarray,mtime, breakpointclass) # 3 dimensional array with class
       t3darrbfamul3<-tasd.bfa.mul(newarray,mtime, breakpoints) #  time can also be attributes
       
       t3darrbfamul[which(is.na(t3darrbfamul1))]<-0    # na to 0
       t3darrbfamul[which(is.na(t3darrbfamul2))]<-0    
       t3darrbfamul[which(is.na(t3darrbfamul3))]<-0    
       
       t3df<-as.data.frame.table (t3darrbfamul)$Freq   # output double 
       t3df2<-as.data.frame.table (t3darrbfamul2)$Freq   
       t3df3<-as.data.frame.table (t3darrbfamul3)$Freq   
       
       list( the_j,the_i,the_t,t3df,t3df2, t3df3  )'),rexectmc)",return=TRUE)       


# rebuild scidb array 

iquery("store(attribute_rename(project(unpack(rexectmc, tmpDim), 
       expr_value_0,expr_value_1,expr_value_20,expr_value_3,expr_value_4,expr_value_5),
       expr_value_0, thei, expr_value_1, thej, expr_value_2, thet, expr_value_3,value
       ,expr_value_4,value2,expr_value_5,value3), rexectmc2)");


#Restore the dimensions (attribute double to int64, then redimension)
iquery("store(project(apply(rexectmc2, i, int64(thei),
       j, int64(thej), t, int64(thet)), i, j, t, value, value2, value3), rexectmc3)");

iquery("store(redimension(rexectmc3, <value:double,value2:double,value3:double> [i=0:150,3,0,j=0:150,3,0,t=0:635,636,0]), rexectmc4)");
#rexec2121re is the 3d scidb array with bfast results 