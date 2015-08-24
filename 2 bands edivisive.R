 
iquery("
       store(
          subarray(
                  repart(
                           MOD09Q1_JUARA,<red:int16,nir:int16,quality:uint16> [col_id=58828:59679,502,0,row_id=48103:49050,502,0,time_id=0:9200,1,0]
                        ),
                          58930,48210,6,59079,48359,643
                   ),
       sub150)",return=TRUE)

iquery("
       store(
       project(
       apply(      
       repart(
    
       subarray(
       repart(
       MOD09Q1_JUARA,
       <red:int16,nir:int16,quality:uint16> [col_id=58828:59679,502,0,row_id=48103:49050,502,0,time_id=0:9200,1,0]),
       
       58930,48210,6,59079,48359,643),
      
       <red:int16,nir:int16,quality:uint16,evi2:double>[col_id=0:149,1,1,row_id=0:149,1,1,time_id=0:637,638,0]),
       the_i,double(col_id),the_j,double(row_id),the_t,double(time_id)),
       the_i, the_j,the_t, red,nir),
       twobands1)")

iquery("dimensions(twobands1)",return=TRUE)
eviinus6<-scidb('twobands1')
ired<-eviinus6[0:149,0:149,]$red[]
inir<-eviinus6[0:149,0:149,]$nir[]
red1<-array(ired,c(150,150,636))
nir1<-array(inir,c(150,150,636))
 
red1a<-aperm (red1,c(2,1,3)) #rotate # new second array
nir1a<-aperm (nir1,c(2,1,3))
 
save(nir1a,file='nir1a.Rdata')
save(red1a,file='red1a.Rdata')


aggrefun <- function(x,dates=a1,aggre='month') 
{
  
  if (aggre  == "month")
  {  
    spt<-zoo(x,dates)
    monmean <- aggregate(spt, as.Date(as.yearmon(dates)), mean)
    
    frequency(monmean)<-12
    na.new <- function(x) ts(na.exclude(x), frequency = 12)
    
    stlmon<-stl(monmean, na.action = na.new, s.window = "per")
    
    datamon <- ts(rowSums(stlmon$time.series)) 
  }
}

redagg<-apply(red1a ,c(1,2), aggrefun)
niragg<-apply(nir1a ,c(1,2), aggrefun)

redagg<-redagg[,,-1]
niragg<-niragg[,,-1]

plot(red8mon[1,1,],typ='l',col='red')
 
red8mon<-aperm(red8agg, c(2,3,1))
nir8mon<-aperm(nir8agg, c(2,3,1))
lines(nir8mon[1,1,])
save(nir8mon,file='nir8mon.Rdata')
save(red8mon,file='red8mon.Rdata')
  
edivisive2bands<-array(,c(150,150,167))
 
# try to use 2 bands: red and nir, deseasonlised woth loess
 for(i in 1:150)
{
  for (j in 1:150)
  {
    
    #fit1<-bfast(fevi3b3,  h=0.15 , season = "harmonic" , max.iter = 1,   level = 0.05, type= "OLS-MOSUM")
 
    minsize=0.15*167
    nir1<-ts(nir8mon[i,j ,],frequency=12)
    red1<-ts(red8mon[i,j ,],frequency=12)
    stlred1<-red1-stl(red1, s.window = "per")$time.series[,"seasonal"]
    stlnir1<-nir1-stl(nir1, s.window = "per")$time.series[,"seasonal"]
    
    plot(red1)
    wm<- as.matrix(na.omit(cbind(stlnir1, stlred1)))
    e.di1<-e.divisive(wm,min.size=minsize) 
    
    if(length( e.di1$estimates)>2)
    {         
      date.edi<-e.di1$estimates[-1]
      date.edi<-date.edi[-length(date.edi)]
      edivisive2bands[i,j,date.edi]<-200
      
      print(paste('si:',i,'sj:',j))
      
      # jpeg(paste('i',i,'j',j,"edivisive vs bfast.jpg"), height=4, width=6, res=800,unit="in")
      plot(nir1 ,ylab="evi2 monthly",typ='l')
      lines(red1,col='red' )
      abline(v=date.edi,col='blue')
      # abline(v= date.trend1/12+1,col='red')
      
      #  legend('bottomleft',c("e.divisive",'bfast'),lty=1,col=c('blue','red'))
      # dev.off()
    }
  }
}
save(edivisive2bands,file="edivisive2bands.Rdata")

