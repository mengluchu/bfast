
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

fevi8agg<-apply(fevi8 ,c(1,2), aggrefun)
length(a1)
length(fevi8[1,1,])
fevi8<-fevi8[,,-1]

fevi8agg<-apply(fevi8 ,c(1,2), aggrefun)
fevi8mon<-aperm(fevi8agg, c(2,3,1))
save(fevi8mon,file='fevi8mon.Rdata')
load('fevi8mon.Rdata')
str(fevi8mon)
 
bfast1<-array(,c(150,150,636))
edivisive1<-array(,c(150,150,636))
 
for(i in 1:150)
{
  for (j in 1:150)
  {
    
    fevi3b3<-ts(fevi8mon[i,j,],frequency=12)
    fit1<-bfast(fevi3b3,  h=0.15 , season = "harmonic" , max.iter = 1,   level = 0.05, type= "OLS-MOSUM")
    minsize=0.15*length(fevi3b3)
 
    f2<-fevi8[i:(i+2),j:(j+2),]
    e.di1<-e.divisive(as.matrix(f2),min.size=minsize) 
   
    
    if(fit1$nobp$Vt==FALSE || length( e.di1$estimates)>2)
    {      
     
      date.trend1<-as.integer(fit1$output[[1]]$Vt.bp)
     
      date.edi<-e.di1$estimates[-1]
      date.edi<-date.edi[-length(date.edi)]
      
      bfast1[i,j, date.trend1]<-100
      edivisive1[i,j,date.edi]<-200
      
      print(paste('si:',i,'sj:',j))
      
     # jpeg(paste('i',i,'j',j,"edivisive vs bfast.jpg"), height=4, width=6, res=800,unit="in")
      plot(fevi3b3,ylab="evi2 monthly")
      abline(v=date.edi/12+1,col='blue')
     # abline(v= date.trend1/12+1,col='red')
      
    #  legend('bottomleft',c("e.divisive",'bfast'),lty=1,col=c('blue','red'))
  # dev.off()
    }
  }
}

save(bfast1,file='bfast1.Rdata')
save(edivisive1,file='edivisive1.Rdata')  

 
# try to use spatial neighbor as variables to estimate breakpoint together. should be comparible with sarbfast
edivisivesp<-array(,c(150,150,167))
for(i in 1:148)
{
  for (j in 1:148)
  {
   
    #fit1<-bfast(fevi3b3,  h=0.15 , season = "harmonic" , max.iter = 1,   level = 0.05, type= "OLS-MOSUM")
    minsize=0.15*167
    
    f2<-fevi8mon[i:(i+2),j:(j+2),]
    wm<- wrap(f2, map=list(3,NA))
    e.di1<-e.divisive(wm,min.size=minsize) 
  
    if(length( e.di1$estimates)>2)
    {         
      date.edi<-e.di1$estimates[-1]
      date.edi<-date.edi[-length(date.edi)]
      edivisivesp[i,j,date.edi]<-200
      
      print(paste('si:',i,'sj:',j))
      
      # jpeg(paste('i',i,'j',j,"edivisive vs bfast.jpg"), height=4, width=6, res=800,unit="in")
      plot(f2[1,1,],ylab="evi2 monthly",typ='l')
      abline(v=date.edi,col='blue')
      # abline(v= date.trend1/12+1,col='red')
      
      #  legend('bottomleft',c("e.divisive",'bfast'),lty=1,col=c('blue','red'))
      # dev.off()
    }
  }
}
save(edivisivesp,file="edivisivesp.Rdata")
