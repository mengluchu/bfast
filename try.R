example(NY_data)
lm0 <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata)
summary(lm0)
lm0w <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata, weights=POP8)
summary(lm0w)
esar0 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
                    listw=listw_NY)
summary(esar0)
system.time(esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                               data=nydata, listw=listw_NY, family="SAR", method="eigen", verbose=TRUE))
summary(esar1f)
system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                               data=nydata, listw=listw_NY, family="SAR", method="Matrix", verbose=TRUE))
summary(esar1M)
library(devtools)
install_github('dutri001/bfastSpatial')

data(tura)
tura
cn<-cell2nb(7,7, type ="queen",torus =FALSE)
 
library(devtools)
install_github('dutri001/bfastSpatial')
# load the package
library(bfastSpatial)

install.packages("C:/Users/m_lu0002/Dropbox/mengluchu/try01.zip", repos = NULL,type='source')
library(try01)

library(xts)
 
 

# Generate dummy time series
from <- as.Date("1950-01-01")
to <- as.Date("1990-12-31")
myDates <- seq.Date(from=from,to=to,by="day")
myTS <- as.xts(runif(length(myDates)),order.by=myDates)

# SPLIT THE TIME SERIES INTO CALENDAR YEARS
myList <- tapply(myTS, format(myDates, "%Y"), c)
plot(myList[[1]])

?tapply
# SPLIT THE TIME SERIES INTO HYDROLOGICAL YEARS

# calculate the number of hydrological years
nHY <- length(split(myTS[.indexmon(myTS) %in% 0:8], f="years"))-1 

# create an empty table , to be populate by a loop
myList <- list()

for ( counter in 1:nHY ){
  
  oct2dec <- split(myTS[.indexmon(myTS) %in% 9:11], f="years")[[counter]]
  jan2sep <- split(myTS[.indexmon(myTS) %in% 0:8], f="years")[[counter + 1]]
  
  myList[[counter]] <- rbind(oct2dec, jan2sep)
  
}

# Access the element of the list using the index, e.g. for plotting:
plot(myList[[1]])
data(tura)
dates
bfm09 <- bfmSpatial(b,dates=dates,start=7,end=10, history=5, order=1,aggre=NULL)

nba<-bfast(ndvi,h=0.15,max.iter=1)
str(nba)
?breakpoints
(Vt ~ breakfactor(bp.Vt)/ti
bn<-breakpoints(ndvi~1)
nba$Mags

l1<-lm(ndvi ~ breakfactor(bn)/time(ndvi))
extract(str(l1))
coefficients(l1)
# extract change
plot(2)
dev.off()
dev.new()
plot(fitted(l1))
 
change09 <- raster(bfm09, 1)
plot(change09)
months <- changeMonth(change09)
# set up labels and colourmap for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)
