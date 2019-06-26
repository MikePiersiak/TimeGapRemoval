### This script is designed to remove time gaps in telemetry data for 
### application of HMMs to time-standardized data. Not removing time 
### gaps leaves large sections of data with linear interpolations,
### which inflates data with values of '1' for angle concentration 
### parameters which largely inhibits model convergence. These sections
### need to be removed to get an accurate biological picture of underlying
### states dictating animal movement



crab.dat <- read.csv("CrabData.csv")
str(crab.dat)

crab.dat$Time <- as.POSIXct(crab.dat$Time, format = "%m/%d/%Y %H:%M", tz="UTC")
crab.dat$release.time <- as.POSIXct(crab.dat$release.time, format = "%m/%d/%Y %H:%M", tz="UTC")
crab.dat$TRANSMITTER <- as.factor(crab.dat$TRANSMITTER)


##change name of column 'TRANSMITTER' to 'ID' for later imputation in momentuHMM (as well as X/Y cols to 'x' and 'y' and col(Time) to 'time)
names(crab.dat)[names(crab.dat) == "TRANSMITTER"] <- "ID"
names(crab.dat)[names(crab.dat) == "X"] <- "x"
names(crab.dat)[names(crab.dat) == "Y"] <- "y"
names(crab.dat)[names(crab.dat) == "Time"] <- "time"

###crawlWrap fails if wrapper function creates track with <100 detection locations
###need to tell R to remove/ignore tags based on total time detected : standardized time interval. 
###if this results in fewer than 100 points, tag must be removed as HMMs cannot be fit to tags w/
###less than 100 detections

###t.int = time interval to regularize telemetry data over (desired time step between detections
###for use via HMMs) UNITS ARE HOURS 
t.int <- 1

##this loop determines the number of points in a predicted track with a pre-stipulated time interval
##between detections (t.int)
nDetections <- NULL

for (i in levels(crab.dat$ID)){
  
  Tag <- subset(crab.dat, ID==i)
  tm.dur <- as.numeric(difftime((max(Tag$time)),(min(Tag$time)),units = "hours"))
  nPoints <- tm.dur/t.int
  sheet <- cbind(i,nPoints)
  nDetections <- rbind(nDetections, sheet)

}

nDetections <- as.data.frame(nDetections)
nDetections$nPoints <- as.integer(nDetections$nPoints)

##add Column of 0's and 1's to distinguish good (100+ detections) and bad (<100 detections) tags
nDetections$Boolean <- ifelse(nDetections$nPoints < 101,0,1)
colnames(nDetections) <- c("ID","nPoints","Boolean")

##Boolean operator--0's indicate bad tags
##1's indicate good tags (100+ detections)

bad.tags <- nDetections[which(nDetections$Boolean==0),]
good.tags <- nDetections[which(nDetections$Boolean==1),]

##Remove Tags w/ sub 100 detections from larger data set

a <- bad.tags$ID
crab.dat <- subset(crab.dat, ! ID %in% a)

##double check to make sure factor levels and assciated data were removed
##should match value from ln 59 (number of good tags)
droplevels(crab.dat$ID)


###DISREGARD BELOW. NOT FINISHED YET.
test <- subset(crab.dat, ID=="31288")
str(test)

head(test)
test$seq <- as.factor(seq(1:nrow(test)))

time.difference <- NULL
time.diff.sheet <- NULL

##may need to tell R to start at Row 2(i.e. avoid using header and spitting out NAs)
for (m in levels(test$seq))
{
  Row <- subset(test, seq==m)
    
  RowAdj <- subset(test, seq == (m+1))
  
  ##probably add an ifelse here telling it to stop at max number of rows in the data sheet (or for that particular tag)
  
  a <- difftime(RowAdj$time, Row$time, tz = "UTC", units="hours")
  time.diff.sheet <- as.data.frame(a)
  time.difference <- rbind(time.difference, time.diff.sheet)
}

