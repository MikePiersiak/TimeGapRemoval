### This script is designed to remove time gaps in telemetry data for 
### application of HMMs to time-standardized data. Not removing time 
### gaps leaves large sections of data with linear interpolations,
### which inflates data with values of '1' for angle concentration 
### parameters which largely inhibits model convergence. These sections
### need to be removed to get an accurate biological picture of underlying
### states dictating animal movement

setwd("/Users/PiersiakM/Desktop/Snow Crab Behavior")

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
###for use via HMMs)
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

head(nDetections)
##add Column of 0's and 1's to distinguish good (100+ detections) and bad (<100 detections) tags
#nDetections$Boolean <- as.integer(nDetections$Boolean)
nDetections$Boolean <- ifelse(nDetections$nPoints < 101,0,1)
colnames(nDetections) <- c("ID","nPoints","Boolean")
str(nDetections)
##Boolean operator--0's indicate bad tags
##1's indicate good tags (100+ detections)

bad.tags <- nDetections[which(nDetections$Boolean==0),]
good.tags <- nDetections[which(nDetections$Boolean==1),]

##Remove Tags w/ sub 100 detections from larger data set
a <- bad.tags$ID
crab.dat.rm <- subset(crab.dat, ! ID %in% a)

##leaves 144 tags w/ >100 detections. This is a prerequisite for application of HMMs
##We're only looking for good tags here so the nested loops below identify the time duration
##between concurrent detections for each individual crab.

crab.dat.rm$ID <- droplevels(crab.dat.rm$ID)

all.tags <- NULL
tag.omnibus <- NULL

for (i in levels(crab.dat.rm$ID))
{
  Tag <- subset(crab.dat.rm, ID==i)
  Tag$seq <- as.factor(seq(1:nrow(Tag)))
  
  tag.omnibus <- NULL
  for (m in levels(Tag$seq)){
    Row <- subset(Tag, seq==m)
    x <- as.numeric(m) + 1
    if (x > max(as.numeric(Tag$seq))) {
      RowAdj <- subset(Tag, seq==m)
    } else {
      RowAdj <- subset(Tag, seq ==(as.numeric(m)+1))
    }
    
    TimeDuration <- difftime(RowAdj$time, Row$time, tz="UTC", units = "hours")
    
    tag.sheet <- cbind(i,m,TimeDuration)
    tag.omnibus <- rbind(tag.omnibus, tag.sheet)
  
  }

  all.tags <- rbind(all.tags, tag.omnibus)
    
}

all.tags <- as.data.frame(all.tags)
head(all.tags)

##don't know why this is in minutes but needed to convert back to hours
all.tags$TimeDuration <- (as.numeric(all.tags$TimeDuration)) / 60

colnames(all.tags) <- c("ID","DetectionNumber","TimeDuration")

##The loop below and subsequent two lines will identify and remove tags exceeding a temporal duration
##stipulated in ln117
new.tags <- NULL

for (l in levels(all.tags$ID)) {
  
  Tag <- subset(all.tags, ID==l)
  
  ##this statement here stipulates the time interval that we want to split our movement tracks around
  ##i.e. how much time has to elapse before we split the track. Here I went with four hours.
  Tag$boolean <- ifelse(Tag$TimeDuration > 4, "1", "0")
  Tag$maxBoolean <- max(Tag$boolean)
  new.tags <- rbind(new.tags, Tag)
  
}

test <- subset(new.tags, maxBoolean=="0")
test$ID <- droplevels(test$ID)
str(test)
