### This script is designed to remove time gaps in telemetry data for 
### application of HMMs to time-standardized data. Not removing time 
### gaps leaves large sections of data with linear interpolations,
### which inflates data with values of '1' for angle concentration 
### parameters which largely inhibits model convergence. These sections
### need to be removed to get an accurate biological picture of underlying
### states dictating animal movement



crab.dat <- read.csv("CrabData.csv")
crab.test <- read.csv("CrabData.test.csv")
str(crab.test)

crab.dat$Time <- as.POSIXct(crab.dat$Time, format = "%m/%d/%Y %H:%M", tz="UTC")
crab.test$Time <- as.POSIXct(crab.test$Time, format = "%m/%d/%Y %H:%M", tz="UTC")

crab.dat$TRANSMITTER <- as.factor(crab.dat$TRANSMITTER)
crab.test$TRANSMITTER <- as.factor(crab.test$TRANSMITTER)

##change name of column 'TRANSMITTER' to 'ID' for later imputation in momentuHMM (as well as X/Y cols to 'x' and 'y' and col(Time) to 'time)
names(crab.dat)[names(crab.dat) == "TRANSMITTER"] <- "ID"
names(crab.dat)[names(crab.dat) == "X"] <- "x"
names(crab.dat)[names(crab.dat) == "Y"] <- "y"
names(crab.dat)[names(crab.dat) == "Time"] <- "time"

names(crab.test)[names(crab.test) == "TRANSMITTER"] <- "ID"
names(crab.test)[names(crab.test) == "X"] <- "x"
names(crab.test)[names(crab.test) == "Y"] <- "y"
names(crab.test)[names(crab.test) == "Time"] <- "time"

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

##doesn't work on initial subset of tags. going tag by tag to ID the bad tag(s) and figure out wtf
##isn't working

require(momentuHMM)
crwOut1 <- crawlWrap(obsData=crab.test, timeStep="hour",theta=c(6.855, -0.007), fixPar=c(NA,NA))

levels(crab.test$ID)                     
tag.remove <- subset(crab.test, ! ID %in% c(#"31289",
                                            #"31290",
                                            #"31292",
                                            #"31293",
                                            #"31295",
                                            #"31296",
                                            #"31298",
                                            #"31299",
                                            #"31300",
                                            "31301"))

###crawlWrap fails if wrapper function creates track with <100 detection locations
###need to tell R to remove/ignore tags based on total time detected : standardized time interval. 
###if this results in fewer than 100 points, tag must be removed as HMMs cannot be fit to tags w/
###less than 100 detections

###need to write a for-loop telling R to remove tags not befitting of HMMs

crwOut1 <- crawlWrap(obsData=a, timeStep="30 mins",theta=c(6.855, -0.007), fixPar=c(NA,NA))
plot(crwOut1)                     

View(crab.test[which(crab.test$ID=="31300"),])
a <- subset(crab.test, ID=="31300")
