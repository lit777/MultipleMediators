#----- Load the main dataset
load("Master.RData")

#----- Required libraries
library(maptools)
library(maps)
library(data.table)

#----- Link monitoring locations to the main dataset
Master_Data <- data.table(Master)
Monitor_Data <- data.table(L_all)
Master_Data <-Master_Data[,list(Facility.ID..ORISPL.,SO2.SC)]
Monitor_Data <- Monitor_Data[,list(Facility.ID..ORISPL.,Longitude.Monitor,Latitude.Monitor,Facility.Longitude=Facility.Longitude.x,Facility.Latitude=Facility.Latitude.x)]
All_Data <- merge(Master_Data, Monitor_Data, by="Facility.ID..ORISPL.")

#----- Power Plants under Z=1 vs Z=0
ppid1 <- with(All_Data, Facility.ID..ORISPL.[which(SO2.SC==1)])
ppid0 <- with(All_Data, Facility.ID..ORISPL.[which(SO2.SC==0)])

#----- Monitors under Z=1 vs Z=0
L <- subset(L, Include.in.Average==1)
monitor1 <- L[Facility.ID..ORISPL. %in% ppid1]
monitor0 <- L[Facility.ID..ORISPL. %in% ppid0]


#----- Plots
pdf("PP_PM.pdf")
par(mfrow=c(2,1))
map('state', mar=c(1,1,0.5,0.1))
with(All_Data, points(Facility.Longitude[which(SO2.SC==1)], Facility.Latitude[which(SO2.SC==1)],col=2,pch=16))
with(monitor1, points(Longitude.Monitor, Latitude.Monitor,col=1,pch=16, cex=0.5))
legend("bottomleft",c("Power Plants (treated)", "Monitors (treated)"), horiz = TRUE, pch = c(16,16), col=c(2,1), cex=0.6)
map('state', mar=c(1,1,0.5,0.1))
with(All_Data, points(Facility.Longitude[which(SO2.SC==0)], Facility.Latitude[which(SO2.SC==0)],col=3,pch=16))
with(monitor0, points(Longitude.Monitor, Latitude.Monitor,col=4,pch=16, cex=0.5))
legend("bottomleft",c("Power Plants (untreated)", "Monitors (untreated)"), horiz = TRUE, pch = c(16,16), col=c(3,4), cex=0.6)
dev.off()


