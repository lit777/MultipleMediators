##---------------------------------------------------------------
## Load data
##---------------------------------------------------------------
load("SurfaceFrame1.RData")
load("SurfaceFrame2.RData")
load("Master.RData")

#------ Set Treatment (TRT), Outcome (OUT), Mediators (M) and Covariates (X)
Data <- Master
OUT <- Data$PM.2.5
TRT <- Data$SO2.SC
M <- cbind(Data$SO2_Annual, Data$NOx_Annual, Data$CO2_Annual)
X <- cbind(Data$S_n_CR, Data$NumNOxControls, Data$Heat_Input/100000, Data$Barometric_Pressure, Data$Temperature,  Data$PctCapacity, Data$sulfur_Content, Data$Phase2_Indicator, Data$Operating_Time/1000)

##---------------------------------------------------------------
## Surface Plots 1
##---------------------------------------------------------------


#----- Extract posterior samples for the SO2 plot
ll <- 1;  uu <- dim(ind1)[1]
d_1 <- apply(result, 2, mean)[ll:uu]
D <- matrix(ncol=30, nrow=30)
index <- 0
for(i in 1:30){
  for(j in 1:30){
    if(mat_1[i,j]==1){index<-1+index; D[i,j]<-d_1[index]}
  }
}

#----- Setting Palette of Colors
nrz <- nrow(D)
ncz <- ncol(D)
jet.colors <- colorRampPalette( c("yellow", "red") )
nbcol <- nrz*ncz
color <- jet.colors(nbcol)
zfacet <- D[-1, -1] + D[-1, -ncz] + D[-nrz, -1] + D[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)

CC<-(C.sample[c(1,4), which((C.sample[1,]) <= max(log(M[,1])) & (C.sample[4,]) <= max(log(M[,1])) & (C.sample[1,]) >= min(log(M[,1])) & (C.sample[4,]) >= min(log(M[,1])))])


# SO2.SD : 1.386
half.SD <- 1.386*0.5
# NOx.SD : 1.084
# half.SD <- 1.084*0.5
# CO2.SD : 0.945
# half.SD <- 0.945*0.5

#----- half SD band
y1 <- x <- y<-  seq(min(log(M[,1])),max(log(M[,1])),length=30)[1:30]
y2 <- x1 <- seq(min(log(M[,1])),max(log(M[,1]))-half.SD/sin(pi/4),length=30)
y3 <- x2 <- seq(min(log(M[,1]))+half.SD/sin(pi/4),max(log(M[,1])),length=30)

#----- Set the bottom
mat_bottom <- matrix(nrow=30,ncol=30)
for(i in 1:30){
  for(j in 1:30){
    if(mat_1[i,j]==1){mat_bottom[i,j] <- -3}
  }
}


#----- Plot 1
# degree <- 0

# png("plot1.png", height=1000, width=1000)
# persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", zlim=c(-3,1.5),cex.lab=1.8, cex.axis=1.5, col=color[facetcol]) -> res1
# persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", cex.lab=1.8, cex.axis=1.5,zlim=c(-3,1.5), col=color[facetcol], main="", cex.main=2)
# par(new=TRUE)
# persp(x,y,mat_bottom,theta= 30+degree, phi=12,zlim=c(-3,1.5), box=FALSE, border=NA, col=rgb(0, 0, 1, 0.1))
# points(trans3d(CC[1,], CC[2,],-3, pmat=res1), col = rgb(0, 0, 1, 0.5), pch = 16)
# lines(trans3d(x,y1,-3, res1), col = "red", lwd = 3)
# lines(trans3d(x1,y2+half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
# lines(trans3d(x2,y3-half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
# dev.off()


#----- Plot 1 (rotated by 105 degrees)
degree <- 105
y1 <- x11<-x[-c(1,2,3,4,5,6,7)]
y2 <- x1 <- seq(min(log(M[,1])),max(log(M[,1]))-half.SD/sin(pi/4),length=30)[-c(1,2,3,4)]
x2<-seq(min(log(M[,1]))+half.SD/sin(pi/4),max(log(M[,1])),length=30)
y3<- x22 <- x2[-c(1,2,3,4,5,6,7)]

png("plot1_rotate.png", height=1000, width=1000)
persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", zlim=c(-3,1.5),cex.lab=1.8, cex.axis=1.5,col=color[facetcol])->res1
persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", cex.lab=1.8, cex.axis=1.5,zlim=c(-3,1.5), col=color[facetcol], main="", cex.main=2)
par(new=TRUE)
persp(x,y,mat_bottom,theta= 30+degree, phi=12,zlim=c(-3,1.5), box=FALSE, border=NA, col=rgb(0, 0, 1, 0.1))
lines(trans3d(x11,y1,-3, res1), col = "red", lwd = 3)
lines(trans3d(x1,y2+half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
lines(trans3d(x22,y3-half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
points(trans3d(CC[1,], CC[2,],-3, pmat=res1), col = rgb(0, 0, 1, 0.5), pch = 16)
dev.off()




##---------------------------------------------------------------
## Surface Plots 2
##---------------------------------------------------------------


#----- Extract posterior samples for the NOx plot
ll <- dim(ind1)[1]+1;  uu <- dim(ind1)[1]+dim(ind2)[1]
d_1 <- apply(result, 2, mean)[ll:uu]
D <- matrix(ncol=30, nrow=30)
index <- 0
for(i in 1:30){
  for(j in 1:30){
    if(mat_2[i,j]==1){index<-1+index; D[i,j]<-d_1[index]}
  }
}

#----- Setting Palette of Colors
nrz <- nrow(D)
ncz <- ncol(D)
jet.colors <- colorRampPalette( c("yellow", "red") )
nbcol <- nrz*ncz
color <- jet.colors(nbcol)
zfacet <- D[-1, -1] + D[-1, -ncz] + D[-nrz, -1] + D[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)

CC<-(C.sample[c(2,5), which((C.sample[2,]) <= max(log(M[,2])) & (C.sample[5,]) <= max(log(M[,2])) & (C.sample[2,]) >= min(log(M[,2])) & (C.sample[5,]) >= min(log(M[,2])))])


# SO2.SD : 1.386
# half.SD <- 1.386*0.5
# NOx.SD : 1.084
half.SD <- 1.084*0.5
# CO2.SD : 0.945
# half.SD <- 0.945*0.5

#----- half SD band
y1 <- x <- y <- seq(min(log(M[,2])),max(log(M[,2])),length=30)[1:30]
y2 <- x1 <- seq(min(log(M[,2])),max(log(M[,2]))-half.SD/sin(pi/4),length=30)
y3 <- x2 <- seq(min(log(M[,2]))+half.SD/sin(pi/4),max(log(M[,2])),length=30)

#----- Set the bottom
mat_bottom <- matrix(nrow=30,ncol=30)
for(i in 1:30){
  for(j in 1:30){
    if(mat_2[i,j]==1){mat_bottom[i,j] <- -3}
  }
}


#----- Plot 2
# degree <- 0

# png("plot2.png", height=1000, width=1000)
# persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", zlim=c(-3,1.5),cex.lab=1.8, cex.axis=1.5, col=color[facetcol]) -> res1
# persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", cex.lab=1.8, cex.axis=1.5,zlim=c(-3,1.5), col=color[facetcol], main="", cex.main=2)
# par(new=TRUE)
# persp(x,y,mat_bottom,theta= 30+degree, phi=12,zlim=c(-3,1.5), box=FALSE, border=NA, col=rgb(0, 0, 1, 0.1))
# points(trans3d(CC[1,], CC[2,],-3, pmat=res1), col = rgb(0, 0, 1, 0.5), pch = 16)
# lines(trans3d(x,y1,-3, res1), col = "red", lwd = 3)
# lines(trans3d(x1,y2+half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
# lines(trans3d(x2,y3-half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
# dev.off()



#----- Plot 2 (rotated by 105 degrees)
degree <- 105
y1 <- x11<-x[-c(1,2,3,4,5,6,7,8,9)]
y2 <- x1 <- seq(min(log(M[,2])),max(log(M[,2]))-half.SD/sin(pi/4),length=30)[-c(1,2,3,4,5,6,7,8,9,10)]
x2<-seq(min(log(M[,2]))+half.SD/sin(pi/4),max(log(M[,2])),length=30)
y3<- x22 <- x2[-c(1,2,3)]

png("plot2_rotate.png", height=1000, width=1000)
persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", zlim=c(-3,1.5),cex.lab=1.8, cex.axis=1.5,col=color[facetcol])->res1
persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", cex.lab=1.8, cex.axis=1.5,zlim=c(-3,1.5), col=color[facetcol], main="", cex.main=2)
par(new=TRUE)
persp(x,y,mat_bottom,theta= 30+degree, phi=12,zlim=c(-3,1.5), box=FALSE, border=NA, col=rgb(0, 0, 1, 0.1))
lines(trans3d(x11,y1,-3, res1), col = "red", lwd = 3)
lines(trans3d(x1,y2+half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
lines(trans3d(x22,y3-half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
points(trans3d(CC[1,], CC[2,],-3, pmat=res1), col = rgb(0, 0, 1, 0.5), pch = 16)
dev.off()





##---------------------------------------------------------------
## Surface Plots 3
##---------------------------------------------------------------


#----- Extract posteriorl samples for the CO2 plot
ll <- dim(ind1)[1]+dim(ind2)[1]+1; uu <- dim(ind1)[1]+dim(ind2)[1]+dim(ind3)[1]
d_1 <- apply(result, 2, mean)[ll:uu]
D <- matrix(ncol=30, nrow=30)
index <- 0
for(i in 1:30){
  for(j in 1:30){
    if(mat_3[i,j]==1){index<-1+index; D[i,j]<-d_1[index]}
  }
}

#----- Setting Palette of Colors
nrz <- nrow(D)
ncz <- ncol(D)
jet.colors <- colorRampPalette( c("yellow", "red") )
nbcol <- nrz*ncz
color <- jet.colors(nbcol)
zfacet <- D[-1, -1] + D[-1, -ncz] + D[-nrz, -1] + D[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)

CC<-(C.sample[c(3,6), which((C.sample[3,]) <= max(log(M[,3])) & (C.sample[6,]) <= max(log(M[,3])) & (C.sample[3,]) >= min(log(M[,3])) & (C.sample[6,]) >= min(log(M[,3])))])


# SO2.SD : 1.386
# half.SD <- 1.386*0.5
# NOx.SD : 1.084
# half.SD <- 1.084*0.5
# CO2.SD : 0.945
half.SD <- 0.945*0.5

#----- half SD band
y1 <- x <- y <- seq(min(log(M[,3])),max(log(M[,3])),length=30)[1:30]
y2 <- x1 <- seq(min(log(M[,3])),max(log(M[,3]))-half.SD/sin(pi/4),length=30)
y3 <- x2 <- seq(min(log(M[,3]))+half.SD/sin(pi/4),max(log(M[,3])),length=30)

#----- Set the bottom
mat_bottom <- matrix(nrow=30,ncol=30)
for(i in 1:30){
  for(j in 1:30){
    if(mat_3[i,j]==1){mat_bottom[i,j] <- -3}
  }
}


#----- Plot 3
# degree <- 0

# png("plot3.png", height=1000, width=1000)
# persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", zlim=c(-3,1.5),cex.lab=1.8, cex.axis=1.5, col=color[facetcol]) -> res1
# persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", cex.lab=1.8, cex.axis=1.5,zlim=c(-3,1.5), col=color[facetcol], main="", cex.main=2)
# par(new=TRUE)
# persp(x,y,mat_bottom,theta= 30+degree, phi=12,zlim=c(-3,1.5), box=FALSE, border=NA, col=rgb(0, 0, 1, 0.1))
# points(trans3d(CC[1,], CC[2,],-3, pmat=res1), col = rgb(0, 0, 1, 0.5), pch = 16)
# lines(trans3d(x,y1,-3, res1), col = "red", lwd = 3)
# lines(trans3d(x1,y2+half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
# lines(trans3d(x2,y3-half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
# dev.off()



#----- Plot 3 (rotated by 105 degrees)
degree <- 105
y1 <- x11<-x[-c(1,2,3,4,5,6,7,8,9,10,11,12,13)]
y2 <- x1 <- seq(min(log(M[,3])),max(log(M[,3]))-half.SD/sin(pi/4),length=30)[-c(1,2,3,4,5,6,7,8,9,10,11)]
x2<-seq(min(log(M[,3]))+half.SD/sin(pi/4),max(log(M[,3])),length=30)
y3<- x22 <- x2[-c(1,2,3,4,5,6,7,8)]

png("plot3_rotate.png", height=1000, width=1000)
persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", zlim=c(-3,1.5),cex.lab=1.8, cex.axis=1.5,col=color[facetcol])->res1
persp(x,y,D, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", cex.lab=1.8, cex.axis=1.5,zlim=c(-3,1.5), col=color[facetcol], main="", cex.main=2)
par(new=TRUE)
persp(x,y,mat_bottom,theta= 30+degree, phi=12,zlim=c(-3,1.5), box=FALSE, border=NA, col=rgb(0, 0, 1, 0.1))
lines(trans3d(x11,y1,-3, res1), col = "red", lwd = 3)
lines(trans3d(x1,y2+half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
lines(trans3d(x22,y3-half.SD/sin(pi/4),-3, res1), col = "red", lty=2, lwd = 3)
points(trans3d(CC[1,], CC[2,],-3, pmat=res1), col = rgb(0, 0, 1, 0.5), pch = 16)
dev.off()










