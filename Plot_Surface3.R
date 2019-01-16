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
M <- (cbind(Data$SO2_Annual/10000, Data$NOx_Annual/1000, Data$CO2_Annual/10000000))
XX <- (cbind(Data$NumNOxControls, log(Data$Heat_Input), Data$Barometric_Pressure/100, Data$Temperature,  Data$PctCapacity/100, Data$Sulfur_Content, log(Data$Operating_Time), log(Data$Heat_Rate)))

##---------------------------------------------------------------
## Surface Plots 1
##---------------------------------------------------------------


#----- Extract posterior samples for the SO2 plot
ll <- 1;  uu <- dim(ind1)[1]
d_1 <- apply(result, 2, function(x) mean(x, na.rm=TRUE))[ll:uu]
D <- matrix(ncol=30, nrow=30)
index <- 0
for(i in 1:30){
  for(j in 1:30){
    if(mat_1[i,j]==1){index<-1+index; D[i,j]<-d_1[index]}
  }
}

DD <- matrix(ncol=30, nrow=30)
for(i in 2:29){
    for(j in 2:29){
        im <- max(i-1,1)
        jm <- max(j-1,1)
        ip <- min(i+1,30)
        jp <- min(j+1,30)

DD[i,j] <- mean(c(D[i,j],D[im,j],D[ip,j],D[i,jm],D[im,jm],D[ip,jm],D[i,jp],D[im,jp],D[ip,jp]),na.rm=TRUE)
    }
}

MIN <- apply(DD, 2, function(x) min(which(!is.na(x))))
MIN <- ifelse(MIN>30, 1, MIN)
MAX <- apply(DD, 2, function(x) max(which(!is.na(x))))
MAX <- ifelse(MIN>30, 1, MAX)
for(i in 1:30){
    DD[MIN[i],i] <- NaN
    DD[MAX[i],i] <- NaN

}
MIN <- apply(DD, 2, function(x) min(which(!is.na(x))))
MIN <- ifelse(MIN>30, 1, MIN)
MAX <- apply(DD, 2, function(x) max(which(!is.na(x))))
MAX <- ifelse(MIN>30, 1, MAX)
for(i in 1:30){
  DD[MIN[i],i] <- NaN
  DD[MAX[i],i] <- NaN
  
}


#----- Setting Palette of Colors
nrz <- nrow(DD)
ncz <- ncol(DD)
jet.colors <- colorRampPalette( c("yellow", "red") )
nbcol <- nrz*ncz
color <- jet.colors(nbcol)
zfacet <- DD[-1, -1] + DD[-1, -ncz] + DD[-nrz, -1] + DD[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)

C.sample <- ME[[1999]]
CC <- log(rbind(C.sample$m0[,1][which(C.sample$m0[,1] <= max((M[,1])) & (C.sample$m1[,1]) <= max((M[,1])) & (C.sample$m0[,1]) >= min((M[,1])) & (C.sample$m1[,1]) >= min((M[,1])))]*10000,
            C.sample$m1[,1][which(C.sample$m0[,1] <= max((M[,1])) & (C.sample$m1[,1]) <= max((M[,1])) & (C.sample$m0[,1]) >= min((M[,1])) & (C.sample$m1[,1]) >= min((M[,1])))]*10000))


half.SD <- 1.19*0.25 # (log scale)

#----- half SD band
y1 <- x <- y<-  seq(min(log(M[,1]*10000)),max(log(M[,1]*10000)),length=30)[1:30]
y2 <- x1 <- seq(min(log(M[,1]*10000)),max(log(M[,1]*10000))-half.SD/sin(pi/4),length=30)
y3 <- x2 <- seq(min(log(M[,1]*10000))+half.SD/sin(pi/4),max(log(M[,1]*10000)),length=30)

#----- Set the bottom
mat_bottom <- matrix(nrow=30,ncol=30)
for(i in 1:30){
  for(j in 1:30){
    if(mat_1[i,j]==1){mat_bottom[i,j] <- -4}
  }
}




#----- Plot 1 (rotated by 105 degrees)
degree <- 105
y1 <- x11<-x
y2 <- x1 <- seq(min(log(M[,1]*10000)),max(log(M[,1]*10000))-half.SD/sin(pi/4),length=30)
x2<-seq(min(log(M[,1]*10000))+half.SD/sin(pi/4),max(log(M[,1]*10000)),length=30)
y3<- x22 <- x2[-c(1)]

png("plot1_rotate.png", height=1000, width=1000)
persp(x,y,DD, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", zlim=c(-4,4.5),cex.lab=1.8, cex.axis=1.5,col=color[facetcol])->res1
persp(x,y,DD, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", cex.lab=1.8, cex.axis=1.5,zlim=c(-4,4.5), col=color[facetcol], main="", cex.main=2)
par(new=TRUE)
persp(x,y,mat_bottom,theta= 30+degree, phi=12,zlim=c(-4,4), box=FALSE, border=NA, col=rgb(0, 0, 1, 0.1))
lines(trans3d(x11,y1,-4, res1), col = "red", lwd = 3)
lines(trans3d(x1,y2+half.SD/sin(pi/4),-4, res1), col = "red", lty=2, lwd = 3)
lines(trans3d(x22,y3-half.SD/sin(pi/4),-4, res1), col = "red", lty=2, lwd = 3)
points(trans3d(CC[1,], CC[2,],-4, pmat=res1), col = rgb(0, 0, 1, 0.5), pch = 16)
dev.off()




##---------------------------------------------------------------
## Surface Plots 2
##---------------------------------------------------------------


#----- Extract posterior samples for the NOx plot
ll <- dim(ind1)[1]+1;  uu <- dim(ind1)[1]+dim(ind2)[1]
d_1 <- apply(result, 2, function(x) mean(x, na.rm=TRUE))[ll:uu]
D <- matrix(ncol=30, nrow=30)
index <- 0
for(i in 1:30){
  for(j in 1:30){
    if(mat_2[i,j]==1){index<-1+index; D[i,j]<-d_1[index]}
  }
}


DD <- matrix(ncol=30, nrow=30)
for(i in 1:30){
    for(j in 1:30){
        im <- max(i-1,1)
        jm <- max(j-1,1)
        ip <- min(i+1,30)
        jp <- min(j+1,30)
        
        DD[i,j] <- mean(c(D[i,j],D[im,j],D[ip,j],D[i,jm],D[im,jm],D[ip,jm],D[i,jp],D[im,jp],D[ip,jp]),na.rm=TRUE)
    }
}

MIN <- apply(DD, 2, function(x) min(which(!is.na(x))))
MAX <- apply(DD, 2, function(x) max(which(!is.na(x))))
for(i in 1:30){
    DD[MIN[i],i] <- NaN
    DD[MAX[i],i] <- NaN
    
}
MIN <- apply(DD, 2, function(x) min(which(!is.na(x))))
MAX <- apply(DD, 2, function(x) max(which(!is.na(x))))
for(i in 1:30){
  DD[MIN[i],i] <- NaN
  DD[MAX[i],i] <- NaN
  
}
MIN <- apply(DD, 2, function(x) min(which(!is.na(x))))
MAX <- apply(DD, 2, function(x) max(which(!is.na(x))))
for(i in 1:30){
  DD[MIN[i],i] <- NaN
  DD[MAX[i],i] <- NaN
  
}

#----- Setting Palette of Colors
nrz <- nrow(DD)
ncz <- ncol(DD)
jet.colors <- colorRampPalette( c("yellow", "red") )
nbcol <- nrz*ncz
color <- jet.colors(nbcol)
zfacet <- DD[-1, -1] + DD[-1, -ncz] + DD[-nrz, -1] + DD[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)

C.sample <- ME[[1999]]
CC <- log(rbind(C.sample$m0[,2][which(C.sample$m0[,2] <= max((M[,2])) & (C.sample$m1[,2]) <= max((M[,2])) & (C.sample$m0[,2]) >= min((M[,2])) & (C.sample$m1[,2]) >= min((M[,2])))]*1000,
C.sample$m1[,2][which(C.sample$m0[,2] <= max((M[,2])) & (C.sample$m1[,2]) <= max((M[,2])) & (C.sample$m0[,2]) >= min((M[,2])) & (C.sample$m1[,2]) >= min((M[,2])))]*1000))

half.SD <- 0.76*0.25

#----- half SD band
y1 <- x <- y <- seq(min(log(M[,2]*1000)),max(log(M[,2]*1000)),length=30)[1:30]
y2 <- x1 <- seq(min(log(M[,2]*1000)),max(log(M[,2]*1000))-half.SD/sin(pi/4),length=30)
y3 <- x2 <- seq(min(log(M[,2]*1000))+half.SD/sin(pi/4),max(log(M[,2]*1000)),length=30)

#----- Set the bottom
mat_bottom <- matrix(nrow=30,ncol=30)
for(i in 1:30){
  for(j in 1:30){
    if(mat_2[i,j]==1){mat_bottom[i,j] <- -4.5}
  }
}



#----- Plot 2 (rotated by 105 degrees)
degree <- 105
y1 <- x11<-x[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]
y2 <- x1 <- seq(min(log(M[,2]*1000)),max(log(M[,2]*1000))-half.SD/sin(pi/4),length=30)[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]
x2<-seq(min(log(M[,2]*1000))+half.SD/sin(pi/4),max(log(M[,2]*1000)),length=30)
y3<- x22 <- x2[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]

png("plot2_rotate.png", height=1000, width=1000)
persp(x,y,DD, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", zlim=c(-4.5,5.5),cex.lab=1.8, cex.axis=1.5,col=color[facetcol])->res1
persp(x,y,DD, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", cex.lab=1.8, cex.axis=1.5,zlim=c(-4.5,5.5), col=color[facetcol], main="", cex.main=2)
par(new=TRUE)
persp(x,y,mat_bottom,theta= 30+degree, phi=12,zlim=c(-4.5,5.5), box=FALSE, border=NA, col=rgb(0, 0, 1, 0.1))
lines(trans3d(x11,y1,-4.5, res1), col = "red", lwd = 3)
lines(trans3d(x1,y2+half.SD/sin(pi/4),-4.5, res1), col = "red", lty=2, lwd = 3)
lines(trans3d(x22,y3-half.SD/sin(pi/4),-4.5, res1), col = "red", lty=2, lwd = 3)
points(trans3d(CC[1,], CC[2,],-4.5, pmat=res1), col = rgb(0, 0, 1, 0.5), pch = 16)
dev.off()





##---------------------------------------------------------------
## Surface Plots 3
##---------------------------------------------------------------


#----- Extract posteriorl samples for the CO2 plot
ll <- dim(ind1)[1]+dim(ind2)[1]+1; uu <- dim(ind1)[1]+dim(ind2)[1]+dim(ind3)[1]
d_1 <- apply(result, 2, function(x) mean(x, na.rm=TRUE))[ll:uu]
D <- matrix(ncol=30, nrow=30)
index <- 0
for(i in 1:30){
  for(j in 1:30){
    if(mat_3[i,j]==1){index<-1+index; D[i,j]<-d_1[index]}
  }
}



DD <- matrix(ncol=30, nrow=30)
for(i in 1:30){
    for(j in 1:30){
        im <- max(i-1,1)
        jm <- max(j-1,1)
        ip <- min(i+1,30)
        jp <- min(j+1,30)
        
        DD[i,j] <- mean(c(D[i,j],D[im,j],D[ip,j],D[i,jm],D[im,jm],D[ip,jm],D[i,jp],D[im,jp],D[ip,jp]),na.rm=TRUE)
    }
}

MIN <- apply(DD, 2, function(x) min(which(!is.na(x))))
MAX <- apply(DD, 2, function(x) max(which(!is.na(x))))
for(i in 1:30){
    DD[MIN[i],i] <- NaN
    DD[MAX[i],i] <- NaN
    
}
MIN <- apply(DD, 2, function(x) min(which(!is.na(x))))
MAX <- apply(DD, 2, function(x) max(which(!is.na(x))))
for(i in 1:30){
  DD[MIN[i],i] <- NaN
  DD[MAX[i],i] <- NaN
  
}

#----- Setting Palette of Colors
nrz <- nrow(DD)
ncz <- ncol(DD)
jet.colors <- colorRampPalette( c("yellow", "red") )
nbcol <- nrz*ncz
color <- jet.colors(nbcol)
zfacet <- DD[-1, -1] + DD[-1, -ncz] + DD[-nrz, -1] + DD[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)

C.sample <- ME[[1999]]
CC <- log(rbind(C.sample$m0[,3][which(C.sample$m0[,3] <= max((M[,3])) & (C.sample$m1[,3]) <= max((M[,3])) & (C.sample$m0[,3]) >= min((M[,3])) & (C.sample$m1[,3]) >= min((M[,3])))]*10000000,
C.sample$m1[,3][which(C.sample$m0[,3] <= max((M[,3])) & (C.sample$m1[,3]) <= max((M[,3])) & (C.sample$m0[,3]) >= min((M[,3])) & (C.sample$m1[,3]) >= min((M[,3])))]*10000000))

half.SD <- 0.68*0.25

#----- half SD band
y1 <- x <- y <- seq(min(log(M[,3]*10000000)),max(log(M[,3]*10000000)),length=30)[1:30]
y2 <- x1 <- seq(min(log(M[,3]*10000000)),max(log(M[,3]*10000000))-half.SD/sin(pi/4),length=30)
y3 <- x2 <- seq(min(log(M[,3]*10000000))+half.SD/sin(pi/4),max(log(M[,3]*10000000)),length=30)

#----- Set the bottom
mat_bottom <- matrix(nrow=30,ncol=30)
for(i in 1:30){
  for(j in 1:30){
    if(mat_3[i,j]==1){mat_bottom[i,j] <- -4.5}
  }
}



#----- Plot 3 (rotated by 105 degrees)
degree <- 105
y1 <- x11<-x[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)]
y2 <- x1 <- seq(min(log(M[,3]*10000000)),max(log(M[,3]*10000000))-half.SD/sin(pi/4),length=30)[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]
x2<-seq(min(log(M[,3]*10000000))+half.SD/sin(pi/4),max(log(M[,3]*10000000)),length=30)
y3<- x22 <- x2[-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)]

png("plot3_rotate.png", height=1000, width=1000)
persp(x,y,DD, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", zlim=c(-4.5,5.5),cex.lab=1.8, cex.axis=1.5,col=color[facetcol])->res1
persp(x,y,DD, theta= 30+degree, phi=12, ticktype="detailed", xlab="M(0)", ylab="M(1)", zlab="Y(1)-Y(0)", cex.lab=1.8, cex.axis=1.5,zlim=c(-4.5,5.5), col=color[facetcol], main="", cex.main=2)
par(new=TRUE)
persp(x,y,mat_bottom,theta= 30+degree, phi=12,zlim=c(-4.5,5.5), box=FALSE, border=NA, col=rgb(0, 0, 1, 0.1))
lines(trans3d(x11,y1,-4.5, res1), col = "red", lwd = 3)
lines(trans3d(x1,y2+half.SD/sin(pi/4),-4.5, res1), col = "red", lty=2, lwd = 3)
lines(trans3d(x22,y3-half.SD/sin(pi/4),-4.5, res1), col = "red", lty=2, lwd = 3)
points(trans3d(CC[1,], CC[2,],-4.5, pmat=res1), col = rgb(0, 0, 1, 0.5), pch = 16)
dev.off()










