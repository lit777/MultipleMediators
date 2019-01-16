ede1 <- ede2 <- ede3 <- ede12 <- ede13 <- ede23 <- ede123 <- list()
eae1 <- eae2 <- eae3 <- eae12 <- eae13 <- eae23 <- eae123 <- list()
eeae1 <- eeae2 <- eeae3 <- eeae12 <- eeae13 <- eeae23 <- eeae123 <- list()
s1 <- s2 <- s3 <- s4 <- s5 <- s6 <- s7 <- s8 <- s9 <- list()

EDE1 <- EDE2 <- EDE3 <- EDE12 <- EDE23 <- EDE13 <- EDE123 <- NULL
EAE1 <- EAE2 <- EAE3 <- EAE12 <- EAE23 <- EAE13 <- EAE123 <- NULL
EEAE1 <- EEAE2 <- EEAE3 <- EEAE12 <- EEAE23 <- EEAE13 <- EEAE123 <- NULL
E1 <- E2 <- E3 <- E4 <- E5 <- E6 <- E7 <- E8 <- E9 <- NULL

c1 <- 0.24*0.25
c2 <- 0.42*0.25
c3 <- 0.02*0.25

for(cc in 1:1500){

c<-cc+500

ede1[[cc]] <- which(abs(ME[[c]]$m1[,1]-ME[[c]]$m0[,1]) < c1)
ede2[[cc]] <- which(abs(ME[[c]]$m1[,2]-ME[[c]]$m0[,2]) < c2)
ede3[[cc]] <- which(abs(ME[[c]]$m1[,3]-ME[[c]]$m0[,3]) < c3)

ede12[[cc]] <- Reduce(function(x,y) intersect(x,y), list(ede1[[cc]],ede2[[cc]]))
ede23[[cc]] <- Reduce(function(x,y) intersect(x,y), list(ede2[[cc]],ede3[[cc]]))
ede13[[cc]] <- Reduce(function(x,y) intersect(x,y), list(ede1[[cc]],ede3[[cc]]))

ede123[[cc]] <- Reduce(function(x,y) intersect(x,y), list(ede1[[cc]],ede2[[cc]],ede3[[cc]]))

eae1[[cc]] <- which((ME[[c]]$m1[,1]-ME[[c]]$m0[,1]) > c1)
eae2[[cc]] <- which((ME[[c]]$m1[,2]-ME[[c]]$m0[,2]) > c2)
eae3[[cc]] <- which((ME[[c]]$m1[,3]-ME[[c]]$m0[,3]) > c3)

eae12[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eae1[[cc]],eae2[[cc]]))
eae23[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eae2[[cc]],eae3[[cc]]))
eae13[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eae1[[cc]],eae3[[cc]]))

eae123[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eae1[[cc]],eae2[[cc]],eae3[[cc]]))

eeae1[[cc]] <- which((ME[[c]]$m1[,1]-ME[[c]]$m0[,1]) < c1)
eeae2[[cc]] <- which((ME[[c]]$m1[,2]-ME[[c]]$m0[,2]) < c2)
eeae3[[cc]] <- which((ME[[c]]$m1[,3]-ME[[c]]$m0[,3]) < c3)

eeae12[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eeae1[[cc]],eeae2[[cc]]))
eeae23[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eeae2[[cc]],eeae3[[cc]]))
eeae13[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eeae1[[cc]],eeae3[[cc]]))

eeae123[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eeae1[[cc]],eeae2[[cc]],eeae3[[cc]]))

EDE1[cc] <- mean(Y_1111[[c]][ede1[[cc]]]-Y_0000[[c]][ede1[[cc]]])
EAE1[cc] <- mean(Y_1111[[c]][eae1[[cc]]]-Y_0000[[c]][eae1[[cc]]])
EEAE1[cc] <- mean(Y_1111[[c]][eeae1[[cc]]]-Y_0000[[c]][eeae1[[cc]]])

EDE2[cc] <- mean(Y_1111[[c]][ede2[[cc]]]-Y_0000[[c]][ede2[[cc]]])
EAE2[cc] <- mean(Y_1111[[c]][eae2[[cc]]]-Y_0000[[c]][eae2[[cc]]])
EEAE2[cc] <- mean(Y_1111[[c]][eeae2[[cc]]]-Y_0000[[c]][eeae2[[cc]]])

EDE3[cc] <- mean(Y_1111[[c]][ede3[[cc]]]-Y_0000[[c]][ede3[[cc]]])
EAE3[cc] <- mean(Y_1111[[c]][eae3[[cc]]]-Y_0000[[c]][eae3[[cc]]])
EEAE3[cc] <- mean(Y_1111[[c]][eeae3[[cc]]]-Y_0000[[c]][eeae3[[cc]]])

EDE13[cc] <- mean(Y_1111[[c]][ede13[[cc]]]-Y_0000[[c]][ede13[[cc]]])
EAE13[cc] <- mean(Y_1111[[c]][eae13[[cc]]]-Y_0000[[c]][eae13[[cc]]])
EEAE13[cc] <- mean(Y_1111[[c]][eeae13[[cc]]]-Y_0000[[c]][eeae13[[cc]]])

EDE23[cc] <- mean(Y_1111[[c]][ede23[[cc]]]-Y_0000[[c]][ede23[[cc]]])
EAE23[cc] <- mean(Y_1111[[c]][eae23[[cc]]]-Y_0000[[c]][eae23[[cc]]])
EEAE23[cc] <- mean(Y_1111[[c]][eeae23[[cc]]]-Y_0000[[c]][eeae23[[cc]]])

EDE12[cc] <- mean(Y_1111[[c]][ede12[[cc]]]-Y_0000[[c]][ede12[[cc]]])
EAE12[cc] <- mean(Y_1111[[c]][eae12[[cc]]]-Y_0000[[c]][eae12[[cc]]])
EEAE12[cc] <- mean(Y_1111[[c]][eeae12[[cc]]]-Y_0000[[c]][eeae12[[cc]]])

EDE123[cc] <- mean(Y_1111[[c]][ede123[[cc]]]-Y_0000[[c]][ede123[[cc]]])
EAE123[cc] <- mean(Y_1111[[c]][eae123[[cc]]]-Y_0000[[c]][eae123[[cc]]])
EEAE123[cc] <- mean(Y_1111[[c]][eeae123[[cc]]]-Y_0000[[c]][eeae123[[cc]]])


s1[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eeae1[[cc]],eeae3[[cc]]))
s2[[cc]] <- Reduce(function(x,y) intersect(x,y), list(ede1[[cc]],eeae3[[cc]]))
s3[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eae1[[cc]],eeae3[[cc]]))

s4[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eeae1[[cc]],ede3[[cc]]))
s5[[cc]] <- Reduce(function(x,y) intersect(x,y), list(ede1[[cc]],ede3[[cc]]))
s6[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eae1[[cc]],ede3[[cc]]))

s7[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eeae1[[cc]],eae3[[cc]]))
s8[[cc]] <- Reduce(function(x,y) intersect(x,y), list(ede1[[cc]],eae3[[cc]]))
s9[[cc]] <- Reduce(function(x,y) intersect(x,y), list(eae1[[cc]],eae3[[cc]]))

E1[cc] <- mean(Y_1111[[c]][s1[[cc]]]-Y_0000[[c]][s1[[cc]]])
E2[cc] <- mean(Y_1111[[c]][s2[[cc]]]-Y_0000[[c]][s2[[cc]]])
E3[cc] <- mean(Y_1111[[c]][s3[[cc]]]-Y_0000[[c]][s3[[cc]]])
E4[cc] <- mean(Y_1111[[c]][s4[[cc]]]-Y_0000[[c]][s4[[cc]]])
E5[cc] <- mean(Y_1111[[c]][s5[[cc]]]-Y_0000[[c]][s5[[cc]]])
E6[cc] <- mean(Y_1111[[c]][s6[[cc]]]-Y_0000[[c]][s6[[cc]]])
E7[cc] <- mean(Y_1111[[c]][s7[[cc]]]-Y_0000[[c]][s7[[cc]]])
E8[cc] <- mean(Y_1111[[c]][s8[[cc]]]-Y_0000[[c]][s8[[cc]]])
E9[cc] <- mean(Y_1111[[c]][s9[[cc]]]-Y_0000[[c]][s9[[cc]]])
}

strata.len <- matrix(nrow=1500, ncol=9)
for(i in 1:1500){
  strata.len[i,1] <- length(s1[[i]])
  strata.len[i,2] <- length(s2[[i]])
  strata.len[i,3] <- length(s3[[i]])
  strata.len[i,4] <- length(s4[[i]])
  strata.len[i,5] <- length(s5[[i]])
  strata.len[i,6] <- length(s6[[i]])
  strata.len[i,7] <- length(s7[[i]])
  strata.len[i,8] <- length(s8[[i]])
  strata.len[i,9] <- length(s9[[i]])
}
result.len <- colMeans(strata.len, na.rm=TRUE)
result <- c(mean(E1),mean(E2),mean(E3,na.rm=TRUE),mean(E4),mean(E5),mean(E6,na.rm=TRUE),mean(E7),mean(E8),mean(E9,na.rm=TRUE))
result.sd <- c(sd(E1),sd(E2),sd(E3,na.rm=TRUE),sd(E4),sd(E5),sd(E6,na.rm=TRUE),sd(E7),sd(E8),sd(E9,na.rm=TRUE))


##-----Posterior mean estimates of principal causal effects

##----- Make transparent color codes
makeTransparent = function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}


pdf("PSplot_substrata_CO2.pdf", width=10, height=5)
plot(result, cex=result.len/15, ylim=c(-2, 1.5), xlim=c(0.5,9.5), col=rep(c(makeTransparent("black",alpha=0.8),makeTransparent("black",alpha=0.6),makeTransparent("black",alpha=0.4)),3),pch=19, ylab="E[Y(1)-Y(0)]",xaxt="n", xlab="")

text(1:9, result,  round(c(result.len[1:3]/sum(result.len[1:3]),result.len[4:6]/sum(result.len[4:6]),result.len[7:9]/sum(result.len[7:9])),2),
     cex=0.65, pos=3,col=rep(c("white","black","black"),3))

text(1:9, result, c("(0.57)","(0.62)","(2.23)","(0.49)","(0.68)","(2.57)","(0.73)","(0.99)","(3.78)"),
     cex=0.65, pos=1,col=rep(c("white","black","black"),3))


xtick<-c(2,5,8)
axis(side=1, at=xtick, labels = FALSE)
text(x=xtick,  par("usr")[3],
     labels = c("CO2 decrease", "CO2 no change", "CO2 increase"), pos = 1, xpd = TRUE)

legend("topleft",col=c(makeTransparent("black",alpha=0.8),makeTransparent("black",alpha=0.6),makeTransparent("black",alpha=0.4)),legend=c("SO2 decrease","SO2 no change","SO2 increase"), pch=19)
dev.off()




