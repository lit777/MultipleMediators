##----- Causal Indirect and Direct Effects
MME <- list()
for(ii in 1:1500){
  i <- ii+500
    MME[[ii]] <-list(TE=mean(unlist(Y_1111[[i]]))-mean(unlist(Y_0000[[i]])), 
                    JNIE=mean(unlist(Y_1111[[i]]))-mean(unlist(Y_1000[[i]])), 
                    NDE=mean(unlist(Y_1000[[i]]))-mean(unlist(Y_0000[[i]])), 
                    NIE1=mean(unlist(Y_1111[[i]]))-mean(unlist(Y_1011[[i]])), 
                    NIE2=mean(unlist(Y_1111[[i]]))-mean(unlist(Y_1101[[i]])), 
                    NIE3=mean(unlist(Y_1111[[i]]))-mean(unlist(Y_1110[[i]])), 
                    NIE12=mean(unlist(Y_1111[[i]]))-mean(unlist(Y_1001[[i]])), 
                    NIE13=mean(unlist(Y_1111[[i]]))-mean(unlist(Y_1010[[i]])), 
                    NIE23=mean(unlist(Y_1111[[i]]))-mean(unlist(Y_1100[[i]])))
}


result <- sapply(1:1500, function(y) unlist(MME[[y]]))
rownames(result)[7:9] <- c("JNIE12","JNIE13","JNIE23")

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


##----- Reshape the data
#install.packages("reshape")
library(reshape)
result <- melt(result)
result$X1 <- factor(result$X1, levels=levels(result$X1)[c(9,1,5,6,7,8,2,3,4)])

pdf("dist_uni.pdf", width=10, height=5)
boxplot(value ~ X1, data = result, lwd = 1,outline=FALSE, ylim=c(-3, 1),whisklty = 0, staplelty = 0)
abline(h=0, col="gray")
stripchart(value ~ X1, vertical = TRUE, data = result, 
           method = "jitter", add = TRUE, pch = 20, col = c(rep(makeTransparent("red",alpha=0.05),3),rep(makeTransparent("blue",alpha=0.05),3),rep(makeTransparent("green",alpha=0.05),3) )       )
dev.off()
