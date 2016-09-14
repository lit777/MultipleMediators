#------ Load Master data
load("Master.RData")

Master <- data.frame(Master) # Convert to data.frame

mean1 <- round(apply(Master[which(Master$SO2.SC >= 0.01), 3:dim(Master)[2]], 2, median), 2)
sd1 <- round(apply(Master[which(Master$SO2.SC >= 0.01), 3:dim(Master)[2]], 2, function(x) quantile(x, c(0.25, 0.75))), 2)

mean0 <- round(apply(Master[which(Master$SO2.SC <= 0.01), 3:dim(Master)[2]], 2, median), 2)
sd0 <- round(apply(Master[which(Master$SO2.SC <= 0.01), 3:dim(Master)[2]], 2, function(x) quantile(x, c(0.25, 0.75))), 2)

new.data <- data.frame(cbind(mean1, t(sd1), mean0, t(sd0)))
names(new.data) <- c("Median (T : 63)","25%","75%","Median (C : 195)","25%","75%") # header

#------ latex table format
# install.packages("xtable")
xtable::xtable(new.data) 
