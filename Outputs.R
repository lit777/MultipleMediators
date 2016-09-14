##----- Causal Effects on the Mediators
load("result_Mediators.RData")

index1 <- seq(1,6000,by=3); index2 <- seq(2,6000,by=3); index3 <- seq(3,6000,by=3)
# Means : SO2, NOx and CO2
(Mean.SO2 <- mean(apply(result, 1, mean)[index1]))
(Mean.NOx <- mean(apply(result, 1, mean)[index2]))
(Mean.CO2 <- mean(apply(result, 1, mean)[index3]))
# 95% C.I.s : SO2, NOx and CO2
(quantile(apply(result, 1, mean)[index1], c(0.025, 0.975)))
(quantile(apply(result, 1, mean)[index2], c(0.025, 0.975)))
(quantile(apply(result, 1, mean)[index3], c(0.025, 0.975)))





##----- Causal Indirect and Direct Effects
load("result_Effects.RData")

# Means : TE, JNIE123, NDE, NIE1, NIE2, NIE3, JNIE12, JNIE23, JNIE13
(Means <- apply(result, 2, mean))
# 95% C.I.s :
apply(result, 2, function(x) quantile(x, c(0.025, 0.975)))





##----- Principal Causal Effects
load("result_PCEffects.RData")

# Means : SO2, NOx, CO2, SO2&NOx, SO2&CO2, NOx&CO2, SO2&NOx&CO2
(Means.EAE1 <- apply(result[,1:7], 2, mean))
(Means.EDE <- apply(result[,8:14], 2, mean))
(Means.EAE2 <- apply(result[,15:21], 2, mean))
# 95% C.I.s :
(apply(result[,1:7], 2, function(x) quantile(x, c(0.025, 0.975))))
(apply(result[,8:14], 2, function(x) quantile(x, c(0.025, 0.975))))
(apply(result[,15:21], 2, function(x) quantile(x, c(0.025, 0.975))))

