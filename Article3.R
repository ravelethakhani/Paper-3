library(MASS)
library(e1071)
library(quantreg)
library(qgam)
library(devtools)
library(gefcom2017)
library(caret)
library(relaimpo)
library(glmnet)
library(evgam)
library(extremefit)
library(scoringRules)
library(fitdistrplus)
library(DataExplorer)
library(ggplot2)
library(gridExtra)
library(data.table)
library(rmarkdown)
library(splines)
library(earth) # fit MARS models
library(evmix)
library(copula)

######################Hourly GHI UNVhr#############################################
attach(UNVimp)
head(UNVimp)
win.graph()

####Summary statistics for hours
summary(GHI)
summary(Temp)
summary(RH)
####Standard deviation
sd(GHI)
sd(Temp)
sd(RH)
####skewness
skewness(GHI)
skewness(Temp)
skewness(RH)
####kurtosis
kurtosis(GHI)
kurtosis(Temp)
kurtosis(RH)
########################################################
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

Temp = na.zero(Temp)
GHI = na.zero(GHI)
RH = na.zero(RH)
GHI <- round(GHI,4)

data1 <- data.frame(cbind(GHI,Temp))
head(data1)

#####ts and density plots for GHI
d1 <- ggplot(UNVimp)+aes(x=t,y=GHI)+
  geom_point(position="jitter")+
  labs(x = "Observation number", title="(a) Plot of GHI",
       y ="Hourly GHI(W/sqre m)")
d2 <- ggplot(data1)+aes(x=GHI)+
  geom_histogram(bins=sqrt(nrow(data1)))+
  labs(x = "Hourly GHI(W/sqre m)", 
       title="(b) Histogram of GHI",y ="Frequency")
d3 <- ggplot(data1)+aes(x=GHI)+geom_boxplot()+
  labs(title="(c) Box plot of GHI",
       x ="Hourly GHI(W/sqre m)") 
#d4 <- ggplot(Hourlydata)+aes(x=t,y=GHI)+
#        geom_point(position="jitter")
plot(d1)
p1 <- grid.arrange(d2,d3,nrow=1,ncol=2)

#----------------------------------
# USING GGPLOT for nonlinear trend
#----------------------------------
library(splines)
qplot(t, GHI, data = UNVimp, geom = c("point", "smooth"), method = "lm",
      formula = y ~ ns(x,159), ylab="Hourly GHI (W/sqre m)",
      xlab="Observation number")
#------------------------------------------

#------------------------------------------------
# MARS- Multivariate Adaptive Regression Splines
#-------------------------------------------------
#-------------------------------
# GHI ~ Temp
#-------------------------------
na.zero <- function (x) {
  x[is.na(x)] <- 0
  return(x)
}

Temp = na.zero(Temp)
GHI = na.zero(GHI)
RH = na.zero(RH)
GHI <- round(GHI,4)

data1 <- data.frame(cbind(GHI,Temp))
head(data1)
?earth

set.seed(4321) #minspan=-2
marsModel1 <- earth(GHI ~ .,minspan=-4, data=data1) # build model
plotmo(marsModel1,col="blue",xlab=" Temperature (deg C)", main="MARS model",
       ylab="Global horizontal irradiance (W/sq m)")
summary(marsModel1, digits =1, style = "pmax")

dat = data.frame(x = Temp, y = GHI)

plot(dat)
curve((4e+02
       + 24 * pmax(0, x -   19) 
       - 17 * pmax(0,   30 - x) 
       -18*pmax(0,x- 30)), #can include a shift factor Mean=375.64
      add=TRUE, col="red",lwd=2) 

#For ggplot
ggplot(dat, aes(x , y)) + 
  geom_point() +
  stat_function(fun = function(x)
    4e+02+24*pmax(0,x-19) 
    -17*pmax(0,30-x) 
    -18*pmax(0,x-30),col="red",lwd=1)

TvGHI <- ggplot(data1)+aes(x=Temp,y=GHI)+ geom_point(position="jitter")+
  labs(x = "Temperature (deg C)", y ="Hourly GHI (W/sqre m)") +
  stat_function(fun = function(x)
    4e+02+24*pmax(0,x-19)-17*pmax(0,30-x)-18*pmax(0,x-30),col="red",lwd=1)
plot(TvGHI)

plot(Temp,GHI)
abline(v=22,lwd=2,col="red")
abline(h=750,lwd=2,col="blue")

#-------------------------------------------------------------------------------------------------------
################# GHI ~ RH
#-------------------------------------------------------------------------------------------------------
data2 <- data.frame(cbind(GHI,RH))
head(data2)

set.seed(1235)
marsModel2 <- earth(GHI ~ ., data=data2) # build model
plotmo(marsModel2,col="blue",xlab=" Relative humidity (%)",
       main="MARS model", ylab="Hourly GHI (W/sq m)")
summary(marsModel2, digits = 2, style = "pmax")

dat = data.frame(x = RH, y = GHI)

plot(dat)
curve((5.4e+02 +1.1 * pmax(0, 48 - x)-12 * pmax(0, x - 48)
       +  6.8 * pmax(0, x - 84)), #can include a shift factor Mean=375.64
      add=TRUE, col="red",lwd=2) 


#For ggplot
ggplot(dat, aes(x , y)) + 
  geom_point() +
  stat_function(fun = function(x)
    5.4e+02 +1.1 * pmax(0, 48 - x)-12 * pmax(0, x - 48)
    +  6.8 * pmax(0, x - 84),col="red",lwd=1)

RHvGHI <- ggplot(data2)+aes(x=RH,y=GHI)+ geom_point(position="jitter")+
  labs(x = "Relative humidity (%)", y ="Hourly GHI (W/sqre m)") +
  stat_function(fun = function(x)
    5.4e+02 +1.1 * pmax(0, 48 - x)-12 * pmax(0, x - 48)
    +  6.8 * pmax(0, x - 84),col="red",lwd=1)
plot(RHvGHI)

plot(RH,GHI)
abline(v=70,lwd=2,col="red")
abline(h=500,lwd=2,col="blue")

#---------------------------------------
# Univariate extreme value theory
#------------------------------------------------------
# Estimating thresholds using parametric mixture models
#------------------------------------------------------
#--------------------------------------
#PARAMETRIC EXTREMAL MIXTURE MODELLING
#---------------------------------------
# TEMP
#--------
library(evmix)

par(mfrow = c(1, 1))
a=Temp
aa = seq(0, 50, 5)
y = dweibullgpd(aa)
hist(a, breaks = 100, freq = FALSE, 
     xlab ="Temperature (deg C)",main="", ylab = "Density",ylim = c(0,0.085))
lines(density(a,adjust=2))  
box()

# The extreme value mixture model with a (truncated) weibull distribution for the bulk and
# GPD upper tail, with bulk model based tail fraction 
# is fitted by default

fit.bulk = fweibullgpd(a)
with(fit.bulk, lines(aa, dweibullgpd(aa,
                                     wshape, wscale, u, sigmau, xi, phiu), col = "red"))
abline(v = fit.bulk$u, col = "red", lty = 2)
fit.bulk

# and parameterised tail fraction requires the 
#option phiu=FALSE to be set:
fit.par = fweibullgpd(a, phiu = FALSE)
with(fit.par, lines(aa, dweibullgpd(aa,  
                                    wshape, wscale, u, sigmau, xi, phiu ), 
                    main=" Histogram of Temperature", col = "blue"))
abline(v = fit.par$u, col = "blue", lty = 2)
legend("topleft", c("True Density", "Bulk Tail Fraction",
                    "Parameterised Tail Fraction"), col=c("black", "red",
                                                          "blue"), lty = 1)
lines(density(a,adjust=2)) 
fit.par

## Diagnostic plots for assessing model fit
par(mfrow=c(1,1))

evmix.diag(fit.bulk)

win.graph()
evmix.diag(fit.par)

#---------------------------------------------------

library(texmex)
library(gridExtra)
palette(c("black","purple","cyan","orange"))
set.seed(20120118)

ggplot(data=data.frame(Temp,index=1:length(Temp)),
       aes(index,Temp))+ geom_point(alpha=0.5,col=4)

########################################################################################################
# RH
#--------
library(evmix)
win.graph()
par(mfrow = c(1, 1))
a=RH
aa = seq(0, 125, 5)
y = dweibullgpd(aa)
hist(a, breaks = 100, freq = FALSE, 
     xlab ="Relative humidity (%)",main="", ylab = "Density",ylim = c(0,0.020))
lines(density(a,adjust=2))  
box()

# The extreme value mixture model with a (truncated) weibull distribution for the bulk and
# GPD upper tail, with bulk model based tail fraction 
# is fitted by default

fit.bulk = fweibullgpd(a)
with(fit.bulk, lines(aa, dweibullgpd(aa,
                                     wshape, wscale, u, sigmau, xi, phiu), col = "red"))
abline(v = fit.bulk$u, col = "red", lty = 2)
fit.bulk

# and parameterised tail fraction requires the 
#option phiu=FALSE to be set:
fit.par = fweibullgpd(a, phiu = FALSE)
with(fit.par, lines(aa, dweibullgpd(aa,  
                                    wshape, wscale, u, sigmau, xi, phiu ), 
                    main=" Histogram of Relative humidity", col = "blue"))
abline(v = fit.par$u, col = "blue", lty = 2)
legend("topleft", c("True Density", "Bulk Tail Fraction",
                    "Parameterised Tail Fraction"), col=c("black", "red",
                                                          "blue"), lty = 1)
lines(density(a,adjust=2)) 
fit.par

## Diagnostic plots for assessing model fit
par(mfrow=c(1,1))

evmix.diag(fit.bulk)

win.graph()
evmix.diag(fit.par)

#---------------------------------------------------

library(texmex)
library(gridExtra)
palette(c("black","purple","cyan","orange"))
set.seed(20120118)

ggplot(data=data.frame(RH,index=1:length(RH)),
       aes(index,RH))+ geom_point(alpha=0.5,col=4)


#--------------------------------------------------
#30.9595096  3.4938320 -0.2627935
#QUANTILE ESTIMATION for RH

tau=31
psi= -0.26279
sigma= 3.494
p <- c(0.1,0.05, 0.03, 0.01,0.001,0.0001,0.00001,
       0.000001,0.0000001,0.00000001,0.000000001,0.0000000001)

x=tau+(sigma)/(psi)*((p)^(-psi) - 1) 
x

T <- c(37.04,38.24,39.01,40.33,42.13,43.11,43.65,43.94,44.10,44.19,44.24, 44.26)

GHI_t = 4e+02+24*pmax(0,T-19)-18*pmax(0,T-30)+375.64 
# 375.64 mean GHI
GHI_t

par(mfrow=c(1,2))
plot(T,type="l", ylab = "Temeperature")
plot(GHI_t,type="l", ylab = "GHI")

#--------------------------------------------------
#89.8992267  8.6447229 -0.8522919
#QUANTILE ESTIMATION for RH
tau=90
psi= -0.26279
sigma= 3.494
p <- c(0.1,0.05, 0.03, 0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001,
       0.000000001,0.0000000001)

x=tau+(sigma)/(psi)*((p)^(-psi) - 1) 
x

T <- c(96.04, 97.24, 98.01, 99.33, 101.13, 102.11, 102.65, 102.94, 103.10, 103.19,
       103.24, 103.26)

GHI_rh = 5.4e+02 -12 * pmax(0, T - 48) +  6.8 * pmax(0, T - 84)+375.64 
# 375.64 mean GHI
GHI_rh

par(mfrow=c(1,2))
plot(T,type="l", ylab = "Relative humidity")
plot(GHI_rh,type="l", ylab = "GHI")


########################################################################################################
#####   Extremal dependence#########################
########################################################################################################
library(qgam)

t <-1:length(GHI)
length(t)
tdata <- cbind(GHI,Temp,RH,t)
data <- data.frame(tdata)
# Calibrate learning rate on a grid
set.seed(5235)
tun <- tuneLearnFast(form=GHI~s(t, bs="ad"), data=data,err = 0.05, qu = 0.9) 
tun

fit1 <-qgam(GHI~s(t, bs="ad"), err = 0.05, qu = 0.9, lsig = tun$lsig, data = data)# insample
summary(fit1, se="boot") #se =" ker" " nid" "boot"
lines(fit1$fit, col="red",lwd=2)
fits <- fitted(fit1)
length(fits)

excess_GHI <- data.frame(GHI-fits)
excess_GHI <- round(excess_GHI,4)
excess_GHI

write.table(excess_GHI,"~/excess_UNV_GHI.txt",sep="\t")
#write.table(dpdfits,"~/its_UNV_GHI.txt",sep="\t")

library("evmix")
GHI. <- excess_GHI[excess_GHI>4.14]
min(GHI.)
max(GHI.)
plot(GHI.)
#########################################################################################
### TEMPERATURE#########################
#dataT <- data.frame(dataT)
# Calibrate learning rate on a grid
set.seed(5235)
tun <- tuneLearnFast(form=Temp~s(t, bs="ad"), data=data,err = 0.05, qu = 0.9) 
tun

fit2 <-qgam(Temp~s(t, bs="ad"), err = 0.05, qu = 0.9, lsig = tun$lsig, data = data)# insample
summary(fit2, se="boot") #se =" ker" " nid" "boot"
lines(fit2$fit, col="red",lwd=2)

fits1 <- fitted(fit2)
length(fits1)
length(Temp)

excess_Temp <- data.frame(Temp-fits1)
excess_Temp <- round(excess_Temp,4)

write.table(excess_Temp,"~/excess_UNV_Temp.txt",sep="\t")

#write.table(dpdfits,"~/its_UNV_GHI.txt",sep="\t")

Temperature <- excess_Temp[excess_Temp>0]

min(Temperature)
max(Temperature)

################################ MODELLING USING TEXMEX ##################################
library(evgam)
library(evd)
library(texmex)
library(evmix)
library(gridExtra)

length(Temperature)
length(GHI.)

TempGHI <- data.frame(cbind(Temperature,GHI.))
TempGHI <- na.omit(TempGHI)

plot(Temperature,GHI., col="blue",xlab="Temperature positive exceedances", 
     ylab="GHI positive exceedances",
     main="Vuwani radiometric station")

#palette(c("black","purple","cyan","orange"))
set.seed(20130618)
summary(TempGHI,digits=2)
GGally::ggpairs(TempGHI)

### examining pairwise extremal dependence using the 
## multivariate conditional Spearman's correlation coefficient
mcsTempGHI <- MCS(TempGHI[, c("Temperature", "GHI.")])

g1 <- ggplot(mcsTempGHI, main="MCS: Temperature and GHI.")
gridExtra::grid.arrange(g1,ncol=1)

## Bootstrapping LEAVE OUT FOR NOW

bootmcsMara <- bootMCS(TempGHI[, c("Temperature", "GHI.")],trace=1000)
g1 <- ggplot(bootmcsMara, main="MCS: Temperature and GHI.")
gridExtra::grid.arrange(g1,ncol=1)

###########################################################
## Conditional multivariate extreme value modelling ######
##########################################################
TempGHI <- data.frame(cbind(Temperature,GHI.))

extremal(Temperature)
extremal(GHI.)

## We need declustering
# Model fitting (marginal modelling)
#mex.Temp <- mex(Temp, mqu=0.75, penalty="none", which="excessThapos")
#mex.Temp

## Alternative approach 
marg <- migpd(TempGHI, mqu=0.90, penalty="none")
marg
mex.Temp <- mexDependence(marg, which = "GHI.")
mex.Temp


## Marginal model diagnostics
g <- ggplot(marg)
do.call("grid.arrange", c(g[[1]], list(ncol=2, nrow=2))) 
do.call("grid.arrange", c(g[[2]], list(ncol=2, nrow=2)))


## Dependence model diagnostics
ggplot(mex.Temp)
mrf <- mexRangeFit(marg, "GHI.", trace=11)
ggplot(mrf)

#start <- coef(mex.Temp$dependence)[1:2,] # alternative starting value
#mrf <- mexRangeFit(marg, "excessThapos", trace=11,start=c(0.1,0.1))
#ggplot(mrf)
### Fitted model parameters
mex.Temp


#########################################################################################
### Relative Humidity#########################
#dataT <- data.frame(dataT)
# Calibrate learning rate on a grid
set.seed(5235)
tun <- tuneLearnFast(form=RH~s(t, bs="ad"), data=data,err = 0.05, qu = 0.9) 
tun

fit3 <-qgam(RH~s(t, bs="ad"), err = 0.05, qu = 0.9, lsig = tun$lsig, data = data)# insample
summary(fit3, se="boot") #se =" ker" " nid" "boot"
lines(fit3$fit, col="red",lwd=2)

fits2 <- fitted(fit3)
length(fits2)
length(RH)

excess_RH <- data.frame(RH-fits2)
excess_RH <- round(excess_RH,4)

write.table(excess_RH,"~/excess_UNV_Temp.txt",sep="\t")

#write.table(dpdfits,"~/its_UNV_GHI.txt",sep="\t")
RH. <- excess_RH[excess_RH>0.271]
min(RH.)
max(RH.)

################################ MODELLING USING TEXMEX ##################################
length(RH.)
length(GHI.)

RHGHI <- data.frame(cbind(RH.,GHI.))


plot(RH.,GHI., col="blue",xlab="RH positive exceedances", 
     ylab="GHI positive exceedances",
     main="Vuwani radiometric station")

#palette(c("black","purple","cyan","orange"))
set.seed(20130618)
summary(RHGHI,digits=2)
GGally::ggpairs(RHGHI)

### examining pairwise extremal dependence using the 
## multivariate conditional Spearman's correlation coefficient
mcsRHGHI <- MCS(RHGHI[, c("RH.", "GHI.")])

g1 <- ggplot(mcsRHGHI, main="MCS: RH and GHI")
gridExtra::grid.arrange(g1,ncol=1)

## Bootstrapping LEAVE OUT FOR NOW

bootmcsMara <- bootMCS(RHGHI[, c("RH.", "GHI.")],trace=1000)
g1 <- ggplot(bootmcsMara, main="MCS: RH and GHI")
gridExtra::grid.arrange(g1,ncol=1)

###########################################################
## Conditional multivariate extreme value modelling ######
##########################################################
RHGHI <- data.frame(cbind(RH.,GHI.))

extremal(RH.)
extremal(GHI.)

## We need declustering
# Model fitting (marginal modelling)

## Alternative approach 
marg <- migpd(RHGHI, mqu=0.90, penalty="none")
marg
mex.RH <- mexDependence(marg, which = "GHI.")
mex.RH


## Marginal model diagnostics
g <- ggplot(marg)
do.call("grid.arrange", c(g[[1]], list(ncol=2, nrow=2))) 
do.call("grid.arrange", c(g[[2]], list(ncol=2, nrow=2)))


## Dependence model diagnostics
ggplot(mex.RH)
mrf <- mexRangeFit(marg, "GHI.", trace=11)
ggplot(mrf)

#start <- coef(mex.Temp$dependence)[1:2,] # alternative starting value
#mrf <- mexRangeFit(marg, "excessThapos", trace=11,start=c(0.1,0.1))
#ggplot(mrf)
### Fitted model parameters
mex.RH


#-----------------------------------
# ARCHIMEDEAN COPULA MODELLING     #
#-----------------------------------
############################################

# Copula package
library(copula)
# Fancy 3D plain scatterplots
library(scatterplot3d)
# ggplot2
library(ggplot2)
# Useful package to set ggplot plots one next to the other
library(grid)
set.seed(235)

#-----------------------------------------------
# Estimate copula parameters for Temp vs GHI   #
#-----------------------------------------------
data1 <- data.frame(cbind(Temp,GHI))
head(data1)

m <- pobs(as.matrix(data1))
head(m)

plot(m, col="blue")

library(psych)
pairs.panels(data1, method="kendall") # data1, m
corr.test(data1,method="kendall")# data1, m "spearman"

#-------------------------------
# Correlation
#------------------------------------
# Measure association using Kendall's Tau and Spearman
cor(data1, method = "kendall") #data1, m

# cor(data1, method = "spearman") # data 1, m


#------------------------------------
# Test for independence
#--------------------------------------



#----------------------------------------
# Clayton
#------------------------------------------
cla_model <- claytonCopula(dim = 2)
fit1 <- fitCopula(cla_model, m, method = 'ml')
summary(fit1)
coef(fit1)
AIC(fit1)
BIC(fit1)
logLik(fit1)

##param 0.7097276 # Check Kendall's tau value for 
# the Clayton copula with theta = 0.7097
tau(claytonCopula(param = 0.7097))
## [1] 0.2619109


#-------------------------------------------------
lower.tail.dep.clayton <- 1/(2^(0.7097))
lower.tail.dep.clayton

upper.tail.dep.clayton = 0

tau.clayton <- (0.7097)/(0.7097+2)
tau.clayton
#---------------------------------------------------

#------------------
# Goodness of fit test (gof)
#-------------------------------------
?gofCopula
gof.clayton <- gofCopula(claytonCopula(dim=2), m, N=27374)
gof.clayton

#--------------------------------------------
# Frank
#----------------------------------------------
fra_model <- frankCopula(dim = 2)
fit2 <- fitCopula(fra_model, m, method = 'ml')
summary(fit2)
coef(fit2)
AIC(fit2)
BIC(fit2)
logLik(fit2)

##param 4.039385 # Check Kendall's tau value for the 
# Frank copula with theta = 4.0394
tau(frankCopula(param = 4.0394))
## [1] 0.3911214
#-------------------------------------
lower.tail.dep.frank = 0

upper.tail.dep.frank = 0

tau.frank <- 0.3911
#-----------------------------------------

gof.frank <- gofCopula(frankCopula(dim=2), m, N=27374)
gof.frank

#---------------------------------------------
# Gumbel
#---------------------------------------------
gum_model <- gumbelCopula(dim = 2)
fit3 <- fitCopula(gum_model, m, method = 'ml')
summary(fit3)
coef(fit3)
AIC(fit3)
BIC(fit3)
logLik(fit3)

##param 1.5099 # Check Kendall's tau value for 
# the Frank copula with theta = 1.5099
tau(gumbelCopula(param = 1.5099))
## [1] 0.3377483
#---------------------------------------------
lower.tail.dep.gumbel = 0

upper.tail.dep.gumbel = 2-(2^(-1.5099))
upper.tail.dep.gumbel

tau.gumbel <- (1.5099-1)/(1.5099)
tau.gumbel
#--------------------------------------------------

gof.gumble <- gofCopula(gumbelCopula(dim=2), data1, N=27374)
gof.gumble

#-------------------------------
# Density and Contour plots
#----------------------------------------------------
clayton <- claytonCopula(dim = 2, param = 0.71)
frank <- frankCopula(dim = 2, param = 4.04)
gumbel <- gumbelCopula(dim = 2, param = 1.51)

par(mfrow = c(1,3))

#--------------------
# Density plot
#-----------------------------------
persp(clayton, dCopula, main ="Clayton copula density")
persp(frank, dCopula, main ="Frank copula density")
persp(gumbel, dCopula, main ="Gumbel copula density")

#----------------------------------
# Contour plot of the densities
#--------------------------------------------
contour(clayton, dCopula, xlim = c(0, 1), ylim=c(0, 1), 
        main = "Contour plot Clayton")
contour(frank, dCopula, xlim = c(0, 1), ylim=c(0, 1), 
        main = "Contour plot Frank")
contour(gumbel, dCopula, xlim = c(0, 1), ylim=c(0, 1),
        main = "Contour plot Gumbel")

#--------------------------------------------------
# MIXTURE OF ARCHIMEDEAN COPULAS
#------------------------------------
mC <- mixCopula(list(frankCopula(4.04, dim = 2), gumbelCopula(1.51, dim = 2)))
mC
stopifnot(dim(mC) == 2)
set.seed(17)
uM <- rCopula(600, mC)
splom2(uM, main = "mixCopula( (clayton, gumbel) )")
d.uM <- dCopula(uM, mC)
p.uM <- pCopula(uM, mC)

fit <- fitCopula(mC, m, method = 'ml')
summary(fit)
coef(fit)
AIC(fit)
BIC(fit)
logLik(fit)

#-------------------------------
# Density and Contour plots
#----------------------------------------------------
par(mfrow = c(1,2))
#--------------------
# Density plot
#-----------------------------------
persp(mC, dCopula, main ="MixCopula density")

#----------------------------------
# Contour plot of the densities
#--------------------------------------------
contour(mC, dCopula, xlim = c(0, 1), ylim=c(0, 1), 
        main = "Contour plot MixCopula")

#-----------------------------------------------#
#-----------------------------------------------#
# Estimate copula parameters for RH vs GHI      #
#-----------------------------------------------#
data2 <- data.frame(cbind(RH,GHI))
head(data2)
m1 <- pobs(as.matrix(data2))
head(m1)

plot(m1,col="blue")

pairs.panels(data2)

library(psych)
pairs.panels(data2, method="kendall") # data1, m
corr.test(data2,method="kendall")# data1, m "spearman"

#-------------------------------
# Correlation
#------------------------------------
# Measure association using Kendall's Tau
cor(data2, method = "kendall")

# cor(data2, method = "spearman")

#----------------------------------
# Clayton
#-----------------------------------------
cla_model <- claytonCopula(dim = 2)
fit4 <- fitCopula(cla_model, m1, method = 'itau')
summary(fit4)
coef(fit4)
AIC(fit4)
BIC(fit4)
logLik(fit4)

##param -0.5749056 # Check Kendall's tau value for 
# the Clayton copula with theta = -0.57
tau(claytonCopula(param = -0.5749))
## [1] -0.4034103

gof.claytonRH <- gofCopula(claytonCopula(dim=2), data2,
                           N=50)
gof.claytonRH

#----------------------------------------
# Frank
#------------------------------------------
fra_model <- frankCopula(dim = 2)
fit5 <- fitCopula(fra_model, m1, method = 'ml')
summary(fit5)
coef(fit5)
AIC(fit5)
BIC(fit5)
logLik(fit5)

##param -4.098635 # Check Kendall's tau value for 
# the Frank copula with theta = -4.0986
tau(frankCopula(param = -4.0986))
## [1] -0.3954416

gof.frankRH <- gofCopula(frankCopula(dim=2), data1, N=50)
gof.frankRH

#------------------------
# Gumbel
#----------------------------------------
gum_model <- gumbelCopula(dim = 2)
fit6 <- fitCopula(gum_model, m1, method = 'mpl')
summary(fit6)
coef(fit6)
AIC(fit6)
BIC(fit6)
logLik(fit6)
## param 1 # Check Kendall's tau value for the Frank 
## copula with theta = 1
tau(gumbelCopula(param = 1))
## [1] 0


gof.gumbleRH <- gofCopula(gumbelCopula(dim=2), data1, N= 50)
gof.gumbleRH


#-------------------------------
# Correlation
#------------------------------------
# Measure association using Kendall's Tau
cor(data2, method = "kendall")

# cor(data2, method = "spearman")

#-------------------------------
# Density and Contour plots
#----------------------------------------------------
clayton <- claytonCopula(dim = 2, param = -0.5749)
frank <- frankCopula(dim = 2, param = -4.0986)
gumbel <- gumbelCopula(dim = 2, param = 1)

par(mfrow = c(1,3))

# Density plot
persp(clayton, dCopula, main ="Clayton copula density")
persp(frank, dCopula, main ="Frank copula density")
persp(gumbel, dCopula, main ="Gumbel copula density")

# Contour plot of the densities
contour(clayton, dCopula, xlim = c(0, 1), ylim=c(0, 1), 
        main = "Contour plot Clayton")
contour(frank, dCopula, xlim = c(0, 1), ylim=c(0, 1), 
        main = "Contour plot Frank")
contour(gumbel, dCopula, xlim = c(0, 1), ylim=c(0, 1),
        main = "Contour plot Gumbel")



