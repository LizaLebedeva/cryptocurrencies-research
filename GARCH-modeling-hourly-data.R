# estimate GARCH DCC - daily prices in USD

library(lubridate)
library(xts)
library(stargazer)
library(mgarchBEKK)
library(rmgarch)
library(rugarch)
library(fGarch)
library(stargazer)
library(MTS)
library(PerformanceAnalytics)
library(quantmod)
library(car)
library(ggfortify)
library(dplyr)
library(tibble)  # for `rownames_to_column` and `column_to_rownames`

# estimate hourly data

#extract data
data <- read.csv('hour_prices_top_20_USD.csv')
head(data)

# delete index
data <- data[,-1]
head(data)

# check for NA
sum(is.na(data))
subset(data, is.na(ETH)) 

# remove NA rows
data <- na.omit(data)
sum(is.na(data))

summary(data)
str(data)
head(data)


# convert to time series
# format = "%Y-%m-%d %H:%M:%S"
# format = "%m/%d/%y %H:%M"
data_ts <- xts(data[,-1], order.by = as.POSIXct(data$timestamp, format = "%m/%d/%y %H:%M"))
head(data_ts)
periodicity(data_ts)


descr_stat <- data.frame(0, 0, 0)
descr_stat
names(descr_stat) <- c('coin', 'mean', 'sd')

# descriptive statistics
for (i in 1:23){
  d <- c(names(data_ts)[i], round(mean(data_ts[,i]), digits = 3), round(sd(data_ts[,i]), digits = 3))
  #print(d)
  descr_stat <- rbind(descr_stat, d)
}

descr_stat <- descr_stat[-1,]
descr_stat



for (i in 1:length(names(data_ts))){
  data_ts[,names(data_ts)[i]]<-diff.xts(log(data_ts[,names(data_ts)[i]]))
}

#remove NA
data_ts <-  na.omit(data_ts) 
head(data_ts)

descr_stat <- data.frame(0, 0, 0)
descr_stat
names(descr_stat) <- c('coin', 'mean', 'sd')

# descriptive statistics
for (i in 1:23){
  d <- c(names(data_ts)[i], round(mean(data_ts[,i]), digits = 3), round(sd(data_ts[,i]), digits = 3))
  #print(d)
  descr_stat <- rbind(descr_stat, d)
}

descr_stat[-1,]


# plot density
autoplot(density(data_ts$BTC), fill = 'green')
autoplot(density(data_ts$ETH), fill = 'green')
autoplot(density(data_ts$ETC), fill = 'green')
autoplot(density(data_ts$BCH), fill = 'green')
MarchTest(data_ts)
# there is conditional heteroscedasticity


# DCCtest(data_ts[,c(1, 2,3,4,5)])
DCCtest(data_ts)
# reject H0: not constant correlation



# normality of error check
jarque.bera.test(data_ts[,1])
jarque.bera.test(data_ts[,2])
# errors are not normally distributed


# with normal distribution of innovations
coins_number = length(names(data_ts))

# define a DCCspec object: 2 stage estimation should usually always use Normal for 1-stage
xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)), 
                   variance.model = list(garchOrder = c(1,1), model = 'sGARCH'), 
                   distribution.model = 'norm')
uspec = multispec(replicate(coins_number, xspec))
spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
spec1a = dccspec(uspec = uspec, dccOrder = c(1, 1), model='aDCC', distribution = 'mvnorm')

# Since multiple DCC models are being estimated on the same dataset with the same first stage dynamics, 
# we can estimate the first stage once and pass it to the dccfit routine (rather than re-estimating it 
# every time dccfit is called)
cl = makePSOCKcluster(coins_number)
cl
multf = multifit(uspec, data_ts, cluster = cl)
multf
# the second stage of the DCC model is estimated:
dcc_norm = dccfit(spec1, data = data_ts, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
dcc_norm
dcca_norm = dccfit(spec1a, data = data_ts, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
dcca_norm

# To fit the DCC (MVT) model in practice, one either assumes a first stage QML, else must jointly 
# estimate in 1 stage the common shape parameter. In the example that follows below, an alternative 
# approach is used to approximate the common shape parameter

# First Estimate a QML first stage model (multf already estimated). Then
# estimate the second stage shape parameter.
spec3 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvt')
dcc_mvt = dccfit(spec3, data = data_ts, fit.control = list(eval.se = FALSE), fit = multf)
dcc_mvt
# obtain the multivariate shape parameter:
mvt.shape = rshape(dcc_mvt)
# Plug that into a fixed first stage model and iterate :
mvt.l = rep(0, 6)
mvt.s = rep(0, 6)
mvt.l[1] = likelihood(dcc_mvt)
mvt.s[1] = mvt.shape
for (i in 1:5) {
  xspec = ugarchspec(mean.model = list(armaOrder = c(1, 1)), variance.model = list(garchOrder = c(1,1), model = 'sGARCH'), distribution.model = 'std', fixed.pars = list(shape = mvt.shape))
  spec3 = dccspec(uspec = multispec(replicate(coins_number, xspec)), dccOrder = c(1,1), distribution = 'mvt')
  dcc_mvt = dccfit(spec3, data = data_ts, solver = 'solnp', fit.control = list(eval.se = FALSE))
  mvt.shape = rshape(dcc_mvt)
  mvt.l[i + 1] = likelihood(dcc_mvt)
  mvt.s[i + 1] = mvt.shape
}


# plot to see that they converge
mvt.l
mvt.s
par(mar=c(5, 4, 4, 6) + 0.1)
plot(mvt.s,type='l', col="black", xlab='Iterations')
par(new=T)
plot(mvt.l,type='l',axes = FALSE, bty = "n", xlab = "", ylab = "",col="red")
axis(side=4, col = 'red', col.axis="red")
mtext("lik", side=4, line=3, col = 'red')
legend("topleft",legend=c("shape","log.lik"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))

# Finally, once more, fixing the second stage shape parameter, and
# evaluating the standard errors
xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)), variance.model = list(garchOrder = c(1,1), model = 'sGARCH'), distribution.model = 'std', fixed.pars = list(shape = mvt.shape))
spec3 = dccspec(uspec = multispec(replicate(coins_number, xspec)), dccOrder = c(1, 1), distribution = 'mvt', fixed.pars = list(shape = mvt.shape))
dcc_mvt = dccfit(spec3, data = data_ts, solver = 'solnp', fit.control = list(eval.se = TRUE), cluster = cl)
dcc_mvt

# for adcc model
xspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(garchOrder = c(1,1), model = "sGARCH"),  distribution.model = "norm")
spec3a  = dccspec(uspec = multispec( replicate(coins_number, xspec) ), dccOrder = c(1,1), distribution = "mvt", model="aDCC")
dcca_mvt = dccfit(spec3a, data = data_ts, fit.control = list(eval.se=FALSE), fit = multf)
dcca_mvt
# obtain the multivariate shape parameter:
mvtx.shape = rshape(dcca_mvt)
# Plug that into a fixed first stage model and iterate :
mvtx.l = rep(0, 6)
mvtx.s = rep(0, 6)
mvtx.l[1] = likelihood(dcca_mvt)
mvtx.s[1] = mvtx.shape
for(i in 1:5){
  xspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(garchOrder = c(1,1), model = "sGARCH"),  distribution.model = "std", fixed.pars = list(shape=mvtx.shape))
  spec3a = dccspec(uspec = multispec( replicate(coins_number, xspec) ), dccOrder = c(1,1), model="aDCC", distribution = "mvt")
  dcca_mvt = dccfit(spec3a, data = data_ts, solver = "solnp", fit.control = list(eval.se=FALSE))
  mvtx.shape = rshape(dcca_mvt)
  mvtx.l[i+1] = likelihood(dcca_mvt)
  mvtx.s[i+1] = mvtx.shape
}

# plot to see that they converge
mvtx.l
mvtx.s
par(mar=c(5, 4, 4, 6) + 0.1)
plot(mvtx.s,type='l', col="black", xlab='Iterations')
par(new=T)
plot(mvtx.l,type='l',axes = FALSE, bty = "n", xlab = "", ylab = "",col="red")
axis(side=4, col = 'red', col.axis="red")
mtext("lik", side=4, line=3, col = 'red')
legend("topleft",legend=c("shape","log.lik"),
       text.col=c("black","red"),pch=c(16,15),col=c("black","red"))


# Finally, once more, fixing the second stage shape parameter, and evaluating the standard errors
xspec = ugarchspec(mean.model = list(armaOrder = c(0,0)), variance.model = list(garchOrder = c(1,1), model = "sGARCH"),  distribution.model = "std", fixed.pars = list(shape=mvtx.shape))
spec3a = dccspec(uspec = multispec( replicate(coins_number, xspec) ), dccOrder = c(1,1), model="aDCC", distribution = "mvt", fixed.pars=list(shape=mvtx.shape))
dcca_mvt = dccfit(spec3a, data = data_ts, solver = "solnp", fit.control = list(eval.se=TRUE), cluster = cl)
dcca_mvt

# it is important to remember to close your clusters using
stopCluster(cl)

# str(dcc_norm)
# names(dcc_norm@model)
# "modelinc"      "modeldesc"     "modeldata"     "varmodel"      "pars"          "start.pars"   
# [7] "fixed.pars"    "maxgarchOrder" "maxdccOrder"   "pos.matrix"    "pidx"          "DCC"          
# [13] "mu"            "residuals"     "sigma"         "mpars"         "ipars"         "midx"         
# [19] "eidx"          "umodel"       
  

# names(dcc_norm@mfit)
# [1] "coef"            "matcoef"         "garchnames"      "dccnames"        "cvar"           
# [6] "scores"          "R"               "H"               "Q"               "stdresid"       
# [11] "llh"             "log.likelihoods" "timer"           "convergence"     "Nbar"           
# [16] "Qbar"            "plik" 

dcc_norm@mfit$garchnames
dcc_norm@model$ipars
dcc_norm@model$mpars

dcc_norm@mfit$matcoef
dcc_norm@mfit$Q

show(dcc_norm)
coef(dcc_norm)
# The fitted conditional GARCH sigma xts object.
sigma(dcc_norm)
plot(sigma(dcc_norm)[,1])
plot(sigma(dcc_norm))
# The fitted conditional mean xts object
fitted(dcc_norm)
# The fitted conditional mean residuals xts object.
residuals(dcc_norm)
plot(residuals(dcc_norm))
plot(dcc_norm)
dev.off()

# plot correlation
plot(timeSeries(rcor(dcc_norm)[1,2,])) #BTC-ETH
plot(timeSeries(rcor(dcc_norm)[2,1,]))#ETH-BTC the same
plot(timeSeries(rcor(dcc_norm)[1,13,]))#ETH-BTC the same

dev.off()
par(mfrow = c(2,2))

plot.xts(sigma(dcc_norm), legend.loc = 'topleft', main = 'Fitted conditional GARCH sigma')
plot.xts(sigma(dcca_norm), legend.loc = 'topleft', main = 'Fitted conditional GARCH sigma')
plot.xts(sigma(dcc_mvt), legend.loc = 'topleft', main = 'Fitted conditional GARCH sigma')
plot.xts(sigma(dcca_mvt), legend.loc = 'topleft', main = 'Fitted conditional GARCH sigma')


plot.xts(residuals(dcc_norm), legend.loc = 'topleft', main = 'Residuals')

# infocriteria
infocriteria(dcc_norm)

likelihood(dcc_mvt)
rshape((dcc_mvt))

# models = c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)
# rm(models)

# likelihoods
for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  print(likelihood(model))
}

# information criteria
for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  print(infocriteria(model))
}

# coefficients
for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  printCoefmat(model@mfit$matcoef, digits=3)
}

# number of significant coefficients
for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  print(sum(model@mfit$matcoef[,4]<0.05))
}

# compare to hourly data
for (model in c(dcc_norm_hourly,dcca_norm_hourly,dcc_mvt_hourly,dcca_mvt_hourly)){
  print(sum(model@mfit$matcoef[,4]<0.05))
}

# look only at betas
dcc_norm@mfit$matcoef[grep("beta1", rownames(model@mfit$matcoef)),]

# look only at aphas
dcc_norm@mfit$matcoef[grep("alpha1", rownames(model@mfit$matcoef)),]

# number of significant alphas
for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  print(sum(model@mfit$matcoef[grep("alpha1", rownames(model@mfit$matcoef)),][,4]<0.05))
}

# number of significant betas
for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  print(sum(model@mfit$matcoef[grep("beta1", rownames(model@mfit$matcoef)),][,4]<0.05))
}

# number of significant omegas
for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  print(sum(model@mfit$matcoef[grep("omega", rownames(model@mfit$matcoef)),][,4]<0.05))
}

# number of significant mus
for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  print(sum(model@mfit$matcoef[grep("mu", rownames(model@mfit$matcoef)),][,4]<0.05))
}

# to see only significant coefficients

for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  print(c(rownames_to_column(as.data.frame(model@mfit$matcoef) %>%
                       rownames_to_column('gene') %>%
                       filter(`Pr(>|t|)`<0.05) %>%
                       column_to_rownames('gene'))$rowname, 'model'))
}

# to see not significant coefficients

for (model in c(dcc_norm,dcca_norm,dcc_mvt,dcca_mvt)){
  print(c(rownames_to_column(as.data.frame(model@mfit$matcoef) %>%
                               rownames_to_column('gene') %>%
                               filter(`Pr(>|t|)`>=0.05) %>%
                               column_to_rownames('gene'))$rowname, 'model'))
}

rownames_to_column(as.data.frame(dcc_norm@mfit$matcoef) %>%
  rownames_to_column('gene') %>%
  filter(`Pr(>|t|)`<0.05) %>%
  column_to_rownames('gene'))$rowname

dcc_norm
dcca_norm
dcc_mvt
dcca_mvt


printCoefmat(dcc_norm@mfit$matcoef)


# export
stargazer(dcc_norm@mfit$matcoef,infocriteria(dcc_norm),
          type="html",
          #summary=FALSE,
          report="vc*",
          out="dcc_norm.html", model.names = TRUE, nobs = TRUE)

stargazer(dcca_norm_hourly@mfit$matcoef,infocriteria(dcca_norm),
          type="html",
          report="vc*",
          out="dcca_norm.html", model.names = TRUE, nobs = TRUE)

stargazer(dcc_mvt_hourly@mfit$matcoef,infocriteria(dcc_mvt),
          type="html",
          report="vc*",
          out="dcc_mvt.html", model.names = TRUE, nobs = TRUE)

stargazer(dcca_mvt_hourly@mfit$matcoef,infocriteria(dcca_mvt),
          type="html",
          report="vc*",
          out="dcca_mvt.html", model.names = TRUE, nobs = TRUE)





stargazer(cbind(dcc_mvt@mfit$matcoef,dcca_mvt@mfit$matcoef[-35,]),
          type="html",
          report="vc*",
          out="dcc2.html", model.names = TRUE, nobs = TRUE)

stargazer(dcc_mvt@mfit$matcoef,dcca_mvt@mfit$matcoef[-35,],
          type="html",
          report="vc*",
          out="dcc2.html", model.names = TRUE, nobs = TRUE)

# cbind(dcc_mvt@mfit$matcoef,dcca_mvt@mfit$matcoef[-35,])

# save 
dcc_norm_hourly <- dcc_norm
dcca_norm_hourly <- dcca_norm
dcc_mvt_hourly <- dcc_mvt
dcca_mvt_hourly <-dcca_mvt
