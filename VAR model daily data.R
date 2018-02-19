# VAR modeling of daily data

rm(list=ls(all=TRUE))
library(lubridate)
library(xts)
library(stargazer)
library(vars)
library(tseries)
library(lmtest)
library(aod)
library(autovarCore)
library(frequencyConnectedness)
library(fastSOM)
library(mgarchBEKK)
library(ggplot2)
library(dplyr)
library(scales)
library(directlabels)


#extract data
daily <- read.csv('prices_top_23_BTC.csv')
View(daily)

# delete index
daily <- daily[,-1]
View(daily)

# check for NA
sum(is.na(daily))
subset(daily, is.na(ETH)) 

# remove NA rows
daily <- na.omit(daily)
sum(is.na(daily))

summary(daily)
str(daily)
head(daily)

# when BTC was max
subset(daily, BTC==max(BTC)) 

# convert to time series
d_prices <- xts(daily[,-1], order.by = as.POSIXct(daily$timestamp, format = "%m/%d/%y %H:%M"))
head(d_prices)

# descriptive statistics
for (i in 1:23){
  print(c(names(d_prices)[i], round(mean(d_prices[,i]), digits = 3), round(sd(d_prices[,i]), digits = 3)))
}

#some functions
ndays(d_prices)
nyears(d_prices)
nhours(d_prices)
periodicity(d_prices)

#plot

plot.xts(d_prices[,c(1:5)],legend.loc = "topright", main = "price")
plot.xts(d_prices[,c(15:18)],legend.loc = "topright", main = "price")


dprices <- daily %>%
  mutate(Date=as.Date(as.POSIXct(daily$timestamp, format = "%m/%d/%y %H:%M")))

plot_d_prices_BTC <- 
  ggplot(dprices, aes(x=Date))+
  geom_line(aes(y=BTC, color='Bitcoin'))+
  geom_line(aes(y=BCH, color='Bitcoin Cash'))

plot_d_prices_BTC

plot_d_prices <- 
  ggplot(dprices, aes(x=Date))+
  #geom_line(aes(y=BTC, color='Bitcoin'))+
  #geom_label(data=last, aes(x = Date, y=BTC, label =
  #                           last$BTC, color='Bitcoin'), show.legend=F)+
  geom_line(aes(y=ETH, color='Ethereum'))+
  #geom_label(data=last, aes(x = Date, y=ETH, label =
  #                           last$ETH, color='Ethereum'), show.legend=F)+
  geom_line(aes(y=LTC, color='Litecoin'))+
  #geom_label(data=last, aes(x = Date, y=LTC, label =
  #                           last$LTC, color='Litecoin'), show.legend=F)+
  geom_line(aes(y=DASH, color='DASH'))+
  #geom_label(data=last, aes(x = Date, y=DASH, label =
  #                           last$DASH, color='DASH'),show.legend=F, nudge_y = 0.3)+
  geom_line(aes(y=XRP, color='Ripple'))+
  #geom_line(aes(y=BCH, color='Bitcoin Cash'))+
  geom_line(aes(y=IOT, color='Miota'))+
  geom_line(aes(y=XMR, color='Monero'))+
  geom_line(aes(y=ETC, color='Etherium Classic'))+
  theme_minimal()+
  #scale_y_continuous(trans=log_trans(), breaks=c(1,10,100,200, 500, 5000,10000))+
  labs(y = "Price, $",
       x = "Time",
       colour = " ")+
  # theme(legend.position = c(0.2, 0.9))+
  scale_x_date(labels = date_format("%d%b"),breaks = date_breaks("3 week"))+
  ggtitle('Prices ')

plot_d_prices
dev.off()
# another way
plot.xts(d_prices[,2:16], legend.loc = "topright", main='Prices')



#Augmented-Dickey-Fuller Unit Root Test

for(i in 1:length(names(d_prices))){
  UR <- summary(ur.df(d_prices[, i], type = "none",lags = 1))
  name = names(d_prices)[i]
  if ((as.numeric(UR@teststat)>as.numeric(UR@cval[3]))==TRUE) {
    print(c(name, 'UR'))
  }else {
    print(c(name, ' no UR'))
  }
}


# all are unit root


#with difference
plot.xts(diff(d_prices[, "ETH"]))

# modify time series to continuously compounded returns
log_d_prices <- daily %>%
  mutate_each(funs(log = (log(.))), -timestamp)

View(daily)
View(d_prices)
head(log_d_prices)

# be careful with column numbers
names(log_d_prices)[c(1, 25:47)]
log_d_prices <- log_d_prices[names(log_d_prices)[c(1, 25:47)]]
names(log_d_prices)

head(log_h_prices)
str(log_h_prices)


# convert to time series
log_d_prices_ts <- xts(log_d_prices[,-1], order.by = as.POSIXct(log_d_prices$timestamp, format = "%m/%d/%y %H:%M"))
head(log_d_prices_ts)

#some functions
ndays(log_d_prices_ts)
nyears(log_d_prices_ts)
nhours(log_d_prices_ts)
periodicity(log_d_prices_ts)


#Augmented-Dickey-Fuller Unit Root Test
summary(ur.df(log_d_prices_ts[, 'BTC_log'], type = "none",lags = 1))


# get difference data
diff(log_d_prices_ts[, "ETH_log"])[-1]
summary(ur.df(diff(log_d_prices_ts[, 1])[-1], type = "none",lags = 1))


for(i in 1:length(names(log_d_prices_ts))){
  UR <- summary(ur.df(diff(log_d_prices_ts[, i])[-1], type = "none",lags = 1))
  name = names(log_d_prices_ts)[i]
  if ((as.numeric(UR@teststat)>as.numeric(UR@cval[3]))==TRUE) {
    print(c(name, 'UR'))
  }else {
    print(c(name, ' no UR'))
  }
}

# every time series in not Unit root


# dataframe with log and difference

log_dif_d_prices_ts <- diff(log_d_prices_ts)
View(log_dif_d_prices_ts)
log_dif_d_prices_ts <- na.omit(log_dif_d_prices_ts)


# change names
names(log_dif_d_prices_ts) <- sub("_.*", "", names(log_dif_d_prices_ts))
names(log_dif_d_prices_ts)

# descriptive statistics
for (i in 1:23){
  print(c(names(log_dif_d_prices_ts)[i], round(mean(log_dif_d_prices_ts[,i]), digits = 3), round(sd(log_dif_d_prices_ts[,i]), digits = 3)))
}


# plot - MAIN
dev.off()
plot.xts(log_dif_d_prices_ts[,c(1, 2, 3, 4, 5, 8, 9, 10)], cex.lab=1.8, cex.axis=1.2,legend.loc = "topright", main = "")
title(main = "Daily returns", font= 14)
legend(pch = 2)

plot.xts(log_dif_d_prices_ts[,c(1, 2, 3, 4, 5, 8)], legend.loc = "topright", main = "Daily returns")

plot.xts(log_dif_d_prices_ts[,1:10], legend.loc = "topright", main = "Log price difference ")


#VAR estimation with lag 1
Var1d <- VAR(log_dif_d_prices_ts, p = 1, type = "both")
# The option type determines whether to include an intercept term, a trend or both in the model

str(Var1d)
Var1d$varresult

# summary results for selected equations
summary(Var1d, equation = "BTC")

# R2
summary(Var1d, equation = "BTC")$varresult$BTC$r.squared

# The F statistic compares the joint effect of all the variables together. 
summary(Var1m, equation = "BTC")$varresult$BTC$fstatistic



# PLOT
# for each equation in a VAR a diagram of fit, a residual plot, 
# the auto-correlation and partial auto-correlation function of the residuals in a single layout
dev.off()
names(d_prices)
# "BTC"   "ETH"   "XRP"   "BCH"   "LTC"   "XEM"   "XLM"   "IOT"   "DASH"  "XMR"   "ETC"  
# "LSK"   "XVG"   "ZEC"   "BCN"   "SC"    "VEN"   "STRAT" "BTS"   "VERI"  "EOS"   "OMG"  
# "DOGE" 
plot(Var1d$varresult$BTC)
plot(Var1d$varresult$ETH)
plot(Var1d$varresult$XRP)
plot(Var1d$varresult$BCH)
plot(Var1d$varresult$LTC)
plot(Var1d$varresult$XEM) #good QQ
plot(Var1d$varresult$XLM)
plot(Var1d$varresult$IOT)
plot(Var1d$varresult$DASH)
plot(Var1d$varresult$XMR) #bad location plot
plot(Var1d$varresult$ETC)
plot(Var1d$varresult$LSK)
plot(Var1d$varresult$DOGE)
plot(Var1d$varresult$BCN)

# Diagram of fit and residuals
plot(Var1d, names = "BTC")
dev.off()
plot(Var1d, names = "ETH")
dev.off()
plot(Var1d, names = "XRP")
dev.off()
plot(Var1d, names = "XEM")
dev.off()
plot(Var1d, names = "XMR")


#export coefficients
Var1d$varresult$BTC

# residuals
residuals(Var1d)
plot(residuals(Var1d))


# for one
residuals(Var1d$varresult$BTC)
plot(residuals(Var1d$varresult$BTC))
plot(residuals(Var1d$varresult$ETH))
plot(residuals(Var1d$varresult$DOGE))
plot(residuals(Var1d$varresult$XMR))
plot(residuals(Var1d$varresult$XEM))

# autocorrelation of residuals
dev.off()
acf(residuals(Var1d)[,1])

# cross-correlation matrix of residuals (change columns - not possible to get all together on one plot)
acf(residuals(Var1d)[,1:5])
acf(residuals(Var1d)[,1:4])
acf(residuals(Var1d)[,5:8])
# The plots along the diagonal are the individual ACFs (autocorrelation functions) for each model’s residuals 
# In addition, we now see the cross-correlation plots of each set of residuals. 
# Ideally, these would also resemble white noise, however we do see remaining cross-correlations


stargazer(Var1d$varresult,
          type="text",
          report="vc*",
          out="Var1d_23BTC.csv",
          # colnames = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = names(log_dif_d_prices_ts))


#to choose best lag
VARselect(log_dif_d_prices_ts, lag.max = 8, type = "both")


criterionD <- VARselect(log_dif_d_prices_ts, lag.max = 5, type = "both")
criterionD
min(criterionD[[2]])

# lag 5 has min AIC


# choose best order of VAR based on AIC
var.aicD <- VAR(log_dif_d_prices_ts,type="both",lag.max=5,ic="AIC")
# names(var.aicD)

# to see that best lag is 5
var.aicD$p

summary(var.aicD)
# best lag is 5 based on AIC


# export
stargazer(var.aicD$varresult,
          type="html",
          report="vc*",
          out="VarAIC_23dBTC.html",
          # colnames = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = names(log_dif_d_prices_ts))


# to compare with lag 2

Var2D <- VAR(log_dif_d_prices_ts, p = 2, type = "both")

# R2
summary(Var2D, equation = "BTC")$varresult$BTC$r.squared

# The F statistic compares the joint effect of all the variables together. 
summary(Var2D, equation = "BTC")$varresult$BTC$fstatistic

# degree of freedom
summary(Var2D, equation = "BTC")$varresult$BTC$df

# residuals
residuals(Var2D)
plot(residuals(Var2D))



dev.off()
# Diagram of fit and residuals
plot(Var2D, names = "BTC")
plot(Var2D, names = "ETH")
plot(var.aicD, names = "ETH")

# export VAR2 results
stargazer(Var2D$varresult,
          type="html",
          report="vc*",
          out="Var2_23dBTC.html",
          # colnames = TRUE,
          dep.var.labels.include = FALSE,
          column.labels = names(log_dif_d_prices_ts))


# Diagnostic testing

# tests residuals for heteroscedasticity - multivariate ARCH-LM test

heter_test_1d <- arch.test(Var1d, lags.multi = 1, multivariate.only = TRUE)
print(heter_test_1d)

heter_test_2d <- arch.test(var.aicD, lags.multi = 2, multivariate.only = TRUE)
print(heter_test_2d)

heter_test_3d <- arch.test(var.aicD, lags.multi = 3, multivariate.only = TRUE)
print(heter_test_3d)

heter_test_5d <- arch.test(var.aicD, lags.multi = 5, multivariate.only = TRUE)
print(heter_test_5)

# singular system for daily data

# small p-value: reject the null hypothesis and hence multivariate time series is heteroskedastic

# About test:
# For q=1, it tests for ARCH effects at lag 1: H0:α1=0H0:α1=0.
# For q=2, it tests for ARCH effects at lags 1 and 2 jointly: H0:α1=α2=0H0:α1=α2=0.
# if p-values decreases as lag becomes larger, That suggests the squared series has 
# significant partial autocorrelations (significant ARCH effects) at higher lag orders as well

# Test shows an ARCH effect at lag 2 and likely higher lag orders, 
# which means  data is conditionally heteroskedastic.

plot(heter_test_1d)

# if figure margins are too large
# dev.off()
# and increase size of plots window


# tests for heteroscedasticity - univariate ARCH-LM test

heter_test_uniD <- arch.test(Var1d, lags.single = 16, lags.multi = 1, multivariate.only = FALSE)
print(heter_test_uniD)

# The Jarque-Bera normality tests for univariate and multivariate series are implemented 
# and applied to the residuals of a VAR(p) as well as separate tests for multivariate 
# skewness and kurtosis

# Jarque–Bera test is a goodness-of-fit test of whether sample data have the skewness and kurtosis matching a normal distribution.
# H0: joint hypothesis of the skewness being zero and the excess kurtosis being zero

# multivariate
normality.test(Var1d)

# large J-B value indicates that errors are not normally distributed.

# univariate
normality.test(Var1d, multivariate.only = FALSE)
names(normality.test(Var1d, multivariate.only = FALSE))
normality.test(Var1d, multivariate.only = FALSE)$jb.uni[1]
str(normality.test(Var1d, multivariate.only = FALSE)$jb.uni[1])
normality.test(Var1d, multivariate.only = FALSE)$jb.uni$BTC$p.value
normality.test(Var1d, multivariate.only = FALSE)$jb.uni[[1]][[3]]

for (i in 1:23){
  print(c(i, normality.test(Var1d, multivariate.only = FALSE)$jb.uni[[i]][[3]]))
}

# Portmanteau Test (asymptotic) - For testing the lack of serial correlation in the residuals of a VAR(p)
serial.test(Var1d)
serial.test(Var1d, type = "ES")

stability(Var1d)
plot(stability(Var1d), nc=2)

autovarCore:::model_is_stable(Var1d)
autovarCore:::assess_portmanteau(Var1d)


# impulse response function

# to work faster set boot=FALSE

ir.1d <- irf(Var1d,impulse="BTC",response="ETH",n.ahead = 20,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.1d)

# The plot gives the response of series 2 for the periods 0 to 20 to a shock in series 1 
# in period 0. The function also automatically calculates so-called bootstrap confidence bands. 
# (Bootstrapping is a common procedure in impulse response analysis. But you ought keep in mind 
#   that it has its drawbacks when you work with structural VAR models.)


# here is better to set ortho = TRUE
ir.1td <- irf(Var1d,impulse="BTC",response="ETH",n.ahead = 20,ortho = TRUE, seed=87, boot = FALSE)
plot(ir.1td)


# for two coins
ir.2d <- irf(Var1d,impulse="BTC",response=c("ETH", "LTC"),n.ahead = 20,ortho = TRUE, seed=87, boot = FALSE)
plot(ir.2d)


# From BTC
# for group of coins short ortho false
ir.2d <- irf(Var1d,impulse="BTC",response=c("BCH", "ETH", "LTC", "XRP", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="BTC",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)


# for group of coins long ortho false
ir.2d <- irf(Var1d,impulse="BTC",response=c("BCH", "ETH", "LTC", "XRP", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="BTC",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)
dev.off()
# for group of coins short ortho true
ir.2d <- irf(Var1d,impulse="BTC",response=c("BCH", "ETH", "LTC", "XRP", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = TRUE, ci = 0.9)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="BTC",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = TRUE, ci = 0.9)
plot(ir.2d)


# for group of coins long ortho true
ir.2d <- irf(Var1d,impulse="BTC",response=c("BCH", "ETH", "LTC", "XRP", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE,cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="BTC",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)


ir.2d <- irf(Var1d,impulse="BTC",response=c("ETH", "LTC", "BCH"),n.ahead = 20,ortho = TRUE, seed=87, boot = FALSE)
plot(ir.2d)


# From ETH
# for group of coins short ortho false
ir.2d <- irf(Var1d,impulse="ETH",response=c("BTC", "BCH", "LTC", "XRP", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="ETH",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)


# for group of coins long ortho false
ir.2d <- irf(Var1d,impulse="ETH",response=c("BCH", "BTC", "LTC", "XRP", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="ETH",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

# for group of coins short ortho true
ir.2d <- irf(Var1d,impulse="ETH",response=c("BCH", "BTC", "LTC", "XRP", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = TRUE, ci = 0.9)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="ETH",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = TRUE, ci = 0.9)
plot(ir.2d)


# for group of coins long ortho true
ir.2d <- irf(Var1d,impulse="ETH",response=c("BCH", "BTC", "LTC", "XRP", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE,cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="ETH",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)


# From XRP
# for group of coins short ortho false
ir.2d <- irf(Var1d,impulse="XRP",response=c("BCH", "ETH", "LTC", "BTC", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="XRP",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)


# for group of coins long ortho false
ir.2d <- irf(Var1d,impulse="XRP",response=c("BCH", "ETH", "LTC", "BTC", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="XRP",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

# for group of coins short ortho true
ir.2d <- irf(Var1d,impulse="XRP",response=c("BCH", "ETH", "LTC", "BTC", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = TRUE, ci = 0.9)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="XRP",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = TRUE, ci = 0.9)
plot(ir.2d)


# for group of coins long ortho true
ir.2d <- irf(Var1d,impulse="XRP",response=c("BCH", "ETH", "LTC", "BTC", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE,cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="XRP",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

# From LTC
# for group of coins short ortho false
ir.2d <- irf(Var1d,impulse="LTC",response=c("BCH", "ETH", "XRP", "BTC", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="LTC",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)


# for group of coins long ortho false
ir.2d <- irf(Var1d,impulse="LTC",response=c("BCH", "ETH", "XRP", "BTC", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="LTC",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

# for group of coins short ortho true
ir.2d <- irf(Var1d,impulse="LTC",response=c("BCH", "ETH", "XRP", "BTC", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = TRUE, ci= 0.9)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="LTC",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE)
plot(ir.2d)


# for group of coins long ortho true
ir.2d <- irf(Var1d,impulse="LTC",response=c("BCH", "ETH", "XRP", "BTC", "XEM", "XLM", "IOT", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE,cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="LTC",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)


# From IOT
# for group of coins short ortho false
ir.2d <- irf(Var1d,impulse="IOT",response=c("BCH", "ETH", "XRP", "BTC", "XEM", "XLM", "LTC", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="IOT",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE)
plot(ir.2d)


# for group of coins long ortho false
ir.2d <- irf(Var1d,impulse="IOT",response=c("BCH", "ETH", "XRP", "BTC", "XEM", "XLM", "LTC", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="IOT",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = FALSE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)

# for group of coins short ortho true
ir.2d <- irf(Var1d,impulse="IOT",response=c( "ETH", "XRP",  "XEM",  "LTC",
                                             "ETC", "OMG"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = TRUE, ci = 0.9)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="IOT",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = TRUE, ci = 0.9)
plot(ir.2d)


# for group of coins long ortho true
ir.2d <- irf(Var1d,impulse="IOT",response=c("BCH", "ETH", "XRP", "BTC", "XEM", "XLM", "LTC", "DASH",
                                            "XMR", "ETC"),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE,cumulative = TRUE)
plot(ir.2d)

ir.2d <- irf(Var1d,impulse="IOT",response=c('LSK','XVG', 'ZEC', 'BCN','SC','VEN',
                                            'STRAT','BTS','VERI', 'EOS','OMG','DOGE'),
             n.ahead = 10,ortho = TRUE, seed=87, boot = FALSE, cumulative = TRUE)
plot(ir.2d)




ir.2d <- irf(Var1d,impulse="BTC",response=c("ETH", "LTC", "BCH"),n.ahead = 20,ortho = TRUE, seed=87, boot = FALSE)
plot(ir.2d)
ir.2d <- irf(Var1d,impulse="BTC",response=c("ETH", "LTC", "BCH"),n.ahead = 20,ortho = TRUE, seed=87, boot = FALSE)
plot(ir.2d)
# long-run effects of a shock
ir.ld <- irf(var.aicD,impulse="BTC",response="ETH",n.ahead = 40,ortho = FALSE,
             cumulative = TRUE,seed=87, boot = FALSE)
plot(ir.ld)

# with ortho
ir.l2d <- irf(var.aicD,impulse="BTC",response="ETH",n.ahead = 40,ortho = TRUE,
              cumulative = TRUE,seed=87, boot=FALSE)
plot(ir.l2d)


# orthogonal cumulative provides more reliable result

# to BTC
# short
ir.2dtoB <- irf(var.aicD,impulse="ETH",response="BTC",n.ahead = 20,ortho = TRUE, seed=87, boot = FALSE)
plot(ir.2dtoB)
# long
ir.ldtoB <- irf(var.aicD,impulse="ETH",response="BTC",n.ahead = 80,ortho = TRUE,
                cumulative = TRUE,seed=87, boot = FALSE)
plot(ir.ldtoB)



# granger causality test for all pairs

# careful with number of coins

jj <- c()
for (i in c(1:length(names(log_dif_d_prices_ts)))) {
  for (j in c(1:length(names(log_dif_d_prices_ts)))) {
    if (i != j) {
      a <- grangertest(log_dif_d_prices_ts[,i], log_dif_d_prices_ts[,j], order=1)[[4]][2]
      if (a < 0.101){
        jj <- c(jj,names(log_dif_d_prices_ts)[j])
        print(c(names(log_dif_d_prices_ts)[i], names(log_dif_d_prices_ts)[j], a)) 
      } 
    }
  }
}

jj
table(jj)


# all pairs Granger test
a <- data.frame(0, 0, 0)
a
names(a) <- c('X', "caused by", 'p')

for (i in 1:length(names(log_dif_d_prices_ts))){
  for (j in 1:length(names(log_dif_d_prices_ts))) {
    if (i != j) {
      g <- grangertest(log_dif_d_prices_ts[,c(j,i)], order=2)[[4]][2]
      c <- c(names(log_dif_d_prices_ts)[i], names(log_dif_d_prices_ts)[j], round(g, digits = 3))
      a <- rbind(a, c)
    }
  }
}


View(a)
a %>%
  filter(p<=0.05)


grangertest(log_dif_d_prices_ts[,22], log_dif_d_prices_ts[,23], order=1)

grangertest(log_dif_d_prices_ts[,2:3], order=1)
grangertest(log_dif_d_prices_ts[,2:3], order=2)

grangertest(log_dif_d_prices_ts[,3], log_dif_d_prices_ts[,2], order=2)
grangertest(log_dif_d_prices_ts[,2], log_dif_d_prices_ts[,3], order=2)

grangertest(log_dif_d_prices_ts[,3], log_dif_d_prices_ts[,2], order=2)[[4]][2]
str(grangertest(log_dif_d_prices_ts[,3], log_dif_d_prices_ts[,2], order=2)[4,3])

grangerD2 <- reshape(a, v.names="p", idvar="X", timevar="caused by", direction="wide")
grangerD2 <- grangerD2 %>%
  filter(X != 0)%>%
  select(-p.0)
grangerD2 #p value: row is caused by columns

View(grangerD2)

# export granger test results
stargazer(grangerD2,
          type="html",
          #report="vc*",
          summary=FALSE,
          out="G2_23dBTC.html",
          colnames = TRUE)


# Granger test for order 1
# all pairs Granger test
a <- data.frame(0, 0, 0)
a
names(a) <- c('X', "caused by", 'p')

for (i in 1:length(names(log_dif_d_prices_ts))){
  for (j in 1:length(names(log_dif_d_prices_ts))) {
    if (i != j) {
      g <- grangertest(log_dif_d_prices_ts[,c(j,i)], order=1)[[4]][2]
      c <- c(names(log_dif_d_prices_ts)[i], names(log_dif_d_prices_ts)[j], round(g, digits = 3))
      a <- rbind(a, c)
    }
  }
}

grangerD1 <- reshape(a, v.names="p", idvar="X", timevar="caused by", direction="wide")
grangerD1 <- grangerD1 %>%
  filter(X != 0)%>%
  select(-p.0)
grangerD1 #p value: row is caused by columns
View(grangerD1)

# export granger test results
stargazer(grangerD1,
          type="html",
          #report="vc*",
          summary=FALSE,
          out="G1_23dBTC.html",
          colnames = TRUE)


# about interpretting grangertest 
# 
# Model 1 is the unrestricted model that includes the Granger-causal terms.
# Model 2 is the restricted model where the Granger-causal terms are omitted.
# The test is a Wald test that assesses whether using the restricted Model 2 in 
# place of Model 1 makes statistical sense (roughly speaking).

# if Pr(>F)<α<α (where αα is your desired level of significance), you reject the null hypothesis of no Granger causality. T
# his indicates that Model 2 is too restrictive as compared with Model 1.
# If the inequality is reversed, you do not reject the null hypothesis as the richer 
# Model 1 is preferred to the restricted Model 2.
# Note: you say we are checking to see if one variable can be used to predict another.
# A more precise statement would be we are checking to see if including xx is useful for
# predicting yy when yy's own history is already being used for prediction. 
# That is, do not miss the fact the xx has to be useful beyond (or extra to) the own history of yy.


# FEVD
# FEVD

# change n.ahead

# lag 1

fevdecomd10 <- round(fevd(Var1d, n.ahead=10), digits = 3)*100
fevdecomd10 <- cbind(names(log_dif_d_prices_ts),fevdecomd10)
View(fevdecomd10)

# lag 2

fevdecomd20 <- round(fevd(Var2D, n.ahead=20), digits = 3)*100
fevdecomd20 <- cbind(names(log_dif_d_prices_ts),fevdecomd20)
View(fevdecomd20)

# why they are the same even if lag is different?
# difference btw n.ahead
(fevd(Var15m,n.ahead=1)-fevd(Var15m,n.ahead=100))*100

# export
stargazer(fevdecomd10,
          type="html",
          report="vc*",
          out="fevddBTC1.html")

stargazer(fevdecomd20,
          type="html",
          report="vc*",
          out="fevddBTC2.html")


dev.off()
plot(fevdecom5m)



fevdecomD <- fevd(Var1d)
fevd(Var1d)

plot(fevdecomD)
str(fevd(Var1d))

fevdecom[[2]]
fevdecom[,1]
names(log_dif_h_prices_ts)

plot(fevd(var.aicD, n.ahead=5))
summary(fevd(var.aicD, n.ahead=2))

