# load libraries

library(fastSOM)
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
library(vars)
library(frequencyConnectedness)
library(BigVAR)
library(fBasics)
library(parallel)
library(readr)
set.seed(1)
library(scales)
library(directlabels)

#################### 
# preprocess data
####################
return_file = 'daily_data_USD.csv'
data <- read.csv(return_file)
head(data)
# delete index
data <- data[,-1]

# change names
names(data) <- sub("_.*", "", names(data))
names(data)

# check for NA
sum(is.na(data))

dim(data)
head(data)

#be aware that symbols are different across sources. For example, Bytom: coinmarketcap.com - BTM, cryptocompare - BTM*, r convert it as 'BTM.'
coins_for_analysis = c("BTC","ETH", "XRP", "BCH", "EOS", "LTC", "XLM", "IOT","XMR","DASH",
                     "XEM", "VEN","ETC", "XVG","OMG",   "LSK",  "ZEC","BCN","BTS",
                     "SC", "STRAT",  "DOGE","VERI","NEO","BTM.","STEEM","DGD","DCR",
                     "WAVES","AE","SNT","KMD","REP","ARDR","ARK","DGB","PIVX","FCT",
                     "MONA","BNB","ZRX","LRC","HSR","WTC","SUB","GAS","FUN","SYS","KNC","GBYTE","SKY","RDD",
                     "NXT", "XZC", "MAID",  "BNT", "STORJ",  "PAY","CND", "EMC","ZEN","NXS",
                     "ICN", "CVC", "VTC",   "MANA", "PART",   "BOS","GAME",  "UBQ",  "GNO",
                     "BLOCK",   "MTL",  "MCO", "SAN", "POE", "BTX","ANT","PLR", "NAV",
                     "BTCD", "ION","ADX","PPC",   "XBY","RLC","XDN","QRL", "EDG", "SLS")
length(coins_for_analysis)

coins_90 = coins_for_analysis
coins_90
coins_40 = coins_90[c(1:40)]
coins_40
coins_23 = coins_90[c(1:23)]
coins_23

#################### 
# start analysis
####################

# define coins, for example, 40
coins = coins_40
coins

start_analysis_return = Sys.time()

#set parameters
window = 80
coins_number = length(coins)
type = sprintf("daily_return_%1$s_rw%2$s", coins_number, window)
n.ahead = 10


#make new data frame
df = subset(data, select=c('timestamp',coins))
summary(df)
str(df)

#convert to time series
df_ts <- xts(df[,-1], order.by = as.POSIXct(df$timestamp, format = "%m/%d/%y %H:%M"))
#"%m/%d/%y %H:%M"  "%Y-%m-%d %H:%M:%S"
head(df_ts)
str(df_ts)

#delete duplicates
df_ts <- make.index.unique(df_ts, drop = TRUE)
periodicity(df_ts)

#convert to return
for (i in 1:length(names(df_ts))){
  df_ts[,names(df_ts)[i]]<-diff.xts(log(df_ts[,names(df_ts)[i]]))
}

#View(df_ts)
tail(df_ts)
str(df_ts)
sum(is.na(df_ts))

#remove NA
df_ts <-  na.omit(df_ts) 
str(df_ts) #also removed some periods

dim(df_ts) 

sum(is.na(df_ts))

head(df_ts)
periodicity(df_ts)

data_for_analysis <- df_ts[-1,]

head(data_for_analysis)
tail(data_for_analysis)
#View(data_for_analysis)
#str(data_sep_march)
sum(is.na(data_for_analysis))


sum(is.nan(data_for_analysis))

periodicity(data_for_analysis)

# plot data
plot.xts(data_for_analysis,legend.loc = "topright")
plot.xts(data_for_analysis[,c(11:20)],legend.loc = "topright")
plot.xts(data_for_analysis[,c(21:40)],legend.loc = "topright")
plot.xts(data_for_analysis[,c(1:10)],legend.loc = "topright", main = 'Daily return for selected coins')
#save figure
dev.copy(png,sprintf("%1$s_plot_10coins.png", type))
dev.off()


#mean and sd
for (i in 1:length(names(data_for_analysis))){
  print(c(names(data_for_analysis)[i], round(mean(data_for_analysis[,i]), digits = 3), round(sd(data_for_analysis[,i]), digits = 3)))
}

#descriptive statistics
basicStats(data_for_analysis[,1])
descriptive = basicStats(data_for_analysis)

#export descriptive statistics 
stargazer(descriptive,
          type="html",
          report="vc*",
          summary = FALSE,
          out=sprintf("%1$s_descriptive.html", type))

#define best lag
VARselect(data_for_analysis, lag.max = 4, type = "both")
#best lag is 1
lag_var = 1

#to ease later on + remove duplicates
data_VAR = make.index.unique(data_for_analysis, drop=TRUE, fromLast=FALSE)
periodicity(data_VAR)
dim(data_VAR)
#244  40

####################################################################################
# VAR-LASSO
####################################################################################

# Perform the estimation
start = Sys.time()
VAR_lasso <- big_var_est(data_VAR)
VAR_lasso

end = Sys.time()
time_lasso = end-start
time_lasso

# Coefficient matrix at end of evaluation period
VAR_lasso@betaPred
#sum of non-zero coef
number_non0coef = sum(VAR_lasso@betaPred != 0)
number_non0coef
# Residuals at end of evaluation period 
VAR_lasso@resids
# Lagged Values at end of evaluation period 
VAR_lasso@Zvals

#plot labmda and save image
plot(VAR_lasso)
SparsityPlot.BigVAR.results(VAR_lasso)
dev.copy(png,sprintf("%1$s_sparsity.png", type))
dev.off()

spillovers_lasso <- spilloverDY12(VAR_lasso, n.ahead = n.ahead, no.corr = F)
spillovers_lasso

#check if all are the same
total_connectendess = sum(to(spillovers_lasso)[[1]])
total_connectendess
sum(from(spillovers_lasso)[[1]])
overall(spillovers_lasso)[[1]]

to = unlist(to(spillovers_lasso))
from = unlist(from(spillovers_lasso))
net = unlist(net(spillovers_lasso))
to_from_net = as.data.frame(cbind(to, from, net))
to_from_net['Label'] = rownames(to_from_net)
to_from_net['Id'] = c(1:length(rownames(to_from_net)))
to_from_net = to_from_net[,c('Id', 'Label', 'to', 'from','net')]
to_from_net
#save to later analyze in gephi
write_csv(to_from_net, sprintf("%1$s_to_from_net_lasso.csv", type))


spillovers_table_lasso <- spillovers_to_dataframe(spillovers_lasso)
#View(spillovers_table_lasso)

#export to html
stargazer(spillovers_table_lasso,
          digits = 2,
          type="html",
          summary = FALSE,
          out=sprintf("%1$s_spillovers_table_lasso.html", type))

write_csv(spillovers_table_lasso, sprintf("%1$s_spillovers_table_lasso.csv", type))

####################################
# rolling window estimation
####################################

#get optimal lambda
optimal_lambda = VAR_lasso@OptimalLambda
optimal_lambda

dim(data_VAR)

start = Sys.time()
sp_lasso_rolling_window <- spilloverRollingDY12(data_VAR, n.ahead = n.ahead, no.corr = F,
                                                func_est = "big_var_est_fixed_lambda", 
                                                params_est = list(),window = window)
end=Sys.time()
time_lasso = end-start
time_lasso

#save object on disk
save(sp_lasso_rolling_window, file = sprintf("%1$s_sp_lasso_rolling_window.RData", type))

#plot and save image
plotOverall(sp_lasso_rolling_window)

# or better plot 

break_point <- sum(is.nan(overall(sp_lasso_rolling_window)[[1]]))+1 #sometimes it can produce NA values because of lack of data
total_points <- length(overall(sp_lasso_rolling_window)[[1]])
dynamic_connectedness <- overall(sp_lasso_rolling_window)[[1]][c(break_point:total_points)]
head(dynamic_connectedness)
names(dynamic_connectedness) <- c('Total connectedness')

plot_dynamic_connectedness <- autoplot.zoo(dynamic_connectedness) +
  geom_line(color='red')+
  labs(x='Time', y='Total connectedness')+
  theme_minimal()+
  #ylim(70, 97)+
  # scale_x_date(breaks = date_breaks("months"),labels = time_format("%b"))
  scale_x_datetime(breaks = date_breaks("2 week"), labels = date_format("%b %d"))+
  ggtitle(sprintf("Total Connectedness %1$s", type))

plot_dynamic_connectedness
# and save
ggsave(sprintf("%1$s_dynamic_connectedness.png", type), plot_dynamic_connectedness)

write.csv(data.frame(dynamic_connectedness), sprintf("%1$s_dynamic_total_connectedness.csv", type))


#plot rolling window
plotFrom(sp_lasso_rolling_window, which = c(1:22))
dev.copy(png,sprintf("%1$s_rw_from.png", type))
dev.off()

# if more coins
plotFrom(sp_lasso_rolling_window, which = c(21:40))
dev.copy(png,sprintf("%1$s_rw_from_part2.png", type))
dev.off()

#plot rolling window
plotTo(sp_lasso_rolling_window, which = c(1:22))
dev.copy(png,sprintf("%1$s_rw_to.png", type))
dev.off()
# if more coins
plotTo(sp_lasso_rolling_window, which = c(21:40))
dev.copy(png,sprintf("%1$s_rw_to_part2.png", type))
dev.off()

plotNet(sp_lasso_rolling_window,which = 1:22)
dev.copy(png,sprintf("%1$s_rw_net.png", type))
dev.off()
# if more coins
plotNet(sp_lasso_rolling_window,which = c(21:40))
dev.copy(png,sprintf("%1$s_rw_net_part2.png", type))
dev.off()

coins_for_graphs = c('BTC', 'ETH', 'BCH', 'IOT', 'XMR', 'XRP')

# to dynamic
to_dynamic = to(sp_lasso_rolling_window)[[1]][c(break_point:total_points)]
plot_to_dynamic = autoplot.zoo(to_dynamic[,coins_for_graphs], facets = NULL) +
  #geom_line()+
  labs(x='Time', y='To others connectedness')+
  theme_minimal()+
  # ylim(0, 4)+
  # scale_x_date(breaks = date_breaks("months"),labels = time_format("%b"))
  scale_x_datetime(breaks = date_breaks("2 week"), labels = date_format("%b %d"))+
  ggtitle(sprintf("To others connectedness %1$s", type))
plot_to_dynamic
ggsave(sprintf("%1$s_dynamic_to.png", type), plot_to_dynamic)



# from dynamic
from_dynamic = from(sp_lasso_rolling_window)[[1]][c(break_point:total_points)]
plot_from_dynamic = autoplot.zoo(from_dynamic[,coins_for_graphs], facets = NULL) +
  #geom_line()+
  labs(x='Time', y='From others connectedness')+
  theme_minimal()+
  # ylim(0, 4)+
  # scale_x_date(breaks = date_breaks("months"),labels = time_format("%b"))
  scale_x_datetime(breaks = date_breaks("2 week"), labels = date_format("%b %d"))+
  ggtitle(sprintf("From others connectedness %1$s", type))
plot_from_dynamic
ggsave(sprintf("%1$s_dynamic_from.png", type), plot_from_dynamic)

# net dynamic
net_dynamic = net(sp_lasso_rolling_window)[[1]][c(break_point:total_points)]
plot_net_dynamic = autoplot.zoo(net_dynamic[,coins_for_graphs], facets = NULL) +
  #geom_line()+
  labs(x='Time', y='Net connectedness')+
  theme_minimal()+
  # ylim(0, 4)+
  # scale_x_date(breaks = date_breaks("months"),labels = time_format("%b"))
  scale_x_datetime(breaks = date_breaks("2 week"), labels = date_format("%b %d"))+
  ggtitle(sprintf("Net connectedness %1$s", type))
plot_net_dynamic
ggsave(sprintf("%1$s_dynamic_net.png", type), plot_net_dynamic)

#plotPairwise(rolling_windows_spillovers)
end_analysis_return = Sys.time()
time_analysis_return = end_analysis_return - start_analysis_return
time_analysis_return

print(total_connectendess)
