# ES and M-C VaR
# by Changjie Ma and Yijun Lou
# Feb. 2019

library(quantmod)
library(PerformanceAnalytics)
library(xts)
library(ggplot2)
library(fOptions)
library(mnormt)

# get data from Yahoo
adj.close <- 6  # 6th field is adjusted close
equity.tickers <- c("^GSPC","^DJI","^IXIC","^RUT")
prices <- getSymbols(equity.tickers[1],from = "2004-01-01", to = "2018-12-31", source="yahoo", auto.assign=FALSE,
                     return.class="xts")[,adj.close]
for (i.tmp in 2:length(equity.tickers)) {
  prices.tmp <- getSymbols(equity.tickers[i.tmp],from = "2004-01-01", to = "2018-12-31", source="yahoo",
                           auto.assign=FALSE, return.class="xts")[,adj.close]
  prices <- cbind(prices, prices.tmp)
}
equity.names <- c("S&P 500","Dow 30","Nasdaq","Russell 2000")
colnames(prices) <- equity.names
returns <- diff(prices)/prices # calculate simple returns
returns.portfolio <- rowMeans(returns)
returns.portfolio <- xts(returns.portfolio,order.by = index(returns))
colnames(returns.portfolio) <- 'returns'

# Calculate the ES using historical simulation

VaR.q1 <- c()
ES.q1 <- c()
start <- length(index(returns.portfolio['/2007-12-31'])) # the number of row for the start day in the 2-year horizon
for (i.tmp in index(returns.portfolio['2008/2009'])){
  returns.portfolio.data <- returns.portfolio[(start-1000):(start-1)]
  return.tailloss <- sort(as.vector(returns.portfolio.data),decreasing = FALSE)
  # VaR.tmp <- VaR(returns.portfolio.data,p=0.95,method='historical')
  VaR.tmp <- return.tailloss[50]
  ES.tmp <- mean(return.tailloss[1:50])
  start <- start+1
  VaR.q1 <- c(VaR.q1,VaR.tmp)
  ES.q1 <- c(ES.q1,ES.tmp)
}
results.q1 <- xts(cbind(VaR.q1,ES.q1),order.by = index(returns['2008/2009']))

ggplot( results.q1*4, aes(index(results.q1)) ) + 
  geom_line( aes( y = VaR.q1 , col='VaR')) +
  geom_line( aes( y = ES.q1 , col='ES' ) ) +
  labs( title = "VaR and ES")+
  xlab('Time')+
  ylab('value (million USD)')

# Compute implied volatilities through Black-Scholes model
# p <- c(100*(18.70+19.70)/2,100*(20.00+21.10)/2,100*(1.80+1.96)/2,100*(2.14+2.37)/2) # observed option price
p <- c((18.70+19.70)/2,(20.00+21.10)/2,(1.80+1.96)/2,(2.14+2.37)/2) # observed option price
S <- c(2564.98,2564.98,23273.96/100,23273.96/100) # current prices
X <- c(2565,2565,233,233) # strike price

r <- 0.0124236 # interest rate
q <- c(0.0164791,0.0164791,0.0236134,0.0236134) # dividands

T <- length(as.Date('2017-10-24'):as.Date('2017-11-16'))/360

vol.1 <- GBSVolatility(p[1],'c',S[1],X[1],T,r,(r-q[1]))
vol.2 <- GBSVolatility(p[2],'p',S[2],X[2],T,r,(r-q[2]))
vol.3 <- GBSVolatility(p[3],'c',S[3],X[3],T,r,(r-q[3]))
vol.4 <- GBSVolatility(p[4],'p',S[4],X[4],T,r,(r-q[4]))

implied.volatilities <- c(vol.1,vol.2,vol.3,vol.4)

# Simulate and compute VaR and ES
BlackScholesCall <- function(S0,K,T,r,d,sigma){
  t <- 0
  d1 <- (log(S0/K)+(r-d+sigma^2/2)*(T-t))/(sigma*sqrt(T-t))
  d2 <- d1-sigma*sqrt(T-t)
  S0*exp(-d*(T-t))*pnorm(d1)-K*exp(-r*(T-t))*pnorm(d2)
}
BlackScholesPut <- function(S0,K,T,r,d,sigma){
  t <- 0
  d1 <- (log(S0/K)+(r-d+sigma^2/2)*(T-t))/(sigma*sqrt(T-t))
  d2 <- d1-sigma*sqrt(T-t)
  K*exp(-r*(T-t))*pnorm(-d2)-S0*exp(-d*(T-t))*pnorm(-d1)
}

n <- 10000 # number of simulation

corr.matrix <- matrix(c(1.00,0.97,-0.80,-0.75,
                        0.97,1.00,-0.75,-0.80,
                        -0.80,-0.75,1.00,0.90,
                        -0.75,-0.80,0.90,1.00),4,4)

sd.vector <- c(0.0105,0.0110,0.1250,0.1150)
mean.variables <- c(0.0001,0.0001,0.0078,0.0066)
m.1 <- matrix(rep(sd.vector,4),4,4)
m.2 <- matrix(rbind(sd.vector,sd.vector,sd.vector,sd.vector),4,4)
sim <- data.frame(rmnorm(n,mean.variables,corr.matrix*m.1*m.2))
colnames(sim) <- c('SPX.log.change','DJX.log.change','SPX.volatility','DJX.volatility')

sim$SPX.volatility[sim$SPX.volatility<0] <- 0
sim$DJX.volatility[sim$DJX.volatility<0] <- 0

price.option.1 <- c()
price.option.2 <- c()
price.option.3 <- c()
price.option.4 <- c()

for (i in 1:n){
  price.option.1 <- c(price.option.1, BlackScholesCall(S[1]*exp(sim[i,1]),X[1],T,r,q[1],sim[i,3]))
  price.option.2 <- c(price.option.2, BlackScholesPut(S[2]*exp(sim[i,1]),X[2],T,r,q[2],sim[i,3]))
  price.option.3 <- c(price.option.3, BlackScholesCall(S[3]*exp(sim[i,2]),X[3],T,r,q[3],sim[i,4]))
  price.option.4 <- c(price.option.4, BlackScholesPut(S[4]*exp(sim[i,2]),X[4],T,r,q[4],sim[i,4]))
}

price.option.all <- data.frame(cbind((price.option.1),
                                     (price.option.2),
                                     (price.option.3),
                                     (price.option.4)))

weight <- matrix(c(-50,-50,600,600),4,1)
price.combined <- data.frame(as.matrix(price.option.all)%*%weight)
price.initial <- c(matrix(p,1,4)%*%weight)

return.combined <- price.combined/price.initial-1

colnames(return.combined) <- 'return'
ggplot(return.combined,aes(index(return.combined)))+
  geom_line(aes(y = return,col='return'))
return.combined.sorted <- sort(return.combined[,1],decreasing=FALSE)
VaR.port.sim <- return.combined.sorted[n*0.05]

ES.port.sim <- mean(return.combined.sorted[1:n*0.05])

# Output
cat('\n\nQ.2',
    '\nFair value of the four options are:',round(p,2),
    '\nImplied volatility of the four options are:',round(implied.volatilities,4),
    '\n\nQ.3',
    '\n%5 VaR of the portfolio return is:',round(VaR.port.sim*100,2),'%',
    '\n%5 VaR of the portfolio value is:',round(VaR.port.sim*price.initial,2),
    '\n\nQ.4',
    '\nExpected shprtfall of the portfolio return is:',round(ES.port.sim*100,2),'%',
    '\nExpected shprtfall of the portfolio value is:',round(ES.port.sim*price.initial,2))
