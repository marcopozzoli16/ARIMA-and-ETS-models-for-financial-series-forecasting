#grid costruzione modelli
library(quantmod)
library(dplyr)
library(tidyr)
library(xts)
library(tseries)
library(forecast)
library(zoo)
library(ggplot2)
library(lubridate)
library(gridExtra)
library(parallel)
library(doParallel)

#Ricostruisco logret per valutazione ----
reconstruct_prices_from_fitted_returns <- function(fitted_returns, start_price) {
  if (length(fitted_returns) == 0) return(numeric(0))
  prices <- numeric(length(fitted_returns))
  prev <- start_price
  for (i in seq_along(fitted_returns)) {
    prices[i] <- prev * exp(fitted_returns[i])
    prev <- prices[i]
  }
  prices
}

forecast_prices_from_returns <- function(fc, last_price) {
  mu <- as.numeric(fc$mean)
  last_price * exp(cumsum(mu))
}

mase_denominator <- function(train_y) {
  if (length(train_y) <= 1) return(NA_real_)
  mean(abs(diff(as.numeric(train_y))), na.rm = TRUE)
}

compute_metrics_price <- function(actual, pred, train_y_for_mase=NULL) {
  ok <- !is.na(actual) & !is.na(pred)
  actual2 <- actual[ok]
  pred2 <- pred[ok]
  if (length(actual2) == 0) {
    return(list(RMSE=NA_real_, MAE=NA_real_, MAPE=NA_real_, MASE=NA_real_))
  }
  rmse <- sqrt(mean((actual2 - pred2)^2))
  mae  <- mean(abs(actual2 - pred2))
  mape <- mean(abs((actual2 - pred2)/actual2)) * 100
  denom <- if (!is.null(train_y_for_mase)) mase_denominator(train_y_for_mase) else mase_denominator(actual2)
  mase <- if (is.na(denom) || denom == 0) NA_real_ else (mae / denom)
  list(RMSE=rmse, MAE=mae, MAPE=mape, MASE=mase)
}

acf1_safe <- function(resid) {
  if (length(resid) < 2) return(NA_real_)
  a <- tryCatch(acf(resid, plot = FALSE, lag.max = 1)$acf[2], error = function(e) NA_real_)
  as.numeric(a)
}

#Gridsearch finale ----
gridsearch <- function(serie, fre, oriz, inizio_df)
{
  metriche <- data.frame(
    Serie = character(),
    NOME  = character(),
    RMSE  = numeric(),
    MAPE  = numeric(),
    MASE  = numeric(),
    MAE   = numeric(),
    AIC   = numeric(),
    BIC   = numeric(),
    ACF   = numeric(),
    Frequenza = numeric(),
    Orizzonte = numeric(),
    PreProc = character(),
    TrainTest = character(),
    InizioTrain = character(),
    stringsAsFactors = FALSE
  )
  
  # carica dati
  df <- read.csv(serie, stringsAsFactors = FALSE)
  df$Date <- as.Date(df$Date)
  df <- df[order(df$Date), ]
  
  # train/test
  df_train <- df[1:(nrow(df)-oriz), ]
  df_test  <- df[(nrow(df)-oriz+1):nrow(df), ]
  df_train <- df_train[df_train$Date >= inizio_df, ]
  
  train_ts <- ts(df_train$Close, frequency = fre)
  lr_ts <- log(train_ts) - log(stats::lag(train_ts))
  clean_ts <- tryCatch(tsclean(train_ts), error = function(e) train_ts)
  
  ###Fitto ----
  
  ets_raw     <- tryCatch(ets(train_ts), error=function(e) NULL)
  ets_boxcox  <- tryCatch(ets(train_ts, lambda="auto", biasadj=TRUE), error=function(e) NULL)
  ets_out     <- tryCatch(ets(clean_ts), error=function(e) NULL)
  
  ets_logret <- tryCatch({
    lr <- lr_ts[-1]
    if (all(is.na(lr))) NULL else ets(ts(lr))
  }, error=function(e) NULL)
  
  arima_raw    <- tryCatch(auto.arima(train_ts), error=function(e) NULL)
  arima_boxcox <- tryCatch(auto.arima(train_ts, lambda="auto", biasadj=TRUE), error=function(e) NULL)
  arima_out    <- tryCatch(auto.arima(clean_ts), error=function(e) NULL)
  
  arima_logret <- tryCatch({
    lr <- lr_ts[-1]
    if (all(is.na(lr))) NULL else auto.arima(ts(lr))
  }, error=function(e) NULL)
  
  ###Valuto train ----
  
  val_train <- function(nome, pp, mod) {
    if (is.null(mod)) return()
    fitted_vals <- tryCatch(as.numeric(fitted(mod)), error=function(e) NA)
    m <- compute_metrics_price(as.numeric(train_ts), fitted_vals, as.numeric(train_ts))
    resid <- tryCatch(residuals(mod), error=function(e) NA)
    acf1 <- if (!is.na(resid[1])) acf1_safe(resid) else NA
    metriche <<- rbind(metriche,
                       data.frame(
                         Serie = serie, NOME = nome,
                         RMSE=m$RMSE, MAPE=m$MAPE, MASE=m$MASE, MAE=m$MAE,
                         AIC = mod$aic, BIC = mod$bic, ACF = acf1,
                         Frequenza = fre, Orizzonte = oriz,
                         PreProc = pp, TrainTest="train",
                         InizioTrain = inizio_df,
                         stringsAsFactors = FALSE))
  }
  
  val_train_lr <- function(nome, pp, mod) {
    if (is.null(mod)) return()
    fitted_ret <- tryCatch(as.numeric(fitted(mod)), error=function(e) NA)
    start_price <- df_train$Close[1]
    recon <- reconstruct_prices_from_fitted_returns(fitted_ret, start_price)
    actual <- df_train$Close[-1]
    m <- compute_metrics_price(actual, recon, df_train$Close)
    resid <- actual - recon
    acf1 <- acf1_safe(resid)
    metriche <<- rbind(metriche,
                       data.frame(
                         Serie = serie, NOME = paste0(nome," (lr)"),
                         RMSE=m$RMSE, MAPE=m$MAPE, MASE=m$MASE, MAE=m$MAE,
                         AIC = mod$aic, BIC = NA, ACF = acf1,
                         Frequenza = fre, Orizzonte = oriz,
                         PreProc = pp, TrainTest="train",
                         InizioTrain = inizio_df,
                         stringsAsFactors = FALSE))
  }
  
  #metriche training
  val_train("ETS","Raw",ets_raw)
  val_train("ETS","BoxCox",ets_boxcox)
  val_train("ETS","Outlier",ets_out)
  
  val_train("ARIMA","Raw",arima_raw)
  val_train("ARIMA","BoxCox",arima_boxcox)
  val_train("ARIMA","Outlier",arima_out)
  
  val_train_lr("ETS","Logret",ets_logret)
  val_train_lr("ARIMA","Logret",arima_logret)
  
  
  ###Valuto test ----
  
  val_test <- function(mod, nome, pp) {
    if (is.null(mod)) return()
    prev <- tryCatch(forecast(mod, h=oriz), error=function(e) NULL)
    if (is.null(prev)) return()
    
    if (pp=="Logret") {
      last_price <- tail(df_train$Close,1)
      pred <- forecast_prices_from_returns(prev, last_price)
    } else {
      pred <- as.numeric(prev$mean)
    }
    
    actual <- df_test$Close
    res <- actual - pred
    den <- mase_denominator(df_train$Close)
    m <- list(
      RMSE = sqrt(mean(res^2)),
      MAE  = mean(abs(res)),
      MAPE = mean(abs(res / actual))*100,
      MASE = ifelse(is.na(den)||den==0, NA, mean(abs(res))/den)
    )
    
    metriche <<- rbind(metriche,
                       data.frame(
                         Serie = serie,
                         NOME = ifelse(pp=="Logret", paste0(nome," (lr)"), nome),
                         RMSE=m$RMSE, MAPE=m$MAPE, MASE=m$MASE, MAE=m$MAE,
                         AIC=0, BIC=0, ACF=0,
                         Frequenza = fre, Orizzonte = oriz,
                         PreProc = pp, TrainTest="test",
                         InizioTrain = inizio_df,
                         stringsAsFactors = FALSE))
  }
  
  val_test(ets_raw,"ETS","Raw")
  val_test(ets_boxcox,"ETS","BoxCox")
  val_test(ets_out,"ETS","Outlier")
  
  val_test(arima_raw,"ARIMA","Raw")
  val_test(arima_boxcox,"ARIMA","BoxCox")
  val_test(arima_out,"ARIMA","Outlier")
  
  val_test(ets_logret,"ETS","Logret")
  val_test(arima_logret,"ARIMA","Logret")
  
  return(metriche)
}

#valori per grid
inizio_df_vec <- c("2015-01-01","2020-01-01","2023-01-01","2025-01-01")
step <- c(1,5,20)
o <- c(15, 30, 45, 60, 120)

#File serie
s <- c(
  "AAPL.csv",
  "Bitcoin.csv",
  "Brent.csv",
  "CNYUSD.csv",
  "Copper.csv",
  "EURGBP.csv",
  "EUROSTOXX50.csv",
  "EURUSD.csv",
  "GOLD.csv",
  "NASDAQ100.csv",
  "NVDA.csv",
  "SOLANA.csv",
  "SP500.csv",
  "TNX.csv",
  "TSLA.csv",
  "Wheat.csv",
  "XLE.csv",
  "DAX.csv",
  "FTSE100.csv",
  "NIKKEI225.csv",
  "SHANGHAI.csv",
  "SENSEX.csv",
  "CAC40.csv",
  "TSX.csv",
  "JPYUSD.csv",
  "GBPUSD.csv",
  "AUDUSD.csv",
  "USDCHF.csv",
  "Ethereum.csv",
  "Ripple.csv",
  "Dogecoin.csv",
  "Polkadot.csv",
  "MSFT.csv",
  "AMZN.csv",
  "GOOGL.csv",
  "META.csv",
  "Silver.csv",
  "Corn.csv",
  "Coffee.csv",
  "XLK.csv",
  "XLF.csv",
  "XLY.csv",
  "XLV.csv",
  "Russell2000.csv",
  "MSCI_World.csv",
  "MSCI_EM.csv",
  "HangSeng.csv",
  "KOSPI.csv",
  "ASX200.csv",
  "Bovespa.csv",
  "FTSE_MIB.csv",
  "IBEX35.csv",
  "OMXS30.csv",
  "AEX.csv",
  "BEL20.csv",
  "Singapore_STI.csv",
  "Taiwan_TAIEX.csv",
  "SouthAfrica_Top40.csv",
  "Mexico_IPC.csv",
  "Turkey_BIST100.csv",
  "US5Y.csv",
  "VXD_DowVol.csv",
  "EVZ_EuroVol.csv",
  "MOVE_BondVol.csv",
  "Europe_StoxxBanks.csv",
  "Europe_StoxxHealth.csv",
  "Japan_TopixETF.csv",
  "China_TechETF.csv",
  "Korea_ETF.csv",
  "India_ETF.csv",
  "Brazil_ETF.csv",
  "Canada_ETF.csv",
  "UK_ETF.csv",
  "SP500_Growth.csv",
  "SP500_EqualWeight.csv",
  "Momentum.csv",
  "LowVol.csv",
  "Quality.csv",
  "HighYield.csv",
  "InvestmentGrade.csv",
  "Infra_Global.csv",
  "CleanEnergy.csv",
  "Platinum.csv",
  "Palladium.csv",
  "Sugar.csv",
  "Cotton.csv",
  "Cocoa.csv",
  "Lumber.csv",
  "HeatingOil.csv",
  "Gasoline.csv",
  "REET_Global_REIT.csv",
  "SCHH_US_REIT.csv",
  "IYR_US_REIT.csv",
  "RWX_International_REIT.csv",
  "NZDUSD.csv",
  "USDCAD.csv",
  "USDCNY.csv",
  "USDJPY.csv",
  "USDSEK.csv",
  "USDMXN.csv",
  "USDZAR.csv",
  "EURJPY.csv",
  "GBPJPY.csv",
  "Litecoin.csv",
  "Avalanche.csv",
  "Chainlink.csv",
  "Stellar.csv",
  "Uniswap.csv",
  "Aptos.csv",
  "Maker.csv",
  "Adobe.csv",
  "Netflix.csv",
  "Cisco.csv",
  "Intel.csv",
  "Salesforce.csv",
  "Oracle.csv",
  "AMD.csv",
  "Qualcomm.csv",
  "IBM.csv",
  "Toyota.csv",
  "Samsung.csv",
  "Nestle.csv",
  "Roche.csv",
  "HSBC.csv",
  "BP.csv",
  "Unilever.csv",
  "Siemens.csv",
  "LVMH.csv",
  "Alibaba.csv"
)


task_grid <- expand.grid(
  serie  = s,
  oriz   = o,
  freq   = step,
  inizio = inizio_df_vec,
  stringsAsFactors = FALSE
)

ncl <- detectCores()-1
cl <- makeCluster(ncl)
registerDoParallel(cl)

results_list <- foreach(k = 1:nrow(task_grid),
                        .packages = c("forecast","tseries","quantmod","dplyr","tidyr","xts","zoo","lubridate")) %dopar% {
                          
                          serie     <- task_grid$serie[k]
                          oriz      <- task_grid$oriz[k]
                          freq      <- task_grid$freq[k]
                          inizio_df <- as.Date(task_grid$inizio[k])
                          
                          tryCatch({
                            gridsearch(serie = serie, fre = freq, oriz = oriz, inizio_df = inizio_df)
                          }, error=function(e) NULL)
                        }

stopCluster(cl)

m_grid <- do.call(rbind, results_list)
m_grid$Serie <- sub("\\.csv$","",m_grid$Serie)

m_grid[m_grid$NOME == "ETS (lr)",]$NOME   <- "ETS"
m_grid[m_grid$NOME == "ARIMA (lr)",]$NOME <- "ARIMA"
m_grid$InizioTrain <- as.character(m_grid$InizioTrain)
m_grid[m_grid$InizioTrain == "2015-01-01",]$InizioTrain <- "2015"
m_grid[m_grid$InizioTrain == "2020-01-01",]$InizioTrain <- "2020"
m_grid[m_grid$InizioTrain == "2023-01-01",]$InizioTrain <- "2023"
m_grid[m_grid$InizioTrain == "2025-01-01",]$InizioTrain <- "2025"

write.csv(m_grid,
          "~/Università/Tesi/Tesi2/MetricheMega.csv",
          row.names = FALSE)
