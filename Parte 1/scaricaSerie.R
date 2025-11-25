library(quantmod)
library(dplyr)
library(tidyr)
library(xts)
library(tseries)
library(forecast)
library(zoo)

scarica_serie <- function(nome, ticker, da = "2015-01-01", a = Sys.Date()) {
  cat("Scaricamento di:", nome, "(", ticker, ")\n")
  
  #serie da Yahoo Finance
  tryCatch({
    dati_xts <- getSymbols(ticker, src = "yahoo", from = da, to = a, auto.assign = FALSE)
    dati_df <- data.frame(Date = index(dati_xts), coredata(Cl(dati_xts)))
    colnames(dati_df) <- c("Date", "Close")
    
    #valori mancanti
    dati_df <- dati_df %>%
      arrange(Date) %>%
      mutate(Close = zoo::na.locf(zoo::na.approx(Close, na.rm = FALSE), na.rm = FALSE)) %>%
      drop_na()
    
    
    #salva su CSV
    write.csv(dati_df, paste0(nome, ".csv"), row.names = FALSE)
    cat("Serie", nome, "salvata con", nrow(dati_df), "osservazioni.\n\n")
    
    return(dati_df)
  }, error = function(e) {
    cat("Errore per", nome, ":", e$message, "\n\n")
    return(NULL)
  })
}

#serie da scaricare
serie <- list(
  "AAPL" = "AAPL",
  "Bitcoin" = "BTC-USD",
  "Brent" = "BZ=F",
  "CNYUSD" = "CNY=X",
  "Copper" = "HG=F",
  "EURGBP" = "EURGBP=X",
  "EUROSTOXX50" = "^STOXX50E",
  "EURUSD" = "EURUSD=X",
  "GOLD" = "GC=F",
  "NASDAQ100" = "^NDX",
  "NVDA" = "NVDA",
  "SOLANA" = "SOL-USD",
  "SP500" = "^GSPC",
  "TNX" = "^TNX",
  "TSLA" = "TSLA",
  "Wheat" = "ZW=F",
  "XLE" = "XLE",
  "DAX" = "^GDAXI",
  "FTSE100" = "^FTSE",
  "NIKKEI225" = "^N225",
  "SHANGHAI" = "000001.SS",
  "SENSEX" = "^BSESN",
  "CAC40" = "^FCHI",
  "TSX" = "^GSPTSE",
  "JPYUSD" = "JPY=X",
  "GBPUSD" = "GBPUSD=X",
  "AUDUSD" = "AUDUSD=X",
  "USDCHF" = "CHF=X",
  "Ethereum" = "ETH-USD",
  "Ripple" = "XRP-USD",
  "Dogecoin" = "DOGE-USD",
  "Polkadot" = "DOT-USD",
  "MSFT" = "MSFT",
  "AMZN" = "AMZN",
  "GOOGL" = "GOOGL",
  "META" = "META",
  "Silver" = "SI=F",
  "Corn" = "ZC=F",
  "Coffee" = "KC=F",
  "XLK" = "XLK",
  "XLF" = "XLF",
  "XLY" = "XLY",
  "XLV" = "XLV",
  "Russell2000" = "^RUT",
  "MSCI_World" = "URTH",
  "MSCI_EM" = "EEM",
  "HangSeng" = "^HSI",
  "KOSPI" = "^KS11",
  "ASX200" = "^AXJO",
  "Bovespa" = "^BVSP",
  "FTSE_MIB" = "FTSEMIB.MI",
  "IBEX35" = "^IBEX",
  "S&P500_Vol" = "^VIX",
  "NZDUSD" = "NZDUSD=X",
  "CADUSD" = "CAD=X",
  "TRYUSD" = "TRY=X",
  "MXNUSD" = "MXN=X",
  "BRLUSD" = "BRL=X",
  "ZARUSD" = "ZAR=X",
  "SEKUSD" = "SEK=X",
  "NOKUSD" = "NOK=X",
  "CHFJPY" = "CHFJPY=X",
  "EURJPY" = "EURJPY=X",
  "Litecoin" = "LTC-USD",
  "Avalanche" = "AVAX-USD",
  "Chainlink" = "LINK-USD",
  "Stellar" = "XLM-USD",
  "Polygon" = "MATIC-USD",
  "Tezos" = "XTZ-USD",
  "TRON" = "TRX-USD",
  "Monero" = "XMR-USD",
  "JPM" = "JPM",
  "BAC" = "BAC",
  "WMT" = "WMT",
  "DIS" = "DIS",
  "NFLX" = "NFLX",
  "INTC" = "INTC",
  "ADBE" = "ADBE",
  "PYPL" = "PYPL",
  "CRM" = "CRM",
  "ORCL" = "ORCL",
  "Platinum" = "PL=F",
  "Palladium" = "PA=F",
  "Soybeans" = "ZS=F",
  "Sugar" = "SB=F",
  "Cotton" = "CT=F",
  "Rice" = "ZR=F",
  "Oats" = "ZO=F",
  "Cocoa" = "CC=F",
  "Ethanol" = "EH=F",
  "XLRE" = "XLRE",
  "XLI" = "XLI",
  "XLU" = "XLU",
  "XLB" = "XLB",
  "XLC" = "XLC",
  "SPY" = "SPY",
  "QQQ" = "QQQ",
  "ShanghaiComp" = "000001.SS",
  "Nifty50" = "^NSEI",
  "SMI" = "^SSMI",
  "OMXS30" = "^OMXS30",
  "AEX" = "^AEX",
  "BEL20" = "^BFX",
  "Singapore_STI" = "^STI",
  "Taiwan_TAIEX" = "^TWII",
  "SouthAfrica_Top40" = "TOP40.JO",
  "Mexico_IPC" = "^MXX",
  "Turkey_BIST100" = "^XU100",
  "US5Y" = "^FVX",
  "VXD_DowVol" = "^VXD",
  "EVZ_EuroVol" = "^V2X",
  "MOVE_BondVol" = "^MOVE",
  "Europe_StoxxBanks" = "SX7E.DE",
  "Europe_StoxxHealth" = "SXDP.DE",
  "Japan_TopixETF" = "1306.T",
  "China_TechETF" = "513050.SS",
  "Korea_ETF" = "091170.KQ",
  "India_ETF" = "INDA",
  "Brazil_ETF" = "EWZ",
  "Canada_ETF" = "EWC",
  "UK_ETF" = "EWU",
  "SP500_Growth" = "SPYG",
  "SP500_EqualWeight" = "RSP",
  "Momentum" = "MTUM",
  "LowVol" = "USMV",
  "Quality" = "QUAL",
  "HighYield" = "HYG",
  "InvestmentGrade" = "LQD",
  "Infra_Global" = "TOLZ",
  "CleanEnergy" = "ICLN",
  "Lumber" = "LBS",
  "HeatingOil" = "HO=F",
  "Gasoline" = "RB=F",
  "REET_Global_REIT" = "REET",
  "SCHH_US_REIT" = "SCHH",
  "IYR_US_REIT" = "IYR",
  "RWX_International_REIT" = "RWX",
  "USDCAD" = "USDCAD=X",
  "USDCNY" = "USDCNY=X",
  "USDJPY" = "JPY=X",
  "USDSEK" = "USDSEK=X",
  "USDMXN" = "USDMXN=X",
  "USDZAR" = "USDZAR=X",
  "GBPJPY" = "GBPJPY=X",
  "Uniswap" = "UNI-USD",
  "Aptos" = "APT-USD",
  "Maker" = "MKR-USD",
  "Adobe" = "ADBE",
  "Netflix" = "NFLX",
  "Cisco" = "CSCO",
  "Intel" = "INTC",
  "Salesforce" = "CRM",
  "Oracle" = "ORCL",
  "AMD" = "AMD",
  "Qualcomm" = "QCOM",
  "IBM" = "IBM",
  "Toyota" = "7203.T",
  "Samsung" = "005930.KS",
  "Nestle" = "NESN.SW",
  "Roche" = "ROG.SW",
  "HSBC" = "HSBA.L",
  "BP" = "BP.L",
  "Unilever" = "ULVR.L",
  "Siemens" = "SIE.DE",
  "LVMH" = "MC.PA",
  "Alibaba" = "BABA"
)



#scarica e salva
dati_lista <- lapply(names(serie), function(nome) {
  ticker <- serie[[nome]]
  scarica_serie(nome, ticker)
})
names(dati_lista) <- names(serie)

