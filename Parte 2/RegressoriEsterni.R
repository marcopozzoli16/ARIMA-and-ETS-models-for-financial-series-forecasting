library(forecast)
library(dplyr)
library(lubridate)
library(parallel)
library(ggplot2)
library(crayon)
library(dunn.test)
library(FSA)

prepare_data <- function(s, r, dist, h) {
  # Leggo serie target
  df <- read.csv(s)
  df$DATE <- as.Date(df$DATE)
  df <- df[1:(nrow(df) - dist), ]
  colnames(df) <- c("Date","Y")
  
  # Aggiungo regressori, solo se esistono
  for (rn in r) {
    t <- read.csv(rn)
    t$DATE <- as.Date(t$DATE)
    t <- t[1:(nrow(t) - dist), ]
    cname <- substring(rn, 1, nchar(rn) - 4)
    colnames(t) <- c("Date", cname)
    # Unisco solo sulle date comuni
    df <- merge(df, t, by = "Date", all = FALSE)
  }
  
  # Se non ci sono regressori, creo frame coerente
  if (length(r) == 0) {
    df <- df %>% select(Date, Y)
  }
  
  # Applico lag di h ai regressori (non alla Y)
  if (ncol(df) > 2) {
    for (j in 3:ncol(df)) {
      df[, j] <- dplyr::lag(df[, j], h)
    }
  }
  
  # Rimuovo righe incomplete
  df <- na.omit(df)
  
  # Divido in train e test (test = ultime h osservazioni)
  n <- nrow(df)
  if (n <= h) stop("Troppo poche osservazioni dopo il lag e dist.")
  dftr <- df[1:(n - h), ]
  dfts <- df[(n - h + 1):n, ]
  
  list(dftr = dftr, dfts = dfts)
}


#fit di tutti i modelli per un dist e un orizzonte h ----

fit_models_one_step <- function(s, r, f, dist, h) {
  data_list <- prepare_data(s, r, dist, h)
  dftr <- data_list$dftr
  dfts <- data_list$dfts
  
  target <- ts(dftr$Y, frequency = f)
  xtr <- if (ncol(dftr) > 2) as.matrix(dftr %>% select(-Y, -Date)) else NULL
  xts <- if (ncol(dfts) > 2) as.matrix(dfts %>% select(-Y, -Date)) else NULL
  y_test <- dfts$Y
  
  results <- list()
  
  # Utility per evitare errori di lunghezza forecast/test
  align_lengths <- function(pred, actual) {
    n <- min(length(pred), length(actual))
    list(pred = pred[1:n], actual = actual[1:n])
  }
  

  ###  ARIMAX ----

  if (!is.null(xtr)) {
    mod_arimax <- auto.arima(target, xreg = xtr, seasonal = TRUE, stepwise = FALSE)
    train_acc <- forecast::accuracy(mod_arimax)[1, "MAPE"]
    fc <- forecast(mod_arimax, xreg = xts, h = h)
    aligned <- align_lengths(fc$mean, y_test)
    test_acc <- forecast::accuracy(aligned$pred, aligned$actual)[1, "MAPE"]
    
    last_pred <- tail(aligned$pred, 1)
    last_real <- tail(aligned$actual, 1)
    err_perc <- abs(last_pred - last_real) / last_real * 100
    
    results[[length(results) + 1]] <- data.frame(
      Modello = "ARIMAX",
      Orizzonte = h,
      MAPE_train = train_acc,
      MAPE_test = test_acc,
      dist = dist,
      ultimaossprevista = last_pred,
      actualultimaoss = last_real,
      errore_ultima_perc = err_perc
    )
  }
  

  ### ARIMA base ----

  mod_arima <- auto.arima(target, seasonal = TRUE, stepwise = FALSE)
  train_acc <- forecast::accuracy(mod_arima)[1, "MAPE"]
  fc <- forecast(mod_arima, h = h)
  aligned <- align_lengths(fc$mean, y_test)
  test_acc <- forecast::accuracy(aligned$pred, aligned$actual)[1, "MAPE"]
  
  last_pred <- tail(aligned$pred, 1)
  last_real <- tail(aligned$actual, 1)
  err_perc <- abs(last_pred - last_real) / last_real * 100
  
  results[[length(results) + 1]] <- data.frame(
    Modello = "ARIMA",
    Orizzonte = h,
    MAPE_train = train_acc,
    MAPE_test = test_acc,
    dist = dist,
    ultimaossprevista = last_pred,
    actualultimaoss = last_real,
    errore_ultima_perc = err_perc
  )
  

  ### ETSX (tslm) ----

  
    form <- as.formula(
      paste("target ~ trend +", paste(colnames(dftr)[3:(2 + ncol(xtr))], collapse = " + "))
    )
    mod_etsx <- tslm(form, data = dftr)
    train_acc <- forecast::accuracy(mod_etsx)[1, "MAPE"]
    
    dfts$target <- NA
    fc <- forecast(mod_etsx, newdata = dfts)
    aligned <- align_lengths(fc$mean, y_test)
    test_acc <- forecast::accuracy(aligned$pred, aligned$actual)[1, "MAPE"]
    
    last_pred <- tail(aligned$pred, 1)
    last_real <- tail(aligned$actual, 1)
    err_perc <- abs(last_pred - last_real) / last_real * 100
    
    results[[length(results) + 1]] <- data.frame(
      Modello = "ETSX",
      Orizzonte = h,
      MAPE_train = train_acc,
      MAPE_test = test_acc,
      dist = dist,
      ultimaossprevista = last_pred,
      actualultimaoss = last_real,
      errore_ultima_perc = err_perc
    )
  

  ### ETS base ----

  mod_ets <- ets(target)
  train_acc <- forecast::accuracy(mod_ets)[1, "MAPE"]
  fc <- forecast(mod_ets, h = h)
  aligned <- align_lengths(fc$mean, y_test)
  test_acc <- forecast::accuracy(aligned$pred, aligned$actual)[1, "MAPE"]
  
  last_pred <- tail(aligned$pred, 1)
  last_real <- tail(aligned$actual, 1)
  err_perc <- abs(last_pred - last_real) / last_real * 100
  
  results[[length(results) + 1]] <- data.frame(
    Modello = "ETS",
    Orizzonte = h,
    MAPE_train = train_acc,
    MAPE_test = test_acc,
    dist = dist,
    ultimaossprevista = last_pred,
    actualultimaoss = last_real,
    errore_ultima_perc = err_perc
  )
  
  bind_rows(results)
}


#Ciclo completo----

run_all_models <- function(s, r, f, dist_max, horizons = c(1, 5, 15, 30, 60, 120)) {
  final_results <- list()
  
  for (dist in seq(1, (dist_max * 9), by = 9)) {
    cat("DIST =", dist, "\n")
    for (h in horizons) {
      cat("  Orizzonte =", h, "\n")
      res <- fit_models_one_step(s, r, f, dist, h)
      final_results[[length(final_results) + 1]] <- res
    }
  }
  
  bind_rows(final_results)
}


# # Lista delle serie e regressori (da personalizzare)
# series_list <- list(
#   list(serie = "SP500.csv", 
#        regressori = c("VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv", "NASDAQ100.csv", 
#                       "GOLD.csv", "RUSSELL2000.csv", "HIGH_YIELD_BOND.csv", 
#                       "WTI_OIL.csv", "MSCI_EM.csv", "IRX.csv")),
#   
#   list(serie = "APPLE.csv", 
#        regressori = c("NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv", 
#                       "DOLLAR_INDEX.csv", "MICROSOFT.csv", "AMAZON.csv", 
#                       "TNX.csv", "CONSUMER_DISCR.csv", "VIX.csv", "TSMC.csv")),
#   
#   list(serie = "ETHEREUM.csv", 
#        regressori = c("BITCOIN.csv", "TETHER.csv", "BNB.csv", "GOLD.csv", 
#                       "NASDAQ100.csv", "DOLLAR_INDEX.csv", "VIX.csv", 
#                       "CARDANO.csv", "SOLANA.csv", "SP500.csv")),
#   
#   list(serie = "BRENT_OIL.csv", 
#        regressori = c("WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv", 
#                       "ENERGY_SECTOR.csv", "SP500.csv", "MSCI_EM.csv", 
#                       "TRANSPORTATION_ETF.csv", "TNX.csv", "COPPER.csv", 
#                       "EUROSTOXX50.csv"))
# )


# 4)Creo risultati -----
  ### SP500 ----
results_SP_1reg <- run_all_models(
  s = "SP500.csv",
  r = c("VIX.csv"),
  f = 1,
  dist_max = 25
)
results_SP_1reg$Serie <- "SP500"
results_SP_1reg$Reg <- c("VIX.csv")
results_SP_1reg$n_reg <- 1

results_SP_2reg <- run_all_models(
  s = "SP500.csv",
  r = c("VIX.csv", "TNX.csv"),
  f = 1,
  dist_max = 25
)
results_SP_2reg$Serie <- "SP500"
results_SP_2reg$Reg <- c("VIX.csv", "TNX.csv")
results_SP_2reg$n_reg <- 2
results_SP_3reg <- run_all_models(
  s = "SP500.csv",
  r = c("VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv"),
  f = 1,
  dist_max = 25
)
results_SP_3reg$Serie <- "SP500"
results_SP_3reg$Reg <- c("VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv")
results_SP_3reg$n_reg <- 3

results_SP_5reg <- run_all_models(
  s = "SP500.csv",
  r = c("VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv", "NASDAQ100.csv", "GOLD.csv"),
  f = 1,
  dist_max = 25
)
results_SP_5reg$Serie <- "SP500"
results_SP_5reg$Reg <- c("VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv", "NASDAQ100.csv", "GOLD.csv")
results_SP_5reg$n_reg <- 5

results_SP_10reg <- run_all_models(
  s = "SP500.csv",
  r = c("VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv", "NASDAQ100.csv", "GOLD.csv",
        "RUSSELL2000.csv", "HIGH_YIELD_BOND.csv", "WTI_OIL.csv",
        "MSCI_EM.csv", "IRX.csv"),
  f = 1,
  dist_max = 25
)
results_SP_10reg$Serie <- "SP500"
results_SP_10reg$Reg <- c("VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv", "NASDAQ100.csv", "GOLD.csv",
                          "RUSSELL2000.csv", "HIGH_YIELD_BOND.csv", "WTI_OIL.csv",
                          "MSCI_EM.csv", "IRX.csv")
results_SP_10reg$n_reg <- 10

### APPLE ----

results_AAPL_1reg <- run_all_models(
  s = "APPLE.csv",
  r = c("NASDAQ100.csv"),
  f = 1,
  dist_max = 25
)
results_AAPL_1reg$Serie <- "AAPL"
results_AAPL_1reg$Reg <- c("NASDAQ100.csv")
results_AAPL_1reg$n_reg <- 1

results_AAPL_2reg <- run_all_models(
  s = "APPLE.csv",
  r = c("NASDAQ100.csv", "TECH_SECTOR.csv"),
  f = 1,
  dist_max = 25
)
results_AAPL_2reg$Serie <- "AAPL"
results_AAPL_2reg$Reg <- c("NASDAQ100.csv", "TECH_SECTOR.csv")
results_AAPL_2reg$n_reg <- 2

results_AAPL_3reg <- run_all_models(
  s = "APPLE.csv",
  r = c("NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv"),
  f = 1,
  dist_max = 25
)
results_AAPL_3reg$Serie <- "AAPL"
results_AAPL_3reg$Reg <- c("NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv")
results_AAPL_3reg$n_reg <- 3

results_AAPL_5reg <- run_all_models(
  s = "APPLE.csv",
  r = c("NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv",
        "DOLLAR_INDEX.csv", "MICROSOFT.csv"),
  f = 1,
  dist_max = 25
)
results_AAPL_5reg$Serie <- "AAPL"
results_AAPL_5reg$Reg <- c("NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv",
                           "DOLLAR_INDEX.csv", "MICROSOFT.csv")
results_AAPL_5reg$n_reg <- 5

results_AAPL_10reg <- run_all_models(
  s = "APPLE.csv",
  r = c("NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv",
        "DOLLAR_INDEX.csv", "MICROSOFT.csv", "AMAZON.csv", "TNX.csv",
        "CONSUMER_DISCR.csv", "VIX.csv", "TSMC.csv"),
  f = 1,
  dist_max = 25
)
results_AAPL_10reg$Serie <- "AAPL"
results_AAPL_10reg$Reg <- c("NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv",
                           "DOLLAR_INDEX.csv", "MICROSOFT.csv", "AMAZON.csv", "TNX.csv",
                           "CONSUMER_DISCR.csv", "VIX.csv", "TSMC.csv")
results_AAPL_10reg$n_reg <- 10

### ETH ----

results_ETH_1reg <- run_all_models(
  s = "ETHEREUM.csv",
  r = c("BITCOIN.csv"),
  f = 1,
  dist_max = 25
)
results_ETH_1reg$Serie <- "ETH"
results_ETH_1reg$Reg <- c("BITCOIN.csv")
results_ETH_1reg$n_reg <- 1

results_ETH_2reg <- run_all_models(
  s = "ETHEREUM.csv",
  r = c("BITCOIN.csv", "TETHER.csv"),
  f = 1,
  dist_max = 25
)
results_ETH_2reg$Serie <- "ETH"
results_ETH_2reg$Reg <- c("BITCOIN.csv", "TETHER.csv")
results_ETH_2reg$n_reg <- 2

results_ETH_3reg <- run_all_models(
  s = "ETHEREUM.csv",
  r = c("BITCOIN.csv", "TETHER.csv", "BNB.csv"),
  f = 1,
  dist_max = 25
)
results_ETH_3reg$Serie <- "ETH"
results_ETH_3reg$Reg <- c("BITCOIN.csv", "TETHER.csv", "BNB.csv")
results_ETH_3reg$n_reg <- 3

results_ETH_5reg <- run_all_models(
  s = "ETHEREUM.csv",
  r = c("BITCOIN.csv", "TETHER.csv", "BNB.csv", "GOLD.csv", "NASDAQ100.csv"),
  f = 1,
  dist_max = 25
)
results_ETH_5reg$Serie <- "ETH"
results_ETH_5reg$Reg <- c("BITCOIN.csv", "TETHER.csv", "BNB.csv", "GOLD.csv", "NASDAQ100.csv")
results_ETH_5reg$n_reg <- 5

results_ETH_10reg <- run_all_models(
  s = "ETHEREUM.csv",
  r = c("BITCOIN.csv", "TETHER.csv", "BNB.csv", "GOLD.csv", "NASDAQ100.csv",
        "DOLLAR_INDEX.csv", "VIX.csv", "CARDANO.csv", "SOLANA.csv", "SP500.csv"),
  f = 1,
  dist_max = 25
)
results_ETH_10reg$Serie <- "ETH"
results_ETH_10reg$Reg <- c("BITCOIN.csv", "TETHER.csv", "BNB.csv", "GOLD.csv", "NASDAQ100.csv",
                          "DOLLAR_INDEX.csv", "VIX.csv", "CARDANO.csv", "SOLANA.csv", "SP500.csv")
results_ETH_10reg$n_reg <- 10

### BRENT ----

results_BR_1reg <- run_all_models(
  s = "BRENT_OIL.csv",
  r = c("WTI_OIL.csv"),
  f = 1,
  dist_max = 25
)
results_BR_1reg$Serie <- "Brent"
results_BR_1reg$Reg <- c("WTI_OIL.csv")
results_BR_1reg$n_reg <- 1

results_BR_2reg <- run_all_models(
  s = "BRENT_OIL.csv",
  r = c("WTI_OIL.csv", "NATURAL_GAS.csv"),
  f = 1,
  dist_max = 25
)
results_BR_2reg$Serie <- "Brent"
results_BR_2reg$Reg <- c("WTI_OIL.csv", "NATURAL_GAS.csv")
results_BR_2reg$n_reg <- 2

results_BR_3reg <- run_all_models(
  s = "BRENT_OIL.csv",
  r = c("WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv"),
  f = 1,
  dist_max = 25
)
results_BR_3reg$Serie <- "Brent"
results_BR_3reg$Reg <- c("WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv")
results_BR_3reg$n_reg <- 3

results_BR_5reg <- run_all_models(
  s = "BRENT_OIL.csv",
  r = c("WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv",
        "ENERGY_SECTOR.csv", "SP500.csv"),
  f = 1,
  dist_max = 25
)
results_BR_5reg$Serie <- "Brent"
results_BR_5reg$Reg <- c("WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv",
                         "ENERGY_SECTOR.csv", "SP500.csv")
results_BR_5reg$n_reg <- 5

results_BR_10reg <- run_all_models(
  s = "BRENT_OIL.csv",
  r = c("WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv",
        "ENERGY_SECTOR.csv", "SP500.csv", "MSCI_EM.csv",
        "TRANSPORTATION_ETF.csv", "TNX.csv", "COPPER.csv",
        "EUROSTOXX50.csv"),
  f = 1,
  dist_max = 25
)
results_BR_10reg$Serie <- "Brent"
results_BR_10reg$Reg <- c("WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv",
                         "ENERGY_SECTOR.csv", "SP500.csv", "MSCI_EM.csv",
                         "TRANSPORTATION_ETF.csv", "TNX.csv", "COPPER.csv",
                         "EUROSTOXX50.csv")
results_BR_10reg$n_reg <- 10

###bind finale ----
results_SP_2reg <- results_SP_2reg[results_SP_2reg$Modello != "ARIMA" & results_SP_2reg$Modello != "ETS",]
results_SP_3reg <- results_SP_3reg[results_SP_3reg$Modello != "ARIMA" & results_SP_3reg$Modello != "ETS",]
results_SP_5reg <- results_SP_5reg[results_SP_5reg$Modello != "ARIMA" & results_SP_5reg$Modello != "ETS",]
results_SP_10reg <- results_SP_10reg[results_SP_10reg$Modello != "ARIMA" & results_SP_10reg$Modello != "ETS",]

results_AAPL_2reg <- results_AAPL_2reg[results_AAPL_2reg$Modello != "ARIMA" & results_AAPL_2reg$Modello != "ETS",]
results_AAPL_3reg <- results_AAPL_3reg[results_AAPL_3reg$Modello != "ARIMA" & results_AAPL_3reg$Modello != "ETS",]
results_AAPL_5reg <- results_AAPL_5reg[results_AAPL_5reg$Modello != "ARIMA" & results_AAPL_5reg$Modello != "ETS",]
results_AAPL_10reg <- results_AAPL_10reg[results_AAPL_10reg$Modello != "ARIMA" & results_AAPL_10reg$Modello != "ETS",]

results_BR_2reg <- results_BR_2reg[results_BR_2reg$Modello != "ARIMA" & results_BR_2reg$Modello != "ETS",]
results_BR_3reg <- results_BR_3reg[results_BR_3reg$Modello != "ARIMA" & results_BR_3reg$Modello != "ETS",]
results_BR_5reg <- results_BR_5reg[results_BR_5reg$Modello != "ARIMA" & results_BR_5reg$Modello != "ETS",]
results_BR_10reg <- results_BR_10reg[results_BR_10reg$Modello != "ARIMA" & results_BR_10reg$Modello != "ETS",]

results_ETH_2reg <- results_ETH_2reg[results_ETH_2reg$Modello != "ARIMA" & results_ETH_2reg$Modello != "ETS",]
results_ETH_3reg <- results_ETH_3reg[results_ETH_3reg$Modello != "ARIMA" & results_ETH_3reg$Modello != "ETS",]
results_ETH_5reg <- results_ETH_5reg[results_ETH_5reg$Modello != "ARIMA" & results_ETH_5reg$Modello != "ETS",]
results_ETH_10reg <- results_ETH_10reg[results_ETH_10reg$Modello != "ARIMA" & results_ETH_10reg$Modello != "ETS",]

all_results <- dplyr::bind_rows(
  results_SP_1reg, results_SP_2reg, results_SP_3reg, results_SP_5reg, results_SP_10reg,
  results_AAPL_1reg, results_AAPL_2reg, results_AAPL_3reg, results_AAPL_5reg, results_AAPL_10reg,
  results_ETH_1reg, results_ETH_2reg, results_ETH_3reg, results_ETH_5reg, results_ETH_10reg,
  results_BR_1reg, results_BR_2reg, results_BR_3reg, results_BR_5reg, results_BR_10reg
)

all_results[all_results$Modello != "ETSX" & all_results$Modello != "ARIMAX",]$n_reg <- 0
View(all_results)


# Interpretazione ----

### PLOT ARIMAs ----
all_results[all_results$Modello == "ETS",]$Modello <- "ETSX"
all_results[all_results$Modello == "ARIMA",]$Modello <- "ARIMAX"
plot_df <- all_results[ all_results$Modello == "ETSX" | all_results$Modello == "ETS",] %>%
  group_by(Serie, Orizzonte, n_reg) %>%
  summarise(MAPE_test_mean = mean(MAPE_test, na.rm = TRUE), .groups = "drop")

plot_df$Orizzonte <- factor(plot_df$Orizzonte)

ggplot(plot_df, aes(x = Orizzonte, y = MAPE_test_mean, 
                    group = as.factor(n_reg), color = as.factor(n_reg))) +
  geom_point(size = 2) +
  geom_line(aes(group = as.factor(n_reg)), linewidth = 1) +
  facet_wrap(~ Serie, scales = "free_y") +
  labs(
    title = "PerformancexOrizzonte - ETS",
    x = "Orizzonte",
    y = "MAPE Test (media)",
    color = "Numero regressori"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

### valutazione stat ----
all_results <- all_results %>%
  mutate(ConReg = ifelse(n_reg > 0, "Con regressori", "Senza regressori"))
rank_df <- all_results %>%
  group_by(Orizzonte, Modello, ConReg) %>%
  summarise(MAPE_mean = mean(MAPE_test, na.rm = TRUE), .groups = "drop")
rank_df <- rank_df %>%
  group_by(Orizzonte) %>%
  mutate(Rank = rank(MAPE_mean, ties.method = "first")) %>%
  arrange(Orizzonte, Rank)

### prova mia ----
incrocio <- function(d, m1, n)
{
  for(m in unique(m1))
  {
    t <- d[m1 == m , ]
    cat(magenta$bold("------------", m, "------------\n"))
    
    tabella_risultati <- data.frame(
      Metriche = c("MAPE_test"),
      Media = c(mean(t$MAPE_test, na.rm = TRUE)),
      Mediana = c(median(t$MAPE_test, na.rm = TRUE)),
      DevStd = c(sd(t$MAPE_test, na.rm = TRUE))
    )
    
    print(tabella_risultati)
  }
  
  print(
    ggplot(d, aes(x=as.factor(m1) , y = MAPE_test)) +
      geom_boxplot(alpha = 0.7) +
      theme_minimal() +
      labs(
        title = paste("Distribuzione MAPE_test -", n),
        y = "MAPE_test",
        x = n,
        fill = "Orizzonte"
      )
  )
}

inc_sig <- function(d, m1, m2)
{
  cat(cyan$bold("Test differenze significative:\n"))
  
  a <- kruskal.test(m2 ~ as.factor(m1), data = d)$p.value
  cat("Test Kruskal,", black$bold("p-value:"), a, "\n")
  if(a <= 0.05){ cat(black$bold("Differenza significativa.\n")) }
  if(a > 0.05){ cat(black$bold("Differenza NON significativa.\n")) }
  
  cat(green$bold("Test differenze significative tra coppie:\n"))
  
  ris <- dunnTest(m2 ~ as.factor(m1), data = d, method = "bonferroni")$res
  ris$sig <- ifelse(ris$P.adj < 0.05, "*", "")
  print(ris)
}

inc_tot <- function(ds, categ, name)
{
  incrocio(ds, categ, name)
  cat(red$bold("SIGNIFICATIVITÀ MAPE_test\n"))
  inc_sig(ds, categ, ds$MAPE_test)
}

for (o in unique(all_results$Orizzonte))
{
  cat(blue$bold("----------Oriz:",o,"------------\n"))
  inc_tot(all_results[all_results$Orizzonte==o & all_results$Modello=="ARIMAX",], all_results[all_results$Orizzonte==o& all_results$Modello=="ARIMAX",]$n_reg, "Numero regressori")
}


#Parte2: ----
#i regressori non contengono info o è un problema dei nostri modelli?
prova_SP_self <- run_all_models(s = "SP500.csv", r = c("SP500.csv"), f=1, dist_max = 5)
mean(prova_SP_self[prova_SP_self$Modello == "ARIMA",]$MAPE_test)
mean(prova_SP_self[prova_SP_self$Modello == "ARIMAX",]$MAPE_test)
mean(prova_SP_self[prova_SP_self$Modello == "ETS",]$MAPE_test)
mean(prova_SP_self[prova_SP_self$Modello == "ETSX",]$MAPE_test)
#l'errore aumenta con il solo regressore!

#proviamo la performance con altri modelli
### REG CON MEZZI PESANTI ----
library(xgboost)
fit_models_pesanti <- function(s, r, f, dist, h) {
    data_list <- prepare_data(s, r, dist, h)
    dftr <- data_list$dftr
    dfts <- data_list$dfts
    
    target <- ts(dftr$Y, frequency = f)
    xtr <- if (ncol(dftr) > 2) as.matrix(dftr %>% select(-Y, -Date)) else NULL
    xts <- if (ncol(dfts) > 2) as.matrix(dfts %>% select(-Y, -Date)) else NULL
    y_train <- dftr$Y
    y_test  <- dfts$Y
    
    results <- list()
    
    # Utility per evitare errori di lunghezza forecast/test
    align_lengths <- function(pred, actual) {
      n <- min(length(pred), length(actual))
      list(pred = pred[1:n], actual = actual[1:n])
    }
    dtrain <- xgb.DMatrix(data = xtr, label = y_train)
    dtest  <- xgb.DMatrix(data = xts)
    #xgb
    params <- list(objective = "reg:squarederror",eta = 0.05,max_depth = 6,subsample = 0.8,colsample_bytree = 0.8)
    xgb_model <- xgb.train(params = params,data = dtrain,nrounds = 400,verbose = 0)
    pred_train <- predict(xgb_model, dtrain)
    pred_test  <- predict(xgb_model, dtest)
    train_acc <- mean(abs(y_train - pred_train) / abs(y_train)) * 100
    test_acc  <- mean(abs(y_test - pred_test) / abs(y_test)) * 100
    
    # fc <- forecast(mod_arimax, xreg = xts, h = h)
    # aligned <- align_lengths(fc$mean, y_test)
    
    last_pred <- pred_test[length(pred_test)]
    last_real <- y_test[length(y_test)]
    err_perc <- abs(last_pred - last_real) / last_real * 100
    
    results <- data.frame(
      Modello = "XGBOOST",
      Orizzonte = h,
      MAPE_train = train_acc,
      MAPE_test = test_acc,
      dist = dist,
      ultimaossprevista = last_pred,
      actualultimaoss = last_real,
      errore_ultima_perc = err_perc,
      n_reg = length(r)-1,
      Serie = substring(s,1,nchar(s)-4)
    )
  }
run_all_models_pesanti <- function(s, r, f, dist_max, horizons = c(1, 5, 15, 30, 60, 120)) {
  final_results <- list()
  
  for (dist in seq(1, (dist_max * 9), by = 9)) {
    cat("DIST =", dist, "\n")
    for (h in horizons) {
      cat("  Orizzonte =", h, "\n")
      res <- fit_models_pesanti(s, r, f, dist, h)
      final_results[[length(final_results) + 1]] <- res
    }
  }
  
  bind_rows(final_results)
}

XG_SP_1reg <- run_all_models_pesanti(
  s = "SP500.csv",
  r = c("SP500.csv","VIX.csv"),
  f = 1,
  dist_max = 10
)

XG_SP_2reg <- run_all_models_pesanti(
  s = "SP500.csv",
  r = c("SP500.csv","VIX.csv", "TNX.csv"),
  f = 1,
  dist_max = 10
)

XG_SP_3reg <- run_all_models_pesanti(
  s = "SP500.csv",
  r = c("SP500.csv","VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv"),
  f = 1,
  dist_max = 10
)

XG_SP_5reg <- run_all_models_pesanti(
  s = "SP500.csv",
  r = c("SP500.csv","VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv", "NASDAQ100.csv", "GOLD.csv"),
  f = 1,
  dist_max = 10
)

XG_SP_10reg <- run_all_models_pesanti(
  s = "SP500.csv",
  r = c("SP500.csv","VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv", "NASDAQ100.csv", "GOLD.csv",
        "RUSSELL2000.csv", "HIGH_YIELD_BOND.csv", "WTI_OIL.csv",
        "MSCI_EM.csv", "IRX.csv"),
  f = 1,
  dist_max = 10
)

### APPLE ----

XG_AAPL_1reg <- run_all_models_pesanti(
  s = "APPLE.csv",
  r = c("APPLE.csv","NASDAQ100.csv"),
  f = 1,
  dist_max = 10
)

XG_AAPL_2reg <- run_all_models_pesanti(
  s = "APPLE.csv",
  r = c("APPLE.csv","NASDAQ100.csv", "TECH_SECTOR.csv"),
  f = 1,
  dist_max = 10
)

XG_AAPL_3reg <- run_all_models_pesanti(
  s = "APPLE.csv",
  r = c("APPLE.csv","NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv"),
  f = 1,
  dist_max = 10
)

XG_AAPL_5reg <- run_all_models_pesanti(
  s = "APPLE.csv",
  r = c("APPLE.csv","NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv",
        "DOLLAR_INDEX.csv", "MICROSOFT.csv"),
  f = 1,
  dist_max = 10
)

XG_AAPL_10reg <- run_all_models_pesanti(
  s = "APPLE.csv",
  r = c("APPLE.csv","NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv",
        "DOLLAR_INDEX.csv", "MICROSOFT.csv", "AMAZON.csv", "TNX.csv",
        "CONSUMER_DISCR.csv", "VIX.csv", "TSMC.csv"),
  f = 1,
  dist_max = 10
)

### ETH ----

XG_ETH_1reg <- run_all_models_pesanti(
  s = "ETHEREUM.csv",
  r = c("ETHEREUM.csv","BITCOIN.csv"),
  f = 1,
  dist_max = 10
)

XG_ETH_2reg <- run_all_models_pesanti(
  s = "ETHEREUM.csv",
  r = c("ETHEREUM.csv","BITCOIN.csv", "TETHER.csv"),
  f = 1,
  dist_max = 10
)

XG_ETH_3reg <- run_all_models_pesanti(
  s = "ETHEREUM.csv",
  r = c("ETHEREUM.csv","BITCOIN.csv", "TETHER.csv", "BNB.csv"),
  f = 1,
  dist_max = 10
)

XG_ETH_5reg <- run_all_models_pesanti(
  s = "ETHEREUM.csv",
  r = c("ETHEREUM.csv","BITCOIN.csv", "TETHER.csv", "BNB.csv", "GOLD.csv", "NASDAQ100.csv"),
  f = 1,
  dist_max = 10
)

XG_ETH_10reg <- run_all_models_pesanti(
  s = "ETHEREUM.csv",
  r = c("ETHEREUM.csv","BITCOIN.csv", "TETHER.csv", "BNB.csv", "GOLD.csv", "NASDAQ100.csv",
        "DOLLAR_INDEX.csv", "VIX.csv", "CARDANO.csv", "SOLANA.csv", "SP500.csv"),
  f = 1,
  dist_max = 10
)

### BRENT ----

XG_BR_1reg <- run_all_models_pesanti(
  s = "BRENT_OIL.csv",
  r = c("BRENT_OIL.csv","WTI_OIL.csv"),
  f = 1,
  dist_max = 10
)

XG_BR_2reg <- run_all_models_pesanti(
  s = "BRENT_OIL.csv",
  r = c("BRENT_OIL.csv","WTI_OIL.csv", "NATURAL_GAS.csv"),
  f = 1,
  dist_max = 10
)

XG_BR_3reg <- run_all_models_pesanti(
  s = "BRENT_OIL.csv",
  r = c("BRENT_OIL.csv","WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv"),
  f = 1,
  dist_max = 10
)

XG_BR_5reg <- run_all_models_pesanti(
  s = "BRENT_OIL.csv",
  r = c("BRENT_OIL.csv","WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv",
        "ENERGY_SECTOR.csv", "SP500.csv"),
  f = 1,
  dist_max = 10
)

XG_BR_10reg <- run_all_models_pesanti(
  s = "BRENT_OIL.csv",
  r = c("BRENT_OIL.csv","WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv",
        "ENERGY_SECTOR.csv", "SP500.csv", "MSCI_EM.csv",
        "TRANSPORTATION_ETF.csv", "TNX.csv", "COPPER.csv",
        "EUROSTOXX50.csv"),
  f = 1,
  dist_max = 10
)
XG_SP_0reg <- run_all_models_pesanti(
  s = "SP500.csv",
  r = c("SP500.csv"),
  f = 1,
  dist_max = 10
)
XG_AAPL_0reg <- run_all_models_pesanti(
  s = "APPLE.csv",
  r = c("APPLE.csv"),
  f = 1,
  dist_max = 10
)
XG_BR_0reg <- run_all_models_pesanti(
  s = "BRENT_OIL.csv",
  r = c("BRENT_OIL.csv"),
  f = 1,
  dist_max = 10
)
XG_ETH_0reg <- run_all_models_pesanti(
  s = "ETHEREUM.csv",
  r = c("ETHEREUM.csv"),
  f = 1,
  dist_max = 10
)
all_XG <- dplyr::bind_rows(
  XG_SP_0reg,XG_AAPL_0reg,XG_ETH_0reg,XG_BR_0reg,XG_SP_1reg, XG_SP_2reg, XG_SP_3reg, XG_SP_5reg, XG_SP_10reg,
  XG_AAPL_1reg, XG_AAPL_2reg, XG_AAPL_3reg, XG_AAPL_5reg, XG_AAPL_10reg,
  XG_ETH_1reg, XG_ETH_2reg, XG_ETH_3reg, XG_ETH_5reg, XG_ETH_10reg,
  XG_BR_1reg, XG_BR_2reg, XG_BR_3reg, XG_BR_5reg, XG_BR_10reg
)

all_results <- all_results %>% select(-Reg,-ConReg)
tot <- rbind(all_results,all_XG)
tot$Modello <- as.factor(tot$Modello)
tot$Serie <- as.factor(tot$Serie)
summary(tot)
tot[tot$Serie == "APPLE",]$Serie <- "AAPL"
tot[tot$Serie == "BRENT_OIL", "Serie"] <- "Brent"
tot[tot$Serie == "ETHEREUM",]$Serie <- "ETH"

#Interpreto
#XGBOOST vs ARIMA-ETS senza regressori

base_models <- tot %>% filter(n_reg == 10)
ggplot(base_models, aes(x = Modello, y = MAPE_test, fill = as.factor(Orizzonte))) +
  geom_boxplot() +
  labs(title = "Confronto MAPE_test tra modelli senza regressori",
       y = "MAPE Test",
       x = "Modello") +
  theme_minimal()


anv <- aov((MAPE_test^(1/6))~factor(Modello) , data = base_models)
res <- residuals(anv)
qqnorm(res)
qqline(res,col=2)
summary(anv)
TukeyHSD(anv)

kruskal.test(MAPE_test ~ factor(Modello), data = base_models)
dunnTest(MAPE_test ~ factor(Modello), data = base_models,method="bonferroni")

#XGB da solo

xgb_results <- tot %>% filter(Modello == "XGBOOST")

# Calcola la media del MAPE_test per Serie, n_reg e Orizzonte
xgb_summary <- xgb_results %>%
  group_by(Serie, n_reg, Orizzonte) %>%
  summarise(MAPE_test_mean = mean(MAPE_test, na.rm = TRUE), .groups = "drop")

ggplot(xgb_summary, aes(x = Orizzonte, y = MAPE_test_mean, color = factor(n_reg), group = n_reg)) +
  geom_line(size = 1.2) +
  facet_wrap(~ Serie, scales = "free_y") +
  labs(
    title = "Prestazioni medie di XGBoost per n_reg e Orizzonte",
    x = "Orizzonte (horizon)",
    y = "MAPE test medio (%)",
    color = "Numero regressori"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )


# best arima/ets vs best xgb
# Filtra i modelli da confrontare
compare_results <- tot %>%
  filter(
    (Modello == "XGBOOST" & n_reg == 10) |
      (Modello %in% c("ARIMAX", "ETSX") & n_reg == 0)
  )

# Calcola la media del MAPE_test per Serie, Modello e Orizzonte
compare_summary <- compare_results %>%
  group_by(Serie, Modello, Orizzonte) %>%
  summarise(MAPE_test_mean = mean(MAPE_test, na.rm = TRUE), .groups = "drop")

ggplot(compare_summary, aes(x = Orizzonte, y = MAPE_test_mean, color = Modello, group = Modello)) +
  geom_line(size = 1.2) +
  facet_wrap(~ Serie, scales = "free_y") +
  labs(
    title = "Best XGB vs ARIMA vs ETS",
    x = "Orizzonte",
    y = "MAPE test medio (%)",
    color = "Modello"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )
