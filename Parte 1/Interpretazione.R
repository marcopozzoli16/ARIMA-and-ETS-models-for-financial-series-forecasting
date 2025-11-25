set.seed(123)
library(dplyr)
library(lubridate)
library(tseries)
library(forecast)
library(tsfeatures)
library(moments)
library(plotly)
library(crayon)
library(dunn.test)
library(FSA)
library(tidyverse)
library(recipes)
library(rsample)
library(parsnip)
library(workflows)
library(yardstick)
library(tune)
library(ranger)
library(xgboost)
library(lme4)
library(vip) 
library(DALEX)  
library(robustbase)
library(tidyr)
library(rpart)
library(rpart.plot)
library(nnet)
library(broom)
library(caret)
library(glmnet)
###Carico dati + gestione -----
df <- read.csv("MetricheComplete.csv")
#train/test split
dftr <- df[df$TrainTest == "train" ,]
dfts <- df[df$TrainTest == "test" ,]

#Tolgo colonne inutili
dfts <- dfts[ , c("Serie","NOME","RMSE","MAPE","MASE","MAE","Frequenza","Orizzonte","PreProc","InizioTrain")]


###ANALISI SERIE - COSTRUZIONE MATRIX ----
analyze_series <- function(d)
{
  d <- d %>% arrange(Date)
  tsd <- ts(d$Close, frequency = 5)
  
  #Trend, ACF, Season
  #interessanti sono: trend, seasonal_strenght (se freq!=1), linearity, e_acf1, spike
  stl <- stl_features(tsd)
  stl
  
  #Outliers
  diff <- tsd - tsclean(tsd)
  nout <- sum(diff != 0)
  
  #Volatilità
  logret <- diff(log(tsd))
  vol30 <- sqrt(var(tail(logret, 30)))
  voltot <- sqrt(var(logret, na.rm = TRUE))
  
  list(
    ACF = stl["e_acf1"] ,
    trend = stl["trend"],
    linearity = stl["linearity"],
    seasonality = stl["seasonal_strength"],
    spikes = stl["spike"],
    n_outliers = nout,
    Vol30 = vol30,
    Volatility = voltot
  )
}
build_feature_matrix <- function(files, path = ".") {
  
  results <- list()
  
  for (file in files) {
    
    d <- read.csv(file.path(path, file))
    d$Date <- as.Date(d$Date)
    
    met <- analyze_series(d)
    
    results[[file]] <- data.frame(
      series = substring(file,1,nchar(file)-4),
      ACF = met$ACF ,
      trend = met$trend,
      linearity = met$linearity,
      seasonality = met$seasonality,
      spikes = met$spikes,
      n_outliers = met$n_outliers,
      Vol30 = met$Vol30,
      Volatility = met$Volatility
    )
  }
  
  final <- do.call(rbind, results)
  rownames(final) <- final$series
  return(final)
}

nomeserie <- unique(dfts$Serie)
nomeserie <- paste(nomeserie,".csv",sep="")

p_serie <- build_feature_matrix(nomeserie)
p_serie <- p_serie %>% select(-series)
p_serie$Serie <- rownames(p_serie)

#Creo un dataset che ha per ogni riga anche le caratteristiche della rispettiva serie
dfT <- dfts %>%
  left_join(p_serie , by = "Serie")
### Continuo costr matrix ----

p_serie$Serie <- rownames(p_serie)

#Creo un dataset che ha per ogni riga anche le caratteristiche della rispettiva serie
dfT <- dfts %>%
  left_join(p_serie , by = "Serie")
best_combo <- dfT %>%
  group_by(Serie) %>%
  slice_min(order_by = MAPE, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(best_label = paste(NOME, PreProc, sep = "__")) %>%
  select(Serie, best_label)

dfT <- dfT %>% left_join(best_combo, by = "Serie")

#creo un mape standardizzato
dfT <- dfT %>% group_by(Serie) %>%
  mutate(MAPE_std = (MAPE - min(MAPE)) / (max(MAPE) - min(MAPE) + 1e-12)) %>%
  ungroup()

dfT <- dfT %>%
  mutate(NOME = as.factor(NOME),
         PreProc = as.factor(PreProc),
         Frequenza = as.factor(Frequenza),
         Orizzonte = as.factor(Orizzonte),
         combo = paste(NOME, PreProc, sep = "__"))

serie_features <- p_serie %>% select(Serie, everything()) %>% select(-Serie)
feat_names <- colnames(serie_features)

feat_names_existing <- feat_names[feat_names %in% colnames(dfT)]
reg_df <- dfT %>%
  mutate(NOME = as.factor(NOME), PreProc = as.factor(PreProc)) %>%
  select(Serie, NOME, PreProc, Frequenza, Orizzonte, MAPE_std, all_of(feat_names_existing))

#inutile per ora
### MODELLO E PREPROC -------


# Funzione che costruisce i dataset 'one-row-per-serie' e allena modelli interpretabili
train_interpretable_for_horizon <- function(h, dfts, p_serie, metric = "MAPE") {
  #1) scelgo la metrica colonna
  metric_col <- ifelse(metric == "MAPE", "MAPE", metric) 
  
  #2) ricavo per ogni (Serie, configurazione) la metricae
  df_h <- dfts %>% filter(Orizzonte == h)
  
  best_per_series <- df_h %>%
    group_by(Serie) %>%
    slice_min(order_by = .data[[metric_col]], n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(Serie, bestModel = NOME, bestPreProc = PreProc)
  
  #3)prendo caratteristiche serie da p_serie (una riga per serie)
  feat <- p_serie %>% select(Serie, everything())
  
  df_clf <- feat %>% inner_join(best_per_series, by = "Serie")
  
  #4)Preparo dati: rimuovi Serie da predictors
  predictors <- df_clf %>% select(-Serie, -bestModel, -bestPreProc)
  
  #5) Train/test split
  set.seed(123)
  folds <- createFolds(df_clf$bestModel, k = 5, list = TRUE, returnTrain = FALSE)
  
  #6A)BestModel con CART
  cart_model <- rpart(bestModel ~ ., data = df_clf %>% select(-bestPreProc, -Serie),
                      method = "class", control = rpart.control(cp = 0.03, minsplit = 10))
  
  #stampo l'albero e le regole
  cat("\n--- CART (bestModel) summary ---\n")
  print(cart_model)
  rpart.plot(cart_model, main = paste0("CART bestModel, horizon ", h))
  
  # estrai regole
  if(requireNamespace("rpart.plot", quietly = TRUE)) {
    cat("\n--- Decision rules (bestModel) ---\n")
    try(rpart.plot::rpart.rules(cart_model), silent = TRUE)
  }
  
  #variable importance
  vi_cart_model <- vip::vi(cart_model)
  cat("\n--- Variable importance (CART bestModel) ---\n")
  print(vi_cart_model)
  
  #preproc
  cart_pp <- rpart(bestPreProc ~ ., data = df_clf %>% select(-bestModel, -Serie),
                   method = "class", control = rpart.control(cp = 0.03, minsplit = 10))
  cat("\n--- CART (bestPreProc) ---\n")
  print(cart_pp)
  rpart.plot(cart_pp, main = paste0("CART bestPreProc, horizon ", h))
  if(requireNamespace("rpart.plot", quietly = TRUE)) {
    cat("\n--- Decision rules (bestPreProc) ---\n")
    try(rpart.plot::rpart.rules(cart_pp), silent = TRUE)
  }
  vi_cart_pp <- vip::vi(cart_pp)
  cat("\n--- Variable importance (CART bestPreProc) ---\n")
  print(vi_cart_pp)
  
  recommend <- function(new_feats) {
    # new_feats: one-row tibble with same colnames as p_serie without Serie
    bm_prob_cart <- predict(cart_model, new_feats, type = "prob")
    pp_prob_cart <- predict(cart_pp, new_feats, type = "prob")
    
    list(
      cart_bestModel = bm_prob_cart,
      cart_bestPreProc = pp_prob_cart,
    )
  }
  
  return(list(
    horizon = h,
    df = df_clf,
    cart_model = cart_model,
    # multinom_model = multinom_model,
    cart_pp = cart_pp,
    # multinom_pp = multinom_pp,
    recommend = recommend
  ))
}

res_h15 <- train_interpretable_for_horizon(15, dfts = dfts, p_serie = p_serie)
res_h30 <- train_interpretable_for_horizon(30, dfts = dfts, p_serie = p_serie)
res_h45 <- train_interpretable_for_horizon(45, dfts = dfts, p_serie = p_serie)
res_h60 <- train_interpretable_for_horizon(60, dfts = dfts, p_serie = p_serie)
res_h120 <- train_interpretable_for_horizon(120, dfts = dfts, p_serie = p_serie)

### SISTEMA ESPERTO MODxPREPROC ----

#dataset
dfrf <- dfts %>%
  group_by(Serie, Orizzonte) %>%
  slice_min(order_by = MAPE, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Serie, Orizzonte, bestModel = NOME, bestPreProc = PreProc)

feats <- p_serie %>% select(Serie, everything())
dfrf <- feats %>% inner_join(dfrf, by = "Serie")

dfrf$bestModel  <- factor(dfrf$bestModel)
dfrf$bestPreProc <- factor(dfrf$bestPreProc)

set.seed(1234)

#F1 per caret
f1 <- function(data, lev = NULL, model = NULL) {
  precision <- posPredValue(data$pred, data$obs, positive = lev[1])
  recall    <- sensitivity(data$pred, data$obs, positive = lev[1])
  F1        <- (2 * precision * recall) / (precision + recall)
  c(F1 = F1)
}

Control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  sampling = "up",
  classProbs = TRUE,
  summaryFunction = f1,
  verboseIter = TRUE
)

gridrf <- expand.grid(
  .mtry = 1:8
)

#RF
rf_robust <- train(
  bestModel ~ ACF + trend + linearity + seasonality +
    spikes + n_outliers + Vol30 + Volatility + Orizzonte,
  data = dfrf,
  method = "rf",
  tuneGrid = gridrf,
  trControl = Control,
  metric = "F1",
  ntree = 2000
)

rf_robust
confusionMatrix(rf_robust)

#pp
Control_pp <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  sampling = "up",
  classProbs = TRUE,
  summaryFunction = multiClassSummary,
  verboseIter = TRUE
)


gridrf_pp <- expand.grid(
  .mtry = 1:8
)
library(MLmetrics)
#RF
rf_preproc <- train(
  bestPreProc ~ ACF + trend + linearity + seasonality +
    spikes + n_outliers + Vol30 + Volatility + Orizzonte,
  data = dfrf,
  method = "rf",
  tuneGrid = gridrf_pp,
  trControl = Control_pp,
  metric = "Kappa",        # più robusta dell'accuratezza per classi sbilanciate
  ntree = 2000
)

rf_preproc
### Frequency -----

ggplot(ris_freq, aes(x = reorder(Serie, p.value), y = p.value)) +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Significatività dell'effetto della Frequenza per Serie",
       x = "Serie",
       y = "p-value") +
  theme_minimal()
dfTclean <- dfT %>%
  group_by(Frequenza) %>%
  mutate(
    Q1 = quantile(MAPE, 0.25, na.rm = TRUE),
    Q3 = quantile(MAPE, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 3 * IQR,
    upper = Q3 + 3 * IQR
  ) %>%
  filter(MAPE >= lower, MAPE <= upper) %>%
  select(-Q1, -Q3, -IQR, -lower, -upper)

anov <- aov(log(MAPE+0.00000001) ~ as.factor(Frequenza), data = dfTclean)
summary(anov)
resF <- residuals(anov)
qqnorm(resF)
qqline(resF,col=2)

kruskal.test(MAPE ~ factor(Frequenza), data = dfTclean)
dunnTest(MAPE ~ factor(Frequenza), data = dfTclean, method="bonferroni")


### InizioTrain ----
df_perf <- dfts %>%
  select(Serie, InizioTrain, MAPE, RMSE, MASE, MAE, Orizzonte)
mean(df_perf[df_perf$InizioTrain == 2015,]$MAPE)
mean(df_perf[df_perf$InizioTrain == 2020,]$MAPE)
mean(df_perf[df_perf$InizioTrain == 2023,]$MAPE)
mean(df_perf[df_perf$InizioTrain == 2025,]$MAPE)
summary(df_perf[df_perf$InizioTrain == 2015,]$MAPE)
summary(df_perf[df_perf$InizioTrain == 2020,]$MAPE)
summary(df_perf[df_perf$InizioTrain == 2023,]$MAPE)
summary(df_perf[df_perf$InizioTrain == 2025,]$MAPE)

ggplot(df_perf, aes(x = factor(InizioTrain), y = MAPE)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.3) +
  theme_minimal() +
  labs(
    x = "Inizio del Training",
    y = "MAPE",
    title = "Andamento della performance al variare di InizioTrain"
  )
library(dplyr)

df_perfclean <- df_perf %>%
  group_by(InizioTrain) %>%
  mutate(
    Q1 = quantile(MAPE, 0.25, na.rm = TRUE),
    Q3 = quantile(MAPE, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower = Q1 - 3 * IQR,
    upper = Q3 + 3 * IQR
  ) %>%
  filter(MAPE >= lower, MAPE <= upper) %>%
  select(-Q1, -Q3, -IQR, -lower, -upper)

#test stat
anova <- aov(MAPE ~ factor(InizioTrain), data = df_perfclean)
summary(anova)
res <-residuals(anova)
qqnorm(res)
qqline(res,col=2)
TukeyHSD(anova)
#NON c'è normalità, non sono risultati attendibili

kruskal.test(MAPE ~ factor(InizioTrain), data = df_perf)
library(FSA)
dunnTest(MAPE ~ factor(InizioTrain), data = df_perf,method="bonferroni")

#levo 2025, vediamo
df_perf2 <- df_perf[df_perf$InizioTrain != 2025 ,]
anova <- aov(MAPE ~ factor(InizioTrain), data = df_perf2)
summary(anova)

TukeyHSD(anova)

kruskal.test(MAPE ~ factor(InizioTrain), data = df_perf2)
library(FSA)
dunnTest(MAPE ~ factor(InizioTrain), data = df_perf2,method="Bonferroni")

#Ho differenze significative per quanto riguarda le prime tre date e 2025.
#2023 è leggermente peggiore di 2020, che è più o meno uguale a 2015.
#ciò significa che mi servono almeno due anni di dati per avere buoni risultati.
library(ggpubr)

# Assumiamo df_perf già disponibile con colonne: Serie, InizioTrain, MAPE

# 1. Boxplot MAPE ~ InizioTrain
ggplot(df_perf, aes(x = factor(InizioTrain), y = MAPE, fill = factor(InizioTrain))) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.2, color = "black") +
  labs(x = "Inizio Train", y = "MAPE", title = "Distribuzione MAPE per Inizio Train") +
  theme_minimal() +
  theme(legend.position = "none")

### perf x orizzonte ----

eps <- 1e-6   # offset piccolo e standard in finanza/time series

dfts_inc <- dfts %>%
  group_by(Serie, NOME) %>%
  arrange(Orizzonte) %>%
  mutate(
    MAPE_log = log(MAPE + eps),
    MAPE_base_log = first(MAPE_log),
    MAPE_log_inc = MAPE_log - MAPE_base_log   # questo è lo “slope" robusto
  ) %>%
  ungroup()

dfts_inc <- dfts_inc %>%
  group_by(Serie) %>%
  mutate(MAPE_inc_centered = MAPE_log_inc - median(MAPE_log_inc)) %>%
  ungroup()

lmrobfit <- lmrob(
  MAPE_log_inc ~ as.factor(Orizzonte) * NOME,
  data = dfts_inc,
  control = lmrob.control(maxit=500, tuning.psi=3.0)
)
summary(lmrobfit)$coefficients


#Parte per grafico

eps <- 1e-6

dfts2 <- dfts %>%
  group_by(Serie, NOME) %>%
  arrange(Orizzonte) %>%
  mutate(
    MAPE_log      = log(MAPE + eps),
    MAPE_log_base = first(MAPE_log),
    MAPE_log_inc  = MAPE_log - MAPE_log_base,  # variazione “robusta”
    
    RMSE_log      = log(RMSE + eps),
    RMSE_log_base = first(RMSE_log),
    RMSE_log_inc  = RMSE_log - RMSE_log_base
  ) %>%
  ungroup()

control <- lmrob.control(maxit = 500, tuning.psi = 3.0)

mod_mape <- lmrob(MAPE_log_inc ~ as.factor(Orizzonte), data=dfts2, control=control)
mod_rmse <- lmrob(RMSE_log_inc ~ as.factor(Orizzonte), data=dfts2, control=control)

summary(mod_mape)$coefficients
summary(mod_rmse)$coefficients

mod_mape_ser <- lmrob(MAPE_log_inc ~ as.factor(Orizzonte) * NOME, data=dfts2, control=control)
summary(mod_mape_ser)$coefficients

ggplot(dfts2, aes(x = as.factor(Orizzonte), y = MAPE_log_inc)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "rlm", se = TRUE, formula = y ~ x) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Trend medio del MAPE al crescere dell’Orizzonte",
    x = "Orizzonte",
    y = "Incremento log(MAPE)"
  )

#vedo i coefficienti
coefs <- coef(mod_mape)

# togliamo l’intercetta
coefs_or <- coefs[-1]

# trasformazione in percentuale
perc_change <- (exp(coefs_or) - 1) * 100

round(perc_change, 2)#è l'aumento in percentuale rispetto all'orizzonte 15.

#Solo ETS
mod_mape <- lmrob(MAPE_log_inc ~ as.factor(Orizzonte), data=dfts2[dfts2$NOME == "ETS",], control=control)
summary(mod_mape)$coefficients
#vedo i coefficienti
coefs <- coef(mod_mape)
# togliamo l’intercetta
coefs_or <- coefs[-1]
# trasformazione in percentuale
perc_change <- (exp(coefs_or) - 1) * 100
round(perc_change, 2)

#Solo ARIMA
mod_mape <- lmrob(MAPE_log_inc ~ as.factor(Orizzonte), data=dfts2[dfts2$NOME == "ARIMA",], control=control)
summary(mod_mape)$coefficients
#vedo i coefficienti
coefs <- coef(mod_mape)
# togliamo l’intercetta
coefs_or <- coefs[-1]
# trasformazione in percentuale
perc_change <- (exp(coefs_or) - 1) * 100
round(perc_change, 2)

#ORIZ*NOME
mod_mape <- lmrob(MAPE_log_inc ~ as.factor(Orizzonte)*NOME, data=dfts2, control=control)
summary(mod_mape)$coefficients
#vedo i coefficienti
coefs <- coef(mod_mape)
# togliamo l’intercetta
coefs_or <- coefs[-1]
# trasformazione in percentuale
perc_change <- (exp(coefs_or) - 1) * 100
round(perc_change, 2)
