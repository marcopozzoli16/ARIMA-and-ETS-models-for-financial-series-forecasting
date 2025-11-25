#allineamento dati
library(dplyr)
library(readr)
library(purrr)

align_series_group <- function(series_entry, input_dir = ".", output_dir = "aligned") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  serie <- series_entry$serie
  regressori <- series_entry$regressori
  
  #Lista completa dei file
  files <- c(serie, regressori)
  
  #Carico tutti i dataset
  data_list <- lapply(files, function(f) {
    df <- read_csv(file.path(input_dir, f), show_col_types = FALSE)
    
    #nomi colonne
    names(df) <- toupper(names(df))
    if (!("DATE" %in% names(df))) stop(paste("Manca colonna DATE in:", f))
    
    df$DATE <- as.Date(df$DATE)
    df <- df[order(df$DATE), ]
    return(df)
  })
  names(data_list) <- files
  
  #Trovo le date comuni
  common_dates <- Reduce(intersect, lapply(data_list, function(df) df$DATE))
  
  #Filtro ciascuna serie sulle date comuni
  aligned_list <- lapply(seq_along(data_list), function(i) {
    df <- data_list[[i]] %>% filter(DATE %in% common_dates)
    return(df)
  })
  names(aligned_list) <- names(data_list)
  
  #Salvo i nuovi file allineati
  for (nm in names(aligned_list)) {
    write_csv(aligned_list[[nm]], file.path(output_dir, nm))
  }
  
  return(invisible(aligned_list))
}

series_list <- list(
  list(serie = "SP500.csv", 
       regressori = c("VIX.csv", "TNX.csv", "DOLLAR_INDEX.csv", "NASDAQ100.csv", 
                      "GOLD.csv", "RUSSELL2000.csv", "HIGH_YIELD_BOND.csv", 
                      "WTI_OIL.csv", "MSCI_EM.csv", "IRX.csv")),
  
  list(serie = "APPLE.csv", 
       regressori = c("NASDAQ100.csv", "TECH_SECTOR.csv", "SEMICONDUCTOR_ETF.csv", 
                      "DOLLAR_INDEX.csv", "MICROSOFT.csv", "AMAZON.csv", 
                      "TNX.csv", "CONSUMER_DISCR.csv", "VIX.csv", "TSMC.csv")),
  
  list(serie = "ETHEREUM.csv", 
       regressori = c("BITCOIN.csv", "TETHER.csv", "BNB.csv", "GOLD.csv", 
                      "NASDAQ100.csv", "DOLLAR_INDEX.csv", "VIX.csv", 
                      "CARDANO.csv", "SOLANA.csv", "SP500.csv")),
  
  list(serie = "BRENT_OIL.csv", 
       regressori = c("WTI_OIL.csv", "NATURAL_GAS.csv", "DOLLAR_INDEX.csv", 
                      "ENERGY_SECTOR.csv", "SP500.csv", "MSCI_EM.csv", 
                      "TRANSPORTATION_ETF.csv", "TNX.csv", "COPPER.csv", 
                      "EUROSTOXX50.csv"))
)

# Allineo tutti i gruppi
aligned_data <- lapply(series_list, align_series_group,
                       input_dir = ".", output_dir = "aligned")
