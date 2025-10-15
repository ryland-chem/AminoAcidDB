#wrote a function to read in the data
#this also folds the data into one large block to queue against, for pubchem, WOS, and WOS_deriv

load_mass_spectra_csvs <- function(folder) {
  
  files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
  
  data_list <- lapply(files, function(file) {
    df <- read.csv(file, stringsAsFactors = FALSE)
    
    if (ncol(df) >= 3) {
      colnames(df)[1:3] <- c("Name", "Adduct", "mz")
      return(df[, c("Name", "Adduct", "mz")])
    } else {
      warning(paste("Skipping:", file, "- fewer than 3 columns"))
      return(NULL)
    }
  })
  
  do.call(rbind, data_list)
  
}