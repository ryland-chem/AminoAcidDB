#readData for taking in data (server side)
#doing one by one for each WOS, WOS deriv, and 
#calling function load_mass_spectra_csvs for all directories

source("load_mass_spectra_csvs.R")

pubchem = load_mass_spectra_csvs("pubchem")
WOS = load_mass_spectra_csvs("LotusHMDB")
#WOS_deriv = load_mass_spectra_csvs("WOS_deriv")

PRTs_pubchem <- read.csv("rts/Pubchem_RTs.csv")
PRTs_pubchem <- data.frame(PRTs_pubchem)

PRTs_WOS <- read.csv("rts/Lotus_HMDB_RTs.csv")
PRTs_WOS <- data.frame(PRTs_WOS)

#PRTs_WOS_deriv <- read.csv("rts/WoS_derivatized_RTs.csv")
#PRTs_WOS_deriv <- data.frame(PRTs_WOS_deriv)

# Map dataset selection (e.g., "1", "2", "3") to actual data
db_map <- list("1" = pubchem, "2" = WOS)