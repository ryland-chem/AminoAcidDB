# AminoacidDB (v1.0)
## Authorship: Pawanjit K Sandhu (1), Ryland Giebelhaus (2), Ryan Hayward (3), Tingting Zhao (4), Alix Tucker (1), Daniel Gaudet (1), Tao Huan (4), Susan J Murch (1)
## (1) Department of Chemistry, University of British Columbia, Okanagan
## (2) Department of Chemistry, University of Victoria
## (3) Supra Research and Development
## (4) Department of Chemistry, University of British Columbia, Vancouver

## Contact: Dr. Ryland T. Giebelhaus, rgiebelhaus@uvic.ca

AminoAcidDB is a tool developed in collaboration between UVic by the Giebelhaus Laboratory, UBC Okanagan by the plantSMART research team, Supra Research and Development, and UBC Vancouver by the Huan laboratory to allow users to process their untargeted metabolomics data to putativley identify all known amino acids.

 ## Recent Updates and Changes (v1.4.1)
* Added a counter in to count number of times the tool has been used.
* Small bug fixes
* Fixed some spelling mistakes within the tool.


## Recent Updates and Changes (v1.4.0)
* Users can now search more than one database at a time.
* Small bug fixes
* More comments in the code for ease of use.
* Seperated M+H and monoisotopic into different databases.

## Current Features
* More than 1500 primary amino acids.
* The AminoacidDB shell to allow users to use their own dataset.

## How to run locally:
1. Download the most recent version of AminoacidDB
2. Make sure that your version of R and R Studio are up to date. To do this run the following:

```
install.packages("installr")
library(installr)
updateR()
```

3. Place all the .csv files and R file in the same directory, set this directory (or file) as your working directory in R Studio.

4. Compile the code in R Studio. If everything is up to date, then the code should compile succesfully and you will find a new window has opened with the application in it. You are now free to use HormonomicsDB on your local device, wthout internet connection!!

## Formatting Data:
* To use the main AminoacidDB tool, data should be formatted like this:

m/z | RT | Sample_1 | Sample_2 | Sample_3 
--- | --- | --- | --- | --- |
233.1 | 2.1 | 0 | 123.1 | 441.2
123.9 | 16.1 | 23441.2 | 222 | 0

* To format a custom database for the Shell, your data should like this in a .csv file.

Compound.Name | Adduct | Mass | RTP
--- | --- | --- | ---
Strictosidine | M+H | 368.1663 | 5.86

* If you don't have predicted retention times for the data, fill the RTP column with 1's.
* Experimental metabolomics data is still formatted the same for the Shell.

## Terms and Agreements
* AminoacidDB was developed for research use only and is not intended for use in diagnostic work. Dispite diligent validation and bug fixing, we are not responsible for any mistakes the application makes in data processing. Considering this, please inform us immediately of any bugs that you encounter. 

* We do not save any data that is uploaded to the server, it is immediatley deleted with every new session that you start.

* Please acknowledge the aforementioned authors in any work where AminoacidDB has been used.

