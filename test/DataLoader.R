# library(data.table)
library(readxl)

## ** load Basel dataset
dtW <- data.table::as.data.table(
  readxl::read_excel("../lmbreak/source/SDI_Psilo_22.11.2023.xlsx")
)
basel <- data.table::melt(dtW, id.vars = c("Study","PatientID"), 
                          variable.name = "time", value.name = "score")
basel$Study <- as.factor(basel$Study)
basel$PatientID <- as.factor(basel$PatientID)
basel$time <- as.numeric(as.character(basel$time))
names(basel)[2] = "ID"
basel <- as.data.frame(basel)

## ** load NRU Psilo dataset (same file, sheet n2)
dtW <- data.table::as.data.table(
  readxl::read_excel("../lmbreak/source/SDI_Psilo_22.11.2023.xlsx", 
             sheet = "Psilo NRU")
)
SDIpsilo.ext <- dtW[,c("CIMBI ID", "Time (min) elapsed since drug administration", "SDI score")]
names(SDIpsilo.ext) <- c("ID", "time", "score")
SDIpsilo.ext$ID = as.factor(SDIpsilo.ext$ID)
SDIpsilo.ext$score = as.numeric(SDIpsilo.ext$score)
rm(dtW)