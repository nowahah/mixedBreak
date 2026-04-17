library(data.table)
library(readxl)


# source("R/print.R")
# source("R/plot.R")
# source("R/subset.R")


## ** load Basel dataset
dtW <- as.data.table(read_excel("../lmbreak/source/SDI_Psilo_22.11.2023.xlsx"))
basel <- melt(dtW, id.vars = c("Study","PatientID"), variable.name = "time", value.name = "score")
basel$Study <- as.factor(basel$Study)
basel$PatientID <- as.factor(basel$PatientID)
basel$time.num <- as.numeric(as.character(basel$time))
names(basel)[2] = "ID"

## ** load NRU Psilo dataset
dtW <- as.data.table(read_excel("../lmbreak/source/SDI_Psilo_22.11.2023.xlsx", sheet = "Psilo NRU"))
SDIpsilo.ext <- dtW[, c("CIMBI ID", "Time (min) elapsed since drug administration", "SDI score")]
names(SDIpsilo.ext) <- c("ID", "time", "score")
SDIpsilo.ext$ID = as.factor(SDIpsilo.ext$ID)
SDIpsilo.ext$score = as.numeric(SDIpsilo.ext$score)
rm(dtW)