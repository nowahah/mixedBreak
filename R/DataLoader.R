library(data.table)
library(readxl)

setwd("/indirect/student/noahdegoulange")

## ** load Basel dataset
dtW <- as.data.table(read_excel("lmbreak/source/SDI_Psilo_22.11.2023.xlsx"))
basel <- melt(dtW, id.vars = c("Study","PatientID"), variable.name = "time", value.name = "score")
basel$Study <- as.factor(basel$Study)
basel$PatientID <- as.factor(basel$PatientID)
basel$time.num <- as.numeric(as.character(basel$time))
names(basel)[2] = "ID"

## ** load NRU Psilo dataset
dtW <- as.data.table(read_excel("lmbreak/source/SDI_Psilo_22.11.2023.xlsx", sheet = "Psilo NRU"))
SDIpsilo <- dtW[, c("CIMBI ID", "Time (min) elapsed since drug administration", "SDI score")]
names(SDIpsilo) <- c("ID", "time", "score")
SDIpsilo$ID = as.factor(SDIpsilo$ID)
SDIpsilo$score = as.numeric(SDIpsilo$score)
rm(dtW)