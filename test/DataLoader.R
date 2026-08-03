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

# Reducing length of leading/ending constant sequences of observations
# up to tolerance `tol` (nothing happens if tol<0)
library(dplyr)
# Trimming quasi-constant sequences (up to `tol`) of length >= `m̀in.length`
# to match length = `min.length-1L`
min.length <- 2L # >= 2L , minimum length of sequence of be trimmed
tol <- 2 # >= 0 , numerical tolerance to declare a sequence as constant
# nothing happens if tol < 0
basel <- basel %>%
  filter(!is.na(score)) %>%
  group_by(ID) %>%
  mutate(
    score.diff = c(diff(score, min.length-1L), rep(NA, min.length-1L)),
    var.high = 1*(abs(score.diff) > tol),
    var.cumsum = cumsum(var.high),
    row.nb = if_else(cumsum(var.high) == 0, 
                     row_number(), 0),
    beg.trail = row.nb == 1:n(), # beginning trailing OK
    rev.high = c(rev(var.high[1:(n()-(min.length-1L))]), 
                 var.high[(n()-min.length+2L):n()]),
    # rev.cumsum = cumsum(rev.high),
    # rev.row.nb = c(
    #   (if_else(cumsum(rev.high) == 0, row_number(), 0)[1:(n()-(min.length-1L))]),
    #   rep(NA, min.length-1L)),
    # rev.row.nb = rev(rev.row.nb > 0),
    end.trail = rev(cumsum(rev.high) == 0),
    type = if_else(beg.trail | end.trail, "trailing", "signal"),
    type = if_else(is.na(type), "signal", type)
  ) %>%
  select(names(basel), type)

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