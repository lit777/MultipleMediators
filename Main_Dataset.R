##------------------------------------------------------------------------------
## Required libraries
##------------------------------------------------------------------------------

#----- Faster data manipulation with keys
library(data.table)

#----- Distance matrix
library(fields)

#----- File name manipulation
library(tools)
library(utils)

#----- Date and time manipulation
library(zoo)

#----- Current Year (Time) and Past Year
TC <- 2005
TP <- TC-1

##------------------------------------------------------------------------------
## Big data download
##------------------------------------------------------------------------------

powerplant_raw <- fread("coal_unit_monthly.csv", header=TRUE)

#------ Coal-burning plants as of Jan. Year TC
powerplant_current <- subset(subset(subset(powerplant_raw, Month==1), Fuel1.IsCoal+Fuel2.IsCoal>0), Year==TC)

#------ Coal-burning plants of Year TP
powerplant_pre <- subset(subset(powerplant_raw, Year==TP), Fuel2.IsCoal+Fuel1.IsCoal>0 | !is.na(Heat.Input..MMBtu.))


##------------------------------------------------------------------------------
## Averages of Unit-level covariates weighted by heat input
##------------------------------------------------------------------------------

#------ Sulfur content weighted by Heat Input
powerplant_sul <- subset(powerplant_pre, !is.na(Sulfur.Content))[ ,list(
        avg.sulfur = mean(Sulfur.Content, na.rm=TRUE),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        sulfur_content = weighted.mean(avg.sulfur, avg.heat)),
        by=c("Facility.ID..ORISPL.")]

#------ S_n_CR weighted by Heat Input
powerplant_S_n_CR <- subset(powerplant_pre, !is.na(S_n_CR))
setkeyv(powerplant_S_n_CR, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_S_n_CR <- unique(powerplant_S_n_CR)[ ,list(
        S_n_CR = mean(S_n_CR, na.rm=TRUE),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        S_n_CR = weighted.mean(S_n_CR, avg.heat)),
        by=c("Facility.ID..ORISPL.")]

#------ LowNOxBurner weighted by Heat Input
powerplant_LowNOxBurner <- subset(powerplant_pre, !is.na(LowNOxBurner))
setkeyv(powerplant_LowNOxBurner, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_LowNOxBurner <- unique(powerplant_LowNOxBurner)[ ,list(
        LowNOxBurner = mean(LowNOxBurner, na.rm=TRUE),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        LowNOxBurner = weighted.mean(LowNOxBurner, avg.heat)),
        by=c("Facility.ID..ORISPL.")]

#------ OtherNOxControl weighted by Heat Input
powerplant_OtherNOxControl <- subset(powerplant_pre, !is.na(OtherNOxControl))
setkeyv(powerplant_OtherNOxControl, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_OtherNOxControl <- unique(powerplant_OtherNOxControl)[ ,list(
        OtherNOxControl = mean(OtherNOxControl, na.rm=TRUE),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        OtherNOxControl = weighted.mean(OtherNOxControl, avg.heat)),
        by=c("Facility.ID..ORISPL.")]

#------ anyNOxControl weighted by Heat Input
powerplant_anyNOxControl <- subset(powerplant_pre[ ,anyNOxControl := ifelse(LowNOxBurner + OtherNOxControl==0, 0, 1)], !is.na(anyNOxControl ))
setkeyv(powerplant_anyNOxControl, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_anyNOxControl <- unique(powerplant_anyNOxControl)[ ,list(
        anyNOxControl = mean(anyNOxControl, na.rm=TRUE),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        anyNOxControl = weighted.mean(anyNOxControl, avg.heat)),
        by=c("Facility.ID..ORISPL.")]

#------ NumNOxControls weighted by Heat Input
powerplant_NumNOxControls <- subset(powerplant_pre, !is.na(NumNOxControls))
setkeyv(powerplant_NumNOxControls, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_NumNOxControls <- unique(powerplant_NumNOxControls)[ ,list(
        NumNOxControls = mean(NumNOxControls, na.rm=TRUE),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        NumNOxControls = weighted.mean(NumNOxControls, avg.heat)),
        by=c("Facility.ID..ORISPL.")]

#------ PM2.5 Scrubber weighted by Heat Input
powerplant_pm <- subset(powerplant_pre, !is.na(Has.PM.Scrub))
setkeyv(powerplant_pm, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_pm <- unique(powerplant_pm)[ ,list(
        Has.PM.Scrub = mean(Has.PM.Scrub, na.rm=TRUE),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        Has.PM.Scrub = weighted.mean(Has.PM.Scrub, avg.heat)),
        by=c("Facility.ID..ORISPL.")]

#------ Phase 2 (indicator) weighted by Heat Input
powerplant_p2 <- subset(powerplant_pre, !is.na(Is.Phase2))
setkeyv(powerplant_p2, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_p2 <- unique(powerplant_p2)[ ,list(
        Is.Phase2 = mean(Is.Phase2, na.rm=TRUE),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        Is.Phase2 = weighted.mean(Is.Phase2, avg.heat)),
        by=c("Facility.ID..ORISPL.")]


#------ Missing Initial Year weighted by Heat Input
powerplant_init <- copy(powerplant_pre)
setkeyv(powerplant_init, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_init <- unique(powerplant_init)[ ,list(
        Init = if(is.na(Initial.Year.of.Operation)) 1 else 0,
        Heat.Input..MMBtu. = Heat.Input..MMBtu.),
        by=c("Facility.ID..ORISPL.", "Unit.ID","Month")]
powerplant_init <- powerplant_init[,list(
        Init = mean(Init),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        Init = weighted.mean(Init, avg.heat)),
        by=c("Facility.ID..ORISPL.")]

#------ Capacity
powerplant_Capacity <- subset(powerplant_pre, !is.na(Capacity))
setkeyv(powerplant_Capacity, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_Capacity <- unique(powerplant_Capacity)[ ,list(
        Capacity = sum(Capacity, na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.")]

#------ PctCapacity
powerplant_PctCapacity <- subset(powerplant_pre, !is.na(PctCapacity))
setkeyv(powerplant_PctCapacity, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_PctCapacity <- unique(powerplant_PctCapacity)[ ,list(
        PctCapacity = sum(PctCapacity, na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.")]

#------ HeatRate
powerplant_HeatRate <- subset(powerplant_pre, !is.na(HeatRate))
setkeyv(powerplant_HeatRate, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_HeatRate <- unique(powerplant_HeatRate)[ ,list(
        HeatRate = sum(HeatRate, na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.")]

#------ Operating time weighted by Heat Input
powerplant_time <- subset(powerplant_pre, !is.na(Operating.Time))
setkeyv(powerplant_time, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
powerplant_time <- unique(powerplant_time)[ ,list(
        Time = sum(Operating.Time, na.rm=TRUE),
        avg.heat = sum(Heat.Input..MMBtu., na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")][ ,list(
        Time = weighted.mean(Time, avg.heat)),
        by=c("Facility.ID..ORISPL.")]

#------ AnySO2Control (indicator) weighted by Heat Input
powerplant_SO2_P <- copy(powerplant_pre)
powerplant_SO2_C <- subset(powerplant_current, !is.na(AnySO2Control))

setkeyv(powerplant_SO2_P, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))
setkeyv(powerplant_SO2_C, c("Facility.ID..ORISPL.", "Unit.ID", "Month"))

powerplant_SO2_P <- unique(powerplant_SO2_P)[ ,list(
        Heat = sum(Heat.Input..MMBtu.,na.rm=TRUE)),
        by=c("Facility.ID..ORISPL.", "Unit.ID")]
powerplant_SO2_C <- unique(powerplant_SO2_C)[ ,list(
        SO2.sc = AnySO2Control),
        by=c("Facility.ID..ORISPL.", "Unit.ID")]

powerplant_SO2 <- subset(
        merge(powerplant_SO2_P, powerplant_SO2_C, by=c("Facility.ID..ORISPL.", "Unit.ID"))[ ,list(
        SO2.Sc  = weighted.mean(SO2.sc , Heat)),
        by=c("Facility.ID..ORISPL.")],
    !is.na(SO2.Sc))


##------------------------------------------------------------------------------
## Treatment Assignment
##------------------------------------------------------------------------------

#------ Cut-off value (treated vs. untreated)
cutoff <- 0.1

#------ Define Treatment
powerplant_SO2 <- powerplant_SO2[ ,list(
        SO2 = ifelse(SO2.Sc >= cutoff, 1, 0)),
        by=c("Facility.ID..ORISPL.")]

#------ Index of Power plants with complete observations
ID <- data.table(Reduce(intersect, list(
        powerplant_SO2[,Facility.ID..ORISPL.],
        powerplant_sul[,Facility.ID..ORISPL.],
        powerplant_S_n_CR[,Facility.ID..ORISPL.],
        powerplant_NumNOxControls[,Facility.ID..ORISPL.],
        powerplant_PctCapacity[,Facility.ID..ORISPL.],
        powerplant_time[,Facility.ID..ORISPL.])))

#------ Subsetting Treatment with complete observations
setkeyv(powerplant_SO2, c("Facility.ID..ORISPL."))
SO2 <- powerplant_SO2[ID, SO2]


##------------------------------------------------------------------------------
## Download and create annual AQS database
##  http://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/download_files.html
##
## year: time span
##------------------------------------------------------------------------------


#------ Download from EPA webpage (Annual)

getAQSData.Annual <- function(year = TP:TC) {
  code <- "annual_all_"
  name <- "all"
  dirdata <- file.path("Data_AQS", name)
  dir.create(dirdata, showWarnings = FALSE)
  files <- paste("http://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/", code,
                 year, ".zip", sep = "")
  for (i in 1:length(files)) {
    #------ Loop stage
    print(year[i])
    url <- files[i]
    file <- basename(url)
    download.file(url, file)
    untar(file, compressed = 'gzip', exdir = dirdata)
  }
  print("Purge downloaded zip files")
  zipfiles <- dir(path = ".",  pattern = "\\.zip$")
  file.remove(zipfiles)
}

#----- Download from EPA webpage (Daily)

getAQSData.Daily <- function(year = TP:TC) {
    code <- "daily_88101_"
    name <- "daily"
    dirdata <- file.path("Data_AQS", name)
    dir.create(dirdata, showWarnings = FALSE)
    files <- paste("http://aqsdr1.epa.gov/aqsweb/aqstmp/airdata/", code,
    year, ".zip", sep = "")
    for (i in 1:length(files)) {
        #----- Loop stage
        print(year[i])
        url <- files[i]
        file <- basename(url)
        download.file(url, file)
        untar(file, compressed = 'gzip', exdir = dirdata)
    }
    print("Purge downloaded zip files")
    zipfiles <- dir(path = ".",  pattern = "\\.zip$")
    file.remove(zipfiles)
}



#------ Load annual data

loadAnnualAverage <- function(year) {
  #----- Loop stage
  print(year)
  code <- "annual_all_"
  name <- "all"
  DO <- fread(file.path("Data_AQS", name, paste0(code, year, ".csv")))
  setnames(DO, make.names(colnames(DO)))
  return(DO)
}

#----- Load daily data

loadDailyAverage <- function(year) {
    #----- Loop stage
    print(year)
    code <- "daily_88101_"
    name <- "daily"
    DO <- fread(file.path("Data_AQS", name, paste0(code, year, ".csv")))
    setnames(DO, make.names(colnames(DO)))
    return(DO)
}



#----- Distance matrix between pollutant monitors and facilities

distmatMonitorsFromFacility <- function(Monitors, Fac = fread("arp_facility.csv")) {
  Fac <- Fac[ , c("Facility.Name",
                 "Facility.ID..ORISPL.",
                 "Facility.Latitude",
                 "Facility.Longitude",
                 "State",
                 "County",
                 "FIPS"),
             with = FALSE]
  setnames(Fac, "State", "Facility.State")
  setnames(Fac, "County", "Facility.County")
  setnames(Fac, "FIPS", "Facility.FIPS")
  
  m <- Monitors
  z <- Fac
  
  setkeyv(m, c("Latitude", "Longitude"))
  mu <- unique(m)
  mu <- cbind(paste("Mon", seq(nrow(mu)), sep = ""), mu)
  setnames(mu, "V1", "Monitor.ID")
  
  setkey(z)
  zu <- unique(z)
  
  distmat <- rdist.earth(mu[, list(Longitude, Latitude)],
                         zu[, list(Facility.Longitude, Facility.Latitude)])
  colnames(distmat) <- 1:ncol(distmat)
  rownames(distmat) <- mu$Monitor.ID
  
  dd <- as.data.frame(cbind(1:nrow(distmat), distmat))
  names(dd)[1] <- "Monitor.ID"
  
  Distmat <- data.table(dd)
  CDistmat <- Distmat[ , list(variable = names(.SD),
                              value = unlist(.SD, use.names = FALSE)),
                      by = "Monitor.ID"]
  setnames(CDistmat, "variable", "Facility.Tmp.ID")
  setnames(CDistmat, "value", "Distance")
  setkey(CDistmat, Distance)
  return(list(CDistmat = CDistmat, mu = mu, zu = zu))
}

#----- All monitors within a given range of the facilities

withinMonitorsFromFacility <- function(CDistmat, withinMiles = 18.641136) {
  mu <- CDistmat$mu
  zu <- CDistmat$zu
  CDistmat <- CDistmat$CDistmat
  Closest.Within <- CDistmat[Distance < withinMiles]
  Closest.Within[, Facility.Tmp.ID := as.integer(Facility.Tmp.ID)]
  M <- cbind(mu[Closest.Within$Monitor.ID],
             zu[Closest.Within$Facility.Tmp.ID],
             Closest.Within)
  setkey(M, Monitor.ID)
  setnames(M, "Longitude", "Longitude.Monitor")
  setnames(M, "Latitude", "Latitude.Monitor")
  M$Monitor.ID <- NULL
  M$Monitor.ID <- NULL
  M$Facility.Tmp.ID <- NULL
  setnames(M, make.names(colnames(M)))
  return(M)
}

#----- Link one month

linkMonitorFacility.Annual <- function(M) {
  #-----  Sort by monitor.  Flag = 1 if for closest facility
  M[, Include.in.Average := 0]
  setkeyv(M, c("Distance"))
  setkeyv(M, c("State.Code", "County.Code", "POC", "Site.Num"))
  M[M[, .I[1], by = key(M)]$V1, Include.in.Average := 1]
  #----- Sort by facility.  Include monitor in average if flag == 1
  setkeyv(M, c("Facility.ID..ORISPL."))
  M[, Average.Measure := mean(Arithmetic.Mean, na.rm = TRUE), by = c("Facility.ID..ORISPL.")]
  M[, Weighted.Measure := sum(Arithmetic.Mean * Include.in.Average) / sum(Include.in.Average),by = c("Facility.ID..ORISPL.")]
  #----- Descriptives
  M[, Number.of.Monitors.Within.Facility.Range := length(Include.in.Average), by = c("Facility.ID..ORISPL.")]
  M[, Number.of.Unique.Monitors.Within.Facility.Range := sum(Include.in.Average), by = c("Facility.ID..ORISPL.")]
  #------ Return sorted dataset
  setkeyv(M, c("Facility.ID..ORISPL."))
  return(M)
}


##------------------------------------------------------------------------------
## 'Main Linkage'
##------------------------------------------------------------------------------

#----- Create directories to store data

dir.create("Data_AQS", showWarnings = FALSE)

#-----   "annual_all_" # Annual Summary Data
#--- Note: AQS data for this paper was downloaded on Dec 2. 2015.  User can re-download AQS data,
#--- but data may be slightly different due to EPA updates to AQS.

# getAQSData.Annual(year = c(TC,TP))


#----- Create annual facility dataset

Fac <- fread("coal_unit_monthly.csv")
Fac <- subset(Fac, Fuel2.IsCoal!=0 | Fuel1.IsCoal!=0 )
setkeyv(Fac, "Facility.ID..ORISPL.")
Fac <- Fac[ID]

#----- Summer Period
# Fac <- subset(Fac, Month %in% c(5,6,7,8,9))
Fac <- subset(Fac, Year==TC)[ ,list(
        Facility.Name = Facility.Name.x,
        Facility.Latitude = Facility.Latitude.x,
        Facility.Longitude = Facility.Longitude.x,
        State = State.x,
        County = County.x,
        FIPS = FIPS,
        Total.Heat.Input..MMBtu. = sum(Heat.Input..MMBtu.,na.rm=TRUE),
        SO2.Annual = sum(SO2..tons., na.rm = TRUE),
        NOx.Annual = sum(NOx..tons., na.rm = TRUE),
        CO2.Annual = sum(CO2..short.tons., na.rm = TRUE)),
        by = c("Month","Year", "Facility.ID..ORISPL.")]

Fac.Annual <- Fac[, list(Facility.Name = Facility.Name,
                         Facility.Latitude = Facility.Latitude,
                         Facility.Longitude = Facility.Longitude,
                         State = State,
                         County = County,
                         FIPS = FIPS,
                         Total.Heat.Input..MMBtu. = mean(Total.Heat.Input..MMBtu.,na.rm=TRUE),
                         SO2.Annual = mean(SO2.Annual, na.rm = TRUE),
                         NOx.Annual = mean(NOx.Annual, na.rm = TRUE),
                         CO2.Annual = mean(CO2.Annual, na.rm = TRUE)),
                  by = c("Year", "Facility.ID..ORISPL.")]
setkeyv(Fac.Annual, c("Year", "Facility.ID..ORISPL."))
Fac.Annual.U <- unique(Fac.Annual)
setkeyv(Fac.Annual.U, c("Year", "Facility.Name"))

#------ PM25 subset

d <- loadAnnualAverage(year = TC)
ds <- subset(subset(d, Parameter.Code == 88101), Observation.Percent >= 67) # PM2.5 FRM Mass (88101) / More than 67% Observation percent

#------ Monitors from facility

# Distance "matrix" - Constant per year
DM <- distmatMonitorsFromFacility(Monitors = ds, Fac = Fac.Annual.U)
# Pairwise "within" (150km)
DMW <- withinMonitorsFromFacility(CDistmat = DM, withinMiles = 93.2057)
# Linked data
L <- linkMonitorFacility.Annual(M = DMW)


#------ Barometric Pressure subset (Year TP)

d <- loadAnnualAverage(year = TP)
ds <- subset(subset(d, Parameter.Code == 68108), Observation.Percent >= 67) # Sample Barometric Pressure


#----- Monitors from facility

# Distance "matrix" - Constant per year
DM <- distmatMonitorsFromFacility(Monitors = ds, Fac = Fac.Annual.U)
# Pairwise "within" (150km)
DMW <- withinMonitorsFromFacility(CDistmat = DM, withinMiles = 93.2057)
# Linked data
Lbaro <- linkMonitorFacility.Annual(M = DMW)


#------ Temperature subset (Year TP)

d <- loadAnnualAverage(year = TP)
ds <- subset(subset(d, Parameter.Code == 68105), Observation.Percent >= 67) # Ambient Temperature

#----- Monitors from facility

# Distance "matrix" - Constant per year
DM <- distmatMonitorsFromFacility(Monitors = ds, Fac = Fac.Annual.U)
# Pairwise "within" (150km)
DMW <- withinMonitorsFromFacility(CDistmat = DM, withinMiles = 93.2057)
# Linked data
Ltemp <- linkMonitorFacility.Annual(M = DMW)


#------ PM2.5 subset (Year TP)

d <- loadAnnualAverage(year = TP)
ds <- subset(subset(d, Parameter.Code == 88101), Observation.Percent >= 67) # PM2.5 FRM Mass (88101)

#----- Monitors from facility

# Distance "matrix" - Constant per year
DM <- distmatMonitorsFromFacility(Monitors = ds, Fac = Fac.Annual.U)
# Pairwise "within"
DMW <- withinMonitorsFromFacility(CDistmat = DM, withinMiles = 93.2057)
# Linked data
Lpm <- linkMonitorFacility.Annual(M = DMW)


##------------------------------------------------------------------------------
## 'Merging Step'
##------------------------------------------------------------------------------

#------ Merge PM2.5, Barometric Pressure, Temperature and Past PM2.5
L <- Reduce(function(x, y) merge(x, y, by = c("State.Code", "County.Code", "Site.Num", "Facility.ID..ORISPL.")),
        list(L, Lbaro[ ,list(State.Code, County.Code,
                             Site.Num, Facility.ID..ORISPL.,
                             baro = Weighted.Measure)],
                Ltemp[ ,list(State.Code, County.Code,
                             Site.Num, Facility.ID..ORISPL.,
                             temp = Weighted.Measure)],
                Lpm[ ,list(State.Code, County.Code,
                           Site.Num, Facility.ID..ORISPL.,
                           prePM2.5 = Weighted.Measure)]))

#------ Unique power plants data
setkey(L,"Facility.ID..ORISPL.")
L_sub <- L[ID, mult="first"][, SO2.SC := SO2]

#------ Discard power plants having no monitors
L_monitor <- subset(L_sub, !is.na(Number.of.Unique.Monitors.Within.Facility.Range) & Number.of.Unique.Monitors.Within.Facility.Range != 0)

#------ Merge L_monitor and all unit-level power plants data
L_all <- Reduce(function(x, y) merge(x, y, by="Facility.ID..ORISPL."),
            list(L_monitor, Fac.Annual.U, powerplant_sul,
                 powerplant_S_n_CR, powerplant_NumNOxControls,
                 powerplant_p2,powerplant_PctCapacity,
                 powerplant_time))

#------ Define variables
L_all_var <- L_all[ ,list(Facility.ID..ORISPL. = Facility.ID..ORISPL.,
                          SO2.SC = SO2.SC,
                          S_n_CR = S_n_CR,
                          NumNOxControls = NumNOxControls,
                          Heat_Input = Total.Heat.Input..MMBtu.,
                          Barometric_Pressure = baro,
                          Temperature = temp,
                          SO2_Annual = SO2.Annual,
                          NOx_Annual = NOx.Annual,
                          CO2_Annual = CO2.Annual,
                          PM.2.5 = Weighted.Measure,
                          Phase2_Indicator = Is.Phase2,
                          PctCapacity = PctCapacity,
                          Sulfur_Content = sulfur_content,
                          Operating_Time=Time)]

L_all_var <- subset(L_all_var, !is.na(Barometric_Pressure) & !is.na(Sulfur_Content))


#------ Master dataset

Master <- subset(L_all_var, SO2_Annual != 0) # no SO2 emission record(s)
index <- setdiff(Master$Facility.ID..ORISPL., c(1040, 2706, 3788))
Master <- subset(Master, Facility.ID..ORISPL. %in% index) # power plants installed scrubbers during 2005


save(Master,L_all,L, file = "Master.RData")


