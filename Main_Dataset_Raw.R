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

# dir.create("Data_AQS", showWarnings = FALSE)

#-----   "annual_all_" # Annual Summary Data
#--- Note: AQS data for this paper was downloaded on Dec 2. 2015.  User can re-download AQS data,
#--- but data may be slightly different due to EPA updates to AQS.

# getAQSData.Annual(year = c(TC,TP))


#----- Create annual facility dataset

Fac <- fread("coal_unit_monthly.csv")
Fac <- subset(Fac, Fuel2.IsCoal!=0 | Fuel1.IsCoal!=0 )
setkeyv(Fac, "Facility.ID..ORISPL.")

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
        by = c("Month","Year", "Facility.ID..ORISPL.","Unit.ID")]

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
                  by = c("Year", "Facility.ID..ORISPL.","Unit.ID")]
setkeyv(Fac.Annual, c("Year", "Facility.ID..ORISPL.","Unit.ID"))
Fac.Annual.U <- unique(Fac.Annual)

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
ds <- subset(subset(d, Parameter.Code == 68105), Observation.Percent >= 67) # Ambient Temperature

#----- Monitors from facility

# Distance "matrix" - Constant per year
DM <- distmatMonitorsFromFacility(Monitors = ds, Fac = Fac.Annual.U)
# Pairwise "within" (150km)
DMW <- withinMonitorsFromFacility(CDistmat = DM, withinMiles = 93.2057)
# Linked data
Ltemp <- linkMonitorFacility.Annual(M = DMW)


#------ PM2.5 subset (Year TP)
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


#------ Merge L and unit-level power plants data
L_all <- merge(L, Fac.Annual.U , by="Facility.ID..ORISPL.", allow.cartesian=TRUE)


save(L_all, file = "RawData.RData")


