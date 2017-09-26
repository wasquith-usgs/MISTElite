MISTE <- new.env()
graphics.off()
library(chron)
suppressPackageStartupMessages(library(mgcv))
library(dataRetrieval)

# GLOBAL VARIABLES
zeroValue <- 0.001
messVerb  <- TRUE
devVerb   <- FALSE
my.purple <- "#9b42f4"
my.red2   <- "#f48241"

# MISTE INTERFACE
misteHA <- function(...) h4(...)
misteHB <- function(...) h5(...)
misteHC <- function(...) h6(...)

# MISTE DATA INTERCONNECTION
misteUVget <- function(siteNumber, pCode="00060", sdate="", edate="", stime="", etime="",
                       message="", storestring="", earlyexit=FALSE, debug=FALSE, ...) {
   sdatetime <- paste0(sdate,"T",stime); sdt <- paste0(sdate," ",stime)
   edatetime <- paste0(edate,"T",etime); edt <- paste0(edate," ",etime)
   if(messVerb) message("MISTElite getting ",message," for\n",
                        "                   ",siteNumber,
                         " (",as.character(format(sdt, usetz=TRUE)), " LOC to ",
                              as.character(format(edt, usetz=TRUE)), " LOC)")
   rawData <- dataRetrieval::readNWISuv(siteNumber, parameterCd=pCode,
                                        startDate=sdatetime, endDate=edatetime, tz="", ...)
   rawData <- dataRetrieval::renameNWISColumns(rawData)
   if(earlyexit) return(rawData) # for debugging.
   assign(paste0("rawData","_",storestring), rawData, envir=MISTE) # Archive data
   sdt <- rawData$dateTime[1]; edt <- rawData$dateTime[length(rawData$dateTime)]
   if(messVerb) message("MISTElite got ",message," for\n",
                        "                   ",siteNumber,
                              " (",as.character(format(sdt, usetz=TRUE))," to ",
                                   as.character(format(edt, usetz=TRUE)),")")
   return(rawData)
}
# TODO, need to research the behavior of the TZ
#UVx <- misteUVget("08167000", sdate="2013-10-01", stime="02:00",edate="2013-10-02", etime="05:00")
#UVy <- misteUVget("08167500", sdate="2013-10-01", stime="02:00",edate="2013-10-02", etime="05:00")

misteUVmerge <- function(x, y, core.origin="2005-01-01", ...) {
   if(messVerb) message("MISTElite merging the data retrievals (predictor / response)")
   # The use of an string of the date-time have the timezone ("UTC") should be superfluous, but WHA
   # is trying to get rid of CST/CDT in the "Product Graphic" and output
   x$dateTimeSTRING <- as.character(format(x$dateTime, usetz=TRUE))
   y$dateTimeSTRING <- as.character(format(y$dateTime, usetz=TRUE))
   UVs <- merge(x, y, by="dateTimeSTRING", all=TRUE)
   #print(head(UVs))
   dt <- unlist(strsplit(as.character(UVs$dateTimeSTRING), "\\s+"))
   #print(head(dt))
   ix <- seq(1,length(dt)-2, by=3)
   ymd <- dt[ix]; hms <- dt[ix+1]; tz <- dt[ix+2]
   POSIXepoch <- UVs$dateTime.x - strptime(core.origin,
                                           format="%Y-%d-%m", tz="UTC")
   cObj <- chron(ymd, hms, format=c(dates="y-m-d", times="h:m:s"))
   CHRONepoch <- cObj - chron(core.origin, "00:00:00",
                              format=c(dates="y-m-d", times="h:m:s"))
   POSIXepoch <- as.integer(POSIXepoch*(24*60*60))
   CHRONepoch <- as.integer(CHRONepoch*(24*60*60))
   POSIXdiff <- diff(POSIXepoch);      CHRONdiff <- diff(CHRONepoch)
   delPosix <- sort(unique(POSIXdiff)); delChron <- sort(unique(CHRONdiff))
   if(messVerb) message("MISTElite time gap analysis (POSIX time math, seconds): ",
           paste(delPosix, collapse=","))
   if(messVerb) message("                             (CHRON time math, seconds): ",
           paste(delChron, collapse=","))
   GlobalDeltaTime <- delPosix[1] # TODO, unverified logic
   if(messVerb) message("                         delta time stamp: ",
           GlobalDeltaTime, " seconds")
   assign("GlobalDeltaTime", GlobalDeltaTime, envir=MISTE)
   XY  <- data.frame(dateTimeSTRING=UVs$dateTimeSTRING,
                     dateTime=UVs$dateTime.x, tz_cd=UVs$tz_cd.x,
                     uvPOSIXepoch=POSIXepoch, chronObject=cObj,
                     uvCHRONepoch=CHRONepoch,
                     predictor=UVs$Flow_Inst.x, response=UVs$Flow_Inst.y)
   XY$dateTimeSTRING <- as.character(XY$dateTimeSTRING)
   assign("datamergeXY", XY, envir=MISTE)
   write.table(XY, file="./vartmp/MISTE_regress_datatable.txt",
                   row.names=FALSE, quote=FALSE, sep="\t")
   return(XY)
}
#misteUVmerge(UVx, UVy)


misteDateRange2Dates <- function(shiny_dates, ...) {
   #print(dates)
   tmpA <- unlist(strsplit(as.character(shiny_dates[1]), "-"))
   tmpB <- unlist(strsplit(as.character(shiny_dates[2]), "-"))
   #print(tmpA); print(tmpB)
   sdate = list(sdate=paste(tmpA, collapse="-"), yyyy=tmpA[1], mm=tmpA[2], dd=tmpA[3])
   edate = list(edate=paste(tmpB, collapse="-"), yyyy=tmpB[1], mm=tmpB[2], dd=tmpB[3])
   return(list(sdate=sdate, edate=edate))
}

misteHHMM2time <- function(hhmm) {
   if(messVerb) message("MISTElite parsing time using misteHHMM2time() on string '", hhmm, "'")

   HHMM <- unlist(strsplit(hhmm, "tz")) # optional trailing timezone offset
   if(length(HHMM) == 1) {
      TZ <- "" # note the leading + sign
   } else {
      TZ <- HHMM[2]
      if(TZ > 0) TZ <- paste0("+",TZ) # convert to a string with leading + sign
   }
   hhmm <- HHMM[1]
   HHMM <- unlist(strsplit(hhmm, "to"))
   if(length(HHMM) == 1) HHMM[2] <- "00:00" # set trailing default 00:00 hours
                                            # for ending time
   if(length(HHMM) == 2) {
      time1 <- unlist(strsplit(HHMM[1], ":"))
      time2 <- unlist(strsplit(HHMM[2], ":"))
      if(any(as.numeric(c(time1[1], time2[1])) <  0) |
         any(as.numeric(c(time1[1], time2[1])) > 24)) {
         warning("erroneous hour, stripping out HH:MM and continuing")
         return(list(timeparsed=FALSE, stime="",   etime="",
                                       shhmmss="", ehhmmss=""))
      }
      if(any(as.numeric(c(time1[2], time2[2])) <  0) |
         any(as.numeric(c(time1[2], time2[2])) > 59)) {
         warning("erroneous minute, stripping out HH:MM and continuing")
         return(list(timeparsed=FALSE, stime="",   etime="",
                                       shhmmss="", ehhmmss=""))
      }
      shhmmss <- paste0(time1[1],":", time1[2],":","00")
      ehhmmss <- paste0(time2[1],":", time2[2],":","00")
      return(list(timeparsed=TRUE,
                  stime=paste0(shhmmss,TZ), shhmmss=shhmmss,
                  etime=paste0(ehhmmss,TZ), ehhmmss=ehhmmss,
                           hh1=time1[1], mm1=time1[2], ss1="00",
                           hh2=time2[1], mm2=time2[2], ss2="00",
                  tz=TZ))
   }
   warning("time range parsing failure, stripping out HH:MM and continuing")
   return(list(timeparsed=FALSE, stime="", etime="", shhmmss="", ehhmmss=""))
}


# MISTE MODEL METHODS
misteMODELpreprocess <- function(xy, input, ...) {
   #print(c(input$site1.lower.cutout, input$site1.upper.cutout))
   #print(c(input$site2.lower.cutout, input$site2.upper.cutout))
   if(messVerb) message("MISTElite preprocessing data")
   if(messVerb) message("           prior to discharge cuts, sample size n=",
           length(xy$response[! is.na(xy$response) & ! is.na(xy$predictor)]))
   xy$predictor[xy$predictor < input$site1.lower.cutout |
                xy$predictor > input$site1.upper.cutout ] <- NA
   if(messVerb) message("              after predictor cuts, sample size n=",
           length(xy$response[! is.na(xy$response) & ! is.na(xy$predictor)]))
   xy$response[xy$response   < input$site2.lower.cutout |
               xy$response   > input$site2.upper.cutout ] <- NA
   if(messVerb) message("              after response cuts,  sample size n=",
           length(xy$response[! is.na(xy$response) & ! is.na(xy$predictor)]))
   corlag <- input$correlation.lag


   seconds.to.hours <- 60*60
   corlag <- seconds.to.hours*corlag
   # print(sum(xy$predictor, na.rm=TRUE))
   # NOTE NOTE NOTE The addition of the corlag to the time link
   A <- data.frame(linkinseconds=xy$uvPOSIXepoch+corlag, xy$predictor,
                   dateTime=xy$dateTime)
   B <- data.frame(linkinseconds=xy$uvPOSIXepoch,        xy$response,
                   dateTime=xy$dateTime)
   AB <- merge(A, B, by="linkinseconds", all=TRUE)
   AB <- AB[complete.cases(AB), ]
   #print(sum(AB$xy.predictor, na.rm=TRUE))
   AB$predictor    <- AB$xy.predictor;  AB$response     <- AB$xy.response
   AB$dateTime     <- AB$dateTime.y;    AB$uvPOSIXepoch <- AB$linkinseconds
   AB$xy.predictor <- AB$xy.response <- NULL
   assign("true.time.matched.XY", xy, envir=MISTE)
   assign("XY", AB, envir=MISTE)
   return(AB)
}


misteMODEL_ols <- function(X,Y, input, ...) {
   isXlog <- input$cb.log.site1
   isYlog <- input$cb.log.site2
   poly.order <- input$poly.order

   # separate calls lm() so log file shows the transforms in the "call" section
   if(isXlog & isYlog) {
      X[X == 0] <- NA; Y[Y == 0] <- NA
      LM <- lm(log10(Y)~poly(log10(X), degree=poly.order, raw=TRUE))
   } else if(isXlog) {
      X[X == 0] <- NA
      LM <- lm(Y~poly(log10(X),        degree=poly.order, raw=TRUE))
   } else if(isYlog) {
      Y[Y == 0] <- NA
      LM <- lm(log10(Y)~poly(X,        degree=poly.order, raw=TRUE))
   } else {
      LM <- lm(Y~poly(X,               degree=poly.order, raw=TRUE))
   }
   return(LM)
}


misteMODEL_gam <- function(X,Y, input, ...) {
   isXlog <- input$cb.log.site1
   isYlog <- input$cb.log.site2

   # separate calls lm() so log file shows the transforms in the "call" section
   if(isXlog & isYlog) {
      X[X == 0] <- NA; Y[Y == 0] <- NA
      LM <- gam(log10(Y)~s(log10(X)), ...)
   } else if(isXlog) {
      X[X == 0] <- NA
      LM <- gam(Y~s(log10(X)), ...)
   } else if(isYlog) {
      Y[Y == 0] <- NA
      LM <- gam(log10(Y)~s(X), ...)
   } else {
      LM <- gam(Y~s(X), ...)
   }
   return(LM)
}


misteMODEL <- function(xy, input, ...) {
   XY <- misteMODELpreprocess(xy, input)
   if(messVerb) message("MISTElite making model ", appendLF=FALSE)
   XY <- XY[complete.cases(XY), ]; X <- XY$predictor; Y <- XY$response
   misteMODELminimum.data <- 10
   if(length(X) <= misteMODELminimum.data) {
      if(messVerb) message("** MISTElite ERROR: insufficient predictor data for MISTElite model")
      return(NA)
   }
   if(length(Y) <= misteMODELminimum.data) {
      if(messVerb) message("** MISTElite ERROR: insufficient response data for MISTElite model")
      return(NA)
   }

   isYlog <- input$cb.log.site2; isGAM <- input$cb.modelgam
   ifelse(isGAM, LM <- misteMODEL_gam(X,Y, input), LM <- misteMODEL_ols(X,Y, input))

   # Difficulty in alternative model construction is the logic flow
   # post model building.  Need a predict.lm(), a summary.lm(), print.lm(),
   # residuals.lm() equivalent functions on the alternative

   swathX <- diff(range(X, na.rm=TRUE)) / 200
   swathX <- sort(unique(c(range(X, na.rm=TRUE),
                         seq(min(X, na.rm=TRUE) / 100,
                             max(X, na.rm=TRUE) * 100, by=swathX))))
   assign("swathX", swathX, envir=MISTE)
   swathY <- predict(LM, newdata=data.frame(X=swathX)) # R lm() or gam() dependency
   if(isYlog) swathY <- 10^swathY
   assign("swathY", swathY, envir=MISTE)
   assign("LM", LM, envir=MISTE)
   sumLM <- summary(LM)
   ifelse(isGAM, assign("RSE", sqrt(sumLM$scale),        envir=MISTE),
                 assign("RSE",      sumLM$sigma,         envir=MISTE))
   ifelse(isGAM, assign("RSQ",      sumLM$r.sq,          envir=MISTE),
                 assign("RSQ",      sumLM$adj.r.squared, envir=MISTE))

   duan <- ifelse(isYlog, mean(10^residuals(LM)), 1) # R lm() or gam() dependency
   assign("duansmearing", duan, envir=MISTE)
   if(messVerb) message("--- done")
   return(LM)
}


# Need to verify that predict.lm() handles the log10 on the X on its own. What
# about the issue of negative predictions if not using log10(Y)?
misteMODELpredict <- function(data, input, ...) {
   isYlog <- input$cb.log.site2; isGAM <- input$cb.modelgam
   seconds.to.hours <- 60*60
   corlag <- seconds.to.hours*input$correlation.lag

   # There are two times dealt with in the next lines. This might be a legacy
   # code issue and a special accommodation for "NWIS" time need not be
   # necessary.
   site2UVpredicted <- data
   dateTime_atSite1 <- site2UVpredicted$dateTime
   dateTime_atSite2 <-                  dateTime_atSite1 + corlag
   site2UVpredicted$dateTime <- dateTime_atSite2
   site2UVpredicted$dateTime_atSite1 <-
                             as.character(format(dateTime_atSite1, usetz=TRUE))
   site2UVpredicted$dateTime_atSite2 <-
                             as.character(format(dateTime_atSite2, usetz=TRUE))

   LM <- get("LM", envir=MISTE)
   X <- data$Flow_Inst
   if(length(X) == 0) {
      if(messVerb) message("** MISTElite ERROR: insufficient data for making predictions")
      return(NA)
   }

   if(isGAM) {
      Y <- predict(LM, newdata=data.frame(X=X), se.fit=TRUE) # R gam() dependency
      Y <- as.data.frame(Y)
      sigma <- sqrt(summary(LM)$scale)
      sqrt1hat <- sqrt(1+(Y$se.fit/sigma)^2) # TODO CHECK ^2!!!!
      xx <- abs(qt((1-input$prediction.interval)/2, 100))*sigma*sqrt1hat
      Y$lwr <- Y$fit - xx; Y$upr <- Y$fit + xx
   } else {
      Y <- predict(LM, newdata=data.frame(X=X),  # R lm() dependency
                   interval="prediction", level=input$prediction.interval)
      Y <- as.data.frame(Y)
   }
   if(isYlog) {
     Y <- 10^Y
     if(input$cb.duansmearing) {
        duan <- get("duansmearing", envir=MISTE)
        Y <- duan*Y
     }
     Y <- log10(Y) # transform back for possible shape transformation
     Y <- as.data.frame(Y) # this need seems spurious but appears necessary to recast
   }

   if(input$cb.shape.transform) { # ShapeTransformation
       if(messVerb) message("        *** triggered with shape transformation ***")
       inQ  <- input$target.entry.discharge
       outQ <- input$target.exit.discharge
       if(isYlog) {
          if(inQ  == 0)  inQ <- zeroValue; inQ  <- log10(inQ)
          if(outQ == 0) outQ <- zeroValue; outQ <- log10(outQ)
       }
       if(is.finite(inQ) & is.finite(outQ)) {
          #if(messVerb) message("      inQ=", inQ, " and outQ=",outQ)
          DTin  <- dateTime_atSite2[1]
          DTout <- dateTime_atSite2[length(dateTime_atSite2)]
          #if(messVerb) message("      DTin=", DTin, " and DTout=",DTout)
          inOffset  <- inQ  - Y$fit[1]
          outOffset <- outQ - Y$fit[length(Y$fit)]
          totOffset <- outOffset - inOffset
          #if(messVerb) message("      inOffset=",inOffset, " and outOffset",outOffset)
          transform.baseline <- as.numeric(DTout - DTin)
          #if(messVerb) message("      transform.baseline=",transform.baseline)
          delTimes <- diff(as.numeric(dateTime_atSite2))
          # if not case as numeric, then would have to code for attr(MISTE$TMP,"units")
          totalTime <- sum(delTimes); fracTime <- cumsum(c(0, delTimes))/totalTime
          reshapeYfit <- Y$fit + totOffset*fracTime + inOffset
          reshapeYlwr <- Y$lwr + totOffset*fracTime + inOffset
          reshapeYupr <- Y$upr + totOffset*fracTime + inOffset
          if(isYlog) {
             reshapeYfit[reshapeYfit <= 0] <- zeroValue
             reshapeYlwr[reshapeYlwr <= 0] <- zeroValue
             reshapeYupr[reshapeYupr <= 0] <- zeroValue
          }
          reshapeY <- data.frame(fit=reshapeYfit, lwr=reshapeYlwr, upr=reshapeYupr)
          Y <- reshapeY
       } else {
          warning(" shape transformation desired, but one of the entry/exit cuts is Inf")
       }
   }

   if(isYlog) Y <- 10^Y
   Y <- as.data.frame(round(Y, digits=2))
   site2UVpredicted$site_no <- rep(input$site2, length(X))
   site2UVpredicted$Flow_Inst     <- Y$fit
   site2UVpredicted$Flow_Inst_lwr <- Y$lwr
   site2UVpredicted$Flow_Inst_upr <- Y$upr
   site2UVpredicted$Flow_Inst_cd  <- rep("MISTE", length(X))

   file <- "./vartmp/MISTE_site2_predict(DECODES).csv"
   write.table(site2UVpredicted, file=file, sep=",",
               row.names=FALSE, quote=FALSE)
   if(messVerb) message("MISTElite making predictions in file    '", file, "'")

   assign("predY",    Y$fit,                         envir=MISTE)
   assign("predYlwr", Y$lwr,                         envir=MISTE)
   assign("predYupr", Y$upr,                         envir=MISTE)
   assign("predX",    X,                             envir=MISTE)
   assign("pred.dateTime_atSite1", dateTime_atSite1, envir=MISTE)
   assign("pred.dateTime_atSite2", dateTime_atSite2, envir=MISTE)
}


"Do_misteModel_and_PredictNeed_overlap" <- function() {
     XY <- get("datamergeXY",         envir=MISTE)
     OT <- get("rawData_outputSite1", envir=MISTE)
     OT$dateTimeSTRING <- as.character(format(OT$dateTime, usetz=TRUE))
     num_overlap <- length(intersect(XY$dateTimeSTRING, OT$dateTimeSTRING))
     if(num_overlap > 0) {
        message("MISTElite EXCEPTION: model time overlaps with prediction ",
                "time, ",num_overlap, " time stamps overlap!")
        return(TRUE)
     }
     return(FALSE)
}


# MISTE TEXTUAL OUTPUT
misteMODELlogsummary <- function(sumLM,
                         logfile="./vartmp/MISTE_model_summary.txt", ...) {
  if(messVerb) message("MISTElite sinking model summary to file '",
           logfile, "'")
  sink(file=logfile)
    print(sumLM) # R lm() dependency
  sink()
}

misteMODELdiagnosticText <- function(input, ...) {
  LM <- get("LM", envir=MISTE)
  sumLM <- summary(LM) # R lm() dependency
  RSE <- get("RSE", envir=MISTE)
  RSQ <- get("RSQ", envir=MISTE)
  misteMODELlogsummary(sumLM)
  duan <- get("duansmearing", envir=MISTE)
  #isGAM  <- input$cb.modelgam
     text <- paste0("MODEL SUMMARY: ",
                    "Rsq=",round(RSQ, digits=3), ", ",
                    "RSE=", round(RSE, digits=4)
                   )
  if(input$cb.duansmearing & input$cb.log.site2) {
    text <- paste0(text, ", Duan bias-correction factor=", round(duan, digits=3))
  } else {
    text <- paste0(text, "") # hook left if extra info is ever desired
  }

  residualtext <- paste0("RESIDUAL SUMMARY: ",
                  paste(round(summary(residuals(LM)), digits=4), collapse=", "),
                     " : (min, 1Q, Median, 3Q, max)")
  return(c(text,residualtext))
}

# MISTE GRAPHICAL OUTPUT
mistePDF <- function(file, input, ...) {
   pdf(file=file, useDingbats=FALSE, width=7.5, height=6.5, ...)
}

miste.mtext <- function(text, cex=0.8, ...) {
   mtext(text, cex=cex, ...)
}

misteLogarithms <- function(x, islog=FALSE, truncate=0.001, ...) {
  if(islog) {
     if(truncate > 0) x[x == 0] <- truncate
     return(log10(x))
  }
  return(x)
}

#X <- misteLogarithms(X, islog=input$cb.log.site1)

misteTimeSeriesPlot_Inputs <- function(input, pdfoutput=FALSE, ...) {
   isXlog <- input$cb.log.site1
   isYlog <- input$cb.log.site2
   isXlog <- isYlog # substitution is made to ensure that the vertical axis
   # data are in transformation consitent with Y, which is what we are trying
   # to predict and not have logged values trying to be plotted on a linear axis

   LC1 <- input$site1.lower.cutout; UC1 <- input$site1.upper.cutout
   LC2 <- input$site2.lower.cutout; UC2 <- input$site2.upper.cutout

   XY <- get("datamergeXY", envir=MISTE)
   D <- XY$dateTime; X <- XY$predictor; Y <- XY$response
   XYm <- get("XY", envir=MISTE)
   Dm <- XYm$dateTime; Xm <- XYm$predictor; Ym <- XYm$response

   xlab="CALENDAR TIME"
   rngD <- range(D, na.rm=TRUE)
   DT1  <- unlist(strsplit(as.character(format(rngD[1], usetz=TRUE)), " "))
   DT2  <- unlist(strsplit(as.character(format(rngD[2], usetz=TRUE)), " "))
   if(length(DT1) == 2) DT1[2:3] <- c("00:00:00", DT1[2])
   if(length(DT2) == 2) DT2[2:3] <- c("00:00:00", DT2[2])
   TZ1  <- DT1[3]; TZ2 <- DT2[3]
   YMD1 <- unlist(strsplit(as.character(DT1[1]), "-"))
   YMD2 <- unlist(strsplit(as.character(DT2[1]), "-"))
   #print(YMD1)
   Y1 <- YMD1[1]; Y2 <- YMD2[1]
   M1 <- YMD1[2]; M2 <- YMD2[2]
   D1 <- YMD1[3]; D2 <- YMD2[3]
   Y1 <- as.numeric(Y1); Y2 <- as.numeric(Y2)
   # Adapt the labeling to ensure user 'knows' what years are being shown
   # on the plot.
   sxlab <- paste0(Y1,"-",M1,"-",D1)
   exlab <- paste0(Y2,"-",M2,"-",D2)
   M1 <- as.numeric(M1); M2 <- as.numeric(M2) # must do this AFTER sxlab, exlab are started!
   if(Y2 - Y1 == 0 & M2 - M1 == 0) {
      sxlab <- paste0(sxlab, " ", DT1[2], " ", TZ1)
      exlab <- paste0(exlab, " ", DT2[2], " ", TZ2)
   }
   xlab <- paste0(xlab, "(", sxlab, " to ", exlab, ")")

   # Setup for grid lines on the first of the month
   dogrid <- ifelse(Y2 - Y1 < 7, TRUE, FALSE)
   YYs <- seq(Y1, Y2, by=1); MMs <- 1:12
   i <- 0; gridline.dates <- gridline.widths <- vector()
   for(yy in YYs) {
     for(mm in MMs) {
       i <- i + 1
       gridline.widths[i] <- ifelse(mm == 1, 1.3, 0.6)
       gridline.dates[i] <- as.Date(paste(yy,"-",mm,"-","01", sep=""))
     }
   }

   ylab <- "STREAMFLOW, IN CFS"
   if(isXlog) {
      X[X == 0] <- zeroValue
      X <- log10(X)
      Xm[Xm == 0] <- zeroValue
      Xm <- log10(Xm)
      LC1 <- log10(LC1); if(! is.finite(LC1)) LC1 <- NA
      UC1 <- log10(UC1); if(! is.finite(UC1)) UC1 <- NA
   }
   if(isYlog) {
      Y[Y == 0] <- zeroValue
      Y <- log10(Y)
      Ym[Ym == 0] <- zeroValue
      Ym <- log10(Ym)
      ylab <- "STREAMFLOW, IN log10(CFS)"
      LC2 <- log10(LC2); if(! is.finite(LC2)) LC2 <- NA
      UC2 <- log10(UC2); if(! is.finite(UC2)) UC2 <- NA
   }
   afunc <- function() {
      plot(c(D,NA,D), c(X,NA,Y), type="n",
           xlab=xlab, ylab=ylab, xaxs="r", las=1, tcl=.7, mgp=c(2,0.5,0))
      if(dogrid) {
        usr <- par()$usr
        for(i in 1:length(gridline.dates)) {
          if(usr[1] == as.numeric(gridline.dates[i]) |
             usr[2] == as.numeric(gridline.dates[i])) next
          lines(rep(gridline.dates[i], 2), usr[3:4],
                                      lty=2, col=8, lwd=gridline.widths[i])
        }
      }

      k <- length(Dm)
      lines(rep(Dm[1], 2), par()$usr[3:4], lty=3, col=my.red2, lwd=0.76)
      lines(rep(Dm[k], 2), par()$usr[3:4], lty=3, col=my.red2, lwd=0.76)
      lines(range(D, na.rm=TRUE), rep(LC1, 2), lwd=1.35, col=my.red2, lty=3)
      lines(range(D, na.rm=TRUE), rep(UC1, 2), lwd=1.35, col=my.red2, lty=4)
      lines(range(D, na.rm=TRUE), rep(LC2, 2), lwd=1.35, col=4, lty=3)
      lines(range(D, na.rm=TRUE), rep(UC2, 2), lwd=1.35, col=4, lty=4)
      lines(D,X,   lwd=1.6, col=2, lty=2)
      lines(D,Y,   lwd=1.6, col=4)
      lines(Dm,Xm, lwd=1.6, col=my.red2)
      miste.mtext("INPUT (site1=red, site2=blue, dash='site1 lagged')")
   }
   if(pdfoutput) {
     file <- "./vartmp/MISTE_timeseries_input.pdf"
     if(messVerb) message("MISTElite making input time series plot '", file, "'")
     mistePDF(file, input); afunc(); dev.off()
   } else {
     afunc()
   }
}


misteTimeSeriesPlot_Product <- function(input, pdfoutput=FALSE, ...) {
   isXlog <- input$cb.log.site1
   isYlog <- input$cb.log.site2
   isXlog <- isYlog # substitution is made to ensure that the vertical axis
   # data are in transformation consitent with Y, which is what we are trying
   # to predict and not have logged values trying to be plotted on a linear axis

   inQ  <- input$target.entry.discharge
   outQ <- input$target.exit.discharge
   if(isYlog) {
      if(inQ  == 0)  inQ <- zeroValue; inQ  <- log10(inQ)
      if(outQ == 0) outQ <- zeroValue; outQ <- log10(outQ)
   }

   LC1 <- input$site1.lower.cutout; UC1 <- input$site1.upper.cutout
   LC2 <- input$site2.lower.cutout; UC2 <- input$site2.upper.cutout

   DTsite1 <- get("pred.dateTime_atSite1", envir=MISTE)
   X       <- get("predX",                 envir=MISTE)
   Y       <- get("predY",                 envir=MISTE)
   Ylo     <- get("predYlwr",              envir=MISTE)
   Yhi     <- get("predYupr",              envir=MISTE)
   DTsite2 <- get("pred.dateTime_atSite2", envir=MISTE)

   BEFORE_UVs <- get("augmentbeforeSite2", envir=MISTE)
   AFTER_UVs  <- get("augmentafterSite2",  envir=MISTE)
   bDT <- BEFORE_UVs$dateTime;  aDT <- AFTER_UVs$dateTime
   bQ  <- BEFORE_UVs$Flow_Inst;  aQ <- AFTER_UVs$Flow_Inst

   #print("misteTimeSeriesPlot_Product() DTsite1")
   #print(head(DTsite1))
   #print("misteTimeSeriesPlot_Product() DTsite2")
   #print(head(DTsite2))

   DTcombo <- DTsite1
   # TODO, this next operation converts to LOC time zone
   #DTcombo <- c(DTsite1, DTsite2) # but succinct logic of desires

   xlab="DATE-TIME and at site2:"
   rngD  <- range(DTsite1, na.rm=TRUE)
   rngD1 <- range(DTsite1, na.rm=TRUE)
   rngD2 <- range(DTsite2, na.rm=TRUE)
   #print(rngD); print(rngD1); print(rngD2) # These show UTC

   # These move to seconds, it does make the plot window work
   xmin <- ifelse(rngD1[1] < rngD2[1], rngD1[1], rngD2[1])
   xmax <- ifelse(rngD1[2] > rngD2[2], rngD1[2], rngD2[2])
   xlim <- c(xmin, xmax)+(c(-1,1)*rep(24*60*60,2)) # one-day left/right buffer

   DT1 <- unlist(strsplit(as.character(format(rngD[1], usetz=TRUE)), " "))
   DT2 <- unlist(strsplit(as.character(format(rngD[2], usetz=TRUE)), " "))
   if(length(DT1) == 2) DT1[2:3] <- c("00:00:00", DT1[2])
   if(length(DT2) == 2) DT2[2:3] <- c("00:00:00", DT2[2])
   TZ1 <- DT1[3]; TZ2 <- DT2[3]
   YMD1 <- unlist(strsplit(as.character(DT1[1]), "-"))
   YMD2 <- unlist(strsplit(as.character(DT2[1]), "-"))
   #print(YMD1)
   Y1 <- YMD1[1]; Y2 <- YMD2[1]
   M1 <- YMD1[2]; M2 <- YMD2[2]
   D1 <- YMD1[3]; D2 <- YMD2[3]
   Y1 <- as.numeric(Y1); Y2 <- as.numeric(Y2)
   # Adapt the labeling to ensure user 'knows' what years are being shown
   # on the plot.
   sxlab <- paste0(Y1,"-",M1,"-",D1)
   exlab <- paste0(Y2,"-",M2,"-",D2)
   M1 <- as.numeric(M1)
   M2 <- as.numeric(M2) # must do this AFTER sxlab, exlab are started!
   if(Y2 - Y1 == 0 & M2 - M1 == 0) {
      sxlab <- paste0(sxlab, " ", DT1[2], " ", TZ1)
      exlab <- paste0(exlab, " ", DT2[2], " ", TZ2)
   }
   xlab <- paste0(xlab, "(", sxlab, " to ", exlab, ")")
   # Setup for grid lines on the first of the month
   dogrid <- ifelse(Y2 - Y1 < 7, TRUE, FALSE)
   YYs <- seq(Y1, Y2, by=1); MMs <- 1:12
   i <- 0; gridline.dates <- gridline.widths <- vector()
   for(yy in YYs) {
     for(mm in MMs) {
       i <- i + 1
       gridline.widths[i] <- ifelse(mm == 1, 1.3, 0.6)
       gridline.dates[i] <- as.Date(paste(yy,"-",mm,"-","01", sep=""))
     }
   }

   ylab <- "STREAMFLOW, IN CFS"
   if(isXlog) {
      X[X == 0] <- zeroValue; X <- log10(X)
      LC1 <- log10(LC1); if(! is.finite(LC1)) LC1 <- NA
      UC1 <- log10(UC1); if(! is.finite(UC1)) UC1 <- NA
   }
   if(isYlog) {
      Y[Y == 0]     <- zeroValue;   Y <- log10(Y  )
      Ylo[Ylo == 0] <- zeroValue; Ylo <- log10(Ylo)
      Yhi[Yhi == 0] <- zeroValue; Yhi <- log10(Yhi)
      bQ[bQ == 0]   <- zeroValue;  bQ <- log10(bQ )
      aQ[aQ == 0]   <- zeroValue;  aQ <- log10(aQ )
      ylab <- "STREAMFLOW, IN log10(CFS)"
      LC2 <- log10(LC2); if(! is.finite(LC2)) LC2 <- NA
      UC2 <- log10(UC2); if(! is.finite(UC2)) UC2 <- NA
   }
   afunc <- function() {
      Haxdat <- c(DTsite1,NA,DTsite2,  NA,DTsite2,   NA,DTsite2)
      if(input$prediction.interval != 0 & input$prediction.interval != 1) {
        Vaxdat <- c(X,NA,Y,  NA,Ylo, NA,Yhi)
      } else {
        Vaxdat <- c(X,NA,Y)
      }
      plot(DTsite2, Y, type="n", xlim=xlim, las=1, tcl=.7, mgp=c(2,0.5,0),
           xlab=xlab, ylab=ylab, xaxs="r", ylim=range(Vaxdat, na.rm=TRUE))
      if(dogrid) {
        usr <- par()$usr
        for(i in 1:length(gridline.dates)) {
          if(usr[1] == as.numeric(gridline.dates[i]) |
             usr[2] == as.numeric(gridline.dates[i])) next
          lines(rep(gridline.dates[i], 2), usr[3:4],
                                      lty=2, col=8, lwd=gridline.widths[i])
        }
      }
      lines(bDT, bQ, lty=4, col=4, lwd=0.76)
      lines(aDT, aQ, lty=4, col=4, lwd=0.76)
      k <- length(DTsite1); j <- length(DTsite2)
      lines(rep(DTsite1[1], 2), par()$usr[3:4], lty=3, col=my.red2, lwd=0.76)
      lines(rep(DTsite1[k], 2), par()$usr[3:4], lty=3, col=my.red2, lwd=0.76)
      lines(rep(DTsite2[1], 2), par()$usr[3:4], lty=3, col=4, lwd=0.76)
      lines(rep(DTsite2[j], 2), par()$usr[3:4], lty=3, col=4, lwd=0.76)
      lines(rngD, rep(LC1, 2), lwd=1.35, col=my.red2, lty=3)
      lines(rngD, rep(UC1, 2), lwd=1.35, col=my.red2, lty=4)
      lines(rngD, rep(LC2, 2), lwd=1.35, col=4, lty=3)
      lines(rngD, rep(UC2, 2), lwd=1.35, col=4, lty=4)
      lines(DTsite1, X,   lwd=1.6, col=my.red2)
      lines(DTsite2, Y,   lwd=1.6, col=4)
      if(input$prediction.interval != 0 & input$prediction.interval != 1) {
         lines(DTsite2, Ylo, lwd=1.4, col=4, lty=2)
         lines(DTsite2, Yhi, lwd=1.4, col=4, lty=2)
      }
      miste.mtext("MISTElite'd (site1=red, site2=blue, dash='predict interval')")
      if(input$cb.shape.transform) { # ShapeTransformation
        points(DTsite2[1],                inQ, pch=8, col="#3262c1", cex=2, lwd=.89)
        points(DTsite2[length(DTsite2)], outQ, pch=8, col="#3262c1", cex=2, lwd=.89)
      }
   }
   if(pdfoutput) {
     file <- "./vartmp/MISTE_timeseries_miste.pdf"
     if(messVerb) message("MISTElite making miste time series plot '", file, "'")
     mistePDF(file, input); afunc(); dev.off()
   } else {
     afunc()
   }
}


misteResidualPlot <- function(input, pdfoutput=FALSE, ...) {
   LM <- get("LM", envir=MISTE)
   isYlog <- input$cb.log.site2

   afunc <- function() {
     rv <- residuals(LM); fv <- predict(LM) # R lm() dependency
     ylim <- max(abs(range(rv, na.rm=TRUE))); ylim <- c(-ylim, ylim)
     xlab <- "FITTED VALUES OF MISTE MODEL, CFS"
     ylab <- "RESIDUALS OF MISTE MODEL, CFS"
     if(isYlog) {
       xlab <- "FITTED VALUES OF MISTE MODEL, log10(CFS)"
       ylab <- "RESIDUALS OF MISTE MODEL, log(CFS)"
     }
     plot(fv, rv, type="n", ylim=ylim, xlab=xlab, ylab=ylab,
          xaxs="i", yaxs="i", las=1, tcl=.7, mgp=c(2,0.5,0))
     lines(range(fv, na.rm=TRUE), rep(0,2), lty=2, lwd=.8)
     points(fv, rv, cex=1.2, lwd=0.71, col=my.purple)
     lines(lowess(fv,rv, f=1/3), lwd=2, col=6)
     miste.mtext("MISTElite RESIDUALS")
   }

   if(pdfoutput) {
     file <- "./vartmp/MISTE_residuals.pdf"
     if(messVerb) message("MISTElite making regress-residual plot  '", file, "'")
     mistePDF(file, input); afunc(); dev.off()
   } else {
     afunc()
   }
}

mistePredictionPlot <- function(input, pdfoutput=FALSE, ...) {
   isXlog <- input$cb.log.site1
   isYlog <- input$cb.log.site2

   XY <- get("XY", envir=MISTE)
   XY <- XY[complete.cases(XY), ]
   X <- XY$predictor; Y <- XY$response
   swathX   <- get("swathX",   envir=MISTE)
   swathY   <- get("swathY",   envir=MISTE)
   predX    <- get("predX",    envir=MISTE)
   predY    <- get("predY",    envir=MISTE)
   predYlwr <- get("predYlwr", envir=MISTE)
   predYupr <- get("predYupr", envir=MISTE)
   xlab <- "SITE NO. 1 (PREDICTOR) STREAMFLOW, CFS"
   ylab <- "SITE NO. 2 (RESPONSE) STREAMFLOW, CFS"
   if(isXlog) {
      X <- log10(X); swathX <- log10(swathX)
      predX <- log10(predX)
      xlab <- "SITE NO. 1 (PREDICTOR) STREAMFLOW, log10(CFS)"
   }
   if(isYlog) {
      Y <- log10(Y); swathY <- log10(swathY)
      predY <- log10(predY)
      predYlwr <- log10(predYlwr); predYupr <- log10(predYupr)
      ylab <- "SITE NO. 2 (RESPONSE) STREAMFLOW, log10(CFS)"
   }
   Xtmp <- c(X, predX);  Ytmp <- c(Y, predY)
   xlim <- range(Xtmp[is.finite(Xtmp)])
   ylim <- range(Ytmp[is.finite(Ytmp)])
   afunc <- function() {
      plot(swathX, swathY, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
           las=1, tcl=.7, mgp=c(2,0.5,0))
      points(X, Y, cex=1.2, lwd=0.71, col=my.purple)
      rug(jitter(X), side=1, ticksize=0.02, col=my.purple)
      rug(jitter(Y), side=2, ticksize=0.02, col=my.purple)
      lines(swathX, swathY, lwd=2.3, col="#48ce53")
      for(i in 1:length(predY)) {
         lines(rep(predX[i], 2), c(predYlwr[i], predYupr[i]),
               col=rgb(0,0.435, 0.254), lwd=0.75)
      }
      points(predX, predY,
             col=rgb(0,0.435, 0.254), pch=21, bg=grey(0.90), cex=1.3, lwd=0.8)
      miste.mtext("MISTElite PREDICTIONS")
   }
   if(pdfoutput) {
     file <- "./vartmp/MISTE_predictions.pdf"
     if(messVerb) message("MISTElite making prediction plot        '", file, "'")
     mistePDF(file, input); afunc(); dev.off()
   } else {
     afunc()
   }
}


misteLagCorrPlot <- function(input, pdfoutput=FALSE, ...) {
   isXlog <- input$cb.log.site1
   isYlog <- input$cb.log.site2

   XY <- get("true.time.matched.XY", envir=MISTE)
   XY <- XY[complete.cases(XY), ]
   X <- XY$predictor; Y <- XY$response
   if(isXlog) {
      X[X == 0] <- zeroValue
      X <- log10(X)
   }
   if(isYlog) {
      Y[Y == 0] <- zeroValue
      Y <- log10(Y)
   }
   Z <- XY$uvPOSIXepoch
   delta.time <- get("GlobalDeltaTime", envir=MISTE) # should be always? in seconds
   count.cor.offset <- input$count.cor.offset * delta.time/(60*60)
   dels <- seq(-count.cor.offset, count.cor.offset, by=1)
   dels <- dels*delta.time
   cors <- vector(mode="numeric", length=length(dels))
   cors <- sapply(1:length(dels), function(i) {
      del <- dels[i]
      A <- data.frame(link=Z+del, X); B <- data.frame(link=Z, Y)
      AB <- merge(A, B, by="link", all=TRUE)
      AB <- AB[complete.cases(AB), ]
      cor <- NULL
      if(length(AB$X) < 3) return(0)
      if(sd(AB$X) == 0 | sd(AB$Y) == 0) return(0)
      try(cor <- cor(AB$X, AB$Y, method="pearson"), silent=TRUE)
      if(is.null(cor)) return(0)
      return(cor)
   })
   max.del <- dels[cors == max(cors, na.rm=TRUE)]
   if(length(max.del) != 1) max.del <- max.del[1]
   seconds.to.hours <- 60*60
   corlag <- seconds.to.hours*input$correlation.lag
   afunc <- function() {
      plot(dels, cors, type="l", xaxs="i", yaxs="i", lwd=1.1, ylim=c(-1,1),
           xlab="LAG IN DELTA SECONDS", ylab="KENDALL CORRELATION",
           las=1, tcl=.7, mgp=c(2,0.5,0))
      lines(rep(0,2), c(-1,1), lty=2)
      lines(par()$usr[1:2], rep(0,2), lty=2)
      lines(rep(max.del,2), c(-1,1), lwd=3,   col=6)
      points(max.del, max(cors),     pch=16,  col=6, cex=1.2)
      lines(rep(corlag, 2), c(-1,1), lwd=1.1, col=rgb(0, 0.435, 0.254))
      mtext(paste0("Lag Maximizes at ", max.del/seconds.to.hours,
                   " hours or ", max.del, " seconds"))
   }
   if(pdfoutput) {
     file <- "./vartmp/MISTE_lagcorr.pdf"
     if(messVerb) message("MISTElite making lag autocorr plot      '", file, "'")
     mistePDF(file, input); afunc(); dev.off()
   } else {
     afunc()
   }
}



misteBoxPlot <- function(input, pdfoutput=FALSE, ...) {
   isXlog <- input$cb.log.site1
   isYlog <- input$cb.log.site2

   XY <- get("XY", envir=MISTE)
   XY <- XY[complete.cases(XY), ]

   Xpred <- get("predX", envir=MISTE)
   Xpred <- Xpred[! is.na(Xpred)] # surely not needed at this point?
   Ypred <- get("predY", envir=MISTE)
   Ypred <- Ypred[! is.na(Ypred)] # surely not needed at this point?

   if(isXlog) XY$predictor <- log10(XY$predictor)
   if(isYlog) XY$response  <- log10(XY$response)
   if(isXlog) Xpred        <- log10(Xpred)
   if(isYlog) Ypred        <- log10(Ypred)

   Xpredictor <- XY$predictor[is.finite(XY$predictor)]
   Yresponse  <- XY$response[is.finite(XY$response)]
   Xpred      <- Xpred[is.finite(Xpred)]
   Ypred      <- Ypred[is.finite(Ypred)]

   Data <- c(Xpredictor, Yresponse, Xpred, Ypred)
   Grp  <- c(rep("Input Predictor",  length(Xpredictor)),
             rep("Input Response",   length(Yresponse)),
             rep("Output Predictor", length(Xpred)),
             rep("Output Response",  length(Ypred)) )
   xlab <- "PREDICTOR SITE NONEXCEEDANCE PROBABILITY"
   ylab <- "RESPONSE SITE NONEXCEEDANCE PROBABILITY"
   afunc <- function() {
      ylab <- ifelse(isYlog, "STREAMFLOW, log10(CFS)", "STREAMFLOW, IN CFS")
      boxplot(Data~as.factor(Grp), ylab=ylab, col=grey(0.85),
                   border=c(2,4,my.red2,rgb(0,0.435, 0.254)),
                   las=1, tcl=.7, mgp=c(2,0.5,0))
      miste.mtext("MISTElite BOXPLOT", cex=0.8)
   }
   if(pdfoutput) {
     file <- "./vartmp/MISTE_boxplot.pdf"
     if(messVerb) message("MISTElite making distribution boxplots  '", file, "'")
     mistePDF(file, input); afunc(); dev.off()
   } else {
     afunc()
   }
}
