library(shiny)
plotrendW <- 800; plotrendH <- 800

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  dataExtract_site1 <- reactive({
     #message("Site1 '", input$site1, "'")
     dataExtract_inputdates <- reactive({input$inputdates})
       misteinputdates <- dataExtract_inputdates()
       misteinputdates <- misteDateRange2Dates(misteinputdates)
     dataExtract_inputtimes <- reactive({input$inputtimes})
       misteinputtimes <- dataExtract_inputtimes()
       misteinputtimes <- misteHHMM2time(misteinputtimes)
     UVs <- misteUVget(input$site1, sdate=misteinputdates$sdate$sdate,
                                    edate=misteinputdates$edate$edate,
                                    stime=misteinputtimes$stime,
                                    etime=misteinputtimes$etime,
                       storestring="inputSite1", message="data for regression")
     tmp <- UVs
     tmp$dateTime <- as.character(format(tmp$dateTime, usetz=TRUE))
     file <- "./vartmp/MISTE_site1data(predictor).txt"
     write.table(tmp, file=file, row.names=FALSE, quote=FALSE, sep="\t")
     if(messVerb) message("MISTElite archiving data in file '", file, "'")
     return(UVs)
  })

  dataExtract_site2 <- reactive({
     #message("Site2 '", input$site2, "'")
     dataExtract_inputdates <- reactive({input$inputdates})
       misteinputdates <- dataExtract_inputdates()
       misteinputdates <- misteDateRange2Dates(misteinputdates)
     dataExtract_inputtimes <- reactive({input$inputtimes})
       misteinputtimes <- dataExtract_inputtimes()
       misteinputtimes <- misteHHMM2time(misteinputtimes)
     UVs <- misteUVget(input$site2, sdate=misteinputdates$sdate$sdate,
                                    edate=misteinputdates$edate$edate,
                                    stime=misteinputtimes$stime,
                                    etime=misteinputtimes$etime,
                       storestring="inputSite2", message="data for regression")
     tmp <- UVs
     tmp$dateTime <- as.character(format(tmp$dateTime, usetz=TRUE))
     file <- "./vartmp/MISTE_site2data(response).txt"
     write.table(tmp, file=file, row.names=FALSE, quote=FALSE, sep="\t")
     if(messVerb) message("MISTElite archiving data in file '", file, "'")
     return(UVs)
  })

  dataExtract_site1_output <- reactive({
     corlag <- input$correlation.lag
     dataExtract_outputdates <- reactive({input$outputdates})
       misteoutputdates <- dataExtract_outputdates()
       misteoutputdates <- misteDateRange2Dates(misteoutputdates)
     dataExtract_outputtimes <- reactive({input$outputtimes})
       misteoutputtimes <- dataExtract_outputtimes()
       misteoutputtimes <- misteHHMM2time(misteoutputtimes)
     if(devVerb) message(" *** Now need to silently +/- lag somehow on site 1 ouput ",
                         "retrieval sequence")
     "dtmassage" <- function() { # Logic sucks but does seem to work, place for code cleanup
        if(devVerb) message(" *** dtmassage() manip dates by corlag=", corlag,
                            " hrs ahead of misteUVget() call ***")
        sdt <- paste0(misteoutputdates$sdate$sdate, " ", misteoutputtimes$shhmmss)
        edt <- paste0(misteoutputdates$edate$edate, " ", misteoutputtimes$ehhmmss)
        sdt <- strptime(sdt, format="%Y-%m-%d %H:%M:%S", tz=misteoutputtimes$tz)
        edt <- strptime(edt, format="%Y-%m-%d %H:%M:%S", tz=misteoutputtimes$tz)
        if(devVerb) message(" *** Preds. wanted at site 2 in real times start: ",
                                          as.character(format(sdt, usetz=TRUE)))
        if(devVerb) message(" ***                                         end: ",
                                          as.character(format(edt, usetz=TRUE)))
        real_sdt <- unlist(strsplit(as.character(sdt), " ", perl=TRUE))
        real_edt <- unlist(strsplit(as.character(edt), " ", perl=TRUE))
        real_sdt[2] <- paste0(real_sdt[2],misteoutputtimes$tz)
        real_edt[2] <- paste0(real_edt[2],misteoutputtimes$tz)
        if(real_sdt[2] == "NA") real_sdt[2] <- paste0("00:00:00",misteoutputtimes$tz)
        if(real_edt[2] == "NA") real_edt[2] <- paste0("00:00:00",misteoutputtimes$tz)

        seconds.to.hours <- 60*60
        proposed_sdt <- sdt - corlag*seconds.to.hours # sign convention is to subtract
        proposed_edt <- edt - corlag*seconds.to.hours # sign convention is to subtract
        if(devVerb) message(" *** But need with corlag offset at site 1 start: ",
                                 as.character(format(proposed_sdt, usetz=TRUE)))
        if(devVerb) message(" ***                                         end: ",
                                 as.character(format(proposed_edt, usetz=TRUE)))
        proposed_sdt <- unlist(strsplit(as.character(proposed_sdt), " ", perl=TRUE))
        proposed_edt <- unlist(strsplit(as.character(proposed_edt), " ", perl=TRUE))
        proposed_sdt[2] <- paste0(proposed_sdt[2],misteoutputtimes$tz)
        proposed_edt[2] <- paste0(proposed_edt[2],misteoutputtimes$tz)

        beforsdt <- sdt - 24*seconds.to.hours; #print(beforsdt)
        afteredt <- edt + 24*seconds.to.hours; #print(afteredt)
        if(devVerb) message(" ***         Before time augmentation start date: ",
                                 as.character(format(beforsdt, usetz=TRUE)))
        if(devVerb) message(" ***            After time augmentation end date: ",
                                 as.character(format(afteredt, usetz=TRUE)))

        beforsdt <- unlist(strsplit(as.character(beforsdt), " ", perl=TRUE))
        afteredt <- unlist(strsplit(as.character(afteredt), " ", perl=TRUE))
        beforsdt[2] <- paste0(beforsdt[2],misteoutputtimes$tz)
        afteredt[2] <- paste0(afteredt[2],misteoutputtimes$tz)
        if(beforsdt[2] == "NA") beforsdt[2] <- paste0("00:00:00",misteoutputtimes$tz)
        if(afteredt[2] == "NA") afteredt[2] <- paste0("00:00:00",misteoutputtimes$tz)

        zz <- list(sdate=   proposed_sdt[1], stime=   proposed_sdt[2],
                   edate=   proposed_edt[1], etime=   proposed_edt[2],
                   real_sdate=  real_sdt[1], real_stime=  real_sdt[2],
                   real_edate=  real_edt[1], real_etime=  real_edt[2],
                   before_sdate=beforsdt[1], before_stime=beforsdt[2],
                   after_edate= afteredt[1], after_etime= afteredt[2])
        return(zz)
     }
     ZZ <- dtmassage(); #print(ZZ)
     UVs <- misteUVget(input$site1, sdate=ZZ$sdate, edate=ZZ$edate,
                                    stime=ZZ$stime, etime=ZZ$etime,
                       storestring="outputSite1", message="data for prediction")
     BEFORE_UVs <- misteUVget(input$site2,
                              sdate=ZZ$before_sdate, edate=ZZ$real_sdate,
                              stime=ZZ$before_stime, etime=ZZ$real_stime,
                              storestring="augmentbeforeSite2",
                              message="before data for graphic augmentation")
     AFTER_UVs  <- misteUVget(input$site2,
                              sdate=ZZ$real_edate, edate=ZZ$after_edate,
                              stime=ZZ$real_etime, etime=ZZ$after_etime,
                              storestring="augmentafterSite2",
                              message="after data for graphic augmentation")
     assign("augmentbeforeSite2", BEFORE_UVs, envir=MISTE)
     assign("augmentafterSite2",  AFTER_UVs,  envir=MISTE)

     tmp <- UVs
     tmp$dateTime <- as.character(format(tmp$dateTime, usetz=TRUE))
     file <- "./vartmp/MISTE_site1data(forprediction).txt"
     write.table(tmp, file=file, row.names=FALSE, quote=FALSE, sep="\t")
     if(messVerb) message("MISTElite archiving data in file        '", file, "'")
     Do_misteModel_and_PredictNeed_overlap()

     tmp <- BEFORE_UVs
     tmp$dateTime <- as.character(format(tmp$dateTime, usetz=TRUE))
     file <- "./vartmp/MISTE_site2data(beforeaugment).txt"
     write.table(tmp, file=file, row.names=FALSE, quote=FALSE, sep="\t")
     if(messVerb) message("MISTElite archiving data in file        '", file, "'")

     tmp <- AFTER_UVs
     tmp$dateTime <- as.character(format(tmp$dateTime, usetz=TRUE))
     file <- "./vartmp/MISTE_site2data(afteraugment).txt"
     write.table(tmp, file=file, row.names=FALSE, quote=FALSE, sep="\t")
     if(messVerb) message("MISTElite archiving data in file        '", file, "'")

     return(UVs)
  })

  output$Plots <- renderPlot({
    message("-------- MISTE beginning")
    site1UV <- dataExtract_site1()
    site2UV <- dataExtract_site2()
    XY <- misteUVmerge(site1UV, site2UV)
    if(length(XY$predictor) <= 2) {
       txt <- "MISTE ERROR: Insufficient Data Provided"
       output$lmdia.text <- renderText(txt)
    } else {
       misteMODEL(XY, input);#print("misteMODEL")
       diagnosticText    <- misteMODELdiagnosticText(input);#print("misteMODELdiagnosticText")
       output$lmdia.text <- renderText(diagnosticText[1])
       output$lmres.text <- renderText(diagnosticText[2])
       site1UV4prediction <- dataExtract_site1_output();#print("dataExtract_site1_output")
       misteMODELpredict(site1UV4prediction, input);#print("misteMODELpredict")

       def.par <- par(no.readonly = TRUE) # save default, for resetting...
       par(mai=c(0.75, 0.75, 0.25, 0.25))
       layout(matrix(c(1,2,3,4,5,6), 3,2, byrow=TRUE))
         misteTimeSeriesPlot_Inputs( input)
         misteTimeSeriesPlot_Product(input)
         misteResidualPlot(          input)
         mistePredictionPlot(        input)
         misteLagCorrPlot(           input)
         misteBoxPlot(               input)
       layout(1) # rarely needed but seems to help if intermediate failure
       par(def.par) # resetting plotting parameters
         misteTimeSeriesPlot_Inputs( input, pdfoutput=TRUE)
         misteTimeSeriesPlot_Product(input, pdfoutput=TRUE)
         mistePredictionPlot(        input, pdfoutput=TRUE)
         misteResidualPlot(          input, pdfoutput=TRUE)
         misteLagCorrPlot(           input, pdfoutput=TRUE)
         misteBoxPlot(               input, pdfoutput=TRUE)
    }
    file <- "./vartmp/MISTE.RData"
    if(messVerb) message("MISTElite saving the MISTE R-envir in   '", file, "'")
    save(MISTE, file=file)
    MISTE <- new.env() # WIPE MISTE
    #message("MISTE wiping the MISTE R-envir in runtime session")
    message("  WA-TODO: archive MISTElite settings and archival report.")
    message("  BB-TODO: GAM prediction limits.")
    message(" WA?-TODO: more introspection of local time zone implicitness.")
    message("WABB-TODO: lag cause need of same data as used to build the model.")
    message("  WA-TODO: code cleaning and use misteLogarithms()")
    message("-------- MISTElite ending"); message("")
  }, width=plotrendW, height=plotrendH)
})

