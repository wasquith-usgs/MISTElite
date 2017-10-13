library(shiny)
source("misteuv.R")

plotW <- "auto"; plotH <- 900
wdW   <-  "8em"; hhdW  <- "18em"; wdID  <- "10em"

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  # Application title
  titlePanel(title=h3("Simple USGS NWIS Missing Record Estimator (MISTElite)"),
             windowTitle="USGS-NWIS MISTElite"),
  navlistPanel("Cntrls", widths = c(1, 10),
  tabPanel("Sites",
  fluidRow(
     column(2,textInput("site1",
                 misteHB("Site 1 (predictor)"),
                         value="08173900", width=wdID)),
     column(2,checkboxInput("cb.log.site1",
                 misteHC("log transform"), value = TRUE)),
     column(2,textInput("site2",
                 misteHB("Site 2 (response)"),
                         value="08175800", width=wdID)),
     column(2,checkboxInput("cb.log.site2",
                 misteHC("log transform"), value = TRUE),
              checkboxInput("cb.shape.transform",
                 misteHC("shape transform to entry/exit targets"), value = FALSE))
  ),
  fluidRow(
     column(3,dateRangeInput("inputdates",
                 misteHC("Dates to build MISTE and optional times"),
                         start="2016-03-06", end="2016-03-11")),
     column(3,dateRangeInput("outputdates",
                 misteHC("Dates to make  MISTE and optional times"),
                         start="2016-03-11", end="2016-03-12")),
     column(3,numericInput("target.entry.discharge",
                 misteHC("Target discharge on epoch entry (cfs)"),
                         value = 5660, step=100, width=wdW))
  ),
  fluidRow(
     column(3, textInput("inputtimes",
               misteHC("HH:MMtoHH:MM [UTC]"),
                       value="12:00to14:00", width=hhdW)),
     column(3, textInput("outputtimes",
               misteHC("HH:MMtoHH:MM [UTC]"),
                       value="02:00to23:00", width=hhdW)),
     column(3,numericInput("target.exit.discharge",
                 misteHC("Target discharge on epoch exit (cfs)"),
                         value = 6760, step=100, width=wdW))
  ) # tabs
  ),
  tabPanel("Cuts",
  fluidRow(
     column(2,numericInput("site1.lower.cutout",
                 misteHC("Site 1 (predictor) min cut (cfs)"),
                         value = 0,   step=100, width=wdW)),
     column(2,numericInput("site1.upper.cutout",
                 misteHC("Site 1 (predictor) max cut (cfs)"),
                         value = 1E6, step=100, width=wdW)),
     column(2,numericInput("site2.lower.cutout",
                 misteHC("Site 2 (response) min cut (cfs)"),
                         value = 0,   step=100, width=wdW)),
     column(2,numericInput("site2.upper.cutout",
                 misteHC("Site 2 (response) max cut (cfs)"),
                         value = 1E6, step=100, width=wdW))
  ) # tabs
  ),
  tabPanel("Model",
  fluidRow(
     column(2, numericInput("poly.order",
                 misteHC("OLS poly order"),
                         value = 1, step=1, min=1, max=10, width=wdW)),
     column(2,checkboxInput("cb.duansmearing",
       misteHC("Apply Duan bias-corr factor, iff log10(y)"),
                         value = TRUE)),
     column(2,numericInput("prediction.interval",
                 misteHC("Predict interval"),
                         value=0.90, min=0, max=1, step=.05, width=wdW))
  ),
  fluidRow(
     column(2,numericInput("correlation.lag",
                 misteHC("Corr lag (hrs)"),
                         value = 19.25,   step=1, min=-21, max=21, width=wdW)),
     column(2,numericInput("count.cor.offset",
                 misteHC("Corr window (+-hrs)"),
                         value=1000, min=0, max=10000, step=1, width=wdW)),
     column(2,checkboxInput("cb.modelgam",
                 misteHC("GAM model"),  value = TRUE),
              numericInput("modelgam.knots",
                 misteHC("Basis dimension"), value=-1, min=-1, max=100, step=1))
  ) # tabs
  ),
  fluidRow(
     column(2, submitButton("MISTE it!")),
     column(4, em("MISTElite 0.001 (wha)")),
     column(2,checkboxInput("cb.makemud", misteHC("MUD output"), value=TRUE))
  ),
  fluidRow(
     column(8, em(paste0("Note: ./vartmp/ has text and pdf files including ",
                         "a .mud file for AQ.")))
  )
  ),
   # Show a plot of the generated distribution
   mainPanel(
      plotOutput("Plots", width=plotW, height=plotH),
      h4("MISTE Diagnostics"),
      h5(textOutput("lmdia.text")),
      h5(textOutput("lmres.text")),
      h5(textOutput("lmepc.text"))
   )
))
