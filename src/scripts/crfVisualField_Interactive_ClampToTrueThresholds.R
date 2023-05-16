crfVisualField_Interactive_ClampToTrueThresholds <- function() {  
  
  rm(list=ls())
  
  # load package
  require(visualFields)
  require(shiny)
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "crfVisualFieldFun.R",sep = ""))
  
  # DEFS FROM OUTSIDE
  normThresholdSD <- 30
  blindSpotThreshold <- -2
  experimentNo <- 1 # only for the age
  
  # load normative data
  data(nvsapdefault)
  normValueParameters <- nvsapdefault$p24d2_sitas$agelm
  patientAge <- vf91016right$sage[[experimentNo]]
  normThresholds <- normValueParameters[,1] + patientAge*normValueParameters[,2]
  normThresholds[which(normThresholds %in% NA)] <- blindSpotThreshold  
  
  # create shiny page
  ui <- fluidPage(titlePanel("Visual Field CRF: Modify the parameters!"),
                  fluidRow(column(2,wellPanel(
                    numericInput(inputId = "L1", "Node potential means: L1", value=round(normThresholds[[1]],0), -5, 40, 1),
                    numericInput(inputId = "L2", "L2", value=round(normThresholds[[2]],0), -5, 40, 1), 
                    numericInput(inputId = "L3", "L3", value=round(normThresholds[[3]],0), -5, 40, 1), 
                    numericInput(inputId = "L4", "L4", value=round(normThresholds[[4]],0), -5, 40, 1), 
                    numericInput(inputId = "L5", "L5", value=round(normThresholds[[5]],0), -5, 40, 1), 
                    numericInput(inputId = "L6", "L6", value=round(normThresholds[[6]],0), -5, 40, 1), 
                    numericInput(inputId = "L7", "L7", value=round(normThresholds[[7]],0), -5, 40, 1), 
                    numericInput(inputId = "L8", "L8", value=round(normThresholds[[8]],0), -5, 40, 1),
                    numericInput(inputId = "L9", "L9", value=round(normThresholds[[9]],0), -5, 40, 1),
                    numericInput(inputId = "L10", "L10", value=round(normThresholds[[10]],0), -5, 40, 1),
                    numericInput(inputId = "L11", "L11", value=round(normThresholds[[11]],0), -5, 40, 1),
                    numericInput(inputId = "L12", "L12", value=round(normThresholds[[12]],0), -5, 40, 1),
                    numericInput(inputId = "L13", "L13", value=round(normThresholds[[13]],0), -5, 40, 1),
                    numericInput(inputId = "L14", "L14", value=round(normThresholds[[14]],0), -5, 40, 1),
                    numericInput(inputId = "L15", "L15", value=round(normThresholds[[15]],0), -5, 40, 1),
                    numericInput(inputId = "L16", "L16", value=round(normThresholds[[16]],0), -5, 40, 1),
                    numericInput(inputId = "L17", "L17", value=round(normThresholds[[17]],0), -5, 40, 1),
                    numericInput(inputId = "L18", "L18", value=round(normThresholds[[18]],0), -5, 40, 1),
                    numericInput(inputId = "L19", "L19", value=round(normThresholds[[19]],0), -5, 40, 1),
                    numericInput(inputId = "L20", "L20", value=round(normThresholds[[20]],0), -5, 40, 1),
                    numericInput(inputId = "L21", "L21", value=round(normThresholds[[21]],0), -5, 40, 1),
                    numericInput(inputId = "L22", "L22", value=round(normThresholds[[22]],0), -5, 40, 1),
                    numericInput(inputId = "L23", "L23", value=round(normThresholds[[23]],0), -5, 40, 1),
                    numericInput(inputId = "L24", "L24", value=round(normThresholds[[24]],0), -5, 40, 1),
                    numericInput(inputId = "L25", "L25", value=round(normThresholds[[25]],0), -5, 40, 1),
                    numericInput(inputId = "L26", "L26", value=round(normThresholds[[26]],0), -5, 40, 1),
                    numericInput(inputId = "L27", "L27", value=round(normThresholds[[27]],0), -5, 40, 1),
                    numericInput(inputId = "L28", "L28", value=round(normThresholds[[28]],0), -5, 40, 1),
                    numericInput(inputId = "L29", "L29", value=round(normThresholds[[29]],0), -5, 40, 1),
                    numericInput(inputId = "L30", "L30", value=round(normThresholds[[30]],0), -5, 40, 1),
                    numericInput(inputId = "L31", "L31", value=round(normThresholds[[31]],0), -5, 40, 1),
                    numericInput(inputId = "L32", "L32", value=round(normThresholds[[32]],0), -5, 40, 1),
                    numericInput(inputId = "L33", "L33", value=round(normThresholds[[33]],0), -5, 40, 1),
                    numericInput(inputId = "L34", "L34", value=round(normThresholds[[34]],0), -5, 40, 1),
                    numericInput(inputId = "L35", "L35", value=round(normThresholds[[35]],0), -5, 40, 1),
                    numericInput(inputId = "L36", "L36", value=round(normThresholds[[36]],0), -5, 40, 1),
                    numericInput(inputId = "L37", "L37", value=round(normThresholds[[37]],0), -5, 40, 1),
                    numericInput(inputId = "L38", "L38", value=round(normThresholds[[38]],0), -5, 40, 1),
                    numericInput(inputId = "L39", "L39", value=round(normThresholds[[39]],0), -5, 40, 1),
                    numericInput(inputId = "L40", "L40", value=round(normThresholds[[40]],0), -5, 40, 1),
                    numericInput(inputId = "L41", "L41", value=round(normThresholds[[41]],0), -5, 40, 1),
                    numericInput(inputId = "L42", "L42", value=round(normThresholds[[42]],0), -5, 40, 1),
                    numericInput(inputId = "L43", "L43", value=round(normThresholds[[43]],0), -5, 40, 1),
                    numericInput(inputId = "L44", "L44", value=round(normThresholds[[44]],0), -5, 40, 1),
                    numericInput(inputId = "L45", "L45", value=round(normThresholds[[45]],0), -5, 40, 1),
                    numericInput(inputId = "L46", "L46", value=round(normThresholds[[46]],0), -5, 40, 1),
                    numericInput(inputId = "L47", "L47", value=round(normThresholds[[47]],0), -5, 40, 1),
                    numericInput(inputId = "L48", "L48", value=round(normThresholds[[48]],0), -5, 40, 1),
                    numericInput(inputId = "L49", "L49", value=round(normThresholds[[49]],0), -5, 40, 1),
                    numericInput(inputId = "L50", "L50", value=round(normThresholds[[50]],0), -5, 40, 1),
                    numericInput(inputId = "L51", "L51", value=round(normThresholds[[51]],0), -5, 40, 1),
                    numericInput(inputId = "L52", "L52", value=round(normThresholds[[52]],0), -5, 40, 1),
                    numericInput(inputId = "L53", "L53", value=round(normThresholds[[53]],0), -5, 40, 1),
                    numericInput(inputId = "L54", "L54", value=round(normThresholds[[54]],0), -5, 40, 1)
                  )),
                  
                  column(2,wellPanel(
                    numericInput(inputId = "sdL1", "Node potential StDevs: L1", value=normThresholdSD, 0, 50, 1),
                    numericInput(inputId = "sdL2", "L2", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL3", "L3", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL4", "L4", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL5", "L5", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL6", "L6", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL7", "L7", value=normThresholdSD, 0, 50, 1),
                    numericInput(inputId = "sdL8", "L8", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL9", "L9", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL10", "L10", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL11", "L11", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL12", "L12", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL13", "L13", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL14", "L14", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL15", "L15", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL16", "L16", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL17", "L17", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL18", "L18", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL19", "L19", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL20", "L20", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL21", "L21", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL22", "L22", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL23", "L23", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL24", "L24", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL25", "L25", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL26", "L26", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL27", "L27", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL28", "L28", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL29", "L29", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL30", "L30", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL31", "L31", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL32", "L32", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL33", "L33", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL34", "L34", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL35", "L35", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL36", "L36", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL37", "L37", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL38", "L38", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL39", "L39", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL40", "L40", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL41", "L41", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL42", "L42", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL43", "L43", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL44", "L44", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL45", "L45", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL46", "L46", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL47", "L47", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL48", "L48", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL49", "L49", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL50", "L50", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL51", "L51", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL52", "L52", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL53", "L53", value=normThresholdSD, 0, 50, 1), 
                    numericInput(inputId = "sdL54", "L54", value=normThresholdSD, 0, 50, 1)
                  )),
                  
                  column(2,wellPanel(
                    numericInput(inputId = "clL1", "Clamp nodes (0=unclamped): L1", value=NA, -5, 40, 1),
                    numericInput(inputId = "clL2", "L2", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL3", "L3", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL4", "L4", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL5", "L5", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL6", "L6", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL7", "L7", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL8", "L8", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL9", "L9", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL10", "L10", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL11", "L11", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL12", "L12", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL13", "L13", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL14", "L14", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL15", "L15", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL16", "L16", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL17", "L17", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL18", "L18", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL19", "L19", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL20", "L20", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL21", "L21", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL22", "L22", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL23", "L23", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL24", "L24", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL25", "L25", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL26", "L26", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL27", "L27", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL28", "L28", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL29", "L29", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL30", "L30", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL31", "L31", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL32", "L32", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL33", "L33", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL34", "L34", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL35", "L35", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL36", "L36", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL37", "L37", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL38", "L38", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL39", "L39", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL40", "L40", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL41", "L41", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL42", "L42", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL43", "L43", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL44", "L44", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL45", "L45", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL46", "L46", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL47", "L47", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL48", "L48", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL49", "L49", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL50", "L50", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL51", "L51", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL52", "L52", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL53", "L53", value=NA, -5, 40, 1), 
                    numericInput(inputId = "clL54", "L54", value=NA, -5, 40, 1)
                  )),
                  
                  column(6,
                         numericInput(inputId = "edgeSD", "Edge StDev", value=15, 0, 50, 1),
                         numericInput(inputId = "maxIter", "Max iterations", value=3, 0, 10000, 10),
                         plotOutput("Vf"))
                  ))
  
  server <- function(input,output) {
    output$Vf <- renderPlot({
      crfVisualField(c(input$L1,input$L2,input$L3,input$L4,input$L5,input$L6,input$L7,input$L8,input$L9,input$L10,
                       input$L11,input$L12,input$L13,input$L14,input$L15,input$L16,input$L17,input$L18,input$L19,
                       input$L20,input$L21,input$L22,input$L23,input$L24,input$L25,input$L26,input$L27,input$L28,
                       input$L29,input$L30,input$L31,input$L32,input$L33,input$L34,input$L35,input$L36,input$L37,
                       input$L38,input$L39,input$L40,input$L41,input$L42,input$L43,input$L44,input$L45,input$L46,
                       input$L47,input$L48,input$L49,input$L50,input$L51,input$L52,input$L53,input$L54), 
                     c(input$sdL1,input$sdL2,input$sdL3,input$sdL4,input$sdL5,input$sdL6,input$sdL7,input$sdL8,input$sdL9,input$sdL10,
                       input$sdL11,input$sdL12,input$sdL13,input$sdL14,input$sdL15,input$sdL16,input$sdL17,input$sdL18,input$sdL19,input$sdL20,
                       input$sdL21,input$sdL22,input$sdL23,input$sdL24,input$sdL25,input$sdL26,input$sdL27,input$sdL28,input$sdL29,input$sdL30,
                       input$sdL31,input$sdL32,input$sdL33,input$sdL34,input$sdL35,input$sdL36,input$sdL37,input$sdL38,input$sdL39,input$sdL40,
                       input$sdL41,input$sdL42,input$sdL43,input$sdL44,input$sdL45,input$sdL46,input$sdL47,input$sdL48,input$sdL49,input$sdL50,
                       input$sdL51,input$sdL52,input$sdL53,input$sdL54),
                     input$edgeSD,
                     c(input$clL1,input$clL2,input$clL3,input$clL4,input$clL5,input$clL6,input$clL7,input$clL8,input$clL9,input$clL10,
                       input$clL11,input$clL12,input$clL13,input$clL14,input$clL15,input$clL16,input$clL17,input$clL18,input$clL19,input$clL20,
                       input$clL21,input$clL22,input$clL23,input$clL24,input$clL25,input$clL26,input$clL27,input$clL28,input$clL29,input$clL30,
                       input$clL31,input$clL32,input$clL33,input$clL34,input$clL35,input$clL36,input$clL37,input$clL38,input$clL39,input$clL40,
                       input$clL41,input$clL42,input$clL43,input$clL44,input$clL45,input$clL46,input$clL47,input$clL48,input$clL49,input$clL50,
                       input$clL51,input$clL52,input$clL53,input$clL54),
                     input$maxIter
                     )
      
    },height=650)
    
  }
  
  shinyApp(ui = ui,server = server)
} 
