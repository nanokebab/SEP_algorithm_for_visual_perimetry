crfVisualField_Interactive_Modify <- function() {  
  
  rm(list=ls())
  
  # load package
  # require(visualFields)
  require(shiny)
  
  #load functions
  funPath <- "src/utils/"
  source(paste(funPath, "crfVisualFieldFunWrapper.R",sep = ""))
  
  # create shiny page
  ui <- fluidPage(titlePanel("Clamp arbitrary locations to true values"),
                  fluidRow(column(12,
                    wellPanel(
                        numericInput(inputId = "experimentNo", label="Visual field no. [1-4863]", 1260, 1, 4863, 1),
                        checkboxGroupInput(inputId = "clamped", label="Clamped nodes (e.g. 3, 6, 10, 11, 15, 22, 26, 29, 33, 39, 43, 48 and 51)",
                                       choices = as.character(seq(1,54,1)),inline=TRUE ),
                        numericInput(inputId = "edgeSD", label="Edge Std:", 5, 1, 50, 1)
                      ))),
                    fluidRow(column(12,plotOutput("Vf")))
                    )
  
  server <- function(input,output) {
    output$Vf <- renderPlot({
      # crfVisualFieldWrapper(input$experimentNo,input$edgeSD,input$clamped)
      crfVisualFieldWrapper(input$experimentNo,input$edgeSD,as.numeric(input$clamped))
    })
    
  }
  
  shinyApp(ui = ui,server = server)
}
