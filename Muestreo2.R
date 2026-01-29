library(shiny)
library(bslib)
library(ggplot2)
library(ggforce)
library(DT)
make_population <-function(n){

  res <-data.frame(
    x=runif(n,0,100),
    y=runif(n,0,100),
    diam = round(runif(n,10,50),1)
    )
  
  res$radio_fijo <- 15
  res$area_fijo <-  pi*(res$radio_fijo^2)/10000
  res$fac_exp_fijo<- 1/res$area_fijo
  
  res$radio_variable <- ifelse(res$diam<15,10,20)
  res$area_variable <-  pi*(res$radio_variable^2)/10000
  res$fac_exp_variable <- 1/res$area_variable
  
  res$radio_relascopio <- res$diam/2
  res$area_relascopio <-  pi*(res$radio_relascopio^2)/10000
  res$fac_exp_relascopio <- 1/res$area_relascopio
  return(res)
}

sampling_points <-function(n){data.frame(x=runif(n,0,100),y=runif(n,0,100))}

forest_data <- make_population(20)
sample <- sampling_points(10)

space <- br()
plot_type <- radioButtons("tipo", "Tipo de parcela:",
                          c("Radio fijo 15 m" = "fijo",
                            "Radios anidados d<15 10m, d>=15 10m " = "variable",
                            "Relascopio BAF=1" = "relascopio"),)
pop_size <- sliderInput("N",
                        "Número de árboles",
                        value = 20,
                        min = 1,
                        max = 500)
samp_size <- sliderInput("n",
                         "Número de parcelas",
                         value = 10,
                         min = 1,
                         max = 50)
reps <-sliderInput("r",
                   "Repeticiones",
                   value = 100,
                   min = 1,
                   max = 200)
reset <- actionButton("reset", "Regenera poblacion")
controls <- list(pop_size,space,
                 samp_size,space,reps,space,reset)


# Define UI for random distribution app ----
ui <- page_navbar(
  
  # App title ----
  title = "Opciones población y muestra",
  
  # Sidebar layout with input and output definitions ----
  sidebar=sidebar(controls,open="always"),
  tabsetPanel(type = "tabs",
           tabPanel("Población",
                    plotOutput("plot_poblacion",width=500,height=500),
                    br(),
                    tableOutput('poblacion')),
           tabPanel("Muestras",
                    fluidRow(splitLayout(
                      style = "border: 1px solid silver:", cellWidths = c(500,500,500),
                      plotOutput("plot_fijo",width=500,height=500),
                      plotOutput("plot_variable",width=500,height=500),
                      plotOutput("plot_relascopio",width=500,height=500)
                    ))),
           tabPanel("Estimaciones", 
                    plot_type,
                    br(),
                    DT::dataTableOutput('muestras'),
                    br(),
                    plotOutput("plot_selected1",width=500,height=500),
                    br(),
                    plotOutput("plot_selected2",width=500,height=500)),
           tabPanel("Distribución muestral",
                    plotOutput("dist_parcela",width=500,height=500),
                    br(),
                    plotOutput("dist_final",width=500,height=500))
  )
                 

)
# Define server logic for random distribution app ----
server <- function(input, output) {
  
  
  gg_plot <- reactive({
    input$reset
    input$N
    reset()
    rect <- data.frame(x=c(0,0,100,100,0),y=c(0,100,100,0,0))
    ggplot(forest_data) +
      geom_polygon(data=rect,aes(x=x,y=y),col="red",fill="darkgreen",alpha=0.1)+
      geom_point(aes(x=x,y=y),col="black",size=2)+
      coord_fixed(ratio = 1) +
      theme_bw(base_size = 16) +
      theme(axis.title = element_blank())
  })

  reset <- reactive({
    input$reset
    input$N
    forest_data <<- make_population(input$N)
  })
  
  
    # Generate an HTML table view of the data ----
  output$poblacion <- renderTable({
    reset()
    forest_data
  })
  
  output$plot_poblacion<-renderPlot({
    gg_plot()+
      geom_label(aes(x=x,y=y-5,label=diam))+
      xlim(c(-20,120)) + ylim(c(-20,120))+
      ggtitle("Radio fijo 15 m")
  })
  
  output$plot_fijo <- renderPlot({
    gg_plot()+
      geom_circle(aes(x0=x,y0=y,r=radio_fijo,fill=radio_fijo),alpha=0.2)+
      xlim(c(-20,120)) + ylim(c(-20,120))+
      guides(fill = guide_colourbar("radio",position = "bottom"))+
      ggtitle("Radio fijo 15m")
  })
  
  output$plot_variable <- renderPlot({
    gg_plot()+
      geom_circle(aes(x0=x,y0=y,r=radio_variable,fill=radio_variable),alpha=0.2)+
      xlim(c(-20,120)) + ylim(c(-20,120))+
      guides(fill = guide_colourbar("radio",position = "bottom"))+
      ggtitle("Radios anidados d<15 cm 10, d>=15 cm 20m")
    
  })
  
  output$plot_relascopio <- renderPlot({
    gg_plot()+
      geom_circle(aes(x0=x,y0=y,r=radio_relascopio,fill=radio_relascopio),alpha=0.2)+
      xlim(c(-20,120)) + ylim(c(-20,120))+
      guides(fill = guide_colourbar("radio",position = "bottom"))+
      ggtitle("Relascopio BAF=1")
  })

  output$plot_selected1 <- renderPlot({
    reset()
    type <- input$tipo
    N <- input$N
    n <- input$n
    plot(forest_data$x,forest_data$y,
         main = paste("Parcelas ", type, sep = ""),
         pch=20)
  })
  
  output$plotaverage<- renderPlot({
    reset()
    type <- input$tipo
    N <- input$N
    n <- input$n
    
    plot(forest_data$x,forest_data$y,
         main = paste("Parcelas ", type, sep = ""),
         pch=20)
  })
  
  # Generate a summary of the data ----
  output$distmuest <- renderPrint({
    summary(forest())
  })
  

  
}

# Create Shiny app ----
shinyApp(ui, server)