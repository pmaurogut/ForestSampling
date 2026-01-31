library(shiny)
library(bslib)
library(ggplot2)
library(ggforce)
library(dplyr)
library(purrr)
library(DT)

make_population <-function(N,L){

  A<-(L*L)/10000
  Npop<-round(N*A)
  res <-data.frame(
    id = 1:Npop,
    x=runif(Npop,0,L),
    y=runif(Npop,0,L),
    diam = round(runif(N,5,40),1)
    )
  
  res$ht <- round(5+(res$diam-5)*0.5,1)
  
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

sampling_points <-function(n,L){data.frame(Parcela=1:n,x=runif(n,0,L),y=runif(n,0,L))}

get_points <- function(population,point,type){
  
  pick <- sqrt((population$x-point$x)^2+(population$y-point$y)^2)<=population[[type]]
  if(all(!pick)){
    return(NULL)
  }else{
    res<-population[pick, ]
    res$Parcela <- point$Parcela
    res$type <- gsub("radio_","",type)
    res$x0 <- point$x
    res$y0 <- point$y
    res$gi_m2 <- (pi*(1/4)*res$diam^2)/10000
    res<-res[,c("Parcela","type","x0","y0","x","y","diam","ht","gi_m2",type)]
    res$area_parecela_ha<-(pi*res[,type]^2)/10000
    res$EXP_FAC <- 1/res$area_parecela_ha
    return(res)
  }
  
}

get_n_points <- function(population,points,type){
  
  res<-list()
  for(i in 1:dim(points)[1]){
    res[[i]]<-get_point(population,points[i,],type)
  }
  do.call(rbind,res)
}

plot_selection <- function(p,selected,type,tree_center=TRUE,all=FALSE){

  title <- switch(type,
    radio_fijo = "Radio fijo 15 m",
    radio_variable = "R anidados d<15 cm 10, d>=15 cm 20m",
    radio_relascopio = "Relascopio BAF=1"
  )
  if(all){
    p <- p + geom_circle(aes(x0=x,y0=y,r=.data[[type]]), fill="grey50",alpha=0.2)
  }
  
  if(!is.null(selected)){
    if(tree_center){
      p <- p  + geom_circle(data=selected,aes(x0=x,y0=y,r=.data[[type]]),fill="purple",alpha=0.2)
    }else{
      selected2 <- selected |> group_by(!! sym(type)) |> filter(row_number()==1) |> ungroup()
      p <- p  + geom_circle(data=selected2,aes(x0=x0,y0=y0,r=.data[[type]]),fill="purple",alpha=0.2)  
    }
    p <- p +geom_point(data=selected,aes(x=x,y=y),shape=20,col="green")
  }
  p <- p + geom_point(data=sample_data[1,],aes(x=x,y=y),shape=13,col="red",size=4)
  p <- p + guides(fill=FALSE)+ggtitle(title)
  p
}

parametros_interes <- function(poblacion, lado){
  A<-lado*lado/10000
  res<-data.frame(
    Area_ha = A,
    Total_N = length(poblacion[[1]]),
    Total_G = (1/10000)*sum(pi*poblacion$diam^2)/4
  )
  res$N <- res$Total_N/A
  res$G <- res$Total_G/A
  res$h_media <- mean(poblacion$ht)
  res$dg <- sqrt(mean(poblacion$diam^2))
  if(res$N<100){
    res$Ho<-10
  }else{
    res$Ho<-10
  }
  names<- colnames(res)
  res <- as.data.frame(t(res))
  colnames(res)<-"Valor"
  res$parametro<-names
  res
}

forest_data <- make_population(20,100)
sample_data <- sampling_points(10,100)


space <- br()
plot_type <- radioButtons("tipo", "Tipo de parcela:",
                          c("Radio fijo 15 m" = "fijo",
                            "Radios anidados d<15 10m, d>=15 10m " = "variable",
                            "Relascopio BAF=1" = "relascopio"),selected = "fijo")

lado <- sliderInput("lado","Lado (m)",value = 100,min = 100,max = 500,step=50)

pop_size <- sliderInput("N","Número de árboles/ha",value = 20,min = 5,max = 500,step=5)

samp_size <- sliderInput("n","Número de parcelas",value = 20,min = 1,max = 50)

reps <-sliderInput("r","Repeticiones",value = 100,min = 1,max = 200)

reset <- actionButton("reset", "Regenera poblacion")

muestra <- actionButton("muestra", "Toma una muestra",color="darkgreen",alpha=0.4)

n_muestras <- actionButton("n_muestras", "Toma n muestras",color="blue",alpha=0.4)

areas_inclusion <- checkboxInput("all_trees","Todas las areas de inclusion",value = FALSE)

controls <- list(lado,space,pop_size,space,
                 samp_size,space,reps,space,areas_inclusion,space,reset,space,muestra,space,n_muestras)


# Define UI for random distribution app ----
ui <- page_navbar(
  
  # App title ----
  title = "Opciones población y muestra",
  
  # Sidebar layout with input and output definitions ----
  sidebar=sidebar(controls,open="always"),
  tabsetPanel(type = "tabs",
           tabPanel("Población y parámetros de interés",
                    fluidRow(splitLayout(
                      style = "border: 1px solid silver:", cellWidths = c(800,500,500,500),
                      plotOutput("plot_poblacion",width=700,height=700),
                      tableOutput('poblacion'),
                      tableOutput('parametros_interes')
                    ))
                  ),
           tabPanel("Parcelas",
                    fluidRow(
                      splitLayout(
                        style = "border: 1px solid silver:", cellWidths = c(500,500,500,200),
                        plotOutput("plot_fijo",width=500,height=500),
                        plotOutput("plot_variable",width=500,height=500),
                        plotOutput("plot_relascopio",width=500,height=500)
                      )
                    ),
                    fluidRow(
                      splitLayout(
                        style = "border: 1px solid silver:", cellWidths = c(500,500,500,200),
                        plotOutput("plot_fijo2",width=500,height=500),
                        plotOutput("plot_variable2",width=500,height=500),
                        plotOutput("plot_relascopio2",width=500,height=500),
                        tableOutput('pop_stats2')
                      )
                    )
                  ),
            # tabPanel("Estimación una parcela",
            #         column(3,
            #               fluidRow(plot_type),
            #               fluidRow(textOutput("Placeholder1"))
            #               ),
            #         column(6,
            #               fluidRow(
            #                 fluidRow(plotOutput("plot_selected1",width=500,height=500)),
            #                 fluidRow(plotOutput("plot_selected2",width=500,height=500))
            #                 )),
            #         column(3,
            #                textOutput("Placeholder2")
            #                 )
            #   ),    
           
           
          
           tabPanel("Estimación una parcela",
                    fluidRow(splitLayout(
                      style = "border: 1px solid silver:", cellWidths = c(500,500,500,200),
                      plot_type,
                      plotOutput("plot_selected1",width=500,height=500),
                      plotOutput("plot_selected2",width=500,height=500),
                      tableOutput("parametros_interes2")
                    )),
                    fluidRow(splitLayout(
                      style = "border: 1px solid silver:", cellWidths = c(800,350,350),
                        tableOutput('muestra'),
                        textOutput("Placeholder"),
                        textOutput("Placeholder2")
                      ))
                  ),
           tabPanel("Estimación múltiples parcelas",
                    fluidRow(splitLayout(
                      style = "border: 1px solid silver:", cellWidths = c(500,500,500),
                      plot_type,
                      plotOutput("plot_selected3",width=500,height=500),
                      plotOutput("plot_selected4",width=500,height=500)
                    )),
                    fluidRow(DT::dataTableOutput('muestras'))
           ),
           tabPanel("Distribución muestral",
                    plotOutput("dist_parcela",width=500,height=500),
                    br(),
                    plotOutput("dist_final",width=500,height=500)
                    )
  )
                 

)
# Define server logic for random distribution app ----
server <- function(input, output) {
  
  
  gg_plot <- reactive({
    input$reset
    input$lado
    input$N
    lado <- input$lado
    reset()
    rect <- data.frame(x=c(0,0,lado,lado,0),y=c(0,lado,lado,0,0))
    ggplot(forest_data) +
      geom_polygon(data=rect,aes(x=x,y=y),col="red",fill="darkgreen",alpha=0.1)+
      geom_circle(aes(x0=x,y0=y,r=diam/20),col="black")+
      xlim(c(-20,input$lado+20)) + ylim(c(-20,input$lado+20))+
      coord_fixed(ratio = 1) +
      labs(x="x(m)",y="y(m)") +
      theme_bw(base_size = 16) +
      theme(axis.title = element_blank())
  })

  reset <- reactive({
    input$reset
    input$N
    input$lado
    forest_data <<- make_population(input$N,input$lado)
    sample_data <<-sampling_points(input$n,input$lado)
  })
  
  reset_sample<-reactive({
    input$reset
    input$n
    input$lado
    input$n_muestras
    input$muestra
    sample_data <<-sampling_points(input$n,input$lado)
  })
  
    # Generate an HTML table view of the data ----
  output$poblacion <- renderTable({
    reset()
    forest_data[,c(1:5)]
  })
  
  output$plot_poblacion<-renderPlot({
    all <- input$all_trees
    reset_sample()
    gg_plot()+
      geom_label(aes(x=x,y=y-3,label=diam),size=3,fill="darkgreen",alpha=0.3)+
      geom_label(aes(x=x,y=y-8,label=ht),size=3,fill="blue",alpha=0.3)
      # ggtitle("Población")
  })
  
  output$parametros_interes <- renderTable({
    reset()
    parametros_interes(forest_data,input$lado)
    })
  
  output$plot_fijo <- renderPlot({
    all <- input$all_trees
    reset_sample()
    selected <- get_points(forest_data,sample_data[1,],"radio_fijo")
    plot_selection(gg_plot(),selected,"radio_fijo",all=all)
  })
  
  output$plot_variable <- renderPlot({
    all <- input$all_trees
    reset_sample()
    selected <- get_points(forest_data,sample_data[1,],"radio_variable")
    plot_selection(gg_plot(),selected,"radio_variable",all=all)
  })
  
  output$plot_relascopio <- renderPlot({
    all <- input$all_trees
    reset_sample()
    selected <- get_points(forest_data,sample_data[1,],"radio_relascopio")
    plot_selection(gg_plot(),selected,"radio_relascopio",all=all)
  })
  
  
  output$plot_fijo2 <- renderPlot({
    reset_sample()
    selected <- get_points(forest_data,sample_data[1,],"radio_fijo")
    plot_selection(gg_plot(),selected,"radio_fijo",tree_center = FALSE)
  })
  
  output$plot_variable2 <- renderPlot({
    reset_sample()
    selected <- get_points(forest_data,sample_data[1,],"radio_variable")
    plot_selection(gg_plot(),selected,"radio_variable",tree_center = FALSE)
    
  })
  
  output$plot_relascopio2 <- renderPlot({
    reset_sample()
    selected <- get_points(forest_data,sample_data[1,],"radio_relascopio")
    plot_selection(gg_plot(),selected,"radio_relascopio",tree_center = FALSE)
  })

  output$plot_selected1 <- renderPlot({
    reset()
    input$tipo
    N <- input$N
    n <- input$n
    reset_sample()
    field <- switch(input$tipo,
                    fijo = "radio_fijo",
                    variable = "radio_variable",
                    relascopio = "radio_relascopio"
    )
    selected <- get_points(forest_data,sample_data[1,],field)
    plot_selection(gg_plot(),selected,field,tree_center = TRUE)
  })

  output$plot_selected2 <- renderPlot({
    reset()
    input$tipo
    N <- input$N
    n <- input$n
    reset_sample()
    field <- switch(input$tipo,
                    fijo = "radio_fijo",
                    variable = "radio_variable",
                    relascopio = "radio_relascopio"
    )
    selected <- get_points(forest_data,sample_data[1,],field)
    plot_selection(gg_plot(),selected,field,tree_center = FALSE)
  })
  
  output$muestra <- renderTable({
    reset()
    input$tipo
    N <- input$N
    n <- input$n
    reset_sample()
    field <- switch(input$tipo,
                    fijo = "radio_fijo",
                    variable = "radio_variable",
                    relascopio = "radio_relascopio"
    )
    get_points(forest_data,sample_data[1,],field)
  })
  output$parametros_interes2<-renderTable({
    reset()
    parametros_interes(forest_data,input$lado)
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