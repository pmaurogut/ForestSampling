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
    diam = round(runif(Npop,5,40),1)
    )
  
  res$ht <- round(5+(res$diam-5)*0.5,1)
  
  res$r_fijo <- 15
  res$area_fijo <-  pi*(res$r_fijo^2)/10000
  res$fac_exp_fijo<- 1/res$area_fijo
  
  res$r_variable <- ifelse(res$diam<15,10,20)
  res$area_variable <-  pi*(res$r_variable^2)/10000
  res$fac_exp_variable <- 1/res$area_variable
  
  res$r_relascopio <- res$diam/2
  res$area_relascopio <-  pi*(res$r_relascopio^2)/10000
  res$fac_exp_relascopio <- 1/res$area_relascopio
  return(res)
}

sampling_points <-function(n,L){
  data.frame(Parc=1:n,x=runif(n,0,L),y=runif(n,0,L))
  }

get_points <- function(population,point,type){
  
  pick <- sqrt((population$x-point$x)^2+(population$y-point$y)^2)<=population[[type]]
  if(all(!pick)){
    return(data.frame(
      type=type,x0=NA,y0=NA,x=NA,y=NA,Parc=point$Parc,
      NA,diam=NA,ht=NA,gi_m2=NA,A_parc_ha=NA,EXP_FAC=NA))
  }else{
    res<-population[pick, ]
    res$Parc<- point$Par
    res$type <- gsub("radio_","",type)
    res$x0 <- point$x
    res$y0 <- point$y
    res$gi_m2 <- (pi*(1/4)*res$diam^2)/10000
    res<-res[,c("type","x0","y0","x","y","Parc","diam","ht","gi_m2",type)]
    res$area_parc_ha<-(pi*res[,type]^2)/10000
    res$EXP_FAC <- 1/res$area_parc_ha
    res<-res[order(res$diam,decreasing = TRUE),]
    colnames(res)<-gsub("area_","A_",colnames(res))
    return(res)
  }
  
}

get_n_points <- function(population,points,type){
  
  res<-list()
  for(i in 1:dim(points)[1]){
    res[[i]]<-get_points(population,points[i,],type)
  }
  map_dfr(res,function(x){x})
}

parametros_interes <- function(poblacion, lado,rotate=TRUE){
  A<-lado*lado/10000
  res<-data.frame(
    Area_ha = A,
    Total_N = length(poblacion[[1]]),
    Total_G = (1/10000)*sum(pi*poblacion$diam^2)/4,
    Total_h = sum(poblacion$ht)
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
  if(rotate){
    names<- colnames(res)
    res <- as.data.frame(t(res))
    colnames(res)<-"Valor"
    res$parametro<-names
    res
  }else{
    res
  }
  
}

estimacion <- function(sample,lado,rotate=TRUE){
  
  if(is.na(sample$diam[1])){
    res <- data.frame(Parc=sample$Parc[1],Total_N=0,Total_G=0,
                      Total_h=0,N=0,G=0,hmedia=0,dg=0,Ho=0)
    if(rotate){
      names <- colnames(res)
      res <- data.frame(t(res))
      colnames(res)[1]<-"Estimacion"
      res$Variable <- names
      return(res[,c("Variable","Estimacion")])
    }else{
      return(res)
    }
    
    
    
  }
  sample <- sample[order(sample$diam,decreasing = TRUE),]
  sample$cum_sum<-cumsum(sample$EXP_FAC)
  
  A <- (lado*lado)/10000
  res<-data.frame(Parc=sample$Parc[1])
  res$Total_N <- sum(sample$EXP_FAC)*A
  res$Total_G <- sum(sample$EXP_FAC*sample$gi_m2)*A
  res$Total_h <- sum(sample$EXP_FAC*sample$ht)*A
  res$N <- res$Total_N/A
  res$G<- res$Total_G/A
  res$hmedia<- res$Total_h/res$Total_N
  res$dg<-sqrt((res$G/res$N)*(4/pi))
  if(sum(sample$EXP_FAC)>=100){
    last <- which(sample$cum_sum>100)[1]
    s2 <- sample[1:last,]
    s2[last,"EXP_FAC"]<-s2[last,]$cum_sum-100
    res$Ho<-sum(s2$EXP_FAC*sample$ht)/sum(s2$EXP_FAC)
    
  }else{
    res$Ho <- sum(sample$EXP_FAC*sample$ht)/sum(sample$EXP_FAC)
  }
  if(rotate){
    names <- colnames(res)
    res <- data.frame(t(res))
    colnames(res)[1]<-"Estimacion"
    res$Variable <- names
    return(res[,c("Variable","Estimacion")])
  }else{
    return(res)
  }
  
}

n_estimaciones<-function(sample,lado,rotate=FALSE){
  map_dfr(group_split(sample,Parc),estimacion,lado=lado,rotate=rotate)
}

plot_selection <- function(p,selected,type,tree_center=TRUE,all=FALSE){
  
  title <- switch(type,
                  r_fijo = "Radio fijo 15 m",
                  r_variable = "R anidados d<15 cm 10, d>=15 cm 20m",
                  r_relascopio = "Relascopio BAF=1"
  )
  if(all){
    p <- p + geom_circle(aes(x0=x,y0=y,r=.data[[type]]), fill="grey50",alpha=0.2)
  }
  
  if(!is.na(selected$diam[1])){
    if(tree_center){
      p <- p  + geom_circle(data=selected,aes(x0=x,y0=y,r=.data[[type]]),fill="purple",alpha=0.2)
    }else{
      selected2 <- selected |> group_by(!! sym(type)) |> filter(row_number()==1) |> ungroup()
      p <- p  + geom_circle(data=selected2,aes(x0=x0,y0=y0,r=.data[[type]]),fill="purple",alpha=0.2)  
    }
    p <- p +geom_point(data=selected,aes(x=x,y=y),shape=20,col="green")
  }
  p <- p + geom_point(data=samp_points[1,],aes(x=x,y=y),shape=13,col="red",size=4)
  p <- p + guides(fill=FALSE)+ggtitle(title)
  p
}

plot_n_selections <- function(p,selected,type,tree_center=TRUE,all=FALSE){
  
  title <- switch(type,
                  r_fijo = "Radio fijo 15 m",
                  r_variable = "R anidados d<15 cm 10, d>=15 cm 20m",
                  r_relascopio = "Relascopio BAF=1"
  )
  if(all){
    p <- p + geom_circle(aes(x0=x,y0=y,r=.data[[type]]), fill="grey50",alpha=0.2)
  }
  
  if(!all(is.na(selected$diam))){
    if(tree_center){
      p <- p  + geom_circle(data=selected,aes(x0=x,y0=y,r=.data[[type]],fill= factor(Parc)),alpha=0.2)
    }else{
      selected2 <- selected |> group_by(Parc,!! sym(type)) |> filter(row_number()==1) |> ungroup()
      p <- p  + geom_circle(data=selected2,aes(x0=x0,y0=y0,r=.data[[type]],fill= factor(Parc)),alpha=0.2)  
    }
    p <- p +geom_point(data=selected,aes(x=x,y=y,col=factor(Parc)),shape=20)
  }
  p <- p + geom_point(data=samp_points,aes(x=x,y=y),shape=13,col="red",size=4)
  p <- p + guides(fill=FALSE,color=FALSE)+ggtitle(title)
  p
}


forest_data <- make_population(20,100)
samp_points <- sampling_points(10,100)
sample <- get_points(forest_data,samp_points[1,],"r_fijo")
sample_n <- get_n_points(forest_data,samp_points,"r_fijo")
est <- estimacion(sample,100,rotate=FALSE)
n_est <- n_estimaciones(sample_n,100,rotate=FALSE)

space <- br()
plot_type <- radioButtons("tipo", "Tipo de Parc:",
                          c("Radio fijo 15 m" = "fijo",
                            "Radios anidados d<15 10m, d>=15 10m " = "variable",
                            "Relascopio BAF=1" = "relascopio"),selected = "fijo")
plot_type2 <- radioButtons("tipo2", "Tipo de Parc:",
                          c("Radio fijo 15 m" = "fijo",
                            "Radios anidados d<15 10m, d>=15 10m " = "variable",
                            "Relascopio BAF=1" = "relascopio"),selected = "fijo")

parametro_interes<- radioButtons("parametro_interes", "Parametro a estimar",
                          c("N(pies/ha)" = "N",
                            "G(m2/ha)" = "variable",
                            "h_media(m)" = "hmedia",
                            "dg(cm)" = "dg",
                            "Ho(m)" = "Ho"
                            ),selected = "G")

lado <- sliderInput("lado","Lado (m)",value = 100,min = 100,max = 500,step=50)

pop_size <- sliderInput("N","Número de árboles/ha",value = 20,min = 5,max = 500,step=5)

samp_size <- sliderInput("n","Número de Parcelas",value = 1,min = 1,max = 50)

reps <-sliderInput("r","Repeticiones",value = 100,min = 1,max = 200)

reset_population <- actionButton("reset_pop", "Regenera poblacion")

muestra <- actionButton("muestra", "Toma una muestra",color="darkgreen",alpha=0.4)

n_muestras <- actionButton("n_muestras", "Toma n muestras",color="blue",alpha=0.4)

areas_inclusion <- checkboxInput("all_trees","Todas las areas de inclusion",value = FALSE)

controls <- list(lado,pop_size,samp_size,reps,space,
                 areas_inclusion,space,reset_population,space,muestra,space,n_muestras)




# Define UI for random distribution app ----
##### UI #####
{
  ui <- page_sidebar(
    
    title = "Muestreo forestal",
    # Sidebar layout with input and output definitions ----
    sidebar=sidebar(title = "Opciones población y muestra",controls,open="always"),
    
    navset_card_underline(
      title = "Ejemplos",
      nav_panel("Población y parámetros de interés",
        
        layout_columns(col_widths=c(6,4,2),
                     card(card_header("Mapa Población"),plotOutput("plot_poblacion",width=800,height=800)),
                     card(card_header("Datos Población"),tableOutput('poblacion')),
                     card(card_header("Parámetros de interés"),tableOutput('tabla_interes1'))
          )
      ),
    
      nav_panel("Selección de muestras",
        layout_columns(col_widths=c(4,4,4,4,4,4),
                     card(plotOutput("plot_fijo",width=500,height=500)),
                     card(plotOutput("plot_variable",width=500,height=500)),
                     card(plotOutput("plot_relascopio",width=500,height=500)),
                     card(plotOutput("plot_fijo2",width=500,height=500)),
                     card(plotOutput("plot_variable2",width=500,height=500)),
                     card(plotOutput("plot_relascopio2",width=500,height=500))
          )
      ),
    
      nav_panel("Estimación una Parcela",
          layout_columns(col_widths=c(5,4,3),
                     card(card_header("Tipo de Parcela"),
                         plot_type,
                         tableOutput('tabla_interes2'),
                         tableOutput('muestra'),
                         
                      ),
                     card(card_header("Selección árboles"),
                         plotOutput("plot_selected1",width=400,height=400),
                         tableOutput("estimacion1")
                      ),
                     card(card_header("Estimaciones"),
                         plotOutput("plot_res1",width=400,height=400),
                         plotOutput("plot_res2",width=400,height=400),
                         tableOutput("tabla_acc")
                      )
            )
        ),
      
      nav_panel("Estimación múltiples Parcelas",
                layout_columns(col_widths=c(5,4,3),
                      card(card_header("Estimación múltiples Parcelass"),
                        plot_type2,
                        tableOutput('tabla_interes3'),
                        tableOutput('n_muestras')
                      ),
                      card(card_header("Selección árboles"),
                         plotOutput("plot_selected_n1",width=500,height=500),
                         plotOutput("tabla_estimacion_n",width=500,height=500)
                        ),
                      card(card_header("Estimaciones"),
                          plotOutput("plot_res_n1",width=500,height=500),
                          plotOutput("plot_res_n2",width=500,height=500),
                          tableOutput("estimacion_n")
                      )
                )
        )
    )
  )
}
###### Define server logic for random distribution app ---- #####
{
  server <- function(input, output) {
    
    
    gg_plot <- reactive({
      input$N
      input$reset_pop
      input$lado
      lado <- input$lado
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
      input$reset_pop
      input$N
      input$lado
      N <- input$N
      lado <- input$lado
      forest_data <<- make_population(N,lado)
      reset_accum()
    })
    
    reset_sample<-reactive({
      input$muestra
      input$tipo
      input$tipo2
      input$n
      input$reset_pop
      input$N
      input$lado
      samp_points <<-sampling_points(input$n,input$lado)
    })
    
    
    reset_accum<-reactive({
      input$reset_pop
      input$lado
      input$N
      input$tipo
      input$tipo1
      field <- switch(input$tipo,
                      fijo = "r_fijo",
                      variable = "r_variable",
                      relascopio = "r_relascopio"
      )
      est<<-estimacion(get_points(forest_data,samp_points[1,],field),input$lado,rotate=FALSE)
    })
    
    add_accum<-reactive({
      input$muestra
      reset_sample()
      field <- switch(input$tipo,
                      fijo = "r_fijo",
                      variable = "r_variable",
                      relascopio = "r_relascopio"
      )

      est<<-rbind(est,estimacion(get_points(forest_data,samp_points[1,],field),input$lado,rotate=FALSE))
    })

    observeEvent(input$muestra,add_accum())

    ##### Population #####
    
    output$poblacion <- renderTable({
      reset()
      forest_data[,c(1:5)]
    })
    
    output$plot_poblacion<-renderPlot({
      reset()
      p <- gg_plot()
      p +
        geom_label(aes(x=x,y=y-3,label=diam),size=3,fill="darkgreen",alpha=0.3)+
        geom_label(aes(x=x,y=y-8,label=ht),size=3,fill="blue",alpha=0.3)+
        ggtitle("Población")
    })
    
    output$tabla_interes1 <- renderTable({
      reset()
      parametros_interes(forest_data,input$lado)
    })
    
    ##### Seleccion #####
    output$plot_fijo <- renderPlot({
      reset()
      reset_sample()
      all <- input$all_trees
      selected <- get_points(forest_data,samp_points[1,],"r_fijo")
      plot_selection(gg_plot(),selected,"r_fijo",all=all)
    })
    
    output$plot_variable <- renderPlot({
      reset()
      reset_sample()
      all <- input$all_trees
      selected <- get_points(forest_data,samp_points[1,],"r_variable")
      plot_selection(gg_plot(),selected,"r_variable",all=all)
    })
    
    output$plot_relascopio <- renderPlot({
      reset()
      reset_sample()
      all <- input$all_trees
      selected <- get_points(forest_data,samp_points[1,],"r_relascopio")
      plot_selection(gg_plot(),selected,"r_relascopio",all=all)
    })
    
    
    output$plot_fijo2 <- renderPlot({
      reset()
      reset_sample()
      selected <- get_points(forest_data,samp_points[1,],"r_fijo")
      plot_selection(gg_plot(),selected,"r_fijo",tree_center = FALSE)
    })
    
    output$plot_variable2 <- renderPlot({
      reset()
      reset_sample()
      selected <- get_points(forest_data,samp_points[1,],"r_variable")
      plot_selection(gg_plot(),selected,"r_variable",tree_center = FALSE)
      
    })
    
    output$plot_relascopio2 <- renderPlot({
      reset()
      reset_sample()
      selected <- get_points(forest_data,samp_points[1,],"r_relascopio")
      plot_selection(gg_plot(),selected,"r_relascopio",tree_center = FALSE)
    })
    
    ##### One plot #####
    output$tabla_interes2<-renderTable({
      reset()
      parametros_interes(forest_data,input$lado)
    })
    output$muestra <- renderTable({
      reset()
      reset_sample()
      field <- switch(input$tipo,
                      fijo = "r_fijo",
                      variable = "r_variable",
                      relascopio = "r_relascopio"
      )
      samp<-get_points(forest_data,samp_points[1,],field)
      samp[,-c(2:5)]
    })
    
    output$plot_selected1 <- renderPlot({
      reset()
      reset_sample()
      field <- switch(input$tipo,
                      fijo = "r_fijo",
                      variable = "r_variable",
                      relascopio = "r_relascopio"
      )
      selected <- get_points(forest_data,samp_points[1,],field)
      plot_selection(gg_plot(),selected,type=field,tree_center = FALSE)
    })
    
    output$estimacion1<-renderTable({
      reset()
      reset_sample()
      field <- switch(input$tipo,
                      fijo = "r_fijo",
                      variable = "r_variable",
                      relascopio = "r_relascopio"
      )
      selected <- get_points(forest_data,samp_points[1,],field)
      estimacion(selected,input$lado)
    })
    # output$plot_selected2 <- renderPlot({
    #   reset()
    #   input$tipo
    #   N <- input$N
    #   n <- input$n
    #   reset_sample()
    #   field <- switch(input$tipo,
    #                   fijo = "r_fijo",
    #                   variable = "r_variable",
    #                   relascopio = "r_relascopio"
    #   )
    #   selected <- get_points(forest_data,samp_points[1,],field)
    #   plot_selection(gg_plot(),selected,type=field,tree_center = FALSE)
    # })
    output$tabla_acc <- renderTable({
      reset()
      reset_sample()
      est
    })
    
    output$plot_res1 <- renderPlot({
      reset()
      reset_sample()
      par_int <- parametros_interes(forest_data,input$lado,rotate=FALSE)
      print(par_int)
      est$Parc <- 1:dim(est)[1]
      print(est)
      ggplot(est,aes(x=G))+geom_histogram()+geom_vline(data=par_int,aes(xintercept=G),col="red")
      
    })

    
    output$plotaverage<- renderPlot({
      type <- input$tipo
      N <- input$N
      n <- input$n
      
      plot(forest_data$x,forest_data$y,
           main = paste("Parcs ", type, sep = ""),
           pch=20)
    })
    ##### n plots #####
    output$n_muestras <- renderTable({
      field <- switch(input$tipo2,
                      fijo = "r_fijo",
                      variable = "r_variable",
                      relascopio = "r_relascopio"
      )
      samp<-get_n_points(forest_data,samp_points,field)
      samp[,-c(2:5)]
    })
    
    output$plot_selected_n1 <- renderPlot({
      field <- switch(input$tipo2,
                      fijo = "r_fijo",
                      variable = "r_variable",
                      relascopio = "r_relascopio"
      )
      selected <- get_n_points(forest_data,samp_points,field)
      plot_n_selections(gg_plot(),selected,type=field,tree_center = TRUE)
    })
    # 
    # output$plot_selected_n2 <- renderPlot({
    #   reset()
    #   input$tipo
    #   N <- input$N
    #   n <- input$n
    #   reset_sample()
    #   field <- switch(input$tipo,
    #                   fijo = "r_fijo",
    #                   variable = "r_variable",
    #                   relascopio = "r_relascopio"
    #   )
    #   selected <- get_n_points(forest_data,samp_points,field)
    #   plot_n_selections(gg_plot(),selected,type=field,tree_center = FALSE)
    # })
    
    
    output$tabla_interes3<-renderTable({
      reset()
      parametros_interes(forest_data,input$lado)
    })
    
    output$plotaverage_n<- renderPlot({
      
      reset()
      type <- input$tipo
      N <- input$N
      n <- input$n
      
      plot(forest_data$x,forest_data$y,
           main = paste("Parcs ", type, sep = ""),
           pch=20)
    })
    
    
    
    # Generate a summary of the data ----
    output$distmuest <- renderPrint({
      summary(forest())
    })
    
    
    
  }
  
}

# Create Shiny app ----
shinyApp(ui, server)