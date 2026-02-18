# initial values

source("utils.R")

L <- 100
N <- 5
K<-1000
init_pop <- make_population(N,L)
init_samp_points <- sampling_points(K,L)
all_trees <- get_all_trees(init_pop,init_samp_points)
est <- n_estimaciones(all_trees,L,rotate=FALSE)
par_int <- parametros_interes(init_pop,L,TRUE)

rm(all_trees)
rm(init_samp_points)

server <- function(input, output, session) {
  # thematic::thematic_shiny()
  ##### reactive values #####
  
  forest <- reactiveVal(init_pop)
  est<-reactiveVal(est)
  pos <- reactiveVal(c(1))
  
  variation<-reactive({
    data_long <- pivot_longer(est()[,c("Type","N","G","V","h_media","dg","Ho")],
                              cols = c("N","G","V","h_media","dg","Ho"),
                              names_to = "parametro",values_to = "estimacion")
    
    data_long|> group_by(parametro,Type)|> 
      summarise(mean=mean(estimacion,na.rm=TRUE),sd=sd(estimacion,na.rm=TRUE)) 
  })
  
  table <- reactive({
    input$muestra
    positions <- pos()
    ests <- est()
    cols <- colnames(ests)
    ests<-ests[ests$Type==input$plot_type1,]
    n <- input$n
    k <- length(positions)
    reps <- k/n
    
    indexes <- match(positions,ests$Parc)
    tabla <- ests[indexes,]
    
    tabla$Parc <- rep(n:1,each=reps)
    tabla$Rep <- rep(reps:1,times=n)
    
    tabla[,c(cols,"Rep")]
  })
  par_int <- reactive({
    parametros_interes(forest(),input$lado,TRUE)
  })
  
  base_plot <- reactive({
    pop_plot(forest(),input$lado)
  })
  
  update_pop <- reactive({
    input$reset_pop
    pop <- make_population(input$N,input$lado)
    samp_points <- sampling_points(K,input$lado)
    
    trees <- get_all_trees(pop,samp_points)
    estimates <- n_estimaciones(trees,input$lado,rotate=FALSE)
    forest(pop)
    est(estimates)
    pos(1:input$n)
    
  })
  
  observeEvent(input$reset_pop,update_pop())
  observeEvent(input$N,update_pop())
  observeEvent(input$lado,update_pop())
  observeEvent(input$n,{
    pos(sample(1:K,input$n,replace = TRUE))
  })
  observeEvent(input$muestra,{
    pos(c(sample(1:K,input$n,replace = TRUE),pos()))
  })
  
  observeEvent(input$samp_dist,{
    new_val <- isolate(input$n)
    new_val <- new_val+1
    new_val <- ifelse(new_val>50,1,new_val)
    updateSelectInput(inputId = "n",selected = new_val)
  })

  
  ##### Population #####
  
  output$poblacion <- renderTable({
    forest()
  })
  
  output$plot_poblacion<-renderPlot({
    p <- base_plot()
    if(input$add_hd){
      p <- p +
        geom_label(aes(x=x,y=y-3,label=paste("d: ",dn)),size=4,fill="darkgreen",alpha=0.3)+
        geom_label(aes(x=x,y=y-8,label=paste("h: ",dn)),size=4,fill="blue",alpha=0.3)
    }
    p + ggtitle("Población")
  })
  
  output$tabla_interes1 <- renderTable({
    par_int()
  })
  
  
  
  ##### Seleccion #####
  output$plot_fijo <- renderPlot({
    trees <- get_trees(forest(),est()[pos()[1],],"r_fijo")
    plot_selection(base_plot(),trees,
                   all=input$all_trees,add_hd=input$add_hd)
  })

  output$plot_variable <- renderPlot({
    trees <- get_trees(forest(),est()[pos()[1],],"r_variable")
    plot_selection(base_plot(),trees,all=input$all_trees,add_hd=input$add_hd)
  })

  output$plot_relascopio <- renderPlot({
    trees <- get_trees(forest(),est()[pos()[1],],"r_relascopio")
    plot_selection(base_plot(),trees,all=input$all_trees,add_hd=input$add_hd)
  })


  output$plot_fijo2 <- renderPlot({
    trees <- get_trees(forest(),est()[pos()[1],],"r_fijo")
    plot_selection(base_plot(),trees,tree_center = FALSE,add_hd=input$add_hd)
  })

  output$plot_variable2 <- renderPlot({
    trees <- get_trees(forest(),est()[pos()[1],],"r_variable")
    plot_selection(base_plot(),trees,tree_center = FALSE,add_hd=input$add_hd)
  })

  output$plot_relascopio2 <- renderPlot({
    trees <- get_trees(forest(),est()[pos()[1],],"r_relascopio")
    plot_selection(base_plot(),trees,tree_center = FALSE,add_hd=input$add_hd)
  })

  # ##### One plot #####
  output$tabla_interes2<-renderTable({
    par_int()
  })
  output$muestra <- renderTable({
    samp<-get_trees(forest(),table()[1,],input$plot_type1)
    samp$Parc <- max(table()$Parc)
    samp
  })

  output$plot_selected1 <- renderPlot({
    trees <- get_trees(forest(),table()[1,],input$plot_type1)
    plot_selection(base_plot(),trees,
                   tree_center = input$centered=="arbol",add_hd=input$add_hd)
  })

  output$estimacion1<-renderTable({
    tabla <- table()
    tabla <- tabla
    tabla$Parc <- tabla$Rep
    tabla
  })

  output$tabla_acc <- renderTable({
    tabla <- table()
    tabla <- tabla[tabla$Parc==1,]
    tabla$Parc <- tabla$Rep
    tabla
    })

  output$plot_res1 <- renderPlot({
    tabla <-table()
    tabla <- tabla[tabla$Parc==1,]
    tabla$Parc <- tabla$Rep
    add_samples_plot(par_int(),tabla)
  })



  # ##### n plots #####

  output$tabla_interes3<-renderTable({
    par_int()
  })
  
  output$plot_selected2<-renderPlot({
    
    selected <- table() |> filter(Rep==max(Rep))
    selected_list <- group_split(selected,Rep,Parc)
    forest_all <- forest()
    type <- input$plot_type1
    selected_trees <- map_dfr(selected_list,function(x,population,type){
      res<- get_trees(population,x,type)
      res$Parc <- x$Parc[1]
      res$Rep <- x$Rep[1]
      res
    },population=forest_all,type=type)
    
    plot_n_selections(base_plot(),selected_trees,
                        tree_center = input$centered=="arbol",
                        all=input$all_trees)
  })
  
  output$n_estimaciones<-renderTable({
    table() |> filter(Rep==max(Rep))
  })
  
  output$plot_res2<- renderPlot({
    add_samples_n_plots(par_int(),table())
  })
  
  # ##### samp dist #####
  
  output$tabla_interes4<-renderTable({
    par_int()
  })
  
  output$plot_selected3<-renderPlot({
    base_plot()
  })
  

  output$var_n<- renderPlot({
    a<-data.frame(n=1:50,id=1)
    var <- variation()
    var <- var[var$Type==input$plot_type1,]
    var <- merge(var,a)
    var$sd_n <- ((var$sd)^2)/var$n
    var$color <- ifelse(var$n==input$n,"red","black")
    
    red <- var[var$color=="red",]
    
    ggplot(var,aes(x=n,y=sd_n))+facet_wrap(.~parametro,scales="free_y")+
      geom_point(aes(color=color))+geom_path() + 
      geom_point(data=red,aes(color=color),pch=20,size=3)+
      xlab("Varianza estimador final")+
      scale_color_manual(values=c("red"="red","black"="black"))+
      guides(color="none")+
      ggtitle("Cambio en la varianza al aumentar n")
  })
  
  output$normal_approx<-renderPlot({
    estimates <- est()
    n <- input$n
    type <- input$plot_type1
    variation <- variation()
    normal_approx(estimates,par_int(),n,type,variation,K)
  })
  
  
}


