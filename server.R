# initial values

source("functions.R")
init_pop <- make_population(20,100)
init_samp_points <- sampling_points(10,100)
trees <- get_n_points(init_pop,init_samp_points,"r_fijo")
est_n <- n_estimaciones(trees,100,rotate=FALSE)
names_est_n <- colnames(est_n)
est_n$Rep <- 1
server <- function(input, output, session) {
  # thematic::thematic_shiny()
  ##### reactive values #####
  
  forest <- reactiveVal(init_pop)
  samp_points<-reactiveVal(init_samp_points)
  
  
  est <- reactiveVal({est_n})
  
  par_int <- reactive({parametros_interes(forest(),input$lado,TRUE)})
  base_plot <- reactive({
    pop_plot(forest(),input$lado)
  })
  
  

  
  update_pop <- reactive({
    new_data <- make_population(input$N,input$lado)
    new_points <- sampling_points(input$n,input$lado)
    print(par_int())
    print(new_points)
    forest(new_data)
    samp_points(new_points)
    update_est()
  })
  
  update_est <- reactive({
    trees <- get_n_points(forest(),samp_points(),input$plot_type1)
    new_estimates <- n_estimaciones(trees,input$lado)
    new_estimates$Rep <- 1
    new_estimates <- new_estimates[,c("Rep",names_est_n)]
    est(new_estimates)
  })
  
  
  #
  
  add_estimate<- reactive({
    new_point <- sampling_points(input$n,input$lado)
    trees <- get_n_points(forest(),samp_points(),input$plot_type1)
    new_estimates <- n_estimaciones(trees,input$lado)
    new_estimates$Rep <- est()$Rep[1]+1
    new_estimates <- new_estimates[,c("Rep",names_est_n)]
    print(new_point)
    print(est())
    old <- est()
    new <- rbind(new_estimates,old)
    est(new)
    samp_points(new_point)
  })
  
  
  observeEvent(input$muestra,{
    add_estimate()})
  observeEvent(input$n_muestras,add_estimate())
  observeEvent(input$n,update_pop())
  observeEvent(input$plot_type1,{
    update_est()
    })
  observeEvent(input$N,update_pop())
  observeEvent(input$lado,update_pop())
  observeEvent(input$reset_pop,update_pop())
  
  ##### Population #####
  
  output$poblacion <- renderTable({
    forest()[,c(1:5)]
  })
  
  output$plot_poblacion<-renderPlot({
    print(input$lado)
    p <- base_plot()
    if(input$add_hd){
      p <- p +
        geom_label(aes(x=x,y=y-3,label=paste("d: ",diam)),size=4,fill="darkgreen",alpha=0.3)+
        geom_label(aes(x=x,y=y-8,label=paste("h: ",diam)),size=4,fill="blue",alpha=0.3)
    }
    p + ggtitle("Población")
  })
  
  output$tabla_interes1 <- renderTable({
    print(par_int())
    par_int()
  })
  
  ##### Seleccion #####
  output$plot_fijo <- renderPlot({
    all <- input$all_trees
    print(samp_points()[1,])
    selected <- get_trees(forest(),samp_points()[1,],"r_fijo")
    plot_selection(base_plot(),selected,samp_points()[1,],"r_fijo",
                   all=all,add_hd=input$add_hd)
  })

  output$plot_variable <- renderPlot({
    all <- input$all_trees
    selected <- get_trees(forest(),samp_points()[1,],"r_variable")
    plot_selection(base_plot(),selected,samp_points()[1,],"r_variable",
                   all=all,add_hd=input$add_hd)
  })

  output$plot_relascopio <- renderPlot({
    all <- input$all_trees
    selected <- get_trees(forest(),samp_points()[1,],"r_relascopio")
    plot_selection(base_plot(),selected,samp_points()[1,],"r_relascopio",
                   all=all,add_hd=input$add_hd)
  })


  output$plot_fijo2 <- renderPlot({
    selected <- get_trees(forest(),samp_points()[1,],"r_fijo")
    plot_selection(base_plot(),selected,samp_points()[1,],"r_fijo",
                   tree_center = FALSE,add_hd=input$add_hd)
  })

  output$plot_variable2 <- renderPlot({
    selected <- get_trees(forest(),samp_points()[1,],"r_variable")
    plot_selection(base_plot(),selected,samp_points()[1,],"r_variable",
                   tree_center = FALSE,add_hd=input$add_hd)

  })

  output$plot_relascopio2 <- renderPlot({
    selected <- get_trees(forest(),samp_points()[1,],"r_relascopio")
    plot_selection(base_plot(),selected,samp_points()[1,],"r_relascopio",
                   tree_center = FALSE,add_hd=input$add_hd)
  })
  # 
  # ##### One plot #####
  output$tabla_interes2<-renderTable({
    par_int()
  })
  output$muestra <- renderTable({
    samp<-get_trees(forest(),samp_points()[1,],input$plot_type1)
    samp[,-c(2:5)]
  })

  output$plot_selected1 <- renderPlot({
    selected <- get_trees(forest(),samp_points()[1,],input$plot_type1)
    plot_selection(base_plot(),selected,samp_points()[1,],type=input$plot_type1,tree_center = FALSE)
  })

  output$estimacion1<-renderTable({
    est()[1,]
  })

  output$tabla_acc <- renderTable({
    est()[est()$Parc==1,]
  })

  output$plot_res1 <- renderPlot({
    first <- est()|> group_by(Rep)|> filter(Parc==1)
    add_samples_plot(par_int(),first)

  })



  # ##### n plots #####

  output$tabla_interes3<-renderTable({
    par_int()
  })
  
  output$plot_selected2<-renderPlot({
      selected <- get_n_points(forest(),samp_points(),input$plot_type2)
      print(samp_points())
      plot_n_selections(base_plot(),selected,samp_points(),type=input$plot_type2,tree_center = FALSE)
  })
  
  output$plot_res2<- renderPlot({
    add_samples_n_plots(par_int(),est())
  })
  
  
}


