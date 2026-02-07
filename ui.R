library(shiny)
library(bslib)
library(ggplot2)
library(ggforce)
library(dplyr)
library(purrr)
library(tidyr)
library(DT)
library(thematic)


plot_type<- radioButtons("plot_type1" , "Tipo de parcela:",
                         c("R fijo 15 m" = "r_fijo",
                           "R. anidados d<15 10m, d>=15 10m" = "r_variable",
                           "Relascopio BAF=1" = "r_relascopio"),selected = "r_fijo")

centrado_arbol <- radioButtons("centered" , "Centrar en",
                              c("Arbol" = "arbol",
                                "Punto" = "punto"),selected = "punto")

space <- br()

lado <- sliderInput("lado","Lado (m)",value = 100,min = 100,max = 500,step=50)

pop_size <- sliderInput("N","Número de árboles/ha",value = 20,min = 5,max = 500,step=5)

samp_size <- sliderInput("n","Número de Parcelas",value = 1,min = 1,max = 50)

reps <-sliderInput("r","Repeticiones",value = 100,min = 1,max = 200)

reset_population <- actionButton("reset_pop", "Regenera poblacion")

muestra <- actionButton("muestra", "Muestrea",color="darkgreen",alpha=0.4)

n_muestras <- actionButton("n_muestras", "Toma n muestras",color="blue",alpha=0.4)

areas_inclusion <- checkboxInput("all_trees","Todas las areas de inclusion",value = FALSE)
add_hd<- checkboxInput("add_hd","Añade altura y diámetro",value = FALSE)



controls <- list(lado,pop_size,samp_size,plot_type,space,
                 centrado_arbol,areas_inclusion,add_hd,space,
                 reset_population,space,muestra)



#### UI ####
ui <- page_navbar(
  
  theme=bs_theme(version=5,preset = "darkly"),
  # tags$head(
  #   tags$style(HTML("
  #     .bslib-page-title {
  #       color: black;
  #       background-color: white;
  #       background: green;
  #       padding: 10px;
  #       bs-emphasis-color: white;
  #     }
  #   "))
  # ),
  title = "Muestreo forestal",
  nav_spacer(),
  sidebar=sidebar(title = "Opciones población y muestra",controls,open="always"),
  
  # navset_card_underline(
  #   title = "Ejemplos",
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
    
    nav_panel("Estimación una parcela",
              layout_columns(col_widths=c(5,7),
                             card(card_header("Parámetros de interés y muestra"),
                                  card(layout_columns(col_widths=c(5,7),
                                                      tableOutput('tabla_interes2'),
                                                      plotOutput("plot_selected1",width=400,height=400)
                                                      )),
                                  card(tableOutput('muestra'),min_height = 450),
                                  
                             ),
                             card(card_header("Estimaciones"),
                                  card(plotOutput("plot_res1",width=950,height=650),min_height = 600),
                                  card(tableOutput("tabla_acc"),min_height = 300)
                             )
              )
    ),
    
    nav_panel("Estimación múltiples parcelas",
              layout_columns(col_widths=c(5,7),
                             card(card_header("Parámetros de interés y muestra"),
                                  card(layout_columns(col_widths=c(5,7),
                                                      tableOutput('tabla_interes3'),
                                                      plotOutput("plot_selected2",width=400,height=400)
                                  )),
                                  card(tableOutput('n_estimaciones'),min_height = 450),
                                  
                             ),
                             card(card_header("Estimación con una parcela vs estimación con n parcelas"),
                                  card(plotOutput("plot_res2",width=950,height=750))
                             )
              )
    ),
    
    nav_panel("Distribución muestral y errores",
              layout_columns(col_widths=c(5,4,3),
                             card(card_header("Tipo de Parcela"),
                                  card(),
                                  card(tableOutput('tabla_interes4')),
                                  card(tableOutput('muestra3')),
                                  
                             ),
                             card(card_header("Selección árboles"),
                                  card(plotOutput("plot_selected3",width=400,height=400)),
                                  card(tableOutput("estimacion3"))
                             ),
                             card(card_header("Estimaciones"),
                                  card(plotOutput("plot_res3",width=400,height=200),min_height = 400),
                                  card(tableOutput("tabla_acc3"),min_height = 400)
                             )
              )
    )
  # )
)