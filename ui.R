library(shiny)
library(bslib)
library(ggplot2)
library(ggforce)
library(dplyr)
library(purrr)
library(tidyr)
library(DT)
library(thematic)
library(htmltools)
library(markdown)

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

samp_dist <- actionButton("samp_dist", "Aumenta n")
reps <- sliderInput("reps","Replicas",value = 3,min = 1,max = 5,step=1)

controls <- list(lado,pop_size,samp_size,plot_type,space,
                 centrado_arbol,areas_inclusion,add_hd,space,
                 reset_population,muestra,samp_dist)



#### UI ####
ui <- page_navbar(
  
  # theme = bs_theme(version = 5, bootswatch = "darkly"),
  title = div(
    img(src = "logo.png", height = "30px"), 
    "Muestreo forestal"
  ),
  
  theme = bs_theme(version = 5, bootswatch = "darkly",
                   "navbar-bg" = "#416e5e",
                   "nav-link-color" = "#60d1b8 !important"
                   )|> 
    bslib::bs_add_rules(
      rules = "
                    .navbar.navbar-default {
                        background-color: #416e5e;
                    }
                    
                    "
    )|>
    bslib::bs_add_rules(".custom-header { background-color: #416e5e ; color: #60d1b8 !important; }"),
  
  
  tags$head(tags$script(
    HTML('
              $(document).ready(function() {
                $(".navbar-brand").replaceWith(
                  $("<a class = \'navbar-brand\' href = \'#\'></a>")
                );
                var containerHeight = $(".navbar .container-fluid").height() + "px";
                $(".navbar-brand")
                  .append(
                    "<img id = \'myImage\' src=\'logo.png\'" +
                    " height = " + containerHeight + ">"
                  );
                });'
    )
  )
  ),
  tags$style(HTML(
    '@media (max-width:992px) { .navbar-brand { padding-top: 0; padding-bottom: 0; }}'
  )
  ),
  
  
  nav_spacer(),
  sidebar=sidebar(title = "Opciones población y muestra",controls,open="always"),
  
  # navset_card_underline(
  #   title = "Ejemplos",
  #### Poblacion ####
  nav_panel("1. Población",
        
        fluidRow({
          layout_columns(col_widths=c(6,4,2),
                         card(card_header("Mapa Población"),plotOutput("plot_poblacion",width=800,height=800)),
                         card(card_header("Datos Población"),tableOutput('poblacion')),
                         card(card_header("Parámetros de interés"),tableOutput('tabla_interes1')),max_height = 800
          )
        }),
        
        fluidRow({   
          card(
            card_header("Explicación",class="custom-header"),
            withMathJax(htmltools::includeMarkdown("help/Poblacion.Rmd")),full_screen = TRUE#, height=200
          )
        })
    ),
  #### Sample selection ####
  nav_panel("2. Selección de muestras",
            fluidRow({
              layout_columns(col_widths=c(4,4,4,4,4,4),
                             card(plotOutput("plot_fijo",width=500,height=500)),
                             card(plotOutput("plot_variable",width=500,height=500)),
                             card(plotOutput("plot_relascopio",width=500,height=500)),
                             card(plotOutput("plot_fijo2",width=500,height=500)),
                             card(plotOutput("plot_variable2",width=500,height=500)),
                             card(plotOutput("plot_relascopio2",width=500,height=500))
              )
            }),
            fluidRow({
              card(
                card_header("Explicación",class="custom-header"),
                withMathJax(htmltools::includeMarkdown("help/Seleccion.Rmd")),full_screen = TRUE#, height=200
              )
            })
    ),
  #### One plot ####
    nav_panel("3. Estimacion 1 parcela",
            fluidRow(
              layout_columns(col_widths=c(5,7),height = 800,
                             card(card_header("Parámetros de interés y muestra"),
                                  card(layout_columns(col_widths=c(5,7),row_heights = 800,
                                                      tableOutput('tabla_interes2'),
                                                      plotOutput("plot_selected1",width=400,height=400)
                                  )),
                                  card(tableOutput('muestra'),min_height = 200),
                                  
                             ),
                             card(card_header("Estimaciones"),
                                  card(plotOutput("plot_res1",width=950,height=650),min_height = 600),
                                  card(tableOutput("tabla_acc"),min_height = 300)
                             )
              )
            ),
            fluidRow(
              card(
                card_header("Explicación",class="custom-header"),
                withMathJax(htmltools::includeMarkdown("help/OnePlot.Rmd")),full_screen = TRUE#, height=200
              )
            )
              
    ),
  #### n plots ####
    nav_panel("4. Estimación con n parcelas",
            fluidRow({
              layout_columns(col_widths=c(5,7),height = 800,
                             card(card_header("Parámetros de interés y muestra"),
                                  card(layout_columns(col_widths=c(5,7),
                                                      tableOutput('tabla_interes3'),
                                                      plotOutput("plot_selected2")
                                  ), min_height=350),
                                  card(tableOutput("n_estimaciones")),
                             ),
                             card(card_header("Estimación con una parcela vs estimación con n parcelas"),
                                  card(plotOutput("plot_res2"))
                             )
              )
            }),
            
            fluidRow(
              card(
                card_header("Explicación",class="custom-header"),
                withMathJax(htmltools::includeMarkdown("help/nPlots.Rmd")),full_screen = TRUE#, height=200
              )
            )
              
    ),
  #### Samp dist ####
  nav_panel("5. Distribución muestral",
          fluidRow({
            layout_columns(col_widths=c(5,7),height = 900,
                           card(
                             card_header("Parámetros de interés y muestra"),
                                card(layout_columns(col_widths=c(5,7),
                                                    tableOutput('tabla_interes4'),
                                                    plotOutput("plot_selected3")
                                ),height=500),
                                card(card_header("Cambio en la desviación típìca al aumentar n"),
                                     plotOutput("var_n")
                                ),
                           ),
                           card(card_header("Aproximación a una normal"),
                                card(plotOutput("normal_approx"))
                           )
            )
          }),
          fluidRow(
            card(
              card_header("Explicación",class="custom-header"),
              withMathJax(htmltools::includeMarkdown("help/Samp_dist.Rmd")),full_screen = TRUE#, height=200
            )
          )
          
            
  ),
  #### IC and error ####
  nav_panel("6. IC y error muestreo",
            fluidRow({
              layout_columns(col_widths=c(5,7),height = 900,
                             card(
                               card_header("Parámetro de interés"),
                               card(height=600,layout_columns(
                                 tableOutput('tabla_interes5'),
                                 selectInput("paramint","Parámetro de interés",
                                             choices=c("N","G","V","h_media","dg","ho"))
                               )
                               ),
                             
                               card(
                                 card_header("Cambio en la desviación típìca aumentar n"),
                                 plotOutput("var_n2")
                               )),
                             card(card_header("Intervalos de confianza"),
                                  plotOutput("intervals"),height=1000
                             )
              )
            }),
            fluidRow(
              card(
                card_header("Explicación",class="custom-header"),
                withMathJax(htmltools::includeMarkdown("help/Interval_error.Rmd")),full_screen = TRUE#, height=200
              )
            )
            
  ),
  #### Sample alloc ####
  nav_panel("7. Cálculo de n",
            fluidRow({
              layout_columns(col_widths=c(5,7),height = 800,
                             card(card_header("Parámetros de interés y muestra"),
                                  card(layout_columns(col_widths=c(5,7),
                                                      tableOutput('tabla_interes6'),
                                                      plotOutput("plot_selected6",width=400,height=400)
                                  )),
                                  card(tableOutput('n_estimaciones3'),min_height = 450),
                                  
                             ),
                             card(card_header("Estimación con una parcela vs estimación con n parcelas"),
                                  card(plotOutput("plot_res4",width=950,height=750))
                             )
              )
            }),
            fluidRow(
              card(
                card_header("Explicación",class="custom-header"),
                withMathJax(htmltools::includeMarkdown("help/Samp_aloc.Rmd")),full_screen = TRUE#, height=200
              )
            )
            
  )
  # )
)
