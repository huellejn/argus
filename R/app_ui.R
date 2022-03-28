#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @importFrom shinythemes shinytheme
#' @import shinyjs
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    
    # UI 
    navbarPage(
      
      title = a("aRgus: multilevel variant visualization & advanced prediction score modeling", href = "https://argus.urz.uni-heidelberg.de/"),
      windowTitle = "aRgus",
      
      theme = shinytheme("flatly"),
      
      tabPanel("Variants", style = "overflow:hidden;",
               fluidRow(
                 
                 sidebarLayout(
                   
                   sidebarPanel(width=2,
                                style = "height:500px;overflow-y: scroll;",
                                
                                # Select gene
                                selectizeInput(
                                  inputId = 'selectgene', label = 'Select gene',
                                  choices = NULL, selected = 1
                                ),
                                
                                # Select plots
                                checkboxGroupInput(
                                  inputId = 'selectplot', label = 'Select plots',
                                  choices = c('Transcript', 'Protein', 'ClinVar', 'gnomAD', 'In silico scores', 'Score statistics'),
                                  selected = c('Transcript', 'Protein', 'ClinVar', 'gnomAD', 'In silico scores', 'Score statistics')
                                  
                                ),
                                
                                # Select Clinvar significance
                                selectInput("selectclinsig", "ClinVar significance", choices = c("pathogenic", "benign"), selected = c("pathogenic", "benign"), multiple = TRUE),
                                
                                # Select scores
                                fluidRow(
                                  column(12, offset = 0, style = 'padding-right:4px', 
                                         shinyjs::disabled(selectizeInput("selectscore", label = "Select scores",
                                                                          choices = NULL,
                                                                          selected = 1,
                                                                          multiple = TRUE, 
                                                                          options = list(
                                                                            maxItems = 3
                                                                          ))
                                         )
                                  )
                                ),
                                
                                tags$hr(style="border-color: #95a5a6;"),
                                
                                h5("Figure download options"),
                                radioButtons(
                                  inputId = 'figureformat', label = 'Download format',
                                  choices = c('PNG', 'SVG'),
                                  selected = 'PNG'
                                ),
                                fluidRow(
                                  column(6, numericInput("figurewidth", "Width", value = 8, min = 0)),
                                  column(6, numericInput("figureheight", "Height", value = 4, min = 0))
                                ),
                                
                                # Download report
                                #br(),
                                
                                #radioButtons(
                                #   inputId = 'reportformat', label = 'Document format',
                                #   choices = c('HTML', 'PDF'),
                                #   selected = 'HTML'
                                # ),
                                
                                #shinyjs::disabled(downloadButton(
                                #   outputId = 'dwnbtn_report', label = 'Generate report'
                                # ))
                                
                   ),
                   
                   mainPanel(width=10,
                             style = "height:500px;overflow-x: scroll;overflow-y: scroll;",
                             
                             conditionalPanel(
                               condition = "input.selectplot.includes('Transcript') && input.selectgene.length > 0", # && $('html').attr('class')!=='shiny-busy'
                               fluidRow(
                                 column(11, shinycssloaders::withSpinner(plotOutput('plot_transcript'))),
                                 column(1, downloadButton('dwnbtn_transcript', '', icon = icon('download') ) )
                               )),                                            
                             conditionalPanel(
                               condition = "input.selectplot.includes('Protein') && input.selectgene.length > 0",
                               fluidRow(
                                 column(11, shinycssloaders::withSpinner(plotOutput('plot_protein')) ),
                                 column(1, downloadButton('dwnbtn_protein', '', icon = icon('download') ) )
                               )),
                             conditionalPanel(
                               condition = "input.selectplot.includes('ClinVar') && input.selectgene.length > 0",
                               fluidRow(
                                 column(11, shinycssloaders::withSpinner(plotOutput('plot_density', height = '250px'))),
                                 column(1, downloadButton('dwnbtn_density', '', icon = icon('download') ))
                               )),
                             conditionalPanel(
                               condition = "input.selectplot.includes('gnomAD') && input.selectgene.length > 0",
                               fluidRow(
                                 column(11, shinycssloaders::withSpinner(plotOutput('plot_gnomad'))),
                                 column(1, downloadButton('dwnbtn_gnomad', '', icon = icon('download') ))
                               )),
                             conditionalPanel(
                               condition = "input.selectplot.includes('In silico scores') && input.selectgene.length > 0",
                               fluidRow(
                                 column(11, shinycssloaders::withSpinner(plotOutput('plot_scores1', height = '250px'))),
                                 column(1, downloadButton('dwnbtn_scores1', '', icon = icon('download') ))
                               )),
                             conditionalPanel(
                               condition = "input.selectplot.includes('In silico scores') && input.selectgene.length > 0 && input.selectscore.length > 1",
                               fluidRow(
                                 column(11, shinycssloaders::withSpinner(plotOutput('plot_scores2', height = '250px'))),
                                 column(1, downloadButton('dwnbtn_scores2', '', icon = icon('download') ))
                               )),
                             conditionalPanel(
                               condition = "input.selectplot.includes('In silico scores') && input.selectgene.length > 0 && input.selectscore.length > 2",
                               fluidRow(
                                 column(11, shinycssloaders::withSpinner(plotOutput('plot_scores3', height = '250px'))),
                                 column(1, downloadButton('dwnbtn_scores3', '', icon = icon('download') ))
                               )),
                             conditionalPanel(
                               condition = "input.selectplot.includes('Score statistics') && input.selectgene.length > 0",
                               fluidRow(
                                 column(11, shinycssloaders::withSpinner(plotOutput('plot_violin'))),
                                 column(1, downloadButton('dwnbtn_violin', '', icon = icon('download') ))
                               ))
                   )
                 )
               ),
               
               fluidRow(
                 column(12, 
                        conditionalPanel(
                          condition = 'input.selectgene.length > 0',
                          br(),
                          tabsetPanel(
                            tabPanel('ClinVar', DT::dataTableOutput('table1')),
                            tabPanel('In silico scores', DT::dataTableOutput('table2'))
                          )
                        )
                 )
               )
               
               ),
      
      navbarMenu("More",
                 tabPanel("About the app", mod_about_ui("about_ui") ),
                 tabPanel("Help", mod_help_ui("help_ui") ),
                 tabPanel("Impressum", mod_impressum_ui("impressum_ui") )
                 )
    )
    
  )
}
#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'argus'
    )
  )
}
