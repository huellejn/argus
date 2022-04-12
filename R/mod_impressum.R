#' impressum UI Function
#'
#' @description Shiny module containing the 'Impressum' section as of the UI.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_impressum_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    h4("Contact"),
    p("Heiko Brennenstuhl, MD, MBA"),
    p("Center for Pediatrics and Adolescent Medicine"),
    p("University Hospital Heidelberg"),
    p("Im Neuenheimer Feld 430"),
    p("D-69120 Heidelberg, Germany"),
    br(),
    p(strong("For individual queries on aRgus, please refer to:")),
    p(a("Heiko.Brennenstuhl@med.uni-heidelberg.de", href="mailto:Heiko.Brennenstuhl@med.uni-heidelberg.de")),
    br(),
    h4("Physical server address"),
    p("IT resources have been kindly provided by heiCLOUD, a service of University Computing Center Heidelberg."),
    br(),
    img(src = "www/URZ_logo.png", height = "100px"),
    h4("heiCLOUD"),
    p(strong("University Computing Centre Heidelberg (URZ)")),
    p("University of Heidelberg"),
    p("Im Neuenheimer Feld 330"),
    p("D-69120 Heidelberg, Germany"),
    p(a("http://www.urz.uni-heidelberg.de", href="http://www.urz.uni-heidelberg.de"))
    
  )
}
    
## To be copied in the UI
# mod_impressum_ui("impressum_ui_1")
    