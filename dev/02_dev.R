# Add a package
usethis::use_package( "jsonlite" )

# Add a module
golem::add_module( name = "variants" )

# Add a function
golem::add_fct( "plot_protein" ) 

golem::add_utils( "utils" )

golem::add_js_file( "script" )
golem::add_js_handler( "handlers" )
golem::add_css_file( "custom" )
usethis::use_data_raw( name = "my_dataset", open = FALSE ) 
usethis::use_test( "app" )
usethis::use_vignette("argus")
devtools::build_vignettes()
usethis::use_coverage()
covrpage::covrpage()
usethis::use_github()
usethis::use_github_action() 
usethis::use_github_action_check_release() 
usethis::use_github_action_check_standard() 
usethis::use_github_action_check_full() 
usethis::use_github_action_pr_commands()
usethis::use_travis() 
usethis::use_travis_badge() 
usethis::use_appveyor() 
usethis::use_appveyor_badge()
usethis::use_circleci()
usethis::use_circleci_badge()
usethis::use_jenkins()
usethis::use_gitlab_ci()
rstudioapi::navigateToFile("dev/03_deploy.R")
