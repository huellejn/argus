#' utils 
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
#' @import RColorBrewer

# Common variables
cols <- c("pathogenic" = "red", "benign" = "blue")
legend_order <- c("pathogenic", "benign")
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))
getPalette2 = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))

# Path to data
path_to_data <- "/Users/huellein/Documents/project/aRgus/app/data"
path_to_db <- "/Users/huellein/Documents/project/aRgus/app/data"

# Read data
## Canonical transcript
load(file.path(path_to_data, "canonical_transcript_dev.RData"))

## ClinVar
load(file.path(path_to_data, "clinvar_jan21_dev.RData"))
clinvar_version <- "January 2021"

## Gene info
load(file.path(path_to_data, "gene_info_dev.RData"))

## dbNSFP scores
load(file.path(path_to_data, "dbNSFP_scores_dev.RData"))

## dbNSFP
file_dbNSFP <- "dbNSFP4.1a_grch38.gz"
dbNSFP_version <- "dbNSFP4.1a_grch38"

## Input values for gene of interest (goi)
val_goi <- canonical_transcript$Gene.symbol

# Define functions ----------------------------
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# x-axis breaks
x_axis_breaks <- function(protein_length) {
  ## The x-axis breaks should start with 1, increase in steps of 50, and end with protein length.
  steps = ifelse(protein_length >= 400, 100, 50)
  x_axis_breaks_cnt <- steps # start of the number series
  x_axis_breaks <- integer() # initiate variable to store the output
  while(x_axis_breaks_cnt < protein_length) { 
    x_axis_breaks = c(x_axis_breaks, x_axis_breaks_cnt)
    x_axis_breaks_cnt = x_axis_breaks_cnt + steps
  }
  # Add protein length. Remove the upper value if it is too close to protein length.
  if( protein_length - x_axis_breaks[[length(x_axis_breaks)]] > 40){
    x_axis_breaks <- c(1, x_axis_breaks, protein_length)
  } else {
    x_axis_breaks <- c(1, x_axis_breaks[1:length(x_axis_breaks)-1], protein_length)
  }
  return(x_axis_breaks)
}
