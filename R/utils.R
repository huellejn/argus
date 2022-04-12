#' utils 
#'
#' @description A collection of utils function
#'
#' @return Returns input data and parameter settings.
#'
#' @noRd
#' @import RColorBrewer

# Common variables
cols <- c("pathogenic" = "red", "benign" = "blue")
legend_order <- c("pathogenic", "benign")
getPalette = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))
getPalette2 = colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))

# Path to data
path_to_data <- "/mnt/data"
path_to_db <- "/mnt/data"

# Read data
## Canonical transcript
load(file.path(path_to_data, "canonical_transcripts.RData"))

## ClinVar
load(file.path(path_to_data, "clinvar.RData"))
clinvar <- clinvar[!clinvar$ClinVar.ID == "CV:15624", ]
#clinvar_version <- gsub(".*clinvar_[].RData", "", clinvar_file)

## Gene info
load(file.path(path_to_data, "gene_info.RData"))

# Protein info
baseurl <- "https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession="

## dbNSFP scores
load(file.path(path_to_data, "gene_info.RData"))

## dbNSFP
dbNSFP_version <- "dbNSFPv4.3a_GRCh38"
file_dbNSFP <- paste(dbNSFP_version, "gz", sep = ".")

## Input values for gene of interest (goi)
val_goi <- canonical_transcript$Gene.symbol

# Settings for violin score plot
## Match bases for genes on reverse strand
reverse_strand <- c(A="T", T="A", C="G", G="C")

violin_color <- setNames(c("#CC6677", "#88CCEE", "#999933","#DDCC77" ), 
                         c("ClinVar_pathogenic", "ClinVar_benign" , "gnomAD", "InSilico"))

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
