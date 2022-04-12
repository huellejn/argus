#' plot_protein 
#'
#' @description Function to create a graphical depiction of protein domains.
#'
#' @return Returns a plot of protein domains.
#'
#' @noRd

plot_protein_domain <- function(protein_id, protein_info_domain, protein_length, color_count, clinvar_goi, toi_nm_short) {
  
  # Suppress 'No visible binding for global variable' message
  type <- begin <- end <- description <- RefSeq.transcript.short <- Label <- AA.position <- NULL
  
  p <- ggplot() + 
    geom_rect(aes(xmin = 1, xmax = protein_length, ymin = -1, ymax = 1), fill = "grey",  size = .25) +
    ylim(-2.5, 10) +
    scale_x_continuous(expand = c(0,0), breaks = x_axis_breaks(protein_length))+
    theme_minimal() +
    labs(title = protein_id, x = "Amino acid position") + 
    theme(axis.text.y = element_blank(),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          plot.margin = margin(.1, .15, .1, .58, unit = "in")
    )
  
  if("CHAIN" %in% protein_info_domain$type) {
    p <- p +
      geom_rect(data = protein_info_domain %>%
                  dplyr::filter(type== "CHAIN") %>% 
                  arrange(desc(begin)), #why arrange desc?
                aes(xmin = begin, xmax = end, ymin = -1, ymax = 1), fill = "grey40", size = .25) + 
      scale_fill_nejm(name = "Chain peptide")
  } 
  
  if("SIGNAL" %in% protein_info_domain$type ){
    p <- p +
      geom_rect(data = protein_info_domain %>%
                  dplyr::filter(type== "SIGNAL") %>%
                  arrange(desc(begin)), aes(xmin = begin, xmax = end, ymin = -1, ymax = 1, fill = description), size = .25) + 
      scale_fill_nejm(name = "Signal peptide")
  } 
  
  if("REGION" %in% protein_info_domain$type ) {
    p <- p +
      new_scale_fill() +
      geom_rect(data=protein_info_domain %>%
                  dplyr::filter(type == "REGION") %>%
                  arrange(desc(begin)), aes(xmin=begin, xmax=end, ymin=-2, ymax=2, fill=description), size=.25) +
      scale_fill_manual(values = getPalette(color_count), name="Region")
  } 
  
  if(nrow(protein_info_domain[! protein_info_domain$type %in% c("CHAIN", "REGION", "SIGNAL"),]) > 0 ) {
    p <- p + 
      new_scale_fill() +
      geom_rect(data=protein_info_domain %>%
                  dplyr::filter(! type %in% c("CHAIN", "REGION", "SIGNAL")) %>%
                  arrange(desc(begin)), aes(xmin=begin, xmax=end, ymin=-.5, ymax=.5, fill=description),size=.25) + 
      scale_fill_manual(values = getPalette2(color_count), name="Domain")
  } 
  
  # Add segments
  if(nrow(clinvar_goi)>0) {
    
    clinvar_canonical_pathogenic <- dplyr::filter(clinvar_goi, 
                                                  RefSeq.transcript.short == toi_nm_short, 
                                                  Label == "pathogenic",
                                                  !is.na(AA.position))
    
    if(nrow(clinvar_canonical_pathogenic)>0) {
      
      p <- p +
        geom_segment(data = clinvar_canonical_pathogenic,
                     aes(x = as.numeric(AA.position), xend = as.numeric(AA.position),
                         y = 2.1, yend = 3),
                     color = "black",
                     size = .4) +
        geom_point(data = clinvar_canonical_pathogenic,
                   aes(x = as.numeric(AA.position), y = 3),
                   shape = 21, fill = "red", alpha = .3)
    }
  }
  
  return(p)
  
}
