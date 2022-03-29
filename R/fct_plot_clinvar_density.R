#' plot_clinvar_density 
#'
#' @description Functions to create a density plot of ClinVar variants.
#'
#' @return Returns a density plot of pathogenic ClinVar variants or an empty plot if no pathogenic entries are available in ClinVar.
#'
#' @noRd

plot_density_clinvar_empty <- function(protein_length) {
  
  p <- ggplot() +
    theme_void() +
    scale_x_continuous(limits = c(1, protein_length), breaks = x_axis_breaks(protein_length))+
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Amino acid position", y = "ClinVar variant density") + 
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(.1, .1, .1, .4, unit = "in"),
      legend.position = "bottom"
    )
  
  return(p)
  
}

plot_density_clinvar <-  function(dat, protein_length, cols_selected) {
  
  p <- ggplot(data = dat) + 
    geom_density(aes(x = AA.position, y = after_stat(ndensity), fill = Label), adjust = 0.1, alpha = .7, linetype = "blank") +  
    geom_segment(aes(x = AA.position, xend = AA.position, y = -.1, yend = 0, color = Label)) +
    geom_hline(yintercept = 0, color = "grey60") +
    scale_x_continuous(expand = c(0,0), limits = c(1,protein_length), breaks = x_axis_breaks(protein_length))+
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual("Clinical significance", values = cols_selected) + 
    scale_color_manual("Clinical significance", values = cols_selected) + 
    labs(x = "Amino acid position", y = "ClinVar variant density") + 
    theme_classic() + 
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(.1, .1, .1, .4, unit = "in")
    )
  
  return(p)
  
} 
