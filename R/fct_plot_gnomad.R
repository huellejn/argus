#' plot_gnomad 
#'
#' @description Function to create a bar plot of gnomAD variants.
#'
#' @return Returns a bar plot of gnomAD variants by amino acid position in the protein.
#'
#' @noRd

plot_gnomad <- function(dat, protein_length, dat_segments) {
  
  # Suppress 'No visible binding for global variable' message
  aa_label <- aa_pos <- aapos <- gnomAD_exomes_AC <- gnomAD_genomes_AC <- Source <- Count <- Sum <- NULL
  
  df_gnomad <- dat %>%
    dplyr::select(aapos, gnomAD_exomes_AC, gnomAD_genomes_AC) %>%
    gather(-aapos, key = "Source", value = "Count") %>%
    group_by(aapos, Source) %>%
    summarise(Sum = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(Sum > 0)
  
  # Common layers
  cm_layers <- list(
    scale_x_continuous(expand = c(0,0), limits = c(1,protein_length), breaks = x_axis_breaks(protein_length)),
    annotation_logticks(base = 10, sides = "l"),
    theme_minimal(),
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(.1, .1, 0, .05, unit = "in")
    )
  )
  
  # Plot gnomAD exomes
  p1 <- ggplot(data = dplyr::filter(df_gnomad, Source == "gnomAD_exomes_AC"), aes(x = aapos, y = Sum)) +
    geom_col(color = "forestgreen", fill = "forestgreen") +
    scale_y_log10() +
    labs(y = "gnomAD exomes count") + 
    cm_layers +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(.1, .1, 0, .05, unit = "in"),
      plot.title = element_text(hjust = .5)
    )
  
  # Plot gnomAD genomes
  p2 <- ggplot(data = dplyr::filter(df_gnomad, Source == "gnomAD_genomes_AC"), aes(x = aapos, y = Sum)) +
    geom_col(color = "deepskyblue", fill = "deepskyblue") +
    scale_y_continuous(trans = reverselog_trans(10), expand = c(0,0)) +
    labs(y = "gnomAD genomes count", x = "Amino acid position") + 
    cm_layers
  
  # Add selected variants
  if(nrow(dat_segments) > 0) {
    
    p1 <- p1 + 
      geom_segment(data = dat_segments, aes(x = aa_pos, xend = aa_pos, y = 0, yend = 1), color = dat_segments$color) +
      geom_text_repel(data = dat_segments, 
                      aes(x = aa_pos, y = 1, 
                          label = aa_label),
                      nudge_y = log10(max(df_gnomad$Sum, na.rm=T)),
                      direction = "x",
                      angle = 90,
                      segment.size = .4,
                      size = 3,
                      segment.linetype = "dotted",
                      max.overlaps = Inf)  
    
    p2 <- p2 + 
      geom_text_repel(data = dat_segments, 
                      aes(x = aa_pos, y = 1,
                          label = aa_label),
                      nudge_y = -log10(max(df_gnomad$Sum, na.rm=T)),
                      direction = "x",
                      angle = 90,
                      segment.size = .4,
                      size = 3,
                      segment.linetype = "dotted",
                      max.overlaps = Inf)  
    
  }
  
  # Combine gnomAD exomes and genomes plots
  p <- cowplot::plot_grid(p1, p2, nrow = 2)
  
  return(p)
  
}
