#' plot_transcript 
#'
#' @description Function to create a graphical depiction of a transcript. The x-axis shows the genomic position according to GRCh38.
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd

plot_transcript <- function(transcript_info, exon_info, toi, toi_nm, toi_nm_short, clinvar_goi, font_size) {
  
  # Suppress 'No visible binding for global variable' message
  tx_seq_start <- tx_seq_end <- exon_seq_start <- exon_seq_end <- RefSeq.transcript.short <- Label <- Chromosome.position <- NULL
  
  # Base plot
  p <- ggplot() + 
    # Transcript
    geom_rect(data = transcript_info, aes(xmin = tx_seq_start, xmax = tx_seq_end, ymin = -1, ymax = 1), fill = "grey") +
    geom_rect(data = exon_info, aes(xmin = exon_seq_start, xmax = exon_seq_end, ymin = -2, ymax = 2),
              fill="aquamarine4", alpha = .8) +
    ylim(-3,10) +
    theme_minimal() +
    labs(title = ifelse(!is.na( toi_nm ), paste( toi, toi_nm, sep = " / "), toi), x = "Genome position" ) + 
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          axis.title.x = element_text(size = font_size),
          axis.text.x = element_text(size = font_size-4),
          plot.title = element_text(size = font_size)
    )
  
  # Add segments for pathogenic variants
  if(nrow(clinvar_goi) > 0) {
    
    clinvar_canonical_pathogenic <- dplyr::filter(clinvar_goi, RefSeq.transcript.short == toi_nm_short, Label == "pathogenic")
    
    if(nrow(clinvar_canonical_pathogenic)>0) {
      
      p <- p +
        geom_segment(data = clinvar_canonical_pathogenic,
                     aes(x = Chromosome.position, xend=Chromosome.position,
                         y = 2.1, yend = 3),
                     color = "black",
                     size = .4) +
        geom_point(data = clinvar_canonical_pathogenic,
                   aes(x = Chromosome.position, y = 3),
                   shape = 21, fill = "red", alpha = .3)
    }
  }
  
  return(p)
  
}
