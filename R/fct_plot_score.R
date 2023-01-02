#' plot_score 
#'
#' @description Function to create a density plot and a heatmap of functional prediction scores.
#'
#' @return Returns a plot of in-silico and conservation prediction scores.
#'
#' @noRd

plot_score <- function(dat, selected_scores, score_index, protein_length, goi, toi, dbNSFP_scores, dat_segments, font_size, font_size_label) {
  
  # Suppress 'No visible binding for global variable' message
  aa_label <- aa_pos <- score_type <- aapos <- score <- NULL
  
  if(length(selected_scores) >= score_index) {
    
    # Get selected score
    selected_score <- selected_scores[score_index]
    
    # Get data
    dat <- dplyr::filter(dat, score_type %in% {{selected_score}} ) 
    
    if(nrow(dat)>0) {
      
      dat_hp <- dat %>%
        arrange(aapos) %>%
        group_by(aapos) %>%
        summarise(mean = mean(score, na.rm = TRUE))
      
      # Plot
      p1 <- ggplot(dat) +
        geom_smooth(aes(x = aapos, y = score), color = "red", na.rm = TRUE) + 
        theme_minimal() + 
        scale_x_continuous(expand = c(0,0), limits = c(1, protein_length)) +
        labs(title = paste0({{selected_score}}, " scores for gene ", goi, " (", toi, ")") ) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = margin(.2, .1, 0, .65, unit = "in"),
              plot.title = element_text(hjust = .5, size = font_size)
        )
      
      # Meta data for the scores
      dbNSFP_score <- dbNSFP_scores[dbNSFP_scores$score_type == selected_score, ]
      # colors
      col_hp_low <- dbNSFP_score$color_min 
      col_hp_high <- dbNSFP_score$color_max
      
      if(is.na(dbNSFP_score$score_min)) {
        lim_hp <- c(
          min(dat_hp$mean),
          max(dat_hp$mean)
        )
      } else {
        lim_hp <- sort(c(
          dbNSFP_score$score_min,
          dbNSFP_score$score_max
        ))
      }
      
      p2 <- ggplot(dat_hp, aes(x = aapos, y = 1, fill = mean)) + 
        geom_tile() + 
        labs(x = "Amino acid position") + 
        theme_minimal() + 
        scale_fill_gradient(low = col_hp_low, high = col_hp_high, name={{selected_score}}, limits = c(lim_hp[1], lim_hp[2]) ) +
        scale_x_continuous(expand = c(0,0), limits = c(1, protein_length), breaks = x_axis_breaks(protein_length)) +
        theme(legend.position = "bottom",
              panel.grid = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              plot.margin = margin(0, .2, .1, .65, unit = "in"),
              axis.title.x = element_text(size = font_size),
              axis.text.x = element_text(size = font_size-4),
              legend.title = element_text(size = font_size),
              legend.text = element_text(size = font_size-8)
        )
      
      # Add selected variants
      if(nrow(dat_segments) > 0) {
        
        p1 <- p1 + 
          geom_segment(data = dat_segments, aes(x = aa_pos, xend = aa_pos, y = -.1, yend = max(dat$score, na.rm=T)/8), color = dat_segments$color) +
          geom_text_repel(data = dat_segments, 
                          aes(x = aa_pos, y = max(dat$score, na.rm=T)/8, 
                              label = aa_label),
                          nudge_y = max(dat$score, na.rm=T),
                          direction = "x",
                          angle = 90,
                          segment.size = .4,
                          size = font_size_label,
                          segment.linetype = "dotted",
                          max.overlaps = Inf)  
        
      }
      
      p <- plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(1,1))
    } else {
      
      # Generate empty plot
      p <- ggplot() +
        theme_void() +
        theme_minimal() +
        scale_x_continuous(expand = c(0,0), limits = c(1, protein_length )) +
        scale_y_continuous(limits = c(0,1)) +
        labs(title = paste0({{selected_score}}, " scores for gene ", goi, " (", toi, ")") ) +
        theme(
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(.1, .1, .1, .65, unit = "in"),
          axis.title.x = element_text(size = font_size),
          axis.text.x = element_text(size = font_size-4)
        )
      } 
    return(p)
  }
  
}
