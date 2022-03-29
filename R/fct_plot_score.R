#' plot_score 
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#'
#' @noRd

plot_score <- function(dat, selected_scores, score_index, protein_length, goi, toi, dbNSFP_scores) {
  
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
              plot.margin = margin(.1, .1, 0, .65, unit = "in"),
              plot.title = element_text(hjust = .5)
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
              plot.margin = margin(0, .1, .1, .65, unit = "in")
        )
      
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
          plot.margin = margin(.1, .1, .1, .65, unit = "in")
        )
      } 
    return(p)
  }
  
}
