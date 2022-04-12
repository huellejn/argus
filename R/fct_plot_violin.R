#' plot_violin 
#'
#' @description Function to create a violin plot showing the distribution of functional prediction scores.
#'
#' @return Returns a plot of the distribution of functional scores and variant frequencies in gnomAD.
#'
#' @noRd

plot_violin <- function(clinvar_goi, transcript_info, tidy_data, selected_scores) {
  
  # Suppress 'No visible binding for global variable' message
  Label <- Chromosome.position <- REF <- ALT <- ClinVar <- score_type <- score <- gnomAD_exomes_AC <- gnomAD_genomes_AC <- gnomAD <- InSilicoScore <- value <- p.adj.signif <- group2 <- NULL
  
  if(length(selected_scores) > 0 ) {
    
    # Settings
    ## Match bases for genes on reverse strand
    reverse_strand <- c(A="T", T="A", C="G", G="C")
    
    violin_color <- setNames(c("#CC6677", "#88CCEE", "#999933","#DDCC77" ), 
                             c("ClinVar_pathogenic", "ClinVar_benign" , "gnomAD", "InSilico"))
    
    # Add ClinVar annotation to in the silico score tidy_data table
    ## Format ClinVar variants
    clinvar_sub_score <- clinvar_goi %>%
      mutate(ClinVar = paste("ClinVar", Label, sep = "_")) %>%
      dplyr::select(pos = Chromosome.position, ref = REF, alt = ALT, ClinVar)
    
    ## Reverse bases if gene is encoded on the reverse strand
    if(transcript_info$seq_strand == -1) {
      clinvar_sub_score$ref <- names(reverse_strand)[match(clinvar_sub_score$ref, reverse_strand)]
      clinvar_sub_score$alt <- names(reverse_strand)[match(clinvar_sub_score$alt, reverse_strand)]
    }
    
    # Anotate ClinVar variants in the InSilico table
    dat <- tidy_data %>%
      dplyr::filter(score_type %in% selected_scores) 
    
    if(nrow(dat)>0) {
      
      dat <- dat %>%
        left_join( clinvar_sub_score, by = c("pos", "ref", "alt") ) %>%
        mutate(
          InSilicoScore = "InSilico",
          gnomAD = ifelse(!is.na(gnomAD_exomes_AC), "gnomAD", ifelse(!is.na(gnomAD_genomes_AC), "gnomAD", NA))
        ) %>%
        dplyr::select(score_type, score, gnomAD, InSilicoScore, ClinVar) %>%
        gather(-score, -score_type, key = "group", value = "value") %>%
        mutate(value = factor(value, levels = c("ClinVar_pathogenic", "ClinVar_benign", "gnomAD", "InSilico"))) %>%
        dplyr::filter(!is.na(value), ! is.na(score)) 
      
      violin_y_lim <- dat %>%
        group_by(score_type) %>%
        summarise(
          max = max(score, na.rm=T),
          max_plus = max + max * 0.1) %>%
        arrange(score_type) 
      
      dat_stat <- dat %>%
        group_by(score_type) %>%
        t_test(score ~ value) %>%
        adjust_pvalue() %>%
        add_significance("p.adj") %>% 
        add_xy_position(x = "value", scales = "free_y", step.increase = .3) %>%
        mutate(p.adj.signif = gsub("\\*\\*\\*\\*", "\\*\\*\\*", p.adj.signif)) %>%
        filter(group2 == "InSilico")
      
      p <- ggplot(dat, aes(x = value, y = score)) +
        geom_violin(aes(fill = value), draw_quantiles = c(.25, .5, .75), alpha = .5 ) +
        ggpubr::stat_pvalue_manual(dat_stat, label = "p.adj.signif", tip.length = 0, hide.ns = T, x = "group1") +
        scale_fill_manual(values = violin_color, name = "Group") +
        facet_wrap(~ score_type, scales = "free_y") +
        theme_bw() +
        theme(
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        labs(y = "Score value") +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) 
      
    } else {
      
      p <- ggplot() + 
        theme_void() + 
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        labs(y = "Score value")
      
    }
    
    return(p)
    
  }
  
}
