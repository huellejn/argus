#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import dplyr
#' @import ggplot2
#' @import drawProteins
#' @import tidyr
#' @import scales
#' @import ggsci
#' @import ggnewscale
#' @import ggrepel
#' @import cowplot
#' @import ggpubr
#' @import grid
#' @import gridExtra
#' @import grDevices
#' @import rstatix
#' @import DT
#' @import shinyjs
#' @import seqminer
#' @import splitstackshape
#' @noRd
app_server <- function( input, output, session ) {
  
  # Input field for selection of gene symbol
  updateSelectizeInput(session, inputId = 'selectgene', 
                       choices = c(val_goi), server = TRUE,
                       options = 
                         list(
                           placeholder = 'HGNC symbol', 
                           maxItems = 5,
                           openOnFocus = FALSE,
                           closeAfterSelect = TRUE,
                           selectOnTab = TRUE
                         ))
  # Update gene input when search term is added to the url
  observe({
    
    query <- parseQueryString(session$clientData$url_search)
    
    if (!is.null(query[['gene']])) {
      updateSelectizeInput(session, "selectgene", selected = query[['gene']])
    }
    
  })
  
  
  # Section 1: Data --------------------------------------------------------
  
  ## Get selected gene and transcript
  goi <- reactive({
    input$selectgene
  })
  
  toi <- reactive({
    canonical_transcript$Ensembl.transcript[canonical_transcript$Gene.symbol == goi()]
  })
  
  toi_nm <- reactive({
    canonical_transcript$RefSeq.transcript[canonical_transcript$Gene.symbol == goi()]
  })
  
  toi_nm_short <- reactive({
    gsub("\\..*", "", toi_nm())
  })

  ## Get transcript, exon and protein information from ensembl and drawProteins
  transcript_info <- reactive({
    dplyr::filter(edb_transcript_info, tx_id == toi())
  })
  
  exon_info <- reactive({
    dplyr::filter(edb_exon_info, tx_id == toi())
  }) 
  
  protein_info <- reactive({
    dplyr::filter(edb_protein_info, tx_id == toi()) %>%
      dplyr::select(-tx_id) %>%
      distinct()
  }) 
  
  protein_info_domain <- reactive({
    
    req(nrow(protein_info()) > 0)
    
    protein_info_domain_json <- drawProteins::get_features(protein_info()$uniprot_id)
    protein_info_domain <- drawProteins::feature_to_dataframe(protein_info_domain_json)
    protein_info_domain <- protein_info_domain %>%
      dplyr::filter(
        type %in% c("SIGNAL", "CHAIN", "DOMAIN", "ZN_FING", "COILED", "TRANSIT", "NP_BIND", "ACT_SITE", "VAR_SEQ", "REPEAT", "REGION"),
        ! description %in% c("NONE", "Disordered")
      )
    return(protein_info_domain)
  })
  
  protein_length <- reactive({
    
    req(nrow(protein_info()) > 0)
    
    max(protein_info_domain()$end, na.rm = TRUE)
  }) 
  
  ## Get clinvar entries for goi
  clinvar_goi <- reactive({
    clinvar %>%
      dplyr::filter(Gene.symbol == goi() ) %>%
      arrange(Chromosome, Chromosome.position)
  })
  
  ## Get dbNSFP entries for goi
  ## sc: selected columns
  sc <- c(1:4,
          12:16,
          20,
          23:25,
          28,
          79, 116, 37, 40, 43, 46, 64, 76, 67, 61, 69, 72, 102, 53, 58, 120, 81, 86, 88, 49, 90, 93, 96, 99, 122, 126, 129, 132, 135, 137, 140, 143, 146, 
          149, 152, 154, 156, 158, 160, 162, 164, 167, 105, 242, 244, 315, 317
          # REVEL, CADD_phred, SIFT_score, SIFT4G_score, Polyphen2_HDIV_score, Polyphen2_HVAR_score, #PROVEAN_score, M-CAP_score, VEST4_score, 
          # FATHMM_score, MetaSVM_score, MetaLR_score, ClinPred_score, MutationTaster_score, MutationAssessor_score, DANN_score, MutPred_score,
          # MVP_score, MPC_score, LRT_score, PrimateAI_score, DEOGEN2_score, BayesDel_addAF_score, BayesDel_noAF_score, fathmm-MKL_coding_score, 
          # fathmm-XF_coding_score, Eigen-raw_coding, Eigen-PC-raw_coding, GenoCanyon_score, integrated_fitCons_score, GM12878_fitCons_score,
          # H1-hESC_fitCons_score, HUVEC_fitCons_score
          # LINSIGHT, GERP++_RS, phyloP100way_vertebrate, phyloP30way_mammalian, phyloP17way_primate, phastCons100way_vertebrate, phastCons30way_mammalian
          # phastCons17way_primate, SiPhy_29way_logOdds, LIST-S2_score, gnomAD_exomes_AC, gnomAD_exomes_AF, gnomAD_genomes_AC, gnomAD_genomes_AF
  )
  
  genotypes <- reactive({
    
    req(goi(), nrow(transcript_info())>0)
    
    genotypes <- tabix.read.table(file.path(path_to_db, file_dbNSFP), transcript_info()$position, col.names = TRUE, stringsAsFactors = FALSE)
    genotypes <- genotypes[grep(toi(), genotypes$Ensembl_transcriptid),sc]
    
    validate(
      need(nrow(genotypes) > 0, "No transcript information available")
    )
    
    genotypes <- cSplit(genotypes, c("aapos", 
                                     "genename", 
                                     "Ensembl_geneid", 
                                     "Ensembl_transcriptid", 
                                     "Ensembl_proteinid", 
                                     "HGVSp_ANNOVAR",
                                     "HGVSc_VEP", 
                                     "HGVSp_VEP", 
                                     "APPRIS", 
                                     "VEP_canonical", 
                                     "SIFT_score", 
                                     "SIFT4G_score", 
                                     "Polyphen2_HDIV_score", 
                                     "Polyphen2_HVAR_score", 
                                     "MutationTaster_score",
                                     "MutationAssessor_score",
                                     "FATHMM_score",
                                     "PROVEAN_score", 
                                     "VEST4_score", 
                                     "MVP_score", 
                                     "MPC_score",
                                     "DEOGEN2_score",
                                     "LIST.S2_score"), ";", "long", type.convert = FALSE)
    
    genotypes[genotypes == "."] <- NA
    
    genotypes <- genotypes %>%
      dplyr::rename(HGVSc = HGVSc_VEP, HGVSp = HGVSp_VEP, pos = pos.1.based.) %>%
      dplyr::filter(
        Ensembl_transcriptid == toi(),
        aapos > 0
      ) %>%
      rename_with(
        ~ gsub('_score', '', .x)
      ) %>%
      rename("GERP++_RS" = "GERP.._RS") %>%
      rename_with(
        ~ gsub('\\.', '-', .x),
      ) %>%
      dplyr::distinct(HGVSp, .keep_all = TRUE) %>%
      mutate(
        aapos = as.integer(aapos),
        VEP_canonical = as.logical(VEP_canonical),
        across(REVEL:gnomAD_genomes_AF, as.numeric)
      )
    
    return(genotypes)
    
  })
  
  ## Format genotypes as tidy data
  tidy_data <- reactive({
    req(goi(), req(genotypes()))
    gather(genotypes(),
           c(REVEL:`LIST-S2`),
           key = "score_type", value = "score")
  }) 
  
  ## Get selected scores
  all_scores <- reactive({
    req(goi(), req(tidy_data()))
    scores <- tidy_data() %>%
      dplyr::filter(!is.na(score), score != "" ) %>%
      distinct(score_type) %>%
      pull()
    return(scores)
  })
  
  # UI elements
  
  ## Update selectscore
  selected_scores <- reactive(input$selectscore)
  
  observeEvent(input$selectgene, {
    all_scores <- all_scores()
    if(!is.null(all_scores)) {
      
      updateSelectizeInput(session, inputId = 'selectscore',
                           choices = all_scores,
                           selected = all_scores[1:2],
                           options = list(
                             maxItems = 3
                           ))
    }
  })
  
  ## Activate selectscore and update button
  observeEvent(input$selectgene, {
    
    if(is.null(input$selectgene) || ! input$selectgene %in% val_goi ) {
      shinyjs::disable("selectscore")
      shinyjs::disable("updatescore")
    } else {
      shinyjs::enable("selectscore")
      shinyjs::enable("updatescore")
    }
  })
  

  # Section 2: Plots --------------------------------------------------------
  colourCount <- reactive({
    req(goi())
    length(unique(protein_info_domain()$description))
  })  
  
  ## Reactive segments
  rct_table1_row <- reactive({input$table1_rows_selected})
  rct_table2_row <- reactive({input$table2_rows_selected})
  
  rct_segments <- reactive({
    #eventReactive(c(input$table1_rows_selected, input$table2_rows_selected), ignoreNULL=TRUE, {
    
    clinvar_goi <- clinvar_goi()
    genotypes <- genotypes()
    toi_nm <- toi_nm()
    
    ## for selected variants
    rct_table1_row <- rct_table1_row()
    rct_table2_row <- rct_table2_row()
    
    ## Format data
    plot_transcript_add_segments <- data.frame(
      cds_pos = integer(), cds_label = character(),
      aa_pos = integer(), aa_label = character(),
      color = character())
    
    if(length(rct_table1_row) > 0) {
      clinvar_goi_selected <- clinvar_goi[rct_table1_row, ]
      clinvar_goi_selected <- clinvar_goi_selected %>%
        mutate(
          cds_label = ifelse(RefSeq.transcript != toi_nm, paste(RefSeq.transcript, CDS.exchange, sep = ": "), CDS.exchange),
          aa_label = ifelse(RefSeq.transcript != toi_nm & !is.na(AA.exchange), paste(RefSeq.transcript, AA.exchange, sep = ": "), AA.exchange),
          color = ifelse(Label == "pathogenic", "red", "blue")) %>%
        dplyr::select(
          cds_pos = Chromosome.position, cds_label,
          aa_pos = AA.position, aa_label,
          color) 
      
      plot_transcript_add_segments <- rbind(plot_transcript_add_segments, clinvar_goi_selected)
    }
    
    if(length(rct_table2_row) > 0) {
      
      genotypes_selected <- genotypes[rct_table2_row, ]
      genotypes_selected <- genotypes_selected %>%
        mutate(
          cds_label = HGVSc,
          color = "grey40") %>%
        dplyr::select(cds_pos = pos, cds_label, 
                      aa_pos = aapos, aa_label = HGVSp,
                      color) 
      
      plot_transcript_add_segments <- rbind(plot_transcript_add_segments, genotypes_selected)
      
    }
    
    return(plot_transcript_add_segments)
    
  })
  

  ## Draw gene plot
  plotTranscript <- reactive({
    
    toi <- toi()
    toi_nm <- toi_nm()
    toi_nm_short <- toi_nm_short()
    transcript_info <- transcript_info()
    exon_info <- exon_info()
    clinvar_goi <- clinvar_goi()
    rct_segments <- rct_segments()
    
    validate(
      need(nrow(transcript_info) > 0, "No transcript information available")
    )
    
    # Base plot
    plot_chromosome <- ggplot() + 
      # Transcript
      geom_rect(data = transcript_info, aes(xmin = tx_seq_start, xmax = tx_seq_end, ymin = -1, ymax = 1), fill = "grey") +
      geom_rect(data = exon_info, aes(xmin = exon_seq_start, xmax = exon_seq_end, ymin = -2, ymax = 2),
                fill="aquamarine4", alpha = .8) +
      ylim(-3,10) +
      theme_minimal() +
      labs(title = ifelse(!is.na( toi_nm ), paste( toi, toi_nm, sep = " / "), toi), x = "Genome position" ) + 
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            panel.grid = element_blank()
            #plot.margin = margin(.1, .1, .1, .2, unit = "in")
      )
    
    # Add segments for pathogenic variants
    if(nrow(clinvar_goi) > 0) {
      
      clinvar_canonical_pathogenic <- dplyr::filter(clinvar_goi, RefSeq.transcript.short == toi_nm_short, Label == "pathogenic")
      
      if(nrow(clinvar_canonical_pathogenic)>0) {
        
        plot_chromosome <- plot_chromosome +
          geom_segment(data = clinvar_canonical_pathogenic,
                       aes(x= Chromosome.position, xend=Chromosome.position,
                           y=2.1, yend=3),
                       color = "black",
                       size = .4) +
          geom_point(data = clinvar_canonical_pathogenic,
                     aes(x = Chromosome.position, y = 3),
                     shape = 21, fill = "red", alpha = .3)
      }
    }
    
    # Add labels
    # Add color for selected rows
    if(nrow(rct_segments) > 0) {
      
      plot_chromosome <- plot_chromosome +
        geom_text_repel(data = rct_segments,
                        mapping = aes(x = cds_pos,
                                      y = 3,
                                      label = cds_label,
                                      color = color),
                        nudge_y = .5,
                        direction = "x",
                        angle = 90,
                        vjust = 1,
                        segment.size = .4,
                        size = 3,
                        segment.linetype = "dotted",
                        max.overlaps = Inf) + 
        geom_segment(data = rct_segments,
                     aes(x = cds_pos, xend = cds_pos,
                         y = 2.1, yend = 3,
                         color = color),
                     size = .4) +
        geom_point(data = rct_segments,
                   aes(x = cds_pos, y = 3, fill = color, color = color),
                   shape = 21, alpha = .3) +
        scale_color_identity() +
        scale_fill_identity() 
      
    }
    
    # Reverse x-scale if gene is located on reverse strand
    if(transcript_info$seq_strand == -1) {
      
      plot_chromosome <- plot_chromosome +
        scale_x_reverse(labels = comma) #position = "top"
      
    } else {
      
      plot_chromosome <- plot_chromosome +
        scale_x_continuous(labels = comma) #position = "top" 
      
    }
    
    return(plot_chromosome)
    
  })
  
  output$plot_transcript <- renderPlot({
    
    req(goi())
    
    plotTranscript()
    
  })
  
  
  ## Draw protein plot
  plotProtein <- reactive({
    req(goi(), nrow(protein_info())>0)
    
    dat <- protein_info_domain()
    protein_length <- protein_length()
    
    plot_protein <- ggplot() + 
      geom_rect(aes(xmin = 1, xmax = protein_length, ymin = -1, ymax = 1), fill = "grey",  size = .25) +
      ylim(-2.5, 10) +
      scale_x_continuous(expand = c(0,0), breaks = x_axis_breaks(protein_length))+
      theme_minimal() +
      labs(title = protein_info()$uniprot_id, x = "Amino acid position") + 
      theme(axis.text.y = element_blank(),
            panel.grid = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "bottom",
            plot.margin = margin(.1, .15, .1, .58, unit = "in")
      )
    
    if("CHAIN" %in% dat$type) {
      plot_protein <- plot_protein +
        geom_rect(data = dat %>%
                    dplyr::filter(type== "CHAIN") %>% 
                    arrange(desc(begin)), #why arrange desc?
                  aes(xmin = begin, xmax = end, ymin = -1, ymax = 1), fill = "grey40", size = .25) + 
        scale_fill_nejm(name = "Chain peptide")
    } 
    
    if("SIGNAL" %in% dat$type ){
      plot_protein <- plot_protein +
        geom_rect(data = dat %>%
                    dplyr::filter(type== "SIGNAL") %>%
                    arrange(desc(begin)), aes(xmin = begin, xmax = end, ymin = -1, ymax = 1, fill = description), size = .25) + 
        scale_fill_nejm(name = "Signal peptide")
    } 
    
    if("REGION" %in% dat$type ) {
      plot_protein <- plot_protein +
        new_scale_fill() +
        geom_rect(data=dat %>%
                    dplyr::filter(type == "REGION") %>%
                    arrange(desc(begin)), aes(xmin=begin, xmax=end, ymin=-2, ymax=2, fill=description), size=.25) +
        scale_fill_manual(values = getPalette(colourCount()), name="Region")
    } 
    
    if(nrow(dat[! dat$type %in% c("CHAIN", "REGION", "SIGNAL"),]) > 0 ) {
      plot_protein <- plot_protein + 
        new_scale_fill() +
        geom_rect(data=dat %>%
                    dplyr::filter(! type %in% c("CHAIN", "REGION", "SIGNAL")) %>%
                    arrange(desc(begin)), aes(xmin=begin, xmax=end, ymin=-.5, ymax=.5, fill=description),size=.25) + 
        scale_fill_manual(values = getPalette2(colourCount()), name="Domain")
    } 
    
    # Add segments
    clinvar_goi <- clinvar_goi()
    
    if(nrow(clinvar_goi)>0) {
      
      clinvar_canonical_pathogenic <- dplyr::filter(clinvar_goi, RefSeq.transcript.short == toi_nm_short(), Label == "pathogenic")
      
      if(nrow(clinvar_canonical_pathogenic)>0) {
        
        plot_protein <- plot_protein +
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
    
    rct_segments <- rct_segments()
    rct_segments <- dplyr::filter(rct_segments, !is.na(aa_pos))
    
    # Add color for selected rows
    if(nrow(rct_segments) > 0) {
      
      plot_protein <- plot_protein + 
        new_scale_fill() +
        scale_fill_manual(values=cols, breaks=legend_order, name = "Clinical significance") +
        geom_text_repel(data = rct_segments,
                        mapping = aes(
                          x = aa_pos, y = 3, 
                          label = aa_label,
                          color = color),
                        nudge_y = 1,
                        direction = "x",
                        angle = 90,
                        vjust = 1,
                        segment.size = .3,
                        size = 3,
                        segment.linetype = "dotted",
                        force = 2,
                        max.overlaps = 200) +
        scale_color_identity() +
        geom_segment(data = rct_segments, 
                     aes(x = aa_pos, xend = aa_pos, y = 2.1, yend = 3, color = color)) +
        geom_point(data = rct_segments,
                   aes(x = aa_pos, y = 3, fill = color, color = color),
                   shape = 21, alpha = .3)
      
    }
    
    return(plot_protein)
    
  })
  
  output$plot_protein <- renderPlot({
    
    plotProtein()
  })
  
  # Density plots ClinVar
  plotDensity <- reactive({
    
    req(goi())
    
    selected_clinsig <- input$selectclinsig #c("benign", "pathogenic")
    protein_length <- protein_length()
    
    plot_density_clinvar_empty <- ggplot() +
      theme_void() +
      scale_x_continuous(limits = c(1, protein_length), breaks = x_axis_breaks(protein_length))+
      scale_y_continuous(expand = c(0,0)) +
      labs(x = "Amino acid position", y = "ClinVar variant density") + 
      theme_classic() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = margin(.1, .15, .1, .4, unit = "in"),
        legend.position = "bottom"
      )
    
    clinvar_goi <- clinvar_goi()
    
    if(nrow(clinvar_goi) > 0 & length(selected_clinsig) > 0) {
      clinvar_goi <- dplyr::filter(clinvar_goi, Label %in% selected_clinsig)
    }
    
    if(nrow(clinvar_goi) > 0 & length(selected_clinsig) > 0) {
      
      cols_selected <- cols[names(cols) %in% selected_clinsig]
      
      plot_density_clinvar <-  ggplot(data = clinvar_goi) + 
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
          plot.margin = margin(.1, .15, .1, .4, unit = "in")
        )
      
    } else {
      plot_density_clinvar <- plot_density_clinvar_empty
    }
    
    return(plot_density_clinvar)
    
  })
  
  output$plot_density <- renderPlot({
    plotDensity()
  })
  
  ## gnomAD
  plotGnomad <- reactive({
    req(goi())
    
    protein_length <- protein_length()
    
    df_gnomad <- genotypes() %>%
      dplyr::select(aapos, gnomAD_exomes_AC, gnomAD_genomes_AC) %>%
      gather(-aapos, key = "Source", value = "Count") %>%
      group_by(aapos, Source) %>%
      summarise(Sum = sum(Count, na.rm = TRUE), .groups = "drop") 
    
    ## Column plot 1
    p_col1 <- ggplot(data = dplyr::filter(df_gnomad, Source == "gnomAD_exomes_AC"), aes(x = aapos, y = Sum)) +
      geom_col(color = "forestgreen", fill = "forestgreen") +
      scale_x_continuous(expand = c(0,0), limits = c(1,protein_length), breaks = x_axis_breaks(protein_length)) +
      scale_y_log10() +
      annotation_logticks(base = 10, sides = "l") +
      theme_minimal() +
      labs(y = "gnomAD exomes count") + 
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        #plot.margin = margin(5.5, 8, 0, 0),
        plot.margin = margin(.1, .15, 0, .05, unit = "in"),
        plot.title = element_text(hjust = .5)
      )
    
    ## Column plot 2
    p_col2 <- ggplot(data = dplyr::filter(df_gnomad, Source == "gnomAD_genomes_AC"), aes(x = aapos, y = Sum)) +
      geom_col(color = "deepskyblue", fill = "deepskyblue") +
      scale_x_continuous(expand = c(0,0), limits = c(1,protein_length), breaks = x_axis_breaks(protein_length)) +
      scale_y_continuous(trans = reverselog_trans(10), expand = c(0,0)) +
      annotation_logticks(base = 10, sides = "l") +
      theme_minimal() +
      labs(y = "gnomAD genomes count", x = "Amino acid position") + 
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, .15, .1, .05, unit = "in")
        #plot.margin = margin(0, 8, 5.5, 0)
      )
    
    ## Combine column plot 1 and 2
    plot_gnomad <- plot_grid(p_col1, p_col2, nrow = 2)
    
    return(plot_gnomad)
    
  })
  output$plot_gnomad <- renderPlot({
    plotGnomad()
  })

  ## scores
  plotScores1 <- reactive({
    
    req(goi())
    
    selected_score <- input$selectscore
    protein_length <- protein_length()
    
    if(length(selected_score) > 0) {
      
      inputs <- selected_score[1]
      
      dat <- tidy_data() %>%
        dplyr::filter(score_type %in% {{inputs}} ) 
      
      dat_hp <- tidy_data() %>%
        dplyr::filter(score_type == {{inputs}}) %>%
        arrange(aapos) %>%
        group_by(aapos) %>%
        summarise(mean = mean(score, na.rm = TRUE))
      
      p1 <- ggplot(dat) +
        geom_smooth(aes(x = aapos, y = score), color = "red") + 
        theme_minimal() + 
        scale_x_continuous(expand = c(0,0), limits = c(1, protein_length)) +
        labs(title = paste0({{inputs}}, " scores for gene ", goi(), " (", toi(), ")") ) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.ticks.y = element_blank(),
              #plot.margin = margin(5.5, 7, 0, 40),
              plot.margin = margin(.1, .15, 0, .65, unit = "in"),
              plot.title = element_text(hjust = .5)
        )
      
      # colors
      col_hp_low <- dbNSFP_scores$color_min[dbNSFP_scores$score_type == inputs] 
      col_hp_high <- dbNSFP_scores$color_max[dbNSFP_scores$score_type == inputs] 
      
      if(is.na(dbNSFP_scores$score_min[dbNSFP_scores$score_type == inputs])) {
        lim_hp <- c(
          min(dat_hp$mean),
          max(dat_hp$mean)
        )
      } else {
        lim_hp <- sort(c(
          dbNSFP_scores$score_min[dbNSFP_scores$score_type == inputs],
          dbNSFP_scores$score_max[dbNSFP_scores$score_type == inputs]
        ))
      }
      
      p2 <- ggplot(dat_hp, aes(x = aapos, y = 1, fill = mean)) + 
        geom_tile() + 
        labs(x = "Amino acid position") + 
        theme_minimal() + 
        scale_fill_gradient(low = col_hp_low, high = col_hp_high, name={{inputs}}, limits = c(lim_hp[1], lim_hp[2]) ) +
        scale_x_continuous(expand = c(0,0), limits = c(1, protein_length), breaks = x_axis_breaks(protein_length)) +
        theme(legend.position = "bottom",
              panel.grid = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              plot.margin = margin(0, .15, .1, .65, unit = "in")
              #plot.margin = margin(0, 7, 5.5, 40)
        )
      
      plot_scores1 <- plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(1,1))
      
      
    } else {
      
      plot_scores1 <- ggplot() +
        theme_void() +
        theme_minimal() +
        scale_x_continuous(expand = c(0,0), limits = c(1, protein_length )) +
        scale_y_continuous(limits = c(0,1)) +
        labs(title = paste0("Scores for gene ", goi(), " (", toi(), ")") ) +
        theme(
          panel.grid.minor = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(.1, .15, .1, .65, unit = "in")
        )
      
    }
    
    return(plot_scores1)
    
  })
  
  output$plot_scores1 <- renderPlot({
    plotScores1()
  })
  
  ## scores
  plotScores2 <- reactive({
    
    req(goi(), req(length(selected_scores())>1))
    
    inputs <- selected_scores()[2]
    protein_length <- protein_length()
    
    dat <- tidy_data() %>%
      dplyr::filter(score_type %in% {{inputs}} ) 
    
    dat_hp <- tidy_data() %>%
      dplyr::filter(score_type == {{inputs}}) %>%
      arrange(aapos) %>%
      group_by(aapos) %>%
      summarise(mean = mean(score, na.rm = TRUE))
    
    p1 <- ggplot(dat) +
      geom_smooth(aes(x = aapos, y = score), color = "red") + 
      theme_minimal() + 
      scale_x_continuous(expand = c(0,0), limits = c(1, protein_length)) +
      labs(title = paste0({{inputs}}, " scores for gene ", goi(), " (", toi(), ")") ) +
      theme(axis.title = element_blank(),
            axis.text = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(.1, .15, 0, .65, unit = "in"),
            plot.title = element_text(hjust = .5)
      )
    
    # colors
    col_hp_low <- dbNSFP_scores$color_min[dbNSFP_scores$score_type == inputs] 
    col_hp_high <- dbNSFP_scores$color_max[dbNSFP_scores$score_type == inputs] 
    
    if(is.na(dbNSFP_scores$score_min[dbNSFP_scores$score_type == inputs])) {
      
      lim_hp <- c(
        min(dat_hp$mean),
        max(dat_hp$mean)
      )
      
    } else {
      
      lim_hp <- sort(c(
        dbNSFP_scores$score_min[dbNSFP_scores$score_type == inputs],
        dbNSFP_scores$score_max[dbNSFP_scores$score_type == inputs]
      ))
      
    }
    
    p2 <- ggplot(dat_hp, aes(x = aapos, y = 1, fill = mean)) + 
      geom_tile() + 
      labs(x = "Amino acid position") + 
      theme_minimal() + 
      scale_fill_gradient(low = col_hp_low, high = col_hp_high, name={{inputs}}, limits = c(lim_hp[1], lim_hp[2]) ) +
      scale_x_continuous(expand = c(0,0), limits = c(1, protein_length), breaks = x_axis_breaks(protein_length)) +
      theme(legend.position = "bottom",
            panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = margin(0, .15, .1, .65, unit = "in"),
            plot.title = element_text(hjust = .5)
      )
    
    p <- plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(1,1))
    return(p)
    
  })
  
  output$plot_scores2 <- renderPlot({
    plotScores2()
  })
  
  ## scores
  plotScores3 <- reactive({
    
    req(goi(), req(length(selected_scores())>2))
    
    if(length(selected_scores())>2) {
      inputs <- selected_scores()[3]
      protein_length <- protein_length()
      
      dat <- tidy_data() %>%
        dplyr::filter(score_type %in% {{inputs}} ) 
      
      dat_hp <- tidy_data() %>%
        dplyr::filter(score_type == {{inputs}}) %>%
        arrange(aapos) %>%
        group_by(aapos) %>%
        summarise(mean = mean(score, na.rm = TRUE))
      
      p1 <- ggplot(dat) +
        geom_smooth(aes(x = aapos, y = score), color = "red") + 
        theme_minimal() + 
        scale_x_continuous(expand = c(0,0), limits = c(1, protein_length)) +
        labs(title = paste0({{inputs}}, " scores for gene ", goi(), " (", toi(), ")") ) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = margin(.1, .15, 0, .65, unit = "in"),
              plot.title = element_text(hjust = .5)
        )
      
      # colors
      col_hp_low <- dbNSFP_scores$color_min[dbNSFP_scores$score_type == inputs] 
      col_hp_high <- dbNSFP_scores$color_max[dbNSFP_scores$score_type == inputs] 
      
      if(is.na(dbNSFP_scores$score_min[dbNSFP_scores$score_type == inputs])) {
        
        lim_hp <- c(
          min(dat_hp$mean),
          max(dat_hp$mean)
        )
        
      } else {
        
        lim_hp <- sort(c(
          dbNSFP_scores$score_min[dbNSFP_scores$score_type == inputs],
          dbNSFP_scores$score_max[dbNSFP_scores$score_type == inputs]
        ))
        
      }
      
      p2 <- ggplot(dat_hp, aes(x = aapos, y = 1, fill = mean)) + 
        geom_tile() + 
        labs(x = "Amino acid position") + 
        theme_minimal() + 
        scale_fill_gradient(low = col_hp_low, high = col_hp_high, name={{inputs}}, limits = c(lim_hp[1], lim_hp[2]) ) +
        scale_x_continuous(expand = c(0,0), limits = c(1, protein_length), breaks = x_axis_breaks(protein_length)) +
        theme(legend.position = "bottom",
              panel.grid = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              plot.margin = margin(0, .15, .1, .65, unit = "in")
        )
      
      p <- plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(1,1))
      
    } else {
      p <- NULL
    }
    return(p)
    
  })
  
  output$plot_scores3 <- renderPlot({
    plotScores3()
  })
  
  # Violin plot
  
  ## Match bases for genes on reverse strand
  reverse_strand <- c(A="T", T="A", C="G", G="C")
  
  violin_color <- setNames(c("#CC6677", "#88CCEE", "#999933","#DDCC77" ), 
                           c("ClinVar_pathogenic", "ClinVar_benign" , "gnomAD", "InSilico"))
  
  plotViolin <- reactive({
    
    req(goi())
    
    clinvar_goi <- clinvar_goi()
    transcript_info <- transcript_info()
    tidy_data <- tidy_data()
    selected_scores <- input$selectscore
    
    p_empty <- ggplot() + 
      theme_void() + 
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      labs(y = "Score value")
    
    if(length(selected_scores) == 0 ) {
      
      p <- p_empty
      
    } else {
      
      # Format ClinVar variants
      clinvar_sub_score <- clinvar_goi %>%
        mutate(ClinVar = paste("ClinVar", Label, sep = "_")) %>%
        dplyr::select(pos = Chromosome.position, ref = REF, alt = ALT, ClinVar)
      
      # Reverse bases if gene is encoded on the reverse strand
      if(transcript_info$seq_strand == -1) {
        clinvar_sub_score$ref <- names(reverse_strand)[match(clinvar_sub_score$ref, reverse_strand)]
        clinvar_sub_score$alt <- names(reverse_strand)[match(clinvar_sub_score$alt, reverse_strand)]
      }
      
      # Anotate ClinVar variants in the InSilico table
      dat <- tidy_data %>%
        dplyr::filter(score_type %in% selected_scores) %>%
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
        stat_pvalue_manual(dat_stat, label = "p.adj.signif", tip.length = 0, hide.ns = T, x = "group1", 
                           #y.position = rep(violin_y_lim$max_plus, each = 3)
        ) +
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
    }
    return(p)
  })
  
  output$plot_violin <- renderPlot({
    req(goi())
    plotViolin()
  })
  
  ## Download buttons
  output$dwnbtn_transcript <- downloadHandler(
    
    filename = function() {
      paste(goi(), "_transcript_", Sys.Date(), switch(input$figureformat, PNG = '.png', SVG = '.svg'), sep = "")
    },
    content = function(file) {
      ggsave(file, plotTranscript(), device = switch(input$figureformat, PNG = 'png', SVG = 'svg'), bg = "white", width = input$figurewidth, height = input$figureheight)
    }
  )
  
  output$dwnbtn_protein <- downloadHandler(
    
    filename = function() {
      paste(goi(), "_protein_", Sys.Date(), switch(input$figureformat, PNG = '.png', SVG = '.svg'), sep = "")
    },
    content = function(file) {
      ggsave(file, plotProtein(), device = switch(input$figureformat, PNG = 'png', SVG = 'svg'), bg = "white", width = input$figurewidth, height = input$figureheight)
    }
  )
  
  output$dwnbtn_gnomad <- downloadHandler(
    
    filename = function() {
      paste(goi(), "_GnomAD_", Sys.Date(), switch(input$figureformat, PNG = '.png', SVG = '.svg'), sep = "")
    },
    content = function(file) {
      ggsave(file, plotGnomad(), device = switch(input$figureformat, PNG = 'png', SVG = 'svg'), bg = "white", width = input$figurewidth, height = input$figureheight)
    }
  )
  
  output$dwnbtn_density <- downloadHandler(
    
    filename = function() {
      paste(goi(), "_ClinVar_", Sys.Date(), switch(input$figureformat, PNG = '.png', SVG = '.svg'), sep = "")
    },
    content = function(file) {
      ggsave(file, plotDensity(), device = switch(input$figureformat, PNG = 'png', SVG = 'svg'), bg = "white", width = input$figurewidth, height = input$figureheight)
    }
  )
  
  output$dwnbtn_scores1 <- downloadHandler(
    
    filename = function() {
      paste(goi(), "_InSilicoScore_", selected_scores()[1], "_", Sys.Date(), switch(input$figureformat, PNG = '.png', SVG = '.svg'), sep = "")
    },
    content = function(file) {
      ggsave(file, plotScores1(), device = switch(input$figureformat, PNG = 'png', SVG = 'svg'), bg = "white", width = input$figurewidth, height = input$figureheight)
    }
  )
  
  output$dwnbtn_scores2 <- downloadHandler(
    
    filename = function() {
      paste(goi(), "_InSilicoScore_", selected_scores()[2], "_", Sys.Date(), switch(input$figureformat, PNG = '.png', SVG = '.svg'), sep = "")
    },
    content = function(file) {
      ggsave(file, plotScores2(), device = switch(input$figureformat, PNG = 'png', SVG = 'svg'), bg = "white", width = input$figurewidth, height = input$figureheight)
    }
  )
  
  output$dwnbtn_scores3 <- downloadHandler(
    
    filename = function() {
      paste(goi(), "_InSilicoScore_", selected_scores()[3], "_", Sys.Date(), switch(input$figureformat, PNG = '.png', SVG = '.svg'), sep = "")
    },
    content = function(file) {
      ggsave(file, plotScores3(), device = switch(input$figureformat, PNG = 'png', SVG = 'svg'), bg = "white", width = input$figurewidth, height = input$figureheight)
    }
  )
  
  output$dwnbtn_violin <- downloadHandler(
    
    filename = function() {
      paste(goi(), "_ScoreStatistics_", Sys.Date(), switch(input$figureformat, PNG = '.png', SVG = '.svg'), sep = "")
    },
    content = function(file) {
      ggsave(file, plotViolin(), device = switch(input$figureformat, PNG = 'png', SVG = 'svg'), bg = "white", width = input$figurewidth, height = input$figureheight)
    }
  )
  
  # Section 3: Tables --------------------------------------------------------
  # ClinVar table
  
  output$table1 <- DT::renderDataTable({
    req(goi())
    
    tbl <- clinvar_goi()
    
    tbl <- tbl %>%
      mutate(
        `Chromosomal position` = paste(Chromosome, Chromosome.position, sep = ":"),
        Type = gsub('single nucleotide variant', 'SNV', Type)
      ) %>%
      dplyr::select(
        `Chromosomal position`, Gene = Gene.symbol, `RefSeq ID` = RefSeq.transcript, 
        HGVSc = CDS.exchange, HGVSp = AA.exchange, 
        `ClinVar ID` = ClinVar.ID, Type, Origin, Significance = Clinical.significance, Consequence,
        Phenotype, Review) %>%
      datatable(
        extensions = 'Buttons',
        filter = 'top',
        rownames = FALSE,
        options = list(
          dom = 'Bftspl',
          buttons = c('copy', 'csv', 'excel'),
          scrollX = TRUE,
          scrollCollapse = TRUE,
          autoWidth = TRUE, 
          lengthMenu = c(25, 50, 100, 200),
          columnDefs = list(
            list(width = '300px', targets = c(8)),
            list(width = '400px', targets = c(8,10)),
            list(width = '150px', targets = c(11)),
            list(visible = FALSE, targets = c(1))
          )
        )
      ) %>%
      formatStyle(0, target = 'row', rowHeight = '90%', lineHeight = '95%', fontSize = '90%') 
    
    return(tbl)
    
  }, server = FALSE)
  
  # Genotypes table
  output$table2 <- DT::renderDataTable({
    req(goi())
    
    genotypes() %>%
      mutate(
        `Chromosomal position` = paste(chr, pos, sep = ":"),
        HGVSc = factor(HGVSc, levels = HGVSc)
      ) %>%
      dplyr::select(
        `Chromosomal position`, 
        HGVSc, HGVSp,
        REVEL:gnomAD_genomes_AF
      ) %>%
      datatable(
        extensions = 'Buttons',
        filter = 'top',
        rownames = FALSE,
        options = list(
          dom = 'Bftspl',
          buttons = c('copy', 'csv', 'excel'),
          scrollX = TRUE,
          scrollCollapse = TRUE,
          lengthMenu = c(25, 50, 100)
        )
      ) %>%
      formatStyle(0, target = 'row', rowHeight = '90%', lineHeight = '95%', fontSize = '90%') %>%
      formatStyle("SIFT", 
                  fontWeight = styleInterval(.05, values = c("bold", "normal") ),
                  backgroundColor = styleInterval(.05, values = c("darksalmon", "white")) ) %>%
      formatStyle("SIFT4G", 
                  fontWeight = styleInterval(.05, values = c("bold", "normal") ),
                  backgroundColor = styleInterval(.05, values = c("darksalmon", "white")) ) %>%
      formatStyle("Polyphen2_HDIV", 
                  fontWeight = styleInterval(.957, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.957, values = c("white", "darksalmon")) ) %>%
      formatStyle("Polyphen2_HVAR", 
                  fontWeight = styleInterval(.909, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.909, values = c("white", "darksalmon")) ) %>%
      formatStyle("LRT", 
                  fontWeight = styleInterval(.001, values = c("bold", "normal") ),
                  backgroundColor = styleInterval(.001, values = c("darksalmon", "white")) ) %>%
      formatStyle("MutationTaster", 
                  fontWeight = styleInterval(.5, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.5, values = c("white", "darksalmon")) ) %>%
      formatStyle("MutationAssessor", 
                  fontWeight = styleInterval(1.9, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(1.9, values = c("white", "darksalmon")) ) %>%
      formatStyle("FATHMM", 
                  fontWeight = styleInterval(1.5, values = c("bold", "normal") ),
                  backgroundColor = styleInterval(1.5, values = c("darksalmon", "white")) ) %>%
      formatStyle("PROVEAN", 
                  fontWeight = styleInterval(-2.5, values = c("bold", "normal") ),
                  backgroundColor = styleInterval(-2.5, values = c("darksalmon", "white")) ) %>%
      formatStyle("VEST4", 
                  fontWeight = styleInterval(.5, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.5, values = c("white", "darksalmon")) ) %>%
      formatStyle("MetaSVM", 
                  fontWeight = styleInterval(0, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(0, values = c("white", "darksalmon")) ) %>%
      formatStyle("MetaLR", 
                  fontWeight = styleInterval(.5, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.5, values = c("white", "darksalmon")) ) %>%
      formatStyle("M-CAP", 
                  fontWeight = styleInterval(.025, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.025, values = c("white", "darksalmon")) ) %>%
      formatStyle("REVEL", 
                  fontWeight = styleInterval(.4, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.4, values = c("white", "darksalmon")) ) %>%
      formatStyle("MutPred", 
                  fontWeight = styleInterval(.611, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.611, values = c("white", "darksalmon")) ) %>%
      formatStyle("MVP", 
                  fontWeight = styleInterval(c(.7, .75), values = c("normal", "italic", "bold") ),
                  backgroundColor = styleInterval(c(.7, .75), values = c("white", "salmon", "darksalmon")) ) %>%
      formatStyle("MPC", 
                  fontWeight = styleInterval(.6, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.6, values = c("white", "darksalmon")) ) %>%
      formatStyle("PrimateAI", 
                  fontWeight = styleInterval(.803, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.803, values = c("white", "darksalmon")) ) %>%
      formatStyle("DEOGEN2", 
                  fontWeight = styleInterval(.5, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.5, values = c("white", "darksalmon")) ) %>%
      formatStyle("BayesDel_addAF", 
                  fontWeight = styleInterval(.0692655, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.0692655, values = c("white", "darksalmon")) ) %>%
      formatStyle("BayesDel_noAF", 
                  fontWeight = styleInterval(-0.0570105, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(-0.0570105, values = c("white", "darksalmon")) ) %>%
      formatStyle("ClinPred", 
                  fontWeight = styleInterval(.5, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.5, values = c("white", "darksalmon")) ) %>%
      formatStyle("LIST-S2", 
                  fontWeight = styleInterval(.85, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.85, values = c("white", "darksalmon")) ) %>%
      formatStyle("CADD_phred", 
                  fontWeight = styleInterval(20, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(20, values = c("white", "darksalmon")) ) %>%
      formatStyle("DANN", 
                  fontWeight = styleInterval(.99, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.99, values = c("white", "darksalmon")) ) %>%
      formatStyle("fathmm-MKL_coding", 
                  fontWeight = styleInterval(.5, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.5, values = c("white", "darksalmon")) ) %>%
      formatStyle("fathmm-XF_coding", 
                  fontWeight = styleInterval(.5, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.5, values = c("white", "darksalmon")) ) %>%
      formatStyle("Eigen-raw_coding", 
                  fontWeight = styleInterval(0, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(0, values = c("white", "darksalmon")) ) %>%
      formatStyle("Eigen-PC-raw_coding", 
                  fontWeight = styleInterval(0, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(0, values = c("white", "darksalmon")) ) %>%
      formatStyle("GenoCanyon", 
                  fontWeight = styleInterval(.999, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.999, values = c("white", "darksalmon")) ) %>%
      formatStyle("integrated_fitCons", 
                  fontWeight = styleInterval(.7, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.7, values = c("white", "darksalmon")) ) %>%
      formatStyle("GM12878_fitCons", 
                  fontWeight = styleInterval(.7, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.7, values = c("white", "darksalmon")) ) %>%
      formatStyle("H1-hESC_fitCons", 
                  fontWeight = styleInterval(.7, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.7, values = c("white", "darksalmon")) ) %>%
      formatStyle("HUVEC_fitCons", 
                  fontWeight = styleInterval(.99, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.99, values = c("white", "darksalmon")) ) %>%
      formatStyle("GERP++_RS", 
                  fontWeight = styleInterval(2, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(2, values = c("white", "darksalmon")) ) %>%
      formatStyle("phyloP100way_vertebrate", 
                  fontWeight = styleInterval(2, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(2, values = c("white", "darksalmon")) ) %>%
      formatStyle("phyloP30way_mammalian", 
                  fontWeight = styleInterval(2, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(2, values = c("white", "darksalmon")) ) %>%
      formatStyle("phyloP17way_primate", 
                  fontWeight = styleInterval(2, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(2, values = c("white", "darksalmon")) ) %>%
      formatStyle("phastCons100way_vertebrate", 
                  fontWeight = styleInterval(.999, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.999, values = c("white", "darksalmon")) ) %>%
      formatStyle("phastCons30way_mammalian", 
                  fontWeight = styleInterval(.999, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.999, values = c("white", "darksalmon")) ) %>%
      formatStyle("phastCons17way_primate", 
                  fontWeight = styleInterval(.999, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(.999, values = c("white", "darksalmon")) ) %>%
      formatStyle("SiPhy_29way_logOdds", 
                  fontWeight = styleInterval(12, values = c("normal", "bold") ),
                  backgroundColor = styleInterval(12, values = c("white", "darksalmon")) )
    
  }, server = FALSE)
  

  
}
