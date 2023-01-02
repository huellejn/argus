#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import dplyr
#' @import ggplot2
#' @import tidyr
#' @import scales
#' @import ggsci
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggrepel geom_text_repel
#' @import cowplot
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom jsonlite fromJSON
#' @import rstatix
#' @importFrom DT renderDataTable datatable formatStyle styleInterval
#' @import shinyjs
#' @importFrom seqminer tabix.read.table
#' @importFrom splitstackshape cSplit
#' @importFrom stats setNames
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
    edb_transcript_info[edb_transcript_info$tx_id == toi(), ]
  })
  
  exon_info <- reactive({
    edb_exon_info[edb_exon_info$tx_id == toi(), ]
  }) 
  
  protein_id <- reactive({
    dplyr::filter(edb_protein_info, tx_id == toi()) %>%
      dplyr::select(-tx_id) %>%
      distinct(uniprot_id) %>%
      pull()
  }) 
  
  protein_info_domain <- reactive({
    
    protein_id <- protein_id()
    
    req(length(protein_id) > 0)
    
    url <- paste0(baseurl, protein_id)
    protein_info_domain_json <- fromJSON(url, flatten=TRUE)
    
    if(length(protein_info_domain_json) > 0) {
      
      protein_info_domain <- protein_info_domain_json[[6]][[1]]
      protein_info_domain <- protein_info_domain %>%
        dplyr::filter(
          type %in% c("SIGNAL", "CHAIN", "DOMAIN", "ZN_FING", "COILED", "TRANSIT", "NP_BIND", "ACT_SITE", "VAR_SEQ", "REPEAT", "REGION"),
          ! description %in% c("Disordered"),
          ! is.na(description)
        ) %>%
        mutate(
          begin = as.numeric(begin),
          end = as.numeric(end)
        )
      
    } else {
      
      protein_info_domain <- NULL
      
    }
    return(protein_info_domain)
  })
  
  protein_length <- reactive({
    
    req(nrow(protein_info_domain()) > 0)
    
    max(protein_info_domain()$end, na.rm = TRUE)
  }) 
  
  ## Get clinvar entries for goi
  clinvar_goi <- reactive({
    clinvar %>%
      dplyr::filter(Gene.symbol == goi() ) %>%
      arrange(Chromosome, Chromosome.position)
  })
  
  ## Get dbNSFP entries for goi
  
  ### Columns to keep
  sc <- c("chr", "pos.1.based.", "ref", "alt",
          "aapos", "genename", "Ensembl_geneid", "Ensembl_transcriptid", "Ensembl_proteinid",   
          "HGVSp_ANNOVAR",
          "HGVSc_VEP", "HGVSp_VEP", "APPRIS",  
          "VEP_canonical",
          "REVEL_score", "CADD_phred", "SIFT_score" , "SIFT4G_score",            
          "Polyphen2_HDIV_score", "Polyphen2_HVAR_score", "PROVEAN_score", "M.CAP_score",             
          "VEST4_score", "FATHMM_score", "MetaSVM_score", "MetaLR_score",            
          "ClinPred_score", "MutationTaster_score", "MutationAssessor_score", "DANN_score",              
          "MutPred_score", "MVP_score", "MPC_score" , "LRT_score",               
          "PrimateAI_score" , "DEOGEN2_score", "BayesDel_addAF_score", "BayesDel_noAF_score",     
          "fathmm.MKL_coding_score", "fathmm.XF_coding_score", "Eigen.raw_coding", "Eigen.PC.raw_coding",     
          "GenoCanyon_score", "integrated_fitCons_score", "GM12878_fitCons_score", "H1.hESC_fitCons_score",   
          "HUVEC_fitCons_score", "GERP.._RS", "phyloP100way_vertebrate", "phyloP30way_mammalian",     
          "phyloP17way_primate" , "phastCons100way_vertebrate", "phastCons30way_mammalian", "phastCons17way_primate",    
          "SiPhy_29way_logOdds", "LIST.S2_score", "gnomAD_exomes_AC", "gnomAD_exomes_AF", "gnomAD_genomes_AC", "gnomAD_genomes_AF"
  )
  
  genotypes <- reactive({
    
    req(goi(), nrow(transcript_info())>0)
    
    # Retrieve scores for selected transcript from dbNSFP database by genomics position (index)
    genotypes <- tabix.read.table(file.path(path_to_db, file_dbNSFP), transcript_info()$position, col.names = TRUE, stringsAsFactors = FALSE)
    
    # Subset by transcript and select columns of interest
    genotypes <- genotypes[grep(toi(), genotypes$Ensembl_transcriptid),sc]
    
    # Produce error if scores are not available for this transcript
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
                                     "REVEL_score", 
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
    df_add_segments <- data.frame(
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
      
      df_add_segments <- rbind(df_add_segments, clinvar_goi_selected)
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
      
      # Combine selected variants from ClinVar and InSilico scores and add ClinVar status (benign, pathogenic, VUS)
      for(i in 1:nrow(genotypes_selected) ) {
        sub <- genotypes_selected[i,]
        if(sub$aa_label %in% clinvar_goi$AA.exchange) {
          label <- clinvar_goi$Label[match(sub$aa_label,clinvar_goi$AA.exchange)]
          genotypes_selected$color[i] <- cols[label]
        }
      }
      
      df_add_segments <- rbind(df_add_segments, genotypes_selected)
      
    }
    
    return(df_add_segments)
    
  })
  
  
  ## Adjust font size when button is clicked
  
  ### Base font size for the plots
  base_font <- eventReactive(input$btn_textsize, {
    input$textsize
  }, ignoreNULL = FALSE)
  
  ### Base font size for the labels (geom_text_repel)
  base_font_label <- eventReactive(input$btn_textsize, {
    input$textsizelabel
  }, ignoreNULL = FALSE)

  
  ## Draw gene plot
  plotTranscript <- reactive({
    
    toi <- toi()
    toi_nm <- toi_nm()
    toi_nm_short <- toi_nm_short()
    transcript_info <- transcript_info()
    exon_info <- exon_info()
    clinvar_goi <- clinvar_goi()
    rct_segments <- rct_segments()
    base_font <- base_font()
    base_font_label <- base_font_label()
    
    validate(
      need(nrow(transcript_info) > 0, "No transcript information available")
    )
    
    # Base plot
    p <- plot_transcript(transcript_info, exon_info, toi, toi_nm, toi_nm_short, clinvar_goi, base_font)
    
    # Add labels and color for selected rows
    if(nrow(rct_segments) > 0) {
      
      p <- p +
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
                        size = base_font_label,
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
      
      p <- p +
        scale_x_reverse(labels = comma)
      
    } else {
      
      p <- p +
        scale_x_continuous(labels = comma)  
      
    }
    
    return(p)
    
  })
  
  output$plot_transcript <- renderPlot({
    
    req(goi())
    
    plotTranscript()
    
  })
  
  
  ## Draw protein plot
  plotProtein <- reactive({
    req(goi(), length(protein_id())>0)
    
    protein_id <- protein_id()
    protein_info_domain <- protein_info_domain()
    protein_length <- protein_length()
    color_count <- colourCount()
    clinvar_goi <- clinvar_goi()
    toi_nm_short <- toi_nm_short()
    base_font <- base_font()
    base_font_label <- base_font_label()
    
    validate(
      need(nrow(protein_info_domain) > 0, "No protein domain information available")
    )
    
    
    p <- plot_protein_domain(protein_id, protein_info_domain, protein_length, color_count, clinvar_goi, toi_nm_short, base_font)
    # Add segments for selected rows
    rct_segments <- rct_segments()
    rct_segments <- dplyr::filter(rct_segments, !is.na(aa_pos))
    
    # Add color for selected rows
    if(nrow(rct_segments) > 0) {
      
      p <- p + 
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
                        size = base_font_label,
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
    
    return(p)
    
  })
  
  output$plot_protein <- renderPlot({
    
    plotProtein()
  })
  
  # Density plots ClinVar
  plotDensity <- reactive({
    
    req(goi())
    
    selected_clinsig <- input$selectclinsig #c("benign", "pathogenic")
    protein_length <- protein_length()
    clinvar_goi <- clinvar_goi()
    rct_segments <- rct_segments()
    base_font <- base_font()
    base_font_label <- base_font_label()
    
    if(nrow(clinvar_goi) > 0 & length(selected_clinsig) > 0) {
      
      clinvar_goi_density <- dplyr::filter(clinvar_goi, Label %in% selected_clinsig)
      
      if( nrow(clinvar_goi_density) > 0) {
        
        cols_selected <- cols[names(cols) %in% selected_clinsig]
        p <- plot_density_clinvar(dat = clinvar_goi_density, protein_length, cols_selected, base_font)
        
  
        # Add labels and color for selected rows
        if(nrow(rct_segments) > 0) {
          
          p <- p +
            geom_segment(data = rct_segments, aes(x = aa_pos, xend = aa_pos, y = -.1, yend = 0), color = rct_segments$color) +
            geom_text_repel(data = rct_segments,
                            mapping = aes(x = aa_pos,
                                          y = 0,
                                          label = aa_label) ,
                            color = "black",
                            nudge_y = .5,
                            direction = "x",
                            angle = 90,
                            segment.size = .4,
                            size = base_font_label,
                            segment.linetype = "dotted",
                            max.overlaps = Inf) 
          
        }

      } else {
        
        p <- plot_density_clinvar_empty(protein_length)
        
      }
    } else {
      
      p <- plot_density_clinvar_empty(protein_length)
      
    }
    
    return(p)
    
  })
  
  output$plot_density <- renderPlot({
    plotDensity()
  })
  
  
  ## gnomAD
  plotGnomad <- reactive({
    req(goi())
    
    protein_length <- protein_length()
    genotypes <- genotypes()
    rct_segments <- rct_segments()
    base_font <- base_font()
    base_font_label <- base_font_label()
    
    p <- plot_gnomad(dat = genotypes, protein_length, dat_segments = rct_segments, base_font, base_font_label)

    return(p)
    
  })
  
  output$plot_gnomad <- renderPlot({
    plotGnomad()
  })

  
  ## scores
  plotScores1 <- reactive({
    
    req(goi())
    
    tidy_data <- tidy_data()
    selected_scores <- input$selectscore
    protein_length <- protein_length()+1
    score_index <- 1
    dat <- tidy_data
    goi <- goi()
    toi <- toi()
    rct_segments <- rct_segments()
    base_font <- base_font()
    base_font_label <- base_font_label()
    
    plot_score(dat = tidy_data, selected_scores, score_index, protein_length, goi, toi, dbNSFP_scores, rct_segments, base_font, base_font_label)
    
  })
  
  output$plot_scores1 <- renderPlot({
    plotScores1()
  })
  
  ## scores
  plotScores2 <- reactive({
    
    req(goi(), req(length(selected_scores())>1))
    
    tidy_data <- tidy_data()
    selected_scores <- input$selectscore
    protein_length <- protein_length()+1
    score_index <- 2
    dat <- tidy_data
    goi <- goi()
    toi <- toi()
    rct_segments <- rct_segments()
    base_font <- base_font()
    base_font_label <- base_font_label()
    
    plot_score(dat = tidy_data, selected_scores, score_index, protein_length, goi, toi, dbNSFP_scores, rct_segments, base_font, base_font_label)
    
  })
  
  output$plot_scores2 <- renderPlot({
    plotScores2()
  })
  
  ## scores
  plotScores3 <- reactive({
    
    req(goi(), req(length(selected_scores())>2))
    
    tidy_data <- tidy_data()
    selected_scores <- input$selectscore
    protein_length <- protein_length()+1
    score_index <- 3
    dat <- tidy_data
    goi <- goi()
    toi <- toi()
    rct_segments <- rct_segments()
    base_font <- base_font()
    base_font_label <- base_font_label()
    
    plot_score(dat = tidy_data, selected_scores, score_index, protein_length, goi, toi, dbNSFP_scores, rct_segments, base_font, base_font_label)
    
  })
  
  output$plot_scores3 <- renderPlot({
    plotScores3()
  })
  
  # Violin plot
  plotViolin <- reactive({
    
    clinvar_goi <- clinvar_goi()
    transcript_info <- transcript_info()
    tidy_data <- tidy_data()
    selected_scores <- input$selectscore
    base_font <- base_font()
    
    plot_violin(clinvar_goi, transcript_info, tidy_data, selected_scores, base_font)
    
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
  
  AA.exchange <- position <- AA.position<- ALT <- CDS.exchange <- `Chromosomal position` <- Chromosome <- Chromosome.position <- ClinVar <- ClinVar.ID <- 
    Clinical.significance <- Consequence <- Ensembl_transcriptid <- Gene.symbol <-
    HGVSc <- HGVSc_VEP <- HGVSp <- HGVSp_VEP <- InSilicoScore <- `LIST-S2` <- Label <- Origin <- 
    Phenotype <- REF <- REVEL <- RefSeq.transcript <- RefSeq.transcript.short <- Review <- 
    Type <- VEP_canonical <- aa_label <- aa_pos <- aapos <- begin <- cds_label <- cds_pos <-chr <- 
    color <- datatable <- description <- end <- exon_seq_end <- exon_seq_start <- 
    formatStyle <- gnomAD <- gnomAD_exomes_AC <- gnomAD_genomes_AC <- 
    gnomAD_genomes_AF <- group2 <- ndensity <- p.adj.signif <- pos <-pos.1.based. <- score <- 
    score_type <- styleInterval <- tx_id <- tx_seq_end <- tx_seq_start <- type <- uniprot_id <- value <- NULL
  
  
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
