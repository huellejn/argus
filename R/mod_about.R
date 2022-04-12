#' about UI Function
#'
#' @description Shiny module containing the 'About the app' section as of the UI.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_about_ui <- function(id){
  ns <- NS(id)
  tagList(
    
    h2("Abstract"),
    p("The widespread use of high-throughput sequencing techniques is leading to identification of a rapidly increasing
      number of potentially disease-associated genes and pathogenic variants. Pathogenicity assessment of new variants
      can be supported by using publicly available databases and scores. However, these data sources may be difficult to exploit.
      Here, we present aRgus", a("https://argus.urz.uni-heidelberg.de/", href = "https://argus.urz.uni-heidelberg.de/"), 
      ", a stand-alone R/shiny web server application for user-friendly compilation and visualization of gene, protein, variant, 
      and functional impact prediction data. Our application provides a lightweight tool to access multilevel data sources 
      (ENSEMBL, dbNSFP, gnomAD, UniProt, as well as Simple ClinVar), and enables visualization of exon-intron structure and UniProt 
      protein domain annotation, together with ClinVar and gnomAD variant data. aRgus automatically determines the canonical 
      transcript based on the user-supplied HGNC gene symbol and gathers all relevant data. The user can choose from a panel of 
      six visualizations: 1.) unspliced transcript plot; 2.) protein plot; 3.) and 4.) the mutational constraint plots of pathogenic 
      and likely pathogenic ClinVar variants, as well as tolerated gnomAD variants, respectively; 5.) a polynomial regression model 
      with position-coded heatmap depiction of all annotated prediction score values; and 6) groupwise statistical comparison of 
      scores as violin plots. An interactive table is available including all ClinVar variants and all annotated non-synonymous 
      single nucleotide variants with color-coded prediction score values. All plots and tables can be exported separately.
      aRgus enables gene- and position-specific prediction score modeling to assess proteins and identification of regions  susceptible 
      to  missense variation up to single amino acid resolution. It is a powerful tool for enhanced variant interpretation."),
    
    p("This website is free and open to all users and there is no login requirement."),
    br(),
    
    h2("Collaborators"),
    img(src = "www/ukhd_logo.png", height = "100px"),
    p(strong("University Hospital Heidelberg")),
    p(strong("Center for Pediatrics and Adolescent Medicine")),
    p(strong("Division of Pediatric Epileptology")),
    p("Im Neuenheimer Feld 430"),
    p("D-69120 Heidelberg, Germany"),
    br(),
    p("Julian Schr\u00f6ter, MD"),
    p("Steffen Syrbe, MD"),
    br(),
    
    p(strong("University Hospital Heidelberg")),
    p(strong("Center for Pediatrics and Adolescent Medicine")),
    p(strong("Division of Neuropediatrics and Metabolic Medicine")),
    p("Im Neuenheimer Feld 430"),
    p("D-69120 Heidelberg, Germany"),
    br(),
    p("Heiko Brennenstuhl, MD, MBA"),
    p("Tal Dattner"),
    p("Dominic Lenz, MD"),
    p("Prof. Stefan K\u00f6lker, MD"),
    p("Thomas Opladen, MD, MHBA"),
    p("Christian Staufner, MD"),
    p("Prof. Georg F. Hoffmann, MD"),
    br(),
    
    p(strong("University Hospital Heidelberg")),
    p(strong("Institute of Human Genetics")),
    p("Im Neuenheimer Feld 366"),
    p("D-69120 Heidelberg, Germany"),
    br(),
    p("Prof. Christian P. Schaaf, MD"),
    br(),
    
    fluidRow(
      column(1, img(src = "www/iwr_logo.png", height = "100px")), column(1, img(src = "www/emcl_logo.png", height = "100px"))
    ),
    p(strong("Interdisciplinary Center for Scientific Computing (IWR)")),
    p(strong("Engineering Mathematics and Computing Lab (EMCL)")),
    p("Im Neuenheimer Feld 205"),
    p("D-69120 Heidelberg, Germany"),
    br(),
    p("Prof. Vincent Heuveline, PhD (Head of EMCL)"),
    p("Alejandra Jayme, MSc"),
    br(),
    
    img(src = "www/dkfz_logo.png", height = "70px"), 
    img(src = "www/NCT_logo.png", height = "80px"),
    p(strong("German Cancer Research Center (DKFZ)")),
    p(strong("National Center for Tumor Diseases (NCT) Heidelberg")),
    p(strong("Molecular Precision Oncology Program")),
    p(strong("Computational Oncology")),
    p("Im Neuenheimer Feld 460"),
    p("D-69120 Heidelberg, Germany"),
    br(),
    p("Daniel H\u00fcbschmann, MD, PhD"),
    p("Jennifer H\u00fcllein, PhD"),
    p("Sebastian Uhrig, PhD"),
    br(),
    
    img(src = "www/uklp_logo.png", height = "100px"),
    p(strong("Heidelberg University Hospital")),
    p(strong("Institute of Human Genetics")),
    p("Im Neuenheimer Feld 366"),
    p("D-69120 Heidelberg, Germany"),
    br(),
    p("Prof. Christian Schaaf, MD"),
    br(),
    
    p(strong("University Medical Center Leipzig")),
    p(strong("Institute of Human Genetics")),
    p("Philipp-Rosenthal-Str. 55, Building W"),
    p("D-04103 Leipzig, Germany"),
    br(),
    p("Bernt Popp, MD"),
    br(),
    
    h2("Funding"),
    img(src = "www/med_fak_hd_logo.png", height = "100px"),
    p(strong("Physician-Scientist-Program")),
    p(strong("Medical Faculty of the University of Heidelberg")),
    p("Heiko Brennenstuhl, MD, MBA"),
    p("Julian Schr\u00f6ter, MD"),
    br(),
    
    img(src = "www/Logo_DHS.gif", height = "100px"),
    p(strong("Dietmar Hopp Foundation")),
    p("Grant 1DH1813319"),
    p("Steffen Syrbe, MD"),
    p("Julian Schr\u00f6ter, MD"),
    br(),
    
    img(src = "www/DFG_logo.png", height = "70px"),
    p(strong("Deutsche Forschungsgemeinschaft (DFG)")),
    p("Grant PO2366/2-1"),
    p("Bernt Popp, MD"),
    br(),
    
    h2("References"),
    h3("Databases"),
    p("Howe KL, Achuthan P, Allen J, Allen J, Alvarez-Jarreta J, Amode MR, Armean IM, Azov AG, Bennett R, Bhai J, 
                               Billis K, Boddu S, Charkhchi M, Cummins C, Da Rin Fioretto L, Davidson C, Dodiya K, El Houdaigui B, Fatima R, 
                               Gall A, Garcia Giron C, Grego T, Guijarro-Clarke C, Haggerty L, Hemrom A, Hourlier T, Izuogu OG, Juettemann T, 
                               Kaikala V, Kay M, Lavidas I, Le T, Lemos D, Gonzalez Martinez J, Marug\u00e1n JC, Maurel T, McMahon AC, Mohanan S, 
                               Moore B, Muffato M, Oheh DN, Paraschas D, Parker A, Parton A, Prosovetskaia I, Sakthivel MP, Salam AIA, Schmitt BM, 
                               Schuilenburg H, Sheppard D, Steed E, Szpak M, Szuba M, Taylor K, Thormann A, Threadgold G, Walts B, Winterbottom A, 
                               Chakiachvili M, Chaubal A, De Silva N, Flint B, Frankish A, Hunt SE, IIsley GR, Langridge N, Loveland JE, Martin FJ, 
                               Mudge JM, Morales J, Perry E, Ruffier M, Tate J, Thybert D, Trevanion SJ, Cunningham F, Yates AD, Zerbino DR, 
                               Flicek P."),
    
    p(strong("Ensembl 2021.")),
    p("Nucleic Acids Res. 2021 Jan 8;49(D1):D884-D891. doi: 10.1093/nar/gkaa942. PMID: 33137190; 
                               PMCID: PMC7778975."),
    br(),
    
    p("Landrum MJ, Chitipiralla S, Brown GR, Chen C, Gu B, Hart J, Hoffman D, Jang W, Kaur K, Liu C, Lyoshin V, 
                               Maddipatla Z, Maiti R, Mitchell J, O'Leary N, Riley GR, Shi W, Zhou G, Schneider V, Maglott D, Holmes JB, Kattman BL."),
    p(strong("ClinVar: improvements to accessing data.")),
    p("Nucleic Acids Res. 2020 Jan 8;48(D1):D835-D844. doi: 10.1093/nar/gkz972. 
                               PMID: 31777943; PMCID: PMC6943040."),
    br(),
    
    p("P\u00e9rez-Palma E, Gramm M, N\u00fcrnberg P, May P, Lal D."),
    p(strong("Simple ClinVar: an interactive web server to explore and retrieve 
                                   gene and disease variants aggregated in ClinVar database.")),
    p("Nucleic Acids Res. 2019 Jul 2;47(W1):W99-W105. 
                                   doi: 10.1093/nar/gkz411. PMID: 31114901; PMCID: PMC6602488."),
    
    br(),
    p("UniProt Consortium."),
    p(strong("UniProt: the universal protein knowledgebase in 2021.")),
    p("Nucleic Acids Res. 2021 Jan 8;49(D1):D480-D489. 
                                   doi: 10.1093/nar/gkaa1100. PMID: 33237286; PMCID: PMC7778908."),
    br(),
    
    p("Karczewski KJ, Francioli LC, Tiao G, Cummings BB, Alf\u00f6ldi J, Wang Q, Collins RL, Laricchia KM, Ganna A, Birnbaum DP, 
                                   Gauthier LD, Brand H, Solomonson M, Watts NA, Rhodes D, Singer-Berk M, England EM, Seaby EG, Kosmicki JA, Walters RK, 
                                   Tashman K, Farjoun Y, Banks E, Poterba T, Wang A, Seed C, Whiffin N, Chong JX, Samocha KE, Pierce-Hoffman E, Zappala Z, 
                                   O'Donnell-Luria AH, Minikel EV, Weisburd B, Lek M, Ware JS, Vittal C, Armean IM, Bergelson L, Cibulskis K, Connolly KM, 
                                   Covarrubias M, Donnelly S, Ferriera S, Gabriel S, Gentry J, Gupta N, Jeandet T, Kaplan D, Llanwarne C, Munshi R, Novod S, 
                                   Petrillo N, Roazen D, Ruano-Rubio V, Saltzman A, Schleicher M, Soto J, Tibbetts K, Tolonen C, Wade G, Talkowski ME; 
                                   Genome Aggregation Database Consortium, Neale BM, Daly MJ, MacArthur DG."),
    p(strong("The mutational constraint spectrum quantified from variation in 141,456 humans.")),
    p("Nature. 2020 May;581(7809):434-443. doi: 10.1038/s41586-020-2308-7. Epub 2020 May 27. 
                                   Erratum in: Nature. 2021 Feb;590(7846):E53. PMID: 32461654; PMCID: PMC7334197."),
    br(),
    
    p("Liu X, Li C, Mou C, Dong Y, Tu Y."),
    p(strong("dbNSFP v4: a comprehensive database of transcript-specific functional predictions 
                                   and annotations for human nonsynonymous and splice-site SNVs.")),
    p("Genome Med. 2020 Dec 2;12(1):103. doi: 10.1186/s13073-020-00803-9. PMID: 33261662; PMCID: PMC7709417."),
    br(),
    
    p("Rodriguez JM, Pozo F, Cerd\u00e1n-V\u00e9lez D, Di Domenico T, V\u00e1zquez J, Tress ML."),
    p(strong("APPRIS: selecting functionally important isoforms.")),
    p("Nucleic Acids Res. 2021 Nov 10:gkab1058. doi: 10.1093/nar/gkab1058. Epub ahead of print. PMID: 34755885."),
    
    
    h3("Software"),
    
    p("RStudio Team (2020)."),
    p(strong("RStudio: Integrated Development for R.")),
    p("RStudio, PBC, Boston, MA URL ", a("http://www.rstudio.com/", href = "http://www.rstudio.com/"), "."),
    br(),
    
    p("Winston Chang, Joe Cheng, JJ Allaire, Carson Sievert, Barret Schloerke, Yihui Xie, Jeff Allen, Jonathan McPherson, 
                                   Alan Dipert and Barbara Borges (2021)."),
    p(strong("shiny: Web Application Framework for R.")),
    p("R package version 1.7.1. ", a("https://CRAN.R-project.org/package=shiny", href = "https://CRAN.R-project.org/package=shiny"))
    
  )
}
    
    
## To be copied in the UI
# mod_about_ui("about_ui_1")
    