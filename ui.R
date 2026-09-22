# UI -----

ui <- page_navbar(
  title = "Shiny AMMOA (Web Ver 2.0)",
  # Choose appearance theme
  theme = bs_theme(bootswatch = "minty"),
  
  # ---- Disclaimer ----
  nav_panel("Welcome",
            id = "panel_about",
            page_fillable(
              layout_column_wrap(
                width = 1,
                heights_equal = "row",
                
                # Web-version notice
                card(
                  card_header(
                    h3("Shiny AMMOA"),
                    h4("(Shiny Aging Murine Multi-Omic Analyzer)")
                  ),
                  card_body(
                    h4("Web Version Notice"),
                    p(
                      "This web version provides a visual demonstration of Shiny AMMOA, ",
                      "with a limited subset of the features available in the full local version. ",
                      "Large datasets, such as bulk RNA-seq data, may take a long time to process ",
                      "or may fail to complete because the web server has limited computing resources. ",
                      "The web version may also become temporarily unavailable after reaching the monthly ",
                      "usage limit of the hosting service."
                    ),
                    p(
                      strong("The web version is intended primarily as a showcase of the user interface "),
                      "to help prospective users evaluate the software before installing the local version ",
                      strong("and is not suitable for intensive analysis. "),
                      "For more flexible or computationally intensive analyses, please install and use the local version."
                    ),
                    p(
                      a(
                        "Local version GitHub repository",
                        href = "https://github.com/M-Ninomiya-Kanda/Shiny_AMMOA_local",
                        target = "_blank",
                        rel = "noopener noreferrer"
                      )
                    ),
                    p(
                      "For details about features that are unavailable in the web version, ",
                      "see the ",
                      strong("User Guide"),
                      " tab."
                    )
                  )
                ),
                
                # About and version information
                card(
                  card_header(h4("About this App")),
                  card_body(
                    h5("Version History"),
                    tags$ul(
                      tags$li("2025-12-23 — Ver. 0.1 (Pre-release Build)"),
                      tags$li("2026-05-18 — Ver. 1.0 (First-Launch Build)"),
                      tags$li("2026-09-18 — Ver. 2.0")
                    ),
                    
                    hr(),
                    
                    h5("Web Version Source Code"),
                    p(
                      a(
                        "GitHub repository",
                        href = "https://github.com/M-Ninomiya-Kanda/Shiny_AMMOA_web",
                        target = "_blank",
                        rel = "noopener noreferrer"
                      )
                    ),
                    
                    hr(),
                    
                    h5("Citation (Preprint)"),
                    tags$ul(
                      tags$li(
                        tags$p(
                          tags$strong(
                            "Shiny AMMOA: an interactive platform for integrative multi-omics analysis of murine aging"
                          ),
                          style = "margin-bottom: 0.25rem;"
                        ),
                        tags$p(
                          "Mayuka Ninomiya Kanda",
                          style = "margin-bottom: 0.25rem;"
                        ),
                        tags$p(
                          "bioRxiv (2026)",
                          tags$br(),
                          "DOI: ",
                          tags$a(
                            href = "https://doi.org/10.64898/2026.05.18.726091",
                            "10.64898/2026.05.18.726091",
                            target = "_blank",
                            rel = "noopener noreferrer"
                          ),
                          style = "margin-bottom: 0;"
                        )
                      )
                    ),
                    
                    hr(),
                    
                    h5("Author"),
                    p("Mayuka Ninomiya Kanda (QB3, University of California, Berkeley)"),
                    p(
                      "Contact: ",
                      a(
                        "m.k.ninomiya@berkeley.edu",
                        href = "mailto:m.k.ninomiya@berkeley.edu"
                      )
                    )
                  )
                )
              )
            )
  ),
  
  # ---- Navigation A: RNA Enrichment ----
  nav_panel("bulk RNA",
            id = "panel_rna_ea",
            useWaiter(),
            
            # sidebar
            layout_sidebar(
              sidebar = sidebar(
                
                # Appearance of Sidebar
                open = "always", # Sidebar cannot be closed
                width = 300,
                
                # Action Button
                card(
                  helpText("After you choose your input, click 'Submit' button."),
                  actionButton(inputId = "submit_rna_ea",
                               label = strong("Submit"))
                ),
                
                # Users' Choice -----
                card(
                  card_header(strong("Choose Your Input")),
                  
                  # Tissue Select
                  card(
                    selectInput(
                      label = strong("Tissue"),
                      inputId = "rna_tissue_ea",
                      choices = choices_tissue$rna,
                      selected = choices_tissue$rna[1]
                    )
                  ),
                  
                  # DESeq2 Design
                  card(
                    radioButtons(
                      inputId = "rna_ea_algorithm",
                      label = strong("Enrichment Analysis Algorithm"),
                      choices = c("ORA" = "ora",
                                  "GSEA" = "gsea"),
                      selected = "ora"
                    ),
                    
                    helpText(
                      "GSEA is available for analysis across all ages only."
                    ),
                    
                    # Database
                    radioButtons(
                      inputId = "rna_database",
                      label = strong("Pathway Database"),
                      choices = c("KEGG", "GO(BP)", "GO(MF)", "GO(CC)"),
                      selected = "KEGG"
                    )
                  ),
                  
                  # DESeq2 Design
                  card(
                    radioButtons(
                      inputId = "rna_group_design",
                      label = strong("Age Groups"),
                      choices = choices_design[["gsea"]],
                      selected = "linear"
                    ),
                    
                    # age groups to be compared (only applicable when design == two_group)
                    conditionalPanel(
                      condition = "input.rna_group_design == 'two_group'",
                      p(strong("Age Group (month)")),
                      fluidRow(
                        column(6, selectInput("age1", "Young", choices = 3, selected = 3)),
                        column(6, selectInput("age2", "Old", choices = 18, selected = 18))
                      )
                    ),
                    
                    # Web-ver only: Help message (Web version uses a fixed, pre-computed age-group comparison)
                    conditionalPanel(
                      condition = "input.rna_group_design == 'two_group'",
                      helpText(
                        "In the web version, RNA-seq comparison is limited to 3-month vs 21-month. Other age-group comparisons are available in the local version."
                      )
                    )
                  ),
                  
                  conditionalPanel(
                    condition = "input.rna_ea_algorithm == 'ora'",
                    # Threshold to Define DEGs (only for ORA)
                    card(
                      p(strong("DEG Threshold")),
                      p("FDR < 0.05 & |Log2FC| > 0.5"),
                      helpText(
                        "Thresholds are fixed in the web app. Use the local app to adjust them."
                      ),
                      radioButtons(
                        inputId = "rna_de_usage",
                        label = strong("DEs Used for Enrichment Analysis"),
                        choices = c("both", "up-regulated", "down-regulated"),
                        selected = "both"
                      )
                    )
                  ),
                  
                )
              ),
              
              # Main Panel -----
              navset_card_pill(
                selected = "rna_enrichment",
                # Tab A1: Volcano Plot
                nav_panel(
                  title = "Volcano Plot",
                  value = "rna_volcano",
                  
                  div(
                    style = "display:flex; flex-direction:column; height:90vh;",
                    uiOutput("rna_vp_tab_title"),
                    # volcano plot
                    div(style = "flex:4; display:flex; gap:20px;",
                        # left:Volcano plot
                        div(
                          style = "flex:3;",
                          girafeOutput("rna_graph_vp", height = "100%", width = "100%")),
                        # right:Selectize Input
                        div(style = "flex:1;",
                            selectizeInput(
                              inputId = "rna_vp_show_labels",
                              label = "Labels to Show:",
                              choices = NULL,
                              multiple = TRUE
                            ),
                            ## show label action button
                            div(
                              style = "margin-top: 10px;",
                              uiOutput("rna_vp_update_button")
                            ),
                            
                        )
                    ),
                    # download button for plot and table
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("rna_download_vp", "Download Plot"),
                      downloadButton("rna_download_de_table", "Download Table")
                    ),
                    # table showing gene name, fold change, etc.
                    div(style = "flex:5; overflow-y:auto;",
                        p(strong("Results Table")),
                        reactableOutput("rna_table_vp", height = "100%")),
                  )
                ),
                
                # Tab A2: Enrichment Analysis
                nav_panel(
                  title = "Pathway Enrichment",
                  value = "rna_enrichment",
                  div(
                    style = "display:flex; flex-direction:column; height:120vh;",
                    uiOutput("rna_ea_tab_title"),
                    # dot plot showing entiched biological terms
                    div(style = "flex:5;", plotOutput("rna_graph_dp", height = "100%")),
                    # download button for plot and table
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("rna_download_dp", "Download Plot"),
                      downloadButton("rna_download_ea_table", "Download Table")
                    ),
                    # table showing all enriched terms
                    div(style = "flex:3; overflow-y:auto;",
                        p(strong("All Enriched Biological Terms")),
                        reactableOutput("rna_table_ea", height = "100%"))
                  )
                ),
                
                # Tab A3: Sample Metadata
                nav_panel(
                  title = "Sample Metadata",
                  value = "rna_metadata",
                  p(
                    strong("Sample Size Analysed"),
                  ),
                  div(
                    tableOutput("rna_metadata_table")
                  ),
                  div(
                    downloadButton("rna_download_meta_table", "Download Table")
                  )
                )
              )
            )
  ),
  
  # ---- Navigation A': Proteome Enrichment ----
  nav_panel("Proteome",
            id = "panel_prot_ea",
            useWaiter(),
            
            # sidebar
            layout_sidebar(
              sidebar = sidebar(
                
                # Appearance of Sidebar
                open = "always", # Sidebar cannot be closed
                width = 300,
                
                # Action Button
                card(
                  helpText("After you choose your input, click 'Submit' button."),
                  actionButton(inputId = "submit_prot_ea",
                               label = strong("Submit"))
                ),
                
                # Users' Choice -----
                card(
                  card_header(strong("Choose Your Input")),
                  
                  # Tissue Select
                  card(
                    selectInput(
                      label = strong("Tissue"),
                      inputId = "prot_tissue_ea",
                      choices = choices_tissue$protein,
                      selected = choices_tissue$protein[1]
                    )
                  ),
                  
                  # DESeq2 Design
                  card(
                    radioButtons(
                      inputId = "prot_ea_algorithm",
                      label = strong("Enrichment Analysis Algorithm"),
                      choices = c("ORA" = "ora",
                                  "GSEA" = "gsea"),
                      selected = "ora"
                    ),
                    
                    # Database
                    radioButtons(
                      inputId = "prot_database",
                      label = strong("Pathway Database"),
                      choices = c("KEGG", "GO(BP)", "GO(MF)", "GO(CC)"),
                      selected = "KEGG"
                    )
                  ),
                  
                  # DESeq2 Design
                  # Comparison Design
                  # Like datasource input, this section doesn't have multiple choices
                  card(
                    p(strong("Age Group")),
                    p("Young: 8 months"),
                    p("Old: 18 months")
                  ),
                  
                  conditionalPanel(
                    condition = "input.prot_ea_algorithm == 'ora'",
                    # Threshold to Define DEGs (only for ORA)
                    card(
                      p(strong("DEP Threshold")),
                      p("FDR < 0.1 & |Log2FC| > 0"),
                      helpText(
                        "Thresholds are fixed in the web app. Use the local app to adjust them."
                      ),
                      radioButtons(
                        inputId = "prot_de_usage",
                        label = strong("DEs Used for Enrichment Analysis"),
                        choices = c("both", "up-regulated", "down-regulated"),
                        selected = "both"
                      )
                    )
                  ),
                  
                )
              ),
              
              # Main Panel -----
              navset_card_pill(
                selected = "prot_enrichment",
                # Tab A1: Volcano Plot
                nav_panel(
                  title = "Volcano Plot",
                  value = "prot_volcano",
                  
                  div(
                    style = "display:flex; flex-direction:column; height:90vh;",
                    uiOutput("prot_vp_tab_title"),
                    # volcano plot
                    div(style = "flex:4; display:flex; gap:20px;",
                        # left:Volcano plot
                        div(
                          style = "flex:3;",
                          girafeOutput("prot_graph_vp", height = "100%", width = "100%")),
                        # right:Selectize Input
                        div(style = "flex:1;",
                            selectizeInput(
                              inputId = "prot_vp_show_labels",
                              label = "Labels to Show:",
                              choices = NULL,
                              multiple = TRUE
                            ),
                            ## show label action button
                            div(
                              style = "margin-top: 10px;",
                              uiOutput("prot_vp_update_button")
                            ),
                            
                        )
                    ),
                    # download button for plot and table
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("prot_download_vp", "Download Plot"),
                      downloadButton("prot_download_de_table", "Download Table")
                    ),
                    # table showing gene name, fold change, etc.
                    div(style = "flex:5; overflow-y:auto;",
                        p(strong("Results Table")),
                        reactableOutput("prot_table_vp", height = "100%")),
                  )
                ),
                
                # Tab A2: Enrichment Analysis
                nav_panel(
                  title = "Pathway Enrichment",
                  value = "prot_enrichment",
                  div(
                    style = "display:flex; flex-direction:column; height:120vh;",
                    uiOutput("prot_ea_tab_title"),
                    # dot plot showing entiched biological terms
                    div(style = "flex:5;", plotOutput("prot_graph_dp", height = "100%")),
                    # download button for plot and table
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("prot_download_dp", "Download Plot"),
                      downloadButton("prot_download_ea_table", "Download Table")
                    ),
                    # table showing all enriched terms
                    div(style = "flex:3; overflow-y:auto;",
                        p(strong("All Enriched Biological Terms")),
                        reactableOutput("prot_table_ea", height = "100%"))
                  )
                ),
                
                # Tab A3: Sample Metadata
                nav_panel(
                  title = "Sample Metadata",
                  value = "prot_metadata",
                  p(
                    strong("Sample Size Analysed"),
                  ),
                  div(
                    tableOutput("prot_metadata_table")
                  ),
                  div(
                    downloadButton("prot_download_meta_table", "Download Table")
                  )
                )
              )
            )
  ),
  
  # ---- Navigation A'' Metabolome Enrichment ----
  nav_panel("Metabolite",
            id = "panel_metab_ea",
            useWaiter(),
            
            # sidebar
            layout_sidebar(
              sidebar = sidebar(
                
                # Appearance of Sidebar
                open = "always", # Sidebar cannot be closed
                width = 300,
                
                # Action Button
                card(
                  helpText("After you choose your input, click 'Submit' button."),
                  actionButton(inputId = "submit_metab_ea",
                               label = strong("Submit")) # strong() -> bold style
                ),
                
                # Users' Choice -----
                card(
                  card_header(strong("Choose Your Input")),
                  
                  # Tissue Select
                  card(
                    selectInput(
                      label = strong("Tissue"),
                      inputId = "met_tissue_ea",
                      choices = c("Brain", "Colon", "Diaphragm", "Eye", "Heart", "Jejunum", "Kidney", "Liver",
                                  "Lung", "Pancreas", "Quadriceps", "Serum", "Skin(Ear)", "Soleus", "Spleen"),
                      selected = "Kidney"
                    )
                  ),
                  
                  card(
                    radioButtons(
                      label = strong("Enrichment Analysis Algorithm"),
                      inputId = "met_ea_algorithm",
                      choices = c("ORA" = "ora"),
                      selected = "ora",
                      inline = FALSE
                    ),
                    
                    # Database
                    radioButtons(
                      inputId = "metab_database",
                      label = strong("Pathway Database"),
                      choices = c("KEGG", "SMPDB"),
                      selected = "KEGG"
                    )
                  ),
                  
                  # Comparison Design
                  # Like datasource input, this section doesn't have multiple choices
                  card(
                    p(strong("Age Group")),
                    p("Young: 15-19 weeks"),
                    p("Old: 90-99 weeks")
                  ),
                  
                  # Threshold to Define DEGs/DEPs
                  card(
                    p(strong("DE Threshold")),
                    p("FDR < 0.1 & |Log2FC| > 0"),
                    helpText(
                      "Thresholds are fixed in the web app. Use the local app to adjust them."
                    ),
                    
                    radioButtons(
                      inputId = "metab_de_usage",
                      label = strong("DEs Used for Enrichment Analysis"),
                      choices = c("both", "up-regulated", "down-regulated"),
                      selected = "both"
                    ),
                  ),
                  
                  # How to set universe (i.e., background metabolites)
                  card(
                    radioButtons(
                      inputId = "metab_universe",
                      label = strong("Background Metabolites"),
                      choiceNames = c("Measured metabolites only",
                                      "All metabolites registered in DB"),
                      choiceValues = c("measured", "all_known"),
                      selected = "measured"
                    ),
                    
                    helpText(
                      "* 'All metabolites registered in DB' uses all KEGG/SMPDB metabolites; ",
                      "it may detect more pathways but can increase false positives."
                    )
                  )
                )
              ),
              
              # Main Panel -----
              navset_card_pill(
                # Tab A1: Volcano Plot
                nav_panel(
                  title = "Metabolite Volcano Plot",
                  
                  div(
                    style = "display:flex; flex-direction:column; height:90vh;",
                    uiOutput("met_vp_tab_title"),
                    # volcano plot
                    div(style = "flex:4; display:flex; gap:20px;",
                        # left：Volcano plot
                        div(
                          style = "flex:3;",
                          girafeOutput("graph_metab_vp", height = "100%", width = "100%")),
                        # right：Selectize Input
                        div(style = "flex:1;",
                            selectizeInput(
                              inputId = "vp_metab_show_labels",
                              label = "Labels to Show:",
                              choices = NULL,
                              multiple = TRUE
                            ),
                            # show label action button placeholder
                            uiOutput("vp_metab_update_button")
                        )
                    ),
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("download_metab_vp", "Download Plot"),
                      downloadButton("download_metab_de_table", "Download Table")
                    ),
                    # table showing gene name, fold change, etc.
                    div(style = "flex:5; overflow-y:auto;",
                        p(strong("Results Table")),
                        reactableOutput("table_metab_vp", height = "100%")),
                  )
                ),
                
                # Tab A2: Enrichment Analysis
                nav_panel(
                  title = "Metabolite Pathway Enrichment",
                  div(
                    style = "display:flex; flex-direction:column; height:120vh;",
                    uiOutput("met_ea_tab_title"),
                    # dot plot showing entiched biological terms
                    div(style = "flex:5;", plotOutput("graph_metab_dp", height = "100%")),
                    # download button for plot and table
                    div(
                      style = "display: flex; justify-content: left; gap: 15px; margin-bottom: 5px;",
                      downloadButton("download_metab_dp", "Download Plot"),
                      downloadButton("download_metab_ea_table", "Download Table")
                    ),
                    # table showing all enriched terms
                    div(style = "flex:3; overflow-y:auto;",
                        p(strong("All Biological Terms")),
                        reactableOutput("table_metab_ea", height = "100%"))
                  )
                ),
                
                # Tab A3: Sample Metadata
                nav_panel(
                  title = "Sample Metadata",
                  p(
                    strong("Sample Size Analysed"),
                  ),
                  div(
                    tableOutput("met_metadata_table")
                  ),
                  div(
                    downloadButton("met_download_meta_table", "Download Table")
                  )
                )
              )
            )
  ),
  
  # ---- Navigation B ----
  nav_panel("Pathway Mapping",
            id = "panel_pv",
            useWaiter(),
            
            # Sidebar ----
            layout_sidebar(
              sidebar = sidebar(
                
                # Appearance of Sidebar
                open = "always", # Sidebar cannot be closed
                width = 300,
                
                # Action Button
                card(
                  helpText("After you choose your input, click 'Submit' button."),
                  
                  actionButton(inputId = "submit_pv",
                               label = strong("Submit"))
                ),
                
                # parameter
                # tissue
                card(
                  p(strong("Sourec of Gene/Protein Nodes")),
                  
                  radioButtons(
                    label = strong("Data Type"),
                    inputId = "source_rect_pv",
                    choiceNames = list(HTML("bulk RNA-seq<br> (Schaum et al 2020)"),
                                       HTML("Proteomics<br> (Keele et al 2023)"),
                                       HTML("Not Show")),
                    choiceValues = list("rna", "protein", "not_show"),
                    selected = "rna",
                    inline = FALSE
                  ),
                  
                  conditionalPanel(
                    condition = "input.source_rect_pv != 'not_show'",
                    selectInput(
                      label = strong("Tissue"),
                      inputId = "tissue_rect_pv",
                      choices = NULL,
                      selected = NULL
                    )
                  ),
                  
                  conditionalPanel(
                    condition = "input.source_rect_pv == 'rna'",
                    radioButtons(
                      inputId = "design_rect_pv",
                      label = strong("Age Groups Used"),
                      choices = c("across all ages" = "all_group",
                                  "two age groups"  = "two_group"),
                      selected = "two_group",
                      inline = TRUE
                    ),
                    helpText("When you choose 'across all ages', compound nodes cannot be shown.")
                  ),
                  
                  conditionalPanel(
                    condition = "input.source_rect_pv == 'rna' && input.design_rect_pv == 'two_group'",
                    p(strong("Age Group (month)")),
                    fluidRow(
                      column(6, selectInput("age1_pv", "Young(Ref)", choices = NULL)),
                      column(6, selectInput("age2_pv", "Old", choices = NULL))
                    ),
                  )
                ),
                
                card(
                  p(strong("Source of Metabolite Nodes")),
                  radioButtons(
                    label = strong("Data Type"),
                    inputId = "source_circ_pv",
                    choiceNames = list(HTML("Metabolites<br> (Jankowski et al 2025)"),
                                       HTML("Not Show")),
                    choiceValues = list("metab", "not_show"),
                    selected = "metab",
                    inline = FALSE
                  ),
                  
                  conditionalPanel(
                    condition = "input.source_circ_pv != 'not_show'",
                    selectInput(
                      label = strong("Tissue"),
                      inputId = "tissue_circ_pv",
                      choices = c("Brain", "Colon", "Diaphragm", "Eye", "Heart", "Jejunum", "Kidney", "Liver",
                                  "Lung", "Pancreas", "Quadriceps", "Serum", "Skin(Ear)", "Soleus", "Spleen"),
                      selected = "Kidney"
                    )
                  ),
                ),
                
                card(
                  # kegg pathway
                  selectInput(
                    label = strong("KEGG Pathway"),
                    inputId = "kegg_id_pv",
                    choices = kegg_ids$label,
                    selected = "mmu00020: Citrate cycle (TCA cycle)"
                  ),
                  # CSS
                  tags$style(HTML(paste(
                    "/* current selection */",
                    "#kegg_pathway + .selectize-control .selectize-input { font-size: 14px; }",
                    "/* drop-down selection list */",
                    "#kegg_pathway + .selectize-control .selectize-dropdown .option { font-size: 14px; }",
                    sep = "\n"
                  ))),
                  p(tags$span(
                    style = "font-size:12px;",
                    "List of KEGG pathways are available at ",
                    tags$a(href = "https://www.kegg.jp/kegg/pathway.html", 
                           "https://www.kegg.jp/kegg/pathway.html", 
                           target = "_blank")
                  )
                  )
                )
              ),
              
              # Navigation B Main Panel -----
              card(
                p(strong("Mapped Pathway")),
                # Legend
                div(style = "height: 150px;",
                    tags$img(
                      src = "pathview_legend_multi.png",
                      style = "height: 80%; width: auto;"
                    )
                ),
                # message
                div(
                  style = "height: 25px;font-size:12px;",
                  p(
                    "KEGG licensing notice: KEGG pathway maps are subject to separate KEGG licensing terms. Non-academic use requires a commercial license. See the ",
                    a(
                      "KEGG Copyright and Disclaimer",
                      href = "https://www.kegg.jp/kegg/legal.html",
                      target = "_blank"
                    ),
                    "."
                  )
                ),
                # pathview image
                div(
                  style="height: 400px; width: 600px; display:flex; align-items:center; justify-content:center;border: 3px solid #ddd;padding: 6px",
                  p(
                    "KEGG diagrams cannot be shown on the Web Ver. due to license restrictions.",
                    br(),
                    "Use the local version if needed.",
                    style="color:#CD0000; font-size:18px"
                  )
                )
              )
            )
  ),
  
  # ---- Navigation C (User Guide)----
  nav_panel("User Guide",
            id = "panel_userguide",
            navset_card_pill(
              nav_panel(
                title = "Pathway Analysis",
                div(
                  style = "width: 200px; margin-bottom: 10px;",
                  downloadButton("download_userguide_analysis", "Download as PDF")
                ),
                div(
                  style = "border: 1px solid #ddd;border-radius: 6px;padding: 8px;background-color: #fafafa;",
                  tags$iframe(
                    src = "userguide_analysis.html",
                    style = "width:100%; height: 85vh; border:none;"
                  )
                )
              ),
              nav_panel(
                title = "Source Data Description",
                div(
                  style = "width: 200px; margin-bottom: 10px;",
                  downloadButton("download_userguide_sourcedata", "Download as PDF")
                ),
                div(
                  style = "border: 1px solid #ddd;border-radius: 6px;padding: 8px;background-color: #fafafa;",
                  tags$iframe(
                    src = "userguide_sourcedata.html",
                    style = "width:100%; height: 85vh; border:none;"
                  )
                )
              )
            )
  )
  
)
