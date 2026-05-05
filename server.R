# Server -----
server <- function(input, output, session) {
  
  # Download Button for Userguides -----
  output$download_userguide_analysis <- downloadHandler(
    filename = "userguide_analysis.pdf",
    content = function(file) {
      file.copy(from = "pdf_userguide_analysis.pdf", to = file, overwrite = TRUE)
    }
  )
  
  output$download_userguide_sourcedata <- downloadHandler(
    filename = "userguide_sourcedata.pdf",
    content = function(file) {
      file.copy(from = "pdf_userguide_sourcedata.pdf", to = file, overwrite = TRUE)
    }
  )
  
  # Dynamic Input Change: Navigation A (RNA/Protein Enrichment Analysis) -----
  observeEvent(input$datasource, {
    # update tissue list in response to data source (RNA / Protein)
    updateSelectInput(
      session,
      "tissue_ea",
      choices = choices_tissue[[input$datasource]],
      selected = choices_tissue[[input$datasource]][1]
    )
    
    # Update available statistical design (only 8month vs 18 month is available for protein)
    updateRadioButtons(
      session,
      "design",
      choices = choices_design[[input$datasource]],
      selected = "two_group"
    )
    
    # initialize age choice
    # Web-ver only: Age group is fixed and only 3-month vs 21-month comparison is available
    if (input$datasource == "protein") {
      updateSelectInput(session, "age1", choices = 8, selected = 8)
      updateSelectInput(session, "age2", choices = 18, selected = 18)
    } else if (input$datasource == "rna") {
      updateSelectInput(session, "age1", choices = 3, selected = 3)
      updateSelectInput(session, "age2", choices = 21, selected = 21)
    }
    
    # run update even in initial state
    ignoreInit = FALSE
  })
  
  # "Old" age has to be older than "Young" age
  # Web-ver only: Age group is fixed and only 3-month vs 21-month comparison is available
  
  # Navigation A and A': Show Volcano Plot label update button When volcano plot is shown -----
  output$vp_update_button <- renderUI({
    req(!is.null(ea_data$gg))
    actionButton(inputId = "submit_label_vp", label = "Update Label")
  })
  
  output$vp_metab_update_button <- renderUI({
    req(!is.null(ea_metab_data$gg_metab))
    actionButton(inputId = "submit_label_vp_metab", label = "Update Label")
  })

  # Dynamic Input Change: Navigation B (Pathview) -----
  observeEvent(input$source_rect_pv, {
    # update tissue list in response to data source (RNA / Protein)
    updateSelectInput(
      session,
      "tissue_rect_pv",
      choices = choices_tissue[[input$source_rect_pv]],
      selected = choices_tissue[[input$source_rect_pv]][1]
    )
    
    # initialize age choice
    # Web-ver only: Age group is fixed and only 3-month vs 21-month comparison is available
    if (input$source_rect_pv == "rna") {
      updateSelectInput(session, "age1_pv", choices = 3, selected = 3)
      updateSelectInput(session, "age2_pv", choices = 21, selected = 21)
    }
    
    # run update even in initial state
    ignoreInit = FALSE
  })
  
  # "Old" age has to be older than "Young" age
  # Web-ver only: Age group is fixed and only 3-month vs 21-month comparison is available
  
  # when source of rectangle is RNA and use all groups, metabolite cannot be shown
  observeEvent(list(input$source_rect_pv, input$design_rect_pv), ignoreInit = TRUE, {
    
    if (input$source_rect_pv == "rna" & input$design_rect_pv == "all_group") {
      updateRadioButtons(session, "source_circ_pv",
                         choiceNames = list(HTML("Not Show")),
                         choiceValues = list("not_show"),
                         selected = "not_show"
      )
    } else {
      updateRadioButtons(session, "source_circ_pv",
                         choiceNames = list(HTML("Metabolites<br> (Jankowski et al 2025)"),
                                            HTML("Not Show")),
                         choiceValues = list("metab", "not_show"),
                         selected = "metab",
                         inline = FALSE
      )
    }
  })
  
  # Define all the reactive values -----
  # Navigation A: RNA/protein enrichment analysis
  ea_data <- reactiveValues(
    dat_vp = NULL,
    gg = NULL,
    dp = NULL,
    dat_ea_show = NULL
  )
  
  # Navigation A': metabolite enrichment analysis
  # define reactive value
  ea_metab_data <- reactiveValues(
    dat_metab_vp = NULL,
    gg_metab = NULL,
    dp_metab = NULL,
    dat_metab_ea_show = NULL
  )
  
  # Navigation B: Pathview
  file_pv <- reactiveVal(NULL)
  file_pv_legend <- reactiveVal(NULL)
  has_gene <- reactiveVal(TRUE)
  has_cpd <- reactiveVal(TRUE)
  img_bin <- reactiveVal(NULL)
  
  # in-house function to reset all the reactive values ----
  reset_ea_data <- function() {
    ea_data$dat_vp <- NULL
    ea_data$gg <- NULL
    ea_data$dp <- NULL
    ea_data$dat_ea_show <- NULL
    
    # Selectize List for shown label in VP
    updateSelectizeInput(
      session,
      "vp_show_labels",
      choices = NULL,
      server = TRUE
    )
    
    gc();gc()
  }
  
  reset_pv_data <- function(){
    file_pv(NULL)
    file_pv_legend(NULL)
    has_gene(NULL)
    has_cpd(NULL)
    img_bin(NULL)
    
    gc();gc()
  }
  
  # Navigation A (RNA/Protein Enrichment Analysis) -----

  observeEvent(input$submit_ea, {
    library(ggiraph)
    
    # show waiter
    w_ea <- Waiter$new(html = spin_3(), color = transparent(0.5))
    w_ea$show()
    
    # ORA -----
    dat_vp <- NULL
    # format fold change and p-value data
    if (input$datasource == "rna"){
      # Web-ver only: load per-computed DEseq2 result instead of calculating DESeq2 on-demand
      if (input$design == "two_group"){
        dat_vp <- readRDS(paste0("data/RNA_3vs21/bulkTMS_DESeq2_res_", input$tissue_ea, "_3vs21_exported_2025_12_23.rds"))
      } else if (input$design == "linear") {
        dat_vp <- readRDS(paste0("data/RNA_linear/bulkTMS_DESeq2_res_", input$tissue_ea, "_linear_exported_2025_12_23.rds"))
      }
      
    } else if (input$datasource == "protein"){
      # for promote data, filter prepossessed data into selected tissue
      dat_vp <- readRDS(paste0("data/Proteome/Proteome_res_", input$tissue_ea, "_2025_12_23.rds"))
    }
    
    # define DEGs/DEPs
    diff_exps <- dat_vp %>%
      drop_na(ENTREZID) %>%
      filter(padj < input$p_thres, abs(log2FoldChange) > input$fc_thres)
    
    # define background genes
    bg_exp <- dat_vp %>%
      filter(!is.na(padj)) %>%
      drop_na(ENTREZID) %>%
      pull(ENTREZID)
    
    if (input$deg_usage == "both") {
      diff_exps <- diff_exps %>% pull(ENTREZID)
    } else if (input$deg_usage == "up-regulated") {
      diff_exps <- diff_exps %>% filter(log2FoldChange > 0) %>% pull(ENTREZID)
    } else if (input$deg_usage == "down-regulated") {
      diff_exps <- diff_exps %>% filter(log2FoldChange < 0) %>% pull(ENTREZID)
    }
    
    # execute ORA
    # prepare null results
    res_ea <- NULL
    dp <- NULL
    dat_ea_show <- NULL
    
    res_ea_list <- my_cluterprofiler_ora(differential_exps = diff_exps,
                                         background_exps = bg_exp,
                                         bioterm_database = input$database)
    
    # Navigation A Rendering -----
    # change reactive value
    ea_data$gg <- my_volcano(input_data = dat_vp,
                             p_threshold = input$p_thres,
                             fc_threshold = input$fc_thres)
    ea_data$dat_vp <- dat_vp
    ea_data$dp <- res_ea_list$grob_dot_plot
    ea_data$dat_ea_show <- res_ea_list$df_ea_show
    
    
    # Output in Navigation A
    # interactive volcano plot
    output$graph_vp <- ggiraph::renderGirafe({
      req(!is.null(ea_data$gg))
      girafe(ggobj = ea_data$gg,
             width_svg = 6,
             height_svg = 4.5,
             options = list(
               opts_sizing(rescale = FALSE),
               opts_toolbar(
                 saveaspng = FALSE,
                 hidden = c("selection", "zoom", "misc")
               )
             )
      )
    })
    
    # Update Selectize List for shown label in VP
    updateSelectizeInput(
      session,
      "vp_show_labels",
      choices = ea_data$dat_vp$symbol,
      selected = NULL,
      server = TRUE
    )
    
    # Download Button for VP
    output$download_vp <- downloadHandler(
      filename = "volcano_plot.png",
      content = function(file) {
        ggsave(file, plot = ea_data$gg, width = 6, height = 4.5, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_de_table <- downloadHandler(
      filename = function() {
        if_else(input$design == "linear",
                paste0("result_", input$tissue_ea, "_", input$datasource, "_", input$design, ".csv"),
                paste0("result_", input$tissue_ea, "_", input$datasource, "_", input$age1, "vs", input$age2, ".csv")) 
      },
      content = function(file) {
        write.csv(ea_data$dat_vp, file, row.names = FALSE)
      }
    )
    
    # Table
    output$table_vp <- renderReactable({
      if (is.null(ea_data$dat_vp) || nrow(ea_data$dat_vp) == 0) {
        reactable(data.frame(Message = "No results to display"))
      } else {
        reactable(
          ea_data$dat_vp,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "12px")),
          # reduce digit to look window tidy
          columns = list(
            log2FoldChange = colDef(name = "covariate-adjusted log2FC", format = colFormat(digits = 5)),
            pvalue = colDef(name = "p-val", format = colFormat(digits = 3)),
            padj = colDef(name = "FDR (BH-adjusted p)", format = colFormat(digits = 3)),
            ENTREZID = colDef(show = FALSE)
          )
        )
      }
    }) 
    
    # Enrichment Dot Plot
    output$graph_dp <- renderPlot(
      {ea_data$dp}
    )
    
    # Download Button for Dot Plot
    output$download_dp <- downloadHandler(
      filename = "enrichment_dot_plot.png",
      content = function(file) {
        ggsave(file, plot = ea_data$dp, width = 9, height = 9, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_ea_table <- downloadHandler(
      filename = "all_enriched_terms.csv",
      content = function(file) {
        write.csv(ea_data$dat_ea_show, file, row.names = FALSE)
      }
    )
    
    output$table_ea <- renderReactable({
      if (is.null(ea_data$dat_ea_show) || nrow(ea_data$dat_ea_show) == 0) {
        reactable(data.frame(Message = "No enrichment results to display"))
      } else {
        reactable(
          ea_data$dat_ea_show,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "10px"))
        )
      }
    })
    
    # hide waiter
    w_ea$hide()
  })
  
  
  # update volcano plot if labels to show is updated -----
  observeEvent(input$submit_label_vp, ignoreInit = TRUE, {
    req(!is.null(ea_data$dat_vp))
    ea_data$gg <-  my_volcano(input_data = ea_data$dat_vp,
                              p_threshold = input$p_thres,
                              fc_threshold = input$fc_thres,
                              tooltip_label = "symbol")
    
    if (!is.null(input$vp_show_labels) && length(input$vp_show_labels) > 0){
      ea_data$gg <- ea_data$gg +
        geom_point(data = ea_data$dat_vp %>% filter(symbol %in% input$vp_show_labels),
                   aes(x = log2FoldChange, y = -log10(padj)),
                   size = 2, shape = 21, color = "black", fill = "magenta") +
        ggrepel::geom_label_repel(data = ea_data$dat_vp %>% filter(symbol %in% input$vp_show_labels),
                                  aes(label = symbol, x = log2FoldChange, y = -log10(padj)),
                                  alpha = 0.8,
                                  color = "magenta",
                                  size = 4,
                                  force = 4,
                                  max.overlaps = 10)
    }
  })
  
  # Navigation A' (Metabolite Enrichment) -----
  observeEvent(input$submit_metab_ea, {
    
    # show waiter
    w_metab_ea <- Waiter$new(html = spin_3(), color = transparent(0.5))
    w_metab_ea$show()
    
    # load library for volcano plot
    library(ggiraph)
    
    # initialize RNA/protein enrichment analysis data to release memory
    reset_ea_data()
    
    # ORA -----
    dat_metab_vp <- dat_metab %>%
      filter(tissue == input$tissue_metab_ea) %>%
      dplyr::select(-tissue) %>%
      left_join(dat_compound_id, by = "compound") %>%
      dplyr::select(compound, log2FC_Aged, pval_Aged, padj_Aged, HMDB, KEGG, PubChem, METLIN, kegg_name) %>%
      dplyr::rename(log2FoldChange = log2FC_Aged,
                    pvalue = pval_Aged,
                    padj = padj_Aged)
    
    # define differentially expressed metabolites and background metabolites and extract compound IDs
    if(input$metab_database == "KEGG"){
      metab_col_pick <- c("KEGG", "kegg_name")
    } else if (input$metab_database == "SMPDB"){
      metab_col_pick <- c("HMDB", "compound")
    }
    
    # extract background compound IDs
    bg_metabs <- dat_metab_vp %>%
      dplyr::select(all_of(metab_col_pick[1])) %>%
      setNames("id") %>%
      drop_na(id) %>%
      pull(id)
    
    # extract differentially regulated compound IDs
    diff_metabs <- dat_metab_vp %>%
      filter(padj < input$metab_p_thres, abs(log2FoldChange) > input$metab_fc_thres) %>%
      dplyr::select(all_of(metab_col_pick[1]), log2FoldChange, padj) %>%
      dplyr::rename(id = all_of(metab_col_pick[1])) %>%
      drop_na(id)
    
    # define background genes
    
    if (input$metab_de_usage == "both") {
      diff_metabs <- diff_metabs %>% pull(id)
    } else if (input$metab_de_usage == "up-regulated") {
      diff_metabs <- diff_metabs %>% filter(log2FoldChange > 0) %>% pull(id)
    } else if (input$metab_de_usage == "down-regulated") {
      diff_metabs <- diff_metabs %>% filter(log2FoldChange < 0) %>% pull(id)
    }
    
    # execute ORA
    # prepare null results
    res_metab_ea <- NULL
    dp_metab <- NULL
    dat_meatab_ea_show <- NULL
    
    # Perform ORA and Create dotplot object
    if(input$metab_universe == "measured"){
      res_metab_ea_list <- my_microbiomeprofiler_ora(differential_metabs = diff_metabs,
                                                     specify_background = TRUE,
                                                     background_metabs = bg_metabs,
                                                     bioterm_database = input$metab_database,
                                                     cid_table = dat_compound_id)
    } else if (input$metab_universe == "all_known") {
      res_metab_ea_list <- my_microbiomeprofiler_ora(differential_metabs = diff_metabs,
                                                     specify_background = FALSE,
                                                     bioterm_database = input$metab_database,
                                                     cid_table = dat_compound_id)
    }
    
    # Navigation A' Rendering -----
    # change reactive value
    ea_metab_data$gg_metab <- my_volcano(input_data = dat_metab_vp,
                                         p_threshold = input$metab_p_thres,
                                         fc_threshold = input$metab_fc_thres,
                                         tooltip_label = "compound")
    ea_metab_data$dat_metab_vp <- dat_metab_vp
    ea_metab_data$dp_metab <- res_metab_ea_list$grob_metab_dot_plot
    ea_metab_data$dat_metab_ea_show <- res_metab_ea_list$df_metab_ea_show
    
    # Output in Navigation A'
    # interactive volcano plot
    output$graph_metab_vp <- ggiraph::renderGirafe({
      girafe(ggobj = ea_metab_data$gg_metab,
             width_svg = 6,
             height_svg = 4.5,
             options = list(
               opts_sizing(rescale = FALSE),
               opts_toolbar(
                 saveaspng = FALSE,
                 hidden = c("selection", "zoom", "misc")
               )
             )
      )
    })
    
    # Update Selectize List for shown label in VP
    updateSelectizeInput(
      session,
      "vp_metab_show_labels",
      choices = ea_metab_data$dat_metab_vp$compound,
      selected = NULL,
      server = TRUE
    )
    
    # Download Button for VP
    output$download_metab_vp <- downloadHandler(
      filename = "volcano_plot.png",
      content = function(file) {
        ggsave(file, plot = ea_metab_data$gg_metab, width = 6, height = 4.5, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_metab_de_table <- downloadHandler(
      filename = function() {
        paste0("result_", input$tissue_metab_ea, "_", input$metab_datasource, ".csv")
      },
      content = function(file) {
        write.csv(ea_metab_data$dat_metab_vp, file, row.names = FALSE)
      }
    )
    
    # Table
    output$table_metab_vp <- renderReactable({
      if (is.null(ea_metab_data$dat_metab_vp) || nrow(ea_metab_data$dat_metab_vp) == 0) {
        reactable(data.frame(Message = "No results to display"))
      } else {
        reactable(
          ea_metab_data$dat_metab_vp,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "12px")),
          # reduce digit to look window tidy
          columns = list(
            log2FoldChange = colDef(name = "covariate-adjusted log2FC", format = colFormat(digits = 5)),
            pvalue = colDef(name = "p-val", format = colFormat(digits = 3)),
            padj = colDef(name = "FDR (BH-adjusted p)", format = colFormat(digits = 3))
          )
        )
      }
    }) 
    
    # Enrichment Dot Plot
    output$graph_metab_dp <- renderPlot(
      {ea_metab_data$dp_metab}
    )
    
    # Download Button for Dot Plot
    output$download_metab_dp <- downloadHandler(
      filename = "enrichment_dot_plot.png",
      content = function(file) {
        ggsave(file, plot = ea_metab_data$dp_metab, width = 9, height = 9, dpi = 300)
      }
    )
    
    # Download Button for Result Table
    output$download_metab_ea_table <- downloadHandler(
      filename = "all_enriched_terms.csv",
      content = function(file) {
        write.csv(ea_metab_data$dat_metab_ea_show, file, row.names = FALSE)
      }
    )
    
    output$table_metab_ea <- renderReactable({
      if (is.null(ea_metab_data$dat_metab_ea_show) || nrow(ea_metab_data$dat_metab_ea_show) == 0) {
        reactable(data.frame(Message = "No enrichment results to display"))
      } else {
        reactable(
          ea_metab_data$dat_metab_ea_show,
          searchable = TRUE,
          theme = reactableTheme(style = list(fontSize = "10px"))
        )
      }
    })
    
    
    # hide waiter
    w_metab_ea$hide()
  })
  
  # update volcano plot if metabolite labels to show is updated -----
  observeEvent(input$submit_label_vp_metab, ignoreInit = TRUE, {
    ea_metab_data$gg_metab <-  my_volcano(input_data = ea_metab_data$dat_metab_vp,
                                          p_threshold = input$metab_p_thres,
                                          fc_threshold = input$metab_fc_thres,
                                          tooltip_label = "compound")
    
    if (!is.null(input$vp_metab_show_labels) && length(input$vp_metab_show_labels) > 0){
      ea_metab_data$gg_metab <- ea_metab_data$gg_metab +
        geom_point(data = ea_metab_data$dat_metab_vp %>% filter(compound %in% input$vp_metab_show_labels),
                   aes(x = log2FoldChange, y = -log10(padj)),
                   size = 2, shape = 21, color = "black", fill = "magenta") +
        ggrepel::geom_label_repel(data = ea_metab_data$dat_metab_vp %>% filter(compound %in% input$vp_metab_show_labels),
                                  aes(label = compound, x = log2FoldChange, y = -log10(padj)),
                                  alpha = 0.6,
                                  color = "magenta",
                                  size = 4,
                                  force = 4,
                                  max.overlaps = 10)
    }
  })
  
  # Navigation B (Pathview)-----
  observeEvent(input$submit_pv, {
    # show waiter
    w_pv <- Waiter$new(html = spin_3(), color = transparent(0.5))
    w_pv$show()
    
    # initialize RNA/protein enrichment analysis data to release memory
    reset_ea_data()
    
    # hide waiter
    w_pv$hide()
  })
  
}