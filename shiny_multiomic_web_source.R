# Load DESeq2 Results ----
my_load_degs <- function(tissue,
                         group_design,
                         age_young = NULL,
                         age_old = NULL) {
  
  dat_ids <- readRDS("data/Transcriptome/gene_symbol_entrezid_table_2026-09-12_exported.rds")
  
  if (group_design == "linear"){
    de_filename <- paste0("data/Transcriptome/DE_res_", group_design,
                          "/DESeq_", tissue, "_linear_2026-09-12_exported.rds")
  } else if (group_design == "two_group"){
    de_filename <- paste0("data/Transcriptome/DE_res_", group_design,
                          "/DESeq_", tissue, "_m", age_young, "_vs_m", age_old, "_2026-09-12_exported.rds")
  }
  
  de_res <- readRDS(de_filename) %>%
    mutate(ENTREZID = as.character(ENTREZID)) %>%
    left_join(dat_ids, by = "ENTREZID") %>%
    dplyr::select(SYMBOL, everything())
  
  return(de_res)
}

# ORA Visualization ----
my_ora_dp <- function(df_ora_res,
                      bioterm_database,
                      tissue_plot_title = NULL){
  
  grob <- NULL
  
  # Plot
  if (is.null(df_ora_res) || nrow(df_ora_res) == 0) {
    # in case no enriched term was found 
    grob <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("No significantly enriched biological terms were found.",
                              "\nConsider loosing the DE thresholds."),
               size = 6, hjust = 0.5, color = "red3") +
      theme_void() +
      theme(text = element_text(family = "sans"))
  } else if (!is.null(df_ora_res) & nrow(df_ora_res) > 0){
    # if <10 enrichment term found
    if (!(any(colnames(df_ora_res) == "cluster"))){
      
      grob <- df_ora_res %>%
        arrange(p.adjust) %>%
        # wrap biological term because some of them are too long to display
        mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
        mutate(Description = factor(Description, levels = rev(.$Description))) %>%
        ggplot(aes(x = Description, y = -log10(p.adjust), size = Count)) +
        geom_point(color = brewer.pal(3, "Set2")[2]) +
        coord_flip() +
        scale_size_continuous(name = "DEs Count",
                              labels = scales::number_format(accuracy = 1),
                              limits = c(1, NA),
                              range = c(3, 8)) +
        xlab(NULL) +
        ylab("-log10(p.adjust)") +
        labs(title = paste0(bioterm_database, " ORA of ", tissue_plot_title)) +
        theme_light() +
        theme(text = element_text(family = "sans", size = 16),
              axis.text.y = element_text(color = "black", size = 11),
              axis.line = element_line(color = "black"))
    } else {
      # Color by cluster
      
      # if more than 25 pathways were found, show top5 terms from each cluster
      
      if (nrow(df_ora_res) > 25){
        grob <- df_ora_res %>%
          group_by(cluster) %>%
          arrange(p.adjust, .by_group = TRUE) %>%
          dplyr::slice_head(n = 5) %>%
          ungroup()
      } else {
        grob <- df_ora_res
      }
      
      grob <- grob %>%
        arrange(cluster, p.adjust) %>%
        # wrap biological term because some of them are too long to display
        mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
        mutate(Description = factor(Description, levels = rev(.$Description))) %>%
        # dot plot
        ggplot(aes(x = Description, y = -log10(pvalue), size = Count, color = cluster)) +
        geom_point() +
        coord_flip() +
        scale_color_brewer(palette = "Set2") +
        scale_size_continuous(name = "DE Count",
                              labels = scales::number_format(accuracy = 1),
                              limits = c(1, NA),
                              range = c(3, 8)) +
        xlab(NULL) +
        ylab("-log10(p.adjust)") +
        labs(title = paste0(bioterm_database, " ORA of ", tissue_plot_title)) +
        theme_light() +
        theme(text = element_text(family = "sans", size = 16),
              axis.text = element_text(color = "black"),
              axis.text.y = element_text(color = "black", size = 11),
              axis.line = element_line(color = "black"))
    }
  }
  
  return(grob)
}

# GSEA Visualization ----
# As GSEA results are pre-computed, only visualize on-demand
my_gsea_dp <- function(gsea_res,
                       bioterm_database,
                       tissue_plot_title = NULL){
  
  grob <- NULL
  
  if (sum(gsea_res$p.adjust < 0.05) <= 10){
    grob <- gsea_res %>%
      filter(p.adjust < 0.05) %>%
      mutate(signed_log10_pavl = -sign(NES)*log10(p.adjust)) %>%
      arrange(desc(signed_log10_pavl)) %>%
      # wrap biological term because some of them are too long to display
      mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
      mutate(Description = factor(Description, levels = rev(.$Description))) %>%
      ggplot(aes(x = Description, y = signed_log10_pavl, size = abs(NES))) +
      geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
      geom_point(color = brewer.pal(3, "Set2")[2]) +
      coord_flip() +
      scale_size_continuous(name = "|NES|",
                            range = c(3, 8)) +
      xlab(NULL) +
      ylab("signed log10(p.adjust)") +
      labs(title = paste0(bioterm_database, " GSEA of ", tissue_plot_title)) +
      theme_light() +
      theme(text = element_text(family = "sans", size = 16),
            axis.text.y = element_text(color = "black", size = 11),
            axis.line = element_line(color = "black"))
    
  } else {
    
    # if more than 25 pathways were found, show top5 terms from each cluster
    
    if (sum(gsea_res$p.adjust < 0.05) > 25){
      grob <- gsea_res %>%
        group_by(cluster) %>%
        arrange(p.adjust, .by_group = TRUE) %>%
        dplyr::slice_head(n = 5) %>%
        ungroup() 
    } else {
      grob <- gsea_res
    }
    
    grob <- grob %>%
      filter(p.adjust < 0.05) %>%
      mutate(signed_log10_pavl = -sign(NES)*log10(p.adjust)) %>%
      arrange(cluster, desc(signed_log10_pavl)) %>%
      # wrap biological term because some of them are too long to display
      mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
      mutate(Description = factor(Description, levels = rev(.$Description))) %>%
      ggplot(aes(x = Description, y = signed_log10_pavl, size = abs(NES))) +
      geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
      geom_point(aes(color = cluster)) +
      scale_color_brewer(palette = "Set2") +
      coord_flip() +
      scale_size_continuous(name = "|NES|",
                            range = c(3, 8)) +
      xlab(NULL) +
      ylab("signed log10(p.adjust)") +
      labs(title = paste0(bioterm_database, " GSEA of ", tissue_plot_title)) +
      theme_light() +
      theme(text = element_text(family = "sans", size = 16),
            axis.text.y = element_text(color = "black", size = 11),
            axis.line = element_line(color = "black"))
  }
  
  return(grob)
}

# volcano plot ----
my_volcano <- function(input_data, p_threshold, fc_threshold, tooltip_label = "SYMBOL"){
  grob <- input_data %>%
    # remove missing value
    drop_na(log2FoldChange, padj) %>%
    # color label
    mutate(change = case_when(padj < p_threshold & log2FoldChange > fc_threshold ~ "up-regulated",
                              padj < p_threshold & log2FoldChange < -fc_threshold ~ "down-regulated",
                              .default = "others")) %>%
    mutate(change = factor(change, levels = c("up-regulated", "down-regulated", "others"))) %>%
    # tooltip label to show
    dplyr::rename(label_hover = all_of(tooltip_label)) %>%
    mutate(tooltip_id = row_number()) %>%
    ggplot() +
    geom_vline(xintercept = fc_threshold, linetype = "dashed", color = "grey20") +
    geom_vline(xintercept = -fc_threshold, linetype = "dashed", color = "grey20") +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "grey20") +
    ggiraph::geom_point_interactive(aes(x = log2FoldChange, y = -log10(padj), color = change,
                                        tooltip = label_hover, data_id = tooltip_id),
                                    hover_nearest = TRUE,
                                    hover_css = "r:2pt;stroke:black;stroke-width:1px;fill:magenta;",
                                    alpha = 0.8,
                                    size = 1) +
    scale_color_manual(name = "Change with Age",
                       values = c("up-regulated" = "orange",
                                  "down-regulated" = "royalblue",
                                  "others" = "grey30")) +
    xlim(c(-max(abs(input_data$log2FoldChange[!is.na(input_data$log2FoldChange) & !is.na(input_data$padj)])),
           max(abs(input_data$log2FoldChange[!is.na(input_data$log2FoldChange) & !is.na(input_data$padj)])))) +
    xlab("covariate-adjusted log2FC") +
    ylab("-Log10 FDR") +
    theme_classic() +
    theme(text = element_text(family = "sans", size = 14),
          axis.text = element_text(color = "black"))
  
  return(grob)
}


# Metabolome ORA result visualization -----
my_metab_ora_dp <- function(df_ora_res,
                            bioterm_database,
                            tissue_plot_title = NULL) {
  grob <- NULL
  
  # Plot
  if (is.null(df_ora_res) || nrow(df_ora_res) == 0) {
    # in case no enriched term was found 
    grob <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("No biological pathways were found.",
                              "\nConsider loosing the DE thresholds."),
               size = 6, hjust = 0.5, color = "red3") +
      theme_void() +
      theme(text = element_text(family = "sans"))
    
  } else {
    
    # if more than 25 pathways are discovered, show top 25 pathways
    if (nrow(df_ora_res) > 25){
      grob <- df_ora_res %>%
        dplyr::slice_head(n = 25)
      
      g_title <- paste0("Top 25 ", bioterm_database, " Pathways in ", tissue_plot_title)
      
    } else {
      grob <- df_ora_res
      g_title <- paste0(bioterm_database, " Pathways in", tissue_plot_title)
    }
    
    grob <- grob %>%
      arrange(p.adjust) %>%
      # wrap biological term because some of them are too long to display
      mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
      mutate(Description = factor(Description, levels = rev(.$Description))) %>%
      ggplot(aes(x = Description, y = -log10(p.adjust), size = Count)) +
      geom_point(aes(fill = p.adjust), shape = 21, color = "grey20") +
      coord_flip() +
      scale_size_continuous(name = "Metabolites Count",
                            labels = scales::number_format(accuracy = 1),
                            limits = c(1, NA),
                            range = c(3, 8)) +
      xlab(NULL) +
      ylab("-p.adjust") +
      labs(title = g_title) +
      theme_light(base_family = "sans") +
      theme(text = element_text(size = 16),
            axis.text = element_text(color = "black"),
            axis.text.y = element_text(color = "black", size = 11),
            axis.line = element_line(color = "black"))
    
    # if non-significant pathway(s) are picked-up to show, show padj=.05 vertical line
    if(max(grob$data$p.adjust > 0.05)){
      grob <- grob +
        scale_fill_viridis(direction = -1, limit = c(min(0.05, min(grob$data$p.adjust)), NA)) +
        geom_hline(yintercept = -log10(0.05), color = "red3", linetype = "dashed") +
        annotate(
          "text",
          x = 0.5,
          y = -log10(0.05),
          label = "p.adjust = 0.05",
          vjust = -1,
          hjust = 1,
          color = "red3",
          size = 4
        )
    } else {
      grob <- grob +
        scale_fill_viridis(direction = -1)
    }
  }
  
  return(grob)
}