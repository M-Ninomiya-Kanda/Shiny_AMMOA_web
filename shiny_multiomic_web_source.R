# DESeq2 preprocess ----
# Dont't use in the web version because it is too much buden for the limited RAM in Shiny.io


# clusterProfiler Over Represenatation Analysis ----

my_cluterprofiler_ora <- function(differential_exps,
                                  background_exps,
                                  bioterm_database){
  requireNamespace("clusterProfiler")
  requireNamespace("enrichplot")
  requireNamespace("org.Mm.eg.db")
  
  res_ora <- NULL
  df_ea_show <- NA
  dot_plot <- NULL
  
  # Perform ORA
  if (length(differential_exps) == 0) {
    # In case no DEs were extracted
    dot_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("No DEs were detected with the current thresholds."),
               size = 6, hjust = 0.5, color = "red3") +
      theme_void() +
      theme(text = element_text(family = "sans"))
    df_ea_show <- data.frame()
    
  } else if (length(differential_exps) > 0){
    if (bioterm_database == "KEGG"){
      res_ora <- clusterProfiler::enrichKEGG(gene = unique(differential_exps),
                                             universe = unique(background_exps),
                                             organism = "mmu",
                                             pvalueCutoff = 0.05,
                                             pAdjustMethod = "BH")
      if (!is.null(res_ora)){
        res_ora <- DOSE::setReadable(res_ora, OrgDb = org.Mm.eg.db::org.Mm.eg.db, keyType = "ENTREZID")
      }
    } else {
      res_ora <- clusterProfiler::enrichGO(gene = unique(differential_exps),
                                           universe = unique(background_exps),
                                           OrgDb = org.Mm.eg.db::org.Mm.eg.db,
                                           ont = str_sub(bioterm_database, 4, -2),
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           readable = TRUE)
    }
    
    # Plot
    if (is.null(res_ora) || (nrow(res_ora@result) > 0 && sum(res_ora@result$p.adjust < 0.05) == 0)) {
      # in case no enriched term was found 
      dot_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No significantly enriched biological terms were found.",
                                "\nConsider loosing the DE thresholds."),
                 size = 6, hjust = 0.5, color = "red3") +
        theme_void() +
        theme(text = element_text(family = "sans"))
      df_ea_show <- data.frame()
    } else if (!is.null(res_ora) & sum(res_ora@result$p.adjust < 0.05) > 0){
      # if <10 enrichment term found
      if (sum(res_ora@result$p.adjust < 0.05) <= 10){
        df_ea_show <- res_ora@result %>%
          filter(p.adjust < 0.05)
        
        dot_plot <- df_ea_show %>%
          arrange(p.adjust) %>%
          # wrap biological term because some of them are too long to display
          mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
          mutate(Description = factor(Description, levels = rev(.$Description))) %>%
          ggplot(aes(x = Description, y = -log10(p.adjust), size = Count)) +
          geom_point(color = brewer.pal(3, "Set2")[2]) +
          coord_flip() +
          scale_size_continuous(name = "DEs Count",
                                labels = scales::number_format(accuracy = 1),
                                limits = c(1, NA)) +
          xlab(NULL) +
          ylab("-log10(p.adjust)") +
          labs(title = paste0("Enrichment Analysis of ", bioterm_database)) +
          theme_light() +
          theme(text = element_text(family = "sans", size = 16),
                axis.text.y = element_text(color = "black", size = 11),
                axis.line = element_line(color = "black"))
      } else {
        # make the cluster of enriched terms
        # get distance matrix
        dist_mat <- stats::as.dist(1 - enrichplot::pairwise_termsim(res_ora, showCategory = nrow(res_ora@result))@termsim)
        
        # hierarchical clustering
        dot_plot <- stats::cutree(hclust(dist_mat, method = "ward.D2"), k = 5) %>%
          # convert cluster information to data.frame
          data.frame() %>%
          setNames("cluster") %>%
          rownames_to_column("Description") %>%
          mutate(cluster = factor(cluster)) %>%
          left_join(res_ora@result, by = "Description")
        
        df_ea_show <- dot_plot
        
        # if more than 25 pathways were found, show top5 terms from each cluster
        
        if (nrow(dot_plot) > 25){
          dot_plot <- dot_plot %>%
            group_by(cluster) %>%
            arrange(p.adjust, .by_group = TRUE) %>%
            slice_head(n = 5) %>%
            ungroup() 
        }
        
        dot_plot <- dot_plot %>%
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
                                limits = c(1, NA)) +
          xlab(NULL) +
          ylab("-log10(p.adjust)") +
          labs(title = paste0(bioterm_database, " Enrichment Analysis")) +
          theme_light() +
          theme(text = element_text(family = "sans", size = 16),
                axis.text = element_text(color = "black"),
                axis.text.y = element_text(color = "black", size = 11),
                axis.line = element_line(color = "black"))
      }
    }
  }
  
  out_res <- list()
  out_res$df_ea_show <- df_ea_show
  out_res$grob_dot_plot <- dot_plot
  
  gc();gc()
  return(out_res)
}


# volcano plot ----
my_volcano <- function(input_data, p_threshold, fc_threshold, tooltip_label = "symbol"){
  if (is.null(input_data)){
    return(NULL)
  }
  
  idx <- !is.na(input_data$log2FoldChange) & !is.na(input_data$padj)
  xlim_vals <- c(-max(abs(input_data$log2FoldChange[idx])), max(abs(input_data$log2FoldChange[idx])))
  
  grob <- input_data %>%
    # remove missing value
    # drop_na(log2FoldChange, padj) %>%  # not necessary because they are alrealy removed in the web (lightweight version)
    # color label
    mutate(change = case_when(padj < p_threshold & log2FoldChange > fc_threshold ~ 0,
                              padj < p_threshold & log2FoldChange < -fc_threshold ~ 1,
                              .default = 2)) %>%
    mutate(change = factor(change,
                           levels = 0:2,
                           labels = c("up-regulated", "down-regulated", "others"))) %>%
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
    xlim(xlim_vals) +
    xlab("covariate-adjusted log2FC") +
    ylab("-Log10 FDR") +
    theme_classic() +
    theme(text = element_text(family = "sans", size = 14),
          axis.text = element_text(color = "black"))
  
  return(grob)
}


# Annotate Compound Name to enrichResult object -----
# As enrichKEGG() and enrichHMDB() cannot convert KDGG/HMDB compound ID to common name, I make in-house function
my_metab_id_convert <- function(id_string, keytype, cid_table){
  # split "/" separated ID list
  tmp1 <- str_split(id_string, "/", simplify = TRUE) %>%
    as.vector() %>%
    data.frame(ID = .)
  
  if(keytype == "KEGG"){
    col_pick <- c("KEGG", "kegg_name")
  } else if (keytype == "HMDB"){
    col_pick <- c("HMDB", "compound")
  }
  
  tmp2 <- cid_table %>%
    dplyr::select(all_of(col_pick)) %>%
    setNames(c("ID", "name")) %>%
    left_join(tmp1, ., by = "ID") %>%
    pull(name) %>%
    str_flatten(collapse = "/")
  
  return(tmp2)
}

my_metab_setReadable <- function(res_object, cid_table){
  # detect which compound ID is used in the enrichment analysis result object
  if(str_sub(res_object@result$ID[1], 1, 3) == "SMP"){
    keytype = "HMDB"
  } else if (str_sub(res_object@result$ID[1], 1, 3) == "map"){
    keytype = "KEGG"
  }
  
  out_res_object <- res_object
  
  if (is(res_object) == "enrichResult"){
    tmp <- out_res_object@result %>%
      mutate(geneID = map_chr(geneID, ~ my_metab_id_convert(.x, keytype = keytype, cid_table = cid_table))) %>%
      # change column name as it is the same name as gene enrichment result(e.g., "GeneRatio") and confusing
      dplyr::rename(MetabRatio = GeneRatio,
                    CompoundName = geneID)
  } else if (is(res_object) == "gseaResult"){
    tmp <- out_res_object@result %>%
      mutate(core_enrichment = map_chr(core_enrichment, ~ my_metab_id_convert(.x, keytype = keytype, cid_table = cid_table)))
  }
  
  out_res_object@result <- tmp
  
  return(out_res_object)
}


# Metabolite Over Representation Analysis with MicrobiomeProfiler----

my_microbiomeprofiler_ora <- function(differential_metabs,
                                      specify_background = TRUE,
                                      background_metabs = NULL,
                                      bioterm_database = c("KEGG", "SMPDB"),
                                      cid_table){
  requireNamespace("clusterProfiler")
  requireNamespace("MicrobiomeProfiler")
  
  res_metab_ora <- NULL
  df_metab_ea_show <- NULL
  metab_dot_plot <- NULL
  
  # Perform ORA
  if (length(differential_metabs) == 0) {
    # In case no DEs were extracted
    metab_dot_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("No differential metabolites were detected with the current thresholds."),
               size = 6, hjust = 0.5, color = "red3") +
      theme_void()
    df_metab_ea_show <- data.frame()
    
  } else if (length(differential_metabs) > 0){
    if (bioterm_database == "KEGG" & specify_background){
      # KEGG Enrichment & Specify background metabolites
      res_metab_ora <- clusterProfiler::enrichKEGG(gene = unique(differential_metabs),
                                                   universe = unique(background_metabs),
                                                   organism = "cpd",
                                                   pvalueCutoff = 1,
                                                   pAdjustMethod = "BH",
                                                   minGSSize = 5)
    } else if (bioterm_database == "KEGG" & !specify_background) {
      # KEGG Enrichment & NOT Specify background metabolites
      res_metab_ora <- clusterProfiler::enrichKEGG(gene = unique(differential_metabs),
                                                   organism = "cpd",
                                                   pvalueCutoff = 0.05,
                                                   pAdjustMethod = "BH",
                                                   minGSSize = 5)
    } else if (bioterm_database == "SMPDB" & specify_background){
      # SMPDB Enrichment & Specify background metabolites
      res_metab_ora <- MicrobiomeProfiler::enrichHMDB(metabo_list = unique(differential_metabs),
                                                      universe = unique(background_metabs),
                                                      pvalueCutoff = 1,
                                                      pAdjustMethod = "BH",
                                                      minGSSize = 5)
    } else if (bioterm_database == "SMPDB" & !specify_background){
      # SMPDB Enrichment & NOT Specify background metabolites
      res_metab_ora <- MicrobiomeProfiler::enrichHMDB(metabo_list = unique(differential_metabs),
                                                      pvalueCutoff = 0.05,
                                                      pAdjustMethod = "BH",
                                                      minGSSize = 5)
    }
    
    # if at least one pathway is discovered, convert compound ID to common name
    if (!is.null(res_metab_ora)){
      res_metab_ora <- my_metab_setReadable(res_object = res_metab_ora,
                                            cid_table = cid_table)
    }
    
    # Plot
    if (is.null(res_metab_ora)) {
      # in case no enriched term was found 
      metab_dot_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No biological pathways were found.",
                                "\nConsider loosing the DE thresholds."),
                 size = 6, hjust = 0.5, color = "red3") +
        theme_void() +
        theme(text = element_text(family = "sans"))
      df_metab_ea_show <- data.frame()
    } else {
      df_metab_ea_show <- res_metab_ora@result
      
      # if more than 25 pathways are discovered, show top 25 pathways
      if (nrow(res_metab_ora@result) > 25){
        metab_dot_plot <- res_metab_ora@result %>%
          dplyr::slice_head(n = 25)
        
        g_title <- paste0("Top 25 ", bioterm_database, " Pathways")
      } else {
        metab_dot_plot <- res_metab_ora@result
        g_title <- paste0(bioterm_database, " Pathways")
      }
      
      metab_dot_plot <- metab_dot_plot %>%
        arrange(p.adjust) %>%
        # wrap biological term because some of them are too long to display
        mutate(Description = str_wrap(Description, width = 60, indent = 0, exdent = 0)) %>%
        mutate(Description = factor(Description, levels = rev(.$Description))) %>%
        ggplot(aes(x = Description, y = -log10(p.adjust), size = Count)) +
        geom_point(aes(fill = p.adjust), shape = 21, color = "grey20") +
        coord_flip() +
        scale_size_continuous(name = "Metabolites Count",
                              labels = scales::number_format(accuracy = 1),
                              limits = c(1, NA)) +
        xlab(NULL) +
        ylab("-p.adjust") +
        labs(title = g_title) +
        theme_light(base_family = "sans") +
        theme(text = element_text(size = 16),
              axis.text = element_text(color = "black"),
              axis.text.y = element_text(color = "black", size = 11),
              axis.line = element_line(color = "black"))
      
      # if non-significant pathway(s) are picked-up to show, show padj=.05 vertical line
      if(max(metab_dot_plot$data$p.adjust > 0.05)){
        metab_dot_plot <- metab_dot_plot +
          scale_fill_viridis(direction = -1, limit = c(min(0.05, min(metab_dot_plot$data$p.adjust)), NA)) +
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
        metab_dot_plot <- metab_dot_plot +
          scale_fill_viridis(direction = -1)
      }
    }
  }
  
  out_res <- list()
  out_res$res_metab_ora <- res_metab_ora
  out_res$df_metab_ea_show <- df_metab_ea_show
  out_res$grob_metab_dot_plot <- metab_dot_plot
  
  return(out_res)
}