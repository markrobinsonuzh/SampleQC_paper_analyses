# qc04_real.R
suppressPackageStartupMessages({
  library('ggbeeswarm')
})

make_all_outliers <- function(sqc_dt, scater_dt, miqc_dt) {
  all_outliers = list(
      sqc_dt[,    .(method = "SampleQC",  sample_id, cell_id, outlier)],
      scater_dt[, .(method = "scater",    sample_id, cell_id, outlier)],
      miqc_dt[,   .(method = "miQC",      sample_id, cell_id, outlier)]
    ) %>% rbindlist %>%
    .[, method  := factor(method, levels = names(method_cols)) ] %>%
    .[, barcode := str_extract(cell_id, '(?<=:)[ATCG]+(?=-1)') ]

  return(all_outliers)
}

make_outlier_qcs <- function(sqc_dt, scater_dt, miqc_dt, qc_all) {
  # define metrics to plot
  to_plot     = c('log_counts', 'log_feats', 'logit_mito', 'splice_ratio')

  # join all together
  outlier_qcs = list(
      sqc_dt[,    .(method = "SampleQC",  sample_id, cell_id, outlier)],
      scater_dt[, .(method = "scater",    sample_id, cell_id, outlier)],
      miqc_dt[,   .(method = "miQC",      sample_id, cell_id, outlier)]
    ) %>% rbindlist %>%
    merge(qc_all[, c('cell_id', to_plot), with = FALSE], by = 'cell_id') %>%
    melt.data.table( measure = to_plot, var = 'qc_var', val = 'qc_val')

  return(outlier_qcs)
}

plot_outlier_qcs <- function(outlier_qcs, sel_samples) {
  plot_dt = copy(outlier_qcs) %>% 
    .[, outlier := ifelse(outlier, "outlier", "ok") %>% factor ]
  g = ggplot(plot_dt[ (sample_id %in% sel_samples) & (outlier == 'outlier') ]) + 
    aes( y = qc_val, x = method, fill = method, 
      colour = method, alpha = outlier, linetype = outlier ) +
    geom_violin( adjust = 0.25, scale = "count" ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    scale_colour_manual( values = method_cols, guide = "none" ) +
    scale_fill_manual( values = method_cols, guide = "none" ) +
    scale_alpha_manual( values = c(ok = 1, outlier = 0.2), guide = 'none' ) +
    scale_linetype_manual( values = c(ok = 'solid', outlier = 'dotted') ) +
    coord_flip() +
    facet_grid( sample_id ~ qc_var, scales = 'free_x' ) +
    # facet_grid( qc_var ~ sample_id, scales = 'free_y' ) +
    theme_bw() + 
    theme( panel.grid = element_blank(), legend.position = 'bottom' ) +
    labs( x = 'outlier method', y = 'QC metric value',
      linetype = 'called as?', alpha = 'called as?' )

  # g_ok = ggplot(plot_dt[ (sample_id %in% sel_samples) & outlier == 'ok' ]) + 
  #   aes( y = qc_val, x = method, fill = method, 
  #     colour = method, alpha = outlier, linetype = outlier ) +
  #   geom_violin( adjust = 0.25, scale = "count" ) +
  #   scale_y_continuous( breaks = pretty_breaks() ) +
  #   scale_colour_manual( values = method_cols, guide = "none" ) +
  #   scale_fill_manual( values = method_cols, guide = "none" ) +
  #   scale_alpha_manual( values = c(ok = 1, outlier = 0.2), guide = 'none' ) +
  #   scale_linetype_manual( values = c(ok = 'solid', outlier = 'dotted') ) +
  #   coord_flip() +
  #   facet_grid( qc_var ~ sample_id, scales = 'free_y' ) +
  #   theme_bw() + 
  #   theme( panel.grid = element_blank(), legend.position = 'bottom' ) +
  #   labs( y = 'outlier method', x = 'QC metric value',
  #     linetype = 'called as?', alpha = 'called as?' )

  # g_out = ggplot(plot_dt[ (sample_id %in% sel_samples) & outlier == 'outlier' ]) + 
  #   aes( y = qc_val, x = method, fill = method, 
  #     colour = method, alpha = outlier, linetype = outlier ) +
  #   geom_violin( adjust = 0.25, scale = "count" ) +
  #   scale_y_continuous( breaks = pretty_breaks() ) +
  #   scale_colour_manual( values = method_cols, guide = "none" ) +
  #   scale_fill_manual( values = method_cols, guide = "none" ) +
  #   scale_alpha_manual( values = c(ok = 1, outlier = 0.2) ) +
  #   scale_linetype_manual( values = c(ok = 'solid', outlier = 'dotted') ) +
  #   coord_flip() +
  #   facet_grid( qc_var ~ sample_id, scales = 'free_y' ) +
  #   theme_bw() + 
  #   theme( panel.grid = element_blank(), legend.position = 'bottom' ) +
  #   labs( x = NULL, y = NULL,
  #     linetype = 'called as?', alpha = 'called as?' )

  # g = g_ok + g_out

  return(g)
}
