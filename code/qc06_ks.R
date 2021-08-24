# qc07_ks.R
suppressPackageStartupMessages({
  library('patchwork')
})

calc_pr_overall <- function(sampleqc_pr_dt, true_outliers) {
  pr_overall    = rbind(scater_pr_dt, sampleqc_pr_dt) %>%
    merge(true_outliers, by=c('sample_id', 'cell_id')) %>%
    setorder('method', 'value') %>%
    .[, .(
      value       = value,
      N_good_reported = cumsum(!outlier),
      N_reported    = seq.int(.N),
      N_good_total  = sum(!outlier)
      ), by=.(method)] %>%
    .[, precision   := N_good_reported / N_reported ] %>%
    .[, recall    := N_good_reported / N_good_total ]
  return(pr_overall)
}

calc_pr_by_sim <- function(truth_dt, fit_vals_dt) {
  pr_by_sim   = merge(truth_dt, fit_vals_dt, 
    by = c('sample_id', 'cell_id', 'k_sim')) %>%
    setorder('k_sim', 'k_fit', 'value') %>%
    .[, .(
      value           = value,
      N_good_reported = cumsum(!outlier),
      N_reported      = seq.int(.N),
      N_good_total    = sum(!outlier)
      ), by = c('k_sim', 'k_fit')] %>%
    .[, precision   := N_good_reported / N_reported ] %>%
    .[, recall    := N_good_reported / N_good_total ] %>%
    .[, sim_label := paste('True = ', k_sim)]

  return(pr_by_sim)
}

plot_sample_fits <- function(fit_list, K_grid, pr_by_sim) {
  # plot every fit
  suppressMessages({
    g_list = lapply(seq_along(fit_list), function(i) {
      qc_obj      = fit_list[[i]]
      if (is.null(qc_obj))
        return(plot_spacer())
      title_str   = paste0(
        'True k = ', str_extract(K_grid$k_sim[[i]], '[0-9]+'),
        ', assumed k = ', K_grid$k_fit[[i]]
        )
      sel_sample  = sample(rownames(qc_obj), 1)
      g = plot_fit_over_biaxials(qc_obj, sel_sample) +
        scale_fill_distiller(palette = 'RdBu', trans = 'log10',
          limits = c(1, 100), guide = FALSE) +
        labs(subtitle = title_str, x = NULL, y = NULL)
    })
  })
  g_fits    = wrap_plots(g_list, ncol = length(unique(K_grid$k_sim)))

  # define cutpoints
  def_cut   = default_cuts[method == 'SampleQC']$def_value
  cuts_dt   = copy(pr_by_sim) %>%
    .[, def_dist  := abs(value - def_cut) ] %>%
    .[, .SD[ def_dist == min(def_dist) ], by = c('k_sim', 'k_fit')] %>%
    .[, sim_label := paste('True = ', k_sim)]

  # make nice cols
  n_fits    = length(unique(K_grid$k_fit))
  fit_cols  = RColorBrewer::brewer.pal(n_fits + 2, 'Greys')[-(1:2)]
  line_vals = c(`TRUE` = 'solid', `FALSE` = 'dashed')

  # plot PR curves
  g_pr    = ggplot(pr_by_sim) + 
    aes(x = recall, y = precision, colour = factor(k_fit), 
      linetype = k_fit == k_sim) +
    geom_path(size = 1) + geom_point(data = cuts_dt, size = 3, alpha = 1) +
    scale_linetype_manual(values = line_vals) +
    scale_colour_manual(values = fit_cols) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    coord_cartesian(ylim = c(0.8, 1)) +
    facet_grid(. ~ sim_label) +
    theme_bw() + theme(legend.position = 'bottom', 
      legend.key.size = unit(0.4, "in")) +
    labs(
      colour = 'Assumed k',
      linetype = 'Assumed k is correct?',
      title = 'B. Precision-recall curves by fit'
      )

  g = g_fits / g_pr + plot_layout(heights = c(4, 1)) + 
    plot_annotation(title = '          A. Fits over randomly selected samples')

  return(g)
}
