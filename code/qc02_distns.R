# qc02_distns.R

suppressPackageStartupMessages({
  library('data.table')
  library('magrittr')
  library('stringr')
  library('forcats')
  library('assertthat')
  library('ggplot2')
  library('patchwork')
  library('scales')
})

plot_marginals <- function(qc_vals_dt) {
  plot_dt   = copy(qc_vals_dt) %>%
    setorderv(c('sample_id', 'qc_var', 'qc_val')) %>%
    .[, emp_p   := frank(qc_val)/(.N + 1), 
      by = c('sample_id', 'qc_var')] %>%
    .[, z       := qnorm(emp_p)]
  stats_dt  = plot_dt[, .(
    median  = median(qc_val),
    mad     = mad(qc_val)
    ), by = c('sample_id', 'qc_var')
    ] %>%
    .[, sample_lab := sprintf("sample: %s", sample_id), by = sample_id ]
  g = ggplot(plot_dt) + 
    aes(x = qc_val, y = z) +
    geom_point(size = 0.1) +
    geom_abline(data = stats_dt, aes(intercept = -median/mad, slope = 1/mad)) +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    facet_grid(sample_lab ~ qc_var, scales = 'free') +
    theme_bw() + 
    labs(
      y     = 'Normal theoretical quantile', 
      x     = 'QC metric value'
      )

  return(g)
}

plot_biaxials_one_sample <- function(qc_dt, qc_names) {
  # make data biaxial
  qc_melt   = qc_dt %>% 
    melt(id = c('sample_id', 'cell_id'), 
      measure = qc_names, variable.name = 'qc_var', value.name = 'qc_val')
  biax_dt   = merge(qc_melt, qc_melt,
    by = c('sample_id', 'cell_id'), suffixes = c('_i', '_j'),
    allow.cartesian = TRUE) %>%
    .[, qc_var_i  := factor(qc_var_i, levels = qc_names)] %>%
    .[, qc_var_j  := factor(qc_var_j, levels = qc_names)] %>%
    .[as.integer(qc_var_i) > as.integer(qc_var_j)]

  g = ggplot(biax_dt) + aes(x = qc_val_i, y = qc_val_j) +
    geom_bin2d() + scale_fill_distiller(palette = 'RdBu', trans = 'log10') +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    facet_grid(qc_var_j ~ qc_var_i, scales = 'free') + theme_bw() +
    labs(x = NULL, y = NULL, fill = 'no. cells')
}
