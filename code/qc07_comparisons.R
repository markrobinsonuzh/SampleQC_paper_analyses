# qc07_comparisons.R

suppressPackageStartupMessages({
  library('data.table')
  library('forcats')
  library('ggplot2')
  library('magrittr')
  library('patchwork')
  library('scales')
  library('stringr')
  library('viridis')

  library('scater')
  library('flexmix')
})

plot_outliers_one_sample <- function(dataset_dt, scater_dt, qc_dt, qc_names, sel_sample) {
  # merge together
  plot_dt = merge(
    dataset_dt[sample_id == s, c('cell_id', qc_names), with = FALSE], 
    scater_dt[sample_id == s, .(cell_id, scater = outlier)], 
    by = 'cell_id') %>%
    merge(miqc_dt[sample_id == s, .(cell_id, miqc = outlier)], by = 'cell_id') %>%
    melt.data.table(measure = qc_names[-1],
      variable.name='qc_x_name', value.name='qc_x') %>%
    melt.data.table(measure = c('scater', 'miqc'),
      variable.name='outlier_method', value.name='outlier')

  # plot
  dens_dt = plot_dt[, .(cell_id, log_counts, qc_x_name, qc_x)] %>% unique
  g_dens  = ggplot(dens_dt) + 
    aes(y = log_counts, x = qc_x) +
    geom_bin2d() +
    scale_x_continuous( breaks=pretty_breaks() ) + 
    scale_y_continuous( breaks=pretty_breaks() ) +
    scale_fill_distiller(palette = 'RdBu', trans = 'log10') +
    facet_grid( . ~ qc_x_name, scales='free_x' ) +
    theme_bw() + 
    labs( y=qc_names[[1]], x='other QC variable', title = sel_sample )

  # plot
  g_outs  = ggplot() + 
    aes(y = log_counts, x = qc_x) +
    geom_point(data = plot_dt[outlier == FALSE], size = 1, colour = 'grey' ) +
    geom_point(data = plot_dt[outlier == TRUE], size = 1, colour = 'red' ) +
    scale_x_continuous( breaks=pretty_breaks() ) + 
    scale_y_continuous( breaks=pretty_breaks() ) +
    facet_grid( outlier_method ~ qc_x_name, scales='free_x' ) +
    theme_bw() + 
    labs( y=qc_names[[1]], x='other QC variable', colour='outlier?' )

  g = g_dens / g_outs + plot_layout(ncol = 1, heights = c(1, 2))

  return(g)
}
