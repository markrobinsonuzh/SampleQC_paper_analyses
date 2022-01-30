
suppressPackageStartupMessages({
  library('SampleQC')

  library('ggplot2')
  library('patchwork')
  library('scales')
  library('stringr')
  library('forcats')
  library('data.table')
  library('magrittr')
  library('assertthat')

  library('Seurat')
  library('scran')
  library('scRNAseq')
})

get_qc_from_cite <- function(sce) {
  # RNA: total, feats, MT 
  .get_qc_from_assay <- function(assay_name) {
    if (assay_name == c('adt')) {
      x   = counts(altExp(sce))
    } else if (assay_name == 'rna') {
      x   = counts(sce)
    }
    return(data.table(
      log_counts  = Matrix::colSums(x) %>% log10,
      log_feats   = Matrix::colSums(x > 0) %>% log10
      ) %>% setnames(names(.), paste0(names(.), '_', assay_name)))
  }

  # put together
  qc_dt   = c('rna', 'adt') %>%
    lapply(.get_qc_from_assay) %>%
    do.call(cbind, .)

  # add mito prop
  mt_pat    = '^MT-'
  sample_pat  = '^[ATCG]+_(.+)$'
  mt_idx    = sce %>% rownames %>% str_detect(mt_pat)
  x       = counts(sce)
  assert_that( sum(mt_idx) == 13 )
  qc_dt[, logit_mito  := qlogis((Matrix::colSums(x[mt_idx, ]) + 1) / (Matrix::colSums(x) + 2 )) ]

  # get cell info
  meta_dt   = colData(sce) %>% as.data.frame %>%
    as.data.table(keep.rownames='cell_id') %>%
    .[, .(cell_id, tenx_lane, batch)] %>%
    .[, sample_id := paste0('batch', batch, '_', tenx_lane) ] %>%
    setcolorder(c('sample_id', 'cell_id'))

  return(cbind(meta_dt, qc_dt))
}

plot_fit_over_biaxials_grid <- function(qc_obj, sel_sample,
  qc_names=NULL, alpha_cut=0.01) {
  # can we do this?
  assert_that(
    SampleQC:::.check_is_qc_obj(qc_obj) == 'fit',
    msg='SampleQC model must be fit (via `fit_sampleqc`) before calling
    this function')

  # expect to use qc_names specified in qc_obj; qc_names parameter here
  # allows zooming in on subset of qc metrics used
  if (is.null(qc_names))
    qc_names  = metadata(qc_obj)$qc_names

  # get data for cells
  points_dt   = copy(colData(qc_obj)[[sel_sample, 'qc_metrics']]) %>%
    .[, cell_id := colData(qc_obj)$cell_id[[sel_sample]] ] %>%
    melt(
      id    = c('cell_id'),
      measure = qc_names,
      variable.name='qc_var', value.name='qc_val'
      )
  points_dt   = merge(
    points_dt[, .(cell_id, var_x=qc_var, val_x=qc_val)],
    points_dt[, .(cell_id, var_y=qc_var, val_y=qc_val)],
    by='cell_id', allow.cartesian=TRUE
    ) %>% .[ as.integer(var_x) > as.integer(var_y) ]

  # extract ellipses
  ellipses_dt = lapply(seq_along(qc_names),
    function(ii) 
      lapply(setdiff(seq_along(qc_names), ii),
        function(jj)
        SampleQC:::.calc_ellipses_dt(
          qc_obj, sel_sample, ii, jj,
          'dummy', alpha=alpha_cut
          ) %>% .[, var_y := qc_names[[ii]] ] %>%
         .[, var_x := qc_names[[jj]] ]
        ) %>% rbindlist
    ) %>% rbindlist %>%
    .[, other_qc := NULL ] %>%
    setnames(c('qc_x', 'qc_y'), c('val_x', 'val_y')) %>% 
    .[, var_x := factor(var_x, levels=qc_names) ] %>%
    .[, var_y := factor(var_y, levels=qc_names) ] %>%
    .[ as.integer(var_x) > as.integer(var_y) ]

  # extract means
  mu_0    = colData(qc_obj)[[sel_sample, 'mu_0']]
  alpha_j   = colData(qc_obj)[[sel_sample, 'alpha_j']]
  beta_k    = colData(qc_obj)[[sel_sample, 'beta_k']]
  sigma_k   = colData(qc_obj)[[sel_sample, 'sigma_k']]
  means_dt  = sweep(beta_k, 2, mu_0 + alpha_j, '+') %>%
    data.table %>%
    set_colnames(qc_names) %>%
    .[, component := 1:.N] %>%
    melt(id='component', variable.name='qc_var', value.name='qc_val')
  means_dt  = merge(
    means_dt[, .(component, var_x=qc_var, val_x=qc_val)],
    means_dt[, .(component, var_y=qc_var, val_y=qc_val)],
    by='component', allow.cartesian=TRUE
    ) %>% .[ as.integer(var_x) > as.integer(var_y) ]

  # # base range on points, not ellipses
  # x_range = c(
  #   floor(min(points_dt$val_x)*2)/2,
  #   ceiling(max(points_dt$val_x)*2)/2
  #   )
  # y_range = c(
  #   floor(min(points_dt$val_y)*2)/2,
  #   ceiling(max(points_dt$val_y)*2)/2
  #   )

  # plot
  g = ggplot() +
    aes(y=val_y, x=val_x) +
    geom_bin2d( data=points_dt ) +
    geom_path( data=ellipses_dt, aes(group=component) ) +
    geom_point( data=means_dt, size=3) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    scale_shape_manual( values=c(1, 16) ) +
    scale_linetype_manual( values=c('dashed', 'solid') ) +
    scale_fill_distiller( palette='RdBu', trans='log10' ) +
    # coord_cartesian( xlim=x_range, ylim=y_range ) +
    facet_grid( var_y ~ var_x, scales='free' ) +
    theme_bw( base_size=14 ) +
    # theme( aspect.ratio=1 ) +
    labs( y=NULL, x=NULL )

  return(g)
}
