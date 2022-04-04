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

  library('scater')
})

# cols from muscat package! https://github.com/HelenaLC/muscat
nice_cols   = c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999"
  )

load_and_prep_sce_file <- function(sce_f, n_hvgs = 1e3) {
  message('loading ', sel_f)
  sce     = sce_f %>% readRDS %>%
    logNormCounts
  message('  finding HVGs')
  var_mod = modelGeneVar(sce)
  hvgs    = getTopHVGs(var_mod, n = n_hvgs)
  message('  running PCA')
  sce     = runPCA(sce, subset_row = hvgs)
  message('  clustering')
  cls     = clusterCells(sce, use.dimred = "PCA")
  sce$cls = cls
  message('  running UMAP')
  sce     = scater::runUMAP(sce, subset_row = hvgs, 
    pca = 50, spread = 1, min_dist = 1)
  message('done!')
  
  return(sce)  
}

make_plot_dt_for_umap <- function(sce, datasets_dt, sel_f, qc_names) {
  # put UMAP and clusters into dt
  plot_dt   = reducedDim(sce, "UMAP") %>% 
    set_colnames(c("UMAP1", "UMAP2")) %>%
    apply(2, rescale, c(0.05, 0.95)) %>%
    as.data.table(keep.rownames = 'cell_id') %>%
    .[, cell_id   := paste0(sel_f, ':', cell_id) ] %>%
    .[, cl_raw    := sce$cls ]

  # add qc metric values and scale them
  plot_dt   = plot_dt %>%
    merge(datasets_dt[ sample_id == sel_f ], by = 'cell_id') %>%
    .[, c('cell_id', 'UMAP1', 'UMAP2', 'cl_raw', qc_names), with = FALSE] %>%
    melt.data.table(measure = qc_names, 
      variable.name = 'qc_var', value.name = 'qc_val') %>%
    .[, qc_scale  := scale(qc_val), by = qc_var] %>%
    .[, qc_scale  := qc_scale %>% pmin(2) %>% pmax(-2), by = qc_var]

  # put clusters in increasing order of mito%
  mito_dt   = plot_dt[ qc_var == 'logit_mito' ] %>%
    .[, .(med_mt = median(qc_val)), by = cl_raw] %>%
    .[ order(med_mt) ] %>% .[, cl := sprintf("K%02d", 1:.N) %>% factor ]
  plot_dt   = merge(plot_dt, mito_dt, by = 'cl_raw') %>% .[, cl_raw := NULL ]

  return(plot_dt)
}

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

plot_umap_clusters <- function(plot_dt) {
  tmp_dt  = unique(plot_dt[, .(UMAP1, UMAP2, cl)])
  g = ggplot(tmp_dt) +
    aes( x = UMAP1, y = UMAP2, colour = cl ) +
    geom_point( size = 0.5 ) + scale_colour_manual( values = nice_cols ) +
    scale_x_continuous( breaks = pretty_breaks(), limits = c(0,1) ) +
    scale_y_continuous( breaks = pretty_breaks(), limits = c(0,1) ) +
    theme_bw() + 
    theme( panel.grid = element_blank(), axis.text = element_blank() ) +
    labs( colour = "cluster" )
  return(g)
}

plot_qc_metrics_over_umap <- function(plot_dt, nrow = 2) {
  qc_list = unique(plot_dt$qc_var)
  g       = lapply(qc_list, .plot_one_qc_metric_over_umap, plot_dt) %>%
    wrap_plots(nrow = nrow)

  return(g)
}

.plot_one_qc_metric_over_umap <- function(this_v, plot_dt) {
  # truncate values
  tmp_dt  = plot_dt[ qc_var == this_v ]
  max_min = median(tmp_dt$qc_val) + 2 * c(-1, 1) * mad(tmp_dt$qc_val)
  tmp_dt  = tmp_dt %>%
    .[, qc_val := qc_val %>% pmax(max_min[[1]]) %>% pmin(max_min[[2]]) ]

  # which direction?
  dir     = 1
  if (this_v == 'logit_mito')
    dir     = -1

  # define nice breaks
  if (this_v %in% c('log_counts', 'log_feats')) {
    brks  = c(1e2, 2e2, 4e2, 1e3, 2e3, 4e3, 1e4, 2e4, 4e4) %>% log10
    labs  = c('100', '200', '400', '1k', '2k', '4k', '10k', '20k', '40k')
  } else if (this_v == 'logit_mito') {
    brks  = c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.25, 0.5, 
      0.75, 0.9, 0.95) %>% qlogis
    labs  = c('0.1%', '0.2%', '0.5%', '1%', '2%', '5%', '10%', '25%', '50%', 
      '75%', '90%', '95%')
  }
  # define nice labels
  labs_ls   = c(
    log_counts  = 'library\nsize', 
    log_feats   = 'no. features', 
    logit_mito  = 'mito.\nproportion'
    )

  # make plot
  g = ggplot(tmp_dt) +
    aes( x = UMAP1, y = UMAP2, colour = qc_val ) +
    geom_point( size = 0.5 ) +
    scale_colour_distiller( palette = "RdBu", direction = dir,
      breaks = brks, labels = labs ) +
    scale_x_continuous( breaks = pretty_breaks(), limits = c(0,1) ) +
    scale_y_continuous( breaks = pretty_breaks(), limits = c(0,1) ) +
    theme_bw() + 
    theme( panel.grid = element_blank(), axis.text = element_blank() ) +
    labs( colour = labs_ls[[this_v]] )

  return(g)
}

plot_violins_qc_metrics <- function(plot_dt, nrow = 1, titles = NULL) {
  qc_list = unique(plot_dt$qc_var)
  if (is.null(titles)) {
    titles  = rep(list(NULL), length(qc_list))
  } else {
    assert_that(length(qc_list) == length(titles))
  }
  names(titles) = qc_list
  g       = lapply(qc_list, .plot_one_qc_violin, plot_dt, titles) %>%
    wrap_plots(nrow = nrow)
}

.plot_one_qc_violin <- function(this_v, plot_dt, titles) {
  tmp_dt  = plot_dt[ qc_var == this_v ]

  # define nice breaks
  if (this_v %in% c('log_counts', 'log_feats')) {
    brks  = c(1e2, 2e2, 4e2, 1e3, 2e3, 4e3, 1e4, 2e4, 4e4) %>% log10
    labs  = c('100', '200', '400', '1k', '2k', '4k', '10k', '20k', '40k')
  } else if (this_v == 'logit_mito') {
    brks  = c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1, 0.25, 0.5, 
      0.75, 0.9, 0.95) %>% qlogis
    labs  = c('0.1%', '0.2%', '0.5%', '1%', '2%', '5%', '10%', '25%', '50%', 
      '75%', '90%', '95%')
  }
  # define nice labels
  labs_ls   = c(
    log_counts  = 'library size', 
    log_feats   = 'no. features', 
    logit_mito  = 'mito. proportion'
    )

  # make plot
  g = ggplot(tmp_dt) +
    aes( x = cl, colour = cl, y = qc_val ) +
    geom_violin() + scale_colour_manual( values = nice_cols, guide = "none" ) +
    scale_y_continuous( breaks = brks, labels = labs ) +
    theme_bw() + 
    theme(
      panel.grid  = element_blank(), 
      axis.text.x = element_text(size = 6)
      ) +
    labs( x = 'cluster', y = labs_ls[[this_v]], title = titles[[this_v]] )
}

make_patchwork_plot <- function(plot_dt) {
  # first list: clusters, then QC vars
  g_cls   = plot_umap_clusters(plot_dt) + 
    guides( colour = "none" ) + ggtitle('A')
    # guides( colour = guide_legend(override.aes = list(size = 2))) +
    # theme(
    #   legend.title = element_text(size = 10),
    #   legend.text = element_text(size = 8)
    #   )
  qc_list = unique(plot_dt$qc_var)
  g_qcs   = lapply(qc_list, .plot_one_qc_metric_over_umap, plot_dt)
  g_qcs   = lapply(seq_along(g_qcs), function(i) g_qcs[[i]] + ggtitle(LETTERS[i + 1]))
  g_top   = c(list(g_cls), g_qcs) %>% wrap_plots(ncol = 2) # & 
    # theme(legend.position = "bottom")

  # then violins
  g_bot   = plot_violins_qc_metrics(plot_dt, nrow = 1, titles = c("E", "", ""))

  # join nicely
  g       = g_top / g_bot + plot_layout(heights = c(8, 3))

  return(g)
}
