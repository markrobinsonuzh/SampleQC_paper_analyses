# SampleQC03_utils.R

suppressPackageStartupMessages({
  library('data.table')
  library('forcats')
  library('assertthat')
  library('ggplot2')
  library('RColorBrewer')
  library('magrittr')
  library('patchwork')
  library('scales')
  library('stringr')
  library('viridis')

  library('scater')
  library('flexmix')
})

method_cols = c(SampleQC = '#F13C20', scater = '#4056A1', miQC = '#D79922')

#' Apply scater exclusion to qc_dt
#'
#' @param qc_dt
#' @param nmads
#' @param qc_names if missing, takes default of c('log_counts', 'log_feats', 'mito_prop')
#' @param log, type whether log should first be taken, and whether to apply exclusion to both low and high tails
#' @importFrom data.table ":="
#' @importFrom scater isOutlier
#' @return data.table with scater exclusion status for each cell
#' @export
calc_scater_dt <- function(qc_dt, nmads=2.5, qc_names, log, type) {
  # handle defaults
  if(missing(qc_names)) {
    qc_names  = c('log_counts', 'log_feats', 'logit_mito')
    log     = c(FALSE, FALSE, FALSE)
    type    = c("both", "both", "higher")
  } else {
    stopifnot(length(log)==length(qc_names))
    stopifnot(length(type)==length(qc_names))
  }

  # initialize output
  scater_dt   = qc_dt[, .(sample_id, cell_id)]

  # calculate outliers for each batch, QC
  for (i in seq_along(qc_names)) {
    if (qc_names[[i]]=='mito_prop') {
      vals  = (10^qc_dt$log_mito - 1) / 10^qc_dt$log_counts
    } else {
      vals  = qc_dt[[qc_names[i]]]
    }
    scater_dt[ , (paste0('o_', qc_names[[i]])) := scater::isOutlier(
      vals, nmads = nmads, type = type[i], 
      log = log[i], batch = qc_dt$sample_id)]
  }

  # calculate outliers for each batch, QC
  scater_dt[, outlier := apply(.SD, 1, any), by=.(sample_id, cell_id) ]

  return(scater_dt)
}

calc_miqc_dt <- function(qc_dt, post_cut=0.75) {
  # get requested features
  data_dt = qc_dt[, .(sample_id, cell_id, 
    detected              = round(10^log_feats - 1),
    subsets_mito_percent  = plogis(logit_mito) * 100
    )]
  samples = unique(data_dt$sample_id)

  # fit miQC mixture model
  miqc_dt = lapply(samples, function(s) {
    # fit miQC model to this sample
    dt    = data_dt[sample_id == s]
    model = flexmix(subsets_mito_percent~detected, data = dt, k = 2)

    if (model@k == 1) {
      message('(only one component identified)', appendLF = FALSE)
      dt[, p_bad := 0]
      dt[, outlier := FALSE]

      return(dt)
    }
    
    # find bad component
    intercept1 = parameters(model, component = 1)[1]
    intercept2 = parameters(model, component = 2)[1]
    if (intercept1 > intercept2) {
      compromised_dist = 1
    } else {
      compromised_dist = 2
    }

    post = posterior(model)
    dt[, p_bad    := post[, compromised_dist]]
    dt[, outlier  := p_bad > post_cut]

    return(dt)
    }) %>% rbindlist

  assert_that(
    all(qc_dt$cell_id == miqc_dt$cell_id),
    all(qc_dt$sample_id == miqc_dt$sample_id)
    )

  return(miqc_dt)
}

calc_true_outliers <- function(sims_list) {
  true_outliers   = data.table(
    cell_id   = sims_list$qc_out$cell_id,
    sample_id = sims_list$samples,
    outlier   = sims_list$outliers==1,
    z         = sims_list$z,
    z_f       = factor(sprintf('type%d', sims_list$z))
    )
  return(true_outliers)  
}

calc_scater_pr_dt <- function(qc_dt, qc_names, log, type) {
  # handle defaults
  if(missing(qc_names)) {
    qc_names  = c('log_counts', 'log_feats', 'logit_mito')
    log     = c(FALSE, FALSE, FALSE)
    type    = c("both", "both", "higher")
  } else {
    stopifnot(length(log)==length(qc_names))
    stopifnot(length(type)==length(qc_names))
  }

  # initialize output
  scater_pr_dt  = copy(qc_dt[, c('sample_id', 'cell_id', qc_names), with=FALSE])

  # calc medians, MADs

  # calculate outliers for each batch, QC
  scater_pr_dt  = lapply(seq_along(qc_names),
    function (i) {
      if (qc_names[[i]]=='mito_prop') {
        vals  = (10^qc_dt$log_mito - 1) / 10^qc_dt$log_counts
      } else {
        vals  = qc_dt[[qc_names[i]]]
      }
      if (log[[i]]) {
        vals  = log10(vals)
      }
      if (type[[i]] == 'higher') {
        vals  = vals
      } else if (type[[i]] == 'lower') {
        vals  = -vals
      }
      nn  = qc_names[[i]]
      tmp_dt  = qc_dt[, .(sample_id, cell_id, qc_metric=nn, val=vals)] %>%
        .[ , med_val  := median(val), by=sample_id ] %>%
        .[ , mad_val  := mad(val), by=sample_id ] %>%
        .[ , val_scale  := (val - med_val) / mad_val ]
      if (type[[i]] == 'both')
        tmp_dt[ , val_scale := abs(val_scale) ]
      return(tmp_dt)
  }) %>% rbindlist %>%
    .[ val_scale < 0, val_scale := 0 ] %>%
    .[, .(value = max(val_scale)), by=.(sample_id, cell_id)] %>%
    .[, method := 'scater']

  return(scater_pr_dt)
}

calc_miqc_pr_dt <- function(miqc_dt) {
  miqc_pr_dt = copy(miqc_dt) %>%
    .[, .(sample_id, cell_id, value = p_bad, method = 'miQC')]

  return(miqc_pr_dt)
}

calc_sampleqc_pr_dt <- function(qc_obj) {
  # calculate minimum maha distance for every cell
  sampleqc_pr_dt  = qc_obj$outlier %>% 
    lapply(function(dt) {
      maha_cols   = str_subset(names(dt), 'maha')
      maha_mins   = dt[, maha_cols, with=FALSE] %>%
        as.matrix %>% apply(1, min)
      tmp_dt    = dt[, .(sample_id, cell_id, value=maha_mins)]
    }) %>% rbindlist %>%
    .[, method := 'SampleQC']

  return(sampleqc_pr_dt)
}

calc_pr_by_sample <- function(scater_pr_dt, miqc_pr_dt, sampleqc_pr_dt, true_outliers) {
  pr_by_sample    = rbind(scater_pr_dt, miqc_pr_dt, sampleqc_pr_dt) %>%
    merge(true_outliers, by=c('sample_id', 'cell_id')) %>%
    setorder('method', 'sample_id', 'value') %>%
    .[, .(
      value       = value,
      N_good_reported = cumsum(!outlier),
      N_reported    = seq.int(.N),
      N_good_total  = sum(!outlier)
      ), by=.(method, sample_id)] %>%
    .[, precision := N_good_reported / N_reported ] %>%
    .[, recall    := N_good_reported / N_good_total ] %>%
    .[, method    := factor(method, levels = names(method_cols))]

  return(pr_by_sample)
}

make_default_cuts <- function(qc_obj) {
  alpha         = 0.01
  default_cuts  = data.table(
    method    = c('SampleQC', 'scater', 'miQC'),
    def_value = c(qchisq(1 - alpha, df=metadata(qc_obj)$D), 2.5, 0.75)
    ) %>%
    .[, method    := factor(method, levels = names(method_cols))]

  return(default_cuts)
}

calc_pr_overall <- function(scater_pr_dt, miqc_pr_dt, sampleqc_pr_dt, true_outliers) {
  pr_overall    = rbind(scater_pr_dt, miqc_pr_dt, sampleqc_pr_dt) %>%
    merge(true_outliers, by=c('sample_id', 'cell_id')) %>%
    setorder('method', 'value') %>%
    .[, .(
      value       = value,
      N_good_reported = cumsum(!outlier),
      N_reported    = seq.int(.N),
      N_good_total  = sum(!outlier)
      ), by=.(method)] %>%
    .[, precision   := N_good_reported / N_reported ] %>%
    .[, recall    := N_good_reported / N_good_total ] %>%
    .[, method    := factor(method, levels = names(method_cols))]

  return(pr_overall)
}

calc_pr_by_modality <- function(scater_pr_dt, miqc_pr_dt, sampleqc_pr_dt, true_outliers) {
  # get modalities
  modals_dt   = true_outliers[, .(
    modality = ifelse(max(z)==1, 'unimodal', 'multimodal') %>% 
      factor(levels = c('unimodal', 'multimodal'))
      ), by = sample_id]
  pr_by_modality  = rbind(scater_pr_dt, miqc_pr_dt, sampleqc_pr_dt) %>%
    merge(true_outliers, by=c('sample_id', 'cell_id')) %>%
    merge(modals_dt, by = 'sample_id') %>%
    setorder('method', 'modality', 'value') %>%
    .[, .(
      value       = value,
      N_good_reported = cumsum(!outlier),
      N_reported    = seq.int(.N),
      N_good_total  = sum(!outlier)
      ), by=.(method, modality)] %>%
    .[, precision   := N_good_reported / N_reported ] %>%
    .[, recall    := N_good_reported / N_good_total ] %>%
    .[, method    := factor(method, levels = names(method_cols))]

  return(pr_by_modality)
}

calc_cutpoints <- function(pr_dt, default_cuts, by_what = c('sample_id', 'method', 'modal')) {
  if (by_what == 'sample_id') {
    by_what   = c('method', 'sample_id')
  } else if (by_what == 'overall') {
    by_what   = c('method')
  } else if (by_what == 'modal') {
    by_what   = c('method', 'modality')
  }
  cutpoints   = merge(pr_dt, default_cuts, by='method') %>%
    .[, def_dist  := abs(value - def_value) ] %>%
    .[, .SD[ def_dist == min(def_dist) ], by = by_what]

  return(cutpoints)
}

calc_perf_by_celltype <- function(qc_obj, scater_dt, miqc_dt, true_outliers) {
  # get QC outliers
  qc_outliers   = get_outliers(qc_obj) %>%
    .[, .(cell_id, sample_id, method = 'SampleQC', 
      o_report = ifelse(outlier, 'Rout', 'Rok'))]

  # put all together
  comparison_dt   = rbind(
    qc_outliers,
    scater_dt[, .(cell_id, sample_id, method = 'scater', 
      o_report = ifelse(outlier, 'Rout', 'Rok'))],
    miqc_dt[, .(cell_id, sample_id, method = 'miQC', 
      o_report = ifelse(outlier, 'Rout', 'Rok'))]
    ) %>% merge(
    true_outliers[, .(cell_id, sample_id, z_f, 
      o_true = ifelse(outlier, 'Tout', 'Tok'))],
      by = c('cell_id', 'sample_id'))

  # which samples have more than one celltype?
  props_dt  = comparison_dt[, .(N_z=.N), by=.(sample_id, z_f)] %>%
    .[, prop_z := N_z/sum(N_z), by=sample_id]

  # calculate performance in terms of detection of good cells
  perf_by_celltype  = comparison_dt[, .N, by=.(method, sample_id, z_f, o_true, o_report)] %>%
    dcast( sample_id + method + z_f ~ o_true + o_report, value.var='N', fill=0 ) %>%
    .[, .(
    N_Tok_Rok   = sum(Tok_Rok), 
    N_Tok_Rout  = sum(Tok_Rout), 
    N_Tout_Rok  = sum(Tout_Rok), 
    N_Tout_Rout = sum(Tout_Rout), 
    N_Tok     = sum(Tok_Rok, Tok_Rout),
    N_Rok     = sum(Tok_Rok, Tout_Rok),
    precision   = Tok_Rok / (Tok_Rok + Tout_Rok),
    recall    = Tok_Rok / (Tok_Rok + Tok_Rout),
    FOR     = Tok_Rout  / (Tok_Rout +  Tout_Rout),
    FPR     = Tout_Rok  / (Tout_Rok  + Tout_Rout)
    ), by=.(sample_id, method, z_f)]
  perf_by_celltype   = merge(perf_by_celltype, props_dt, by=c('sample_id', 'z_f')) %>%
    .[, prop_is_1   := ifelse(prop_z == 1, 'one celltype', 'multiple types') ] %>%
    .[, prop_is_1   := fct_relevel(prop_is_1, 'one celltype') ] %>%
    .[, method      := factor(method, levels = names(method_cols))]

  return(perf_by_celltype)
}

plot_maha_dists <- function(qc_obj, s) {
  # get fitted params and data
  mu_0    = qc_obj$mu_0[[s]]
  alpha_j   = qc_obj$alpha_j[[s]]
  beta_ks   = qc_obj$beta_k[[s]]
  sigma_ks  = qc_obj$sigma_k[[s]]
  K     = colData(qc_obj)[s, ]$K
  x     = qc_obj$qc_metrics[[s]]

  # calculate mahalanobis distances for each cell from each cluster
  mu_j    = mu_0 + alpha_j
  n     = nrow(x)
  mahas   = vapply(seq_len(K),
    function(k) mahalanobis(x, mu_j + beta_ks[k,], sigma_ks[,,k]),
    numeric(n)
    )
  maha_min  = apply(mahas, 1, min)
  maha_max  = 20
  maha_dt   = data.table(cell_id = qc_obj$cell_id[[s]], maha_min = maha_min, x)
  maha_dt[, maha_trunc := pmin(maha_min, maha_max)]

  # unpack, arrange
  qc_names  = metadata(qc_obj)$qc_names
  qc_1      = qc_names[[1]]
  qc_not_1  = qc_names[-1]

  # join together
  plot_dt   = maha_dt %>%
  melt(
    id  = c('cell_id', 'maha_min', 'maha_trunc', qc_1), 
    measure = qc_not_1,
    variable.name='qc_x_name', value.name='qc_x'
    ) %>%
  setnames(qc_1, 'qc_y') %>%
  setorder('maha_min')

  # plot
  g = ggplot(plot_dt) + 
  aes(y = qc_y, x = qc_x, colour = maha_trunc) +
  geom_point(size = 1) + scale_colour_viridis(limits = c(0, 20)) +
  scale_x_continuous(breaks = pretty_breaks()) + 
  scale_y_continuous(breaks = pretty_breaks()) +
  facet_grid(. ~ qc_x_name, scales = 'free_x') +
  theme_bw() + 
  labs(y = qc_names[[1]], x = 'other QC variable', 
    colour = 'Maha.\ndist. to\nclosest\ncomponent')

  return(g)
}

plot_sampleqc_vs_others_over_biaxials <- function(qc_obj, true_outliers, scater_dt, miqc_dt, s) {
  g_true  = plot_outliers(qc_obj, s, true_outliers) + labs( title = 'A. True outliers' )
  g_fit   = plot_fit_over_biaxials(qc_obj, s) + labs( title = 'B. SampleQC fits' )
  g_maha  = plot_maha_dists(qc_obj, s) + labs( title = 'C. Mahalanobis distances' )
  g_qc    = plot_outliers(qc_obj, s) + labs( title = 'D. SampleQC outliers' )
  g_scat  = plot_outliers(qc_obj, s, scater_dt) + labs( title = 'E. scater outliers' )
  g_miqc  = plot_outliers(qc_obj, s, miqc_dt) + labs( title = 'F. miQC outliers' )
  g       = list(g_true, g_fit, g_maha, g_qc, g_scat, g_miqc) %>% 
    wrap_plots(ncol = 2, byrow = FALSE)

  return(g)
}

calc_mahas_dt <- function(qc_obj) {
  n_dims    = metadata(qc_obj)$D
  alpha     = metadata(qc_obj)$fit_params$alpha
  mahas_dt  = colData(qc_obj)$outlier %>% 
    lapply(function(dt) 
      dt[, .(min_maha = apply(.SD, 1, min), n_comps = length(.SD)), 
        by= .(cell_id, sample_id, outlier)]) %>% 
    rbindlist %>%
    setorder('sample_id', 'min_maha') %>%
    .[, q         := seq.int(.N) / .N, by = sample_id] %>%
    .[, chi_exp   := dchisq(q, df = n_dims) ] %>%
    .[, maha_cut  := qchisq(1 - alpha, df = n_dims) ]
}

plot_maha_qqs <- function(mahas_dt, sel_samples, alpha = 0.01) {
  g = ggplot(mahas_dt[ sample_id %in% sel_samples ]) +
    aes( x = sqrt(min_maha), y = chi_exp ) +
    geom_point( size = 1, alpha = 0.5 ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    scale_x_continuous( breaks = pretty_breaks() ) +
    coord_cartesian( xlim = c(0, 4) ) +
    facet_wrap( ~ sample_id, scales = 'free_x' ) +
    theme_bw() +
    labs(
      x = 'minimum Mahalanobis distance', 
      y = 'expected Chi-squared distance'
      )
}

perf_plot_fn <- function(perf_list) {
  lapply(perf_list,
    function(p) {
      dt_ok   = perf_dt[perf_metric == p] %>%
         setorder(-'min_p_jk')
      g =  ggplot(dt_ok) + 
        aes( x=scater, y=SampleQC, size=N_Tout, fill=min_p_jk ) + 
        geom_point( shape=21, colour='black' ) +
        geom_abline( intercept=0, slope=1, colour='grey' ) +
        expand_limits( x=0, y=0 ) +
        scale_x_continuous( breaks=pretty_breaks() ) +
        scale_y_continuous( breaks=pretty_breaks() ) +
        # coord_fixed( ratio=1 ) +
        theme_bw() + labs( title=p )
      if (p == 'FPR') {
        g = g + coord_cartesian( xlim=c(0,0.2), ylim=c(0,0.2)) +
          scale_fill_distiller( palette='Reds', breaks=pretty_breaks() ) + 
          scale_size( guide=FALSE ) +
          labs(
            fill  = 'minimum\ncelltype\npropn.'
            )
      } else if (p == 'FOR') {
        g = g + coord_cartesian( xlim=c(0,0.2), ylim=c(0,0.2)) +
          scale_fill_distiller( palette='Reds', breaks=pretty_breaks(), guide=FALSE ) + 
          scale_size(  ) +
          labs(
            # size  = '# cells\nin sample'
            size  = '# true\noutliers\nin sample'
            )
      } else {
        g = g + scale_fill_distiller( palette='Reds', breaks=pretty_breaks() ) + 
          scale_size( trans='log10' ) +
          guides( fill=FALSE, size=FALSE )
      }
      return(g)
    })
}

plot_pr_overall <- function(pr_by_sample, pr_overall, cuts_overall, sel_samples) {
  # plot all PR curves
  g = ggplot() +
    aes( x=recall, y=precision, colour=method ) +
    geom_path(
      data  = pr_by_sample, 
      aes(group=interaction(sample_id, method)),
      size=1, alpha=0.2
      ) +
    geom_path(
      data  = pr_overall,
      size=2, alpha=1
      ) +
    geom_point(
      data  = cuts_overall, 
      size=4, alpha=1 ) +
    scale_colour_manual( values=method_cols ) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    coord_cartesian( ylim=c(0.8, 1) ) +
    theme_bw() + theme( legend.position='bottom' )

  return(g)
}

plot_pr_by_modality <- function(pr_by_modality, cuts_by_modal) {
  # plot all PR curves
  g = ggplot(pr_by_modality) +
    aes(x = recall, y = precision, colour = method, linetype = modality) +
    geom_path(size=1, alpha=1) +
    geom_point(data  = cuts_by_modal, size=4, alpha=1) +
    scale_colour_manual( values=method_cols ) +
    scale_linetype_manual(values = c(unimodal = 'dashed', multimodal = 'solid')) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    coord_cartesian( ylim=c(0.8, 1) ) +
    theme_bw() + theme( legend.position='bottom', 
      legend.key.width = unit(0.5,"in") )

  return(g)
}

plot_pr_by_sample <- function(pr_by_sample, cuts_by_sample, sel_samples) { 
  g = ggplot(pr_by_sample[ sample_id %in% sel_samples ]) +
    aes( x=recall, y=precision, colour=method ) +
    geom_path() +
    geom_point( data=cuts_by_sample[ sample_id %in% sel_samples ], size=2 ) +
    # scale_colour_brewer( palette='Set1' ) +
    scale_colour_manual( values=method_cols ) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    coord_cartesian( ylim=c(0.8, 1) ) +
    facet_wrap( ~ sample_id, nrow=4 ) +
    theme_bw()
  return(g)
}

plot_bars_of_lost_cells_by_type <- function(perf_by_celltype, what = c('abs', 'prop')) {
  n_cols    = length(unique(perf_by_celltype$z_f))
  bar_cols  = brewer.pal(n_cols + 1, 'Greys')[-1]
  what      = match.arg(what)
  plot_dt   = copy(perf_by_celltype) %>% 
    .[, delta_ok  := N_Tok - N_Tok_Rok] %>%
    .[, delta_abs := abs(delta_ok)]
  y_var     = c(abs = 'delta_ok', prop = 'delta_abs')[[what]]
  g = ggplot(plot_dt) +
    aes_string(x = 'sample_id', y = y_var, fill = 'z_f' )
  if (what == 'abs') {
    g = g + 
      geom_col(colour = 'black') +
      scale_y_continuous( breaks=seq(0, 2000, by = 200) )
      expand_limits(y = 600)
    y_lab   = 'good cells falsely reported as outliers'
  } else if (what == 'prop') {
    g = g + 
      geom_col(colour = 'black', position = 'fill') +
      scale_y_continuous(breaks= pretty_breaks())    
    y_lab   = 'proportion of misreported good cells'
  }
  g = g +
    scale_fill_manual(values = bar_cols) +
    facet_grid( method ~ prop_is_1, scales='free', space='free' ) +
    # coord_flip() +
    theme_bw() + theme(
      # axis.text.x   = element_text(angle=90, vjust=0.5, size=6),
      axis.text.x   = element_blank(),
      legend.position = 'bottom'
      ) +
    labs(x = 'sample', y = y_lab, fill = 'celltype')

  return(g)
}

.not_used_calcs <- function() {
  # get outliers
  qc_outliers   = get_outliers(qc_out_obj) %>%
    .[, .(cell_id, sample_id, method='SampleQC', o_report=ifelse(outlier, 'Rout', 'Rok'))]

  # join w scater
  outliers_dt   = rbind(
    qc_outliers,
    scater_out[, .(cell_id, sample_id, method='scater', o_report=ifelse(outlier, 'Rout', 'Rok'))]
    )

  # add min p_jk proportion for each sample (?)
  samples_by_g  = split(sims_list$samples, sims_list$groups) %>%
    lapply(unique) %>% 
    lapply(sort)
  p_jks_dt    = names(sims_list$group_sims) %>%
    lapply(function(n) data.table(
      sample_id   = samples_by_g[[n]],
      # N       = table(sims_list$samples[sims_list$groups == n ]) %>% as.vector,
      min_p_jk  = apply(sims_list$group_sims[[n]]$p_jk, 1, min)
      )) %>% rbindlist

  # horrible
  comparison_dt   = merge(
    true_outliers[, .(cell_id, sample_id, z_f, o_true=ifelse(outlier, 'Tout', 'Tok'))], 
    outliers_dt,
    by=c('cell_id', 'sample_id')
    ) 

  # calculate precision, recall etc
  pr_dt   = comparison_dt %>%
    .[, .N, by=.(sample_id, o_true, method, o_report)] %>%
    dcast( sample_id + method ~ o_true + o_report, value.var='N', fill=0 ) %>%
    .[, .(
    N_Tout    = sum(Tout_Rout, Tout_Rok),
    N_Rout    = sum(Tout_Rout, Tok_Rout),
    precision   = Tout_Rout / (Tout_Rout + Tok_Rout),
    FOR     = Tout_Rok  / (Tout_Rok +  Tok_Rok),
    recall    = Tout_Rout / (Tout_Rout + Tout_Rok),
    FPR     = Tok_Rout  / (Tok_Rout  + Tok_Rok)
    ), by=.(sample_id, method)] %>%
    merge(p_jks_dt, by='sample_id')

  perf_dt = melt(pr_dt, id=c('sample_id', 'method', 'N_Tout', 'N_Rout', 'min_p_jk'), 
    variable.name='perf_metric') %>%
    dcast( sample_id + perf_metric + N_Tout + min_p_jk ~ method, value.name='value')

  # calculate sample-level performance
  sample_perf   = comparison_dt[, .N, by=.(method, sample_id, z_f, o_true, o_report)] %>%
    dcast( sample_id + method + z_f ~ o_true + o_report, value.var='N', fill=0 ) %>%
    .[, .(
    N_Tok_Rok   = sum(Tok_Rok), 
    N_Tok_Rout  = sum(Tok_Rout), 
    N_Tout_Rok  = sum(Tout_Rok), 
    N_Tout_Rout = sum(Tout_Rout), 
    N_Tout    = sum(Tout_Rout, Tout_Rok),
    N_Rout    = sum(Tout_Rout, Tok_Rout),
    precision   = Tout_Rout / (Tout_Rout + Tok_Rout),
    FOR     = Tout_Rok  / (Tout_Rok +  Tok_Rok),
    recall    = Tout_Rout / (Tout_Rout + Tout_Rok),
    FPR     = Tok_Rout  / (Tok_Rout  + Tok_Rok)
    ), by=.(sample_id, method, z_f)]
  props_dt  = comparison_dt[, .(N_z=.N), by=.(sample_id, z_f)] %>%
    .[, prop_z := N_z/sum(N_z), by=sample_id]
  sample_perf   = merge(sample_perf, props_dt, by=c('sample_id', 'z_f')) %>%
    .[, prop_is_1   := ifelse(prop_z == 1, 'one celltype', 'multiple types') ] %>%
    .[, prop_is_1   := fct_relevel(prop_is_1, 'one celltype') ]

  # try sample-level performance with task reversed (i.e. trying to detect good cells, not outliers)
  sample_perf2  = comparison_dt[, .N, by=.(method, sample_id, z_f, o_true, o_report)] %>%
    dcast( sample_id + method + z_f ~ o_true + o_report, value.var='N', fill=0 ) %>%
    .[, .(
    N_Tok_Rok   = sum(Tok_Rok), 
    N_Tok_Rout  = sum(Tok_Rout), 
    N_Tout_Rok  = sum(Tout_Rok), 
    N_Tout_Rout = sum(Tout_Rout), 
    N_Tok     = sum(Tok_Rok, Tok_Rout),
    N_Rok     = sum(Tok_Rok, Tout_Rok),
    precision   = Tok_Rok / (Tok_Rok + Tout_Rok),
    recall    = Tok_Rok / (Tok_Rok + Tok_Rout),
    FOR     = Tok_Rout  / (Tok_Rout +  Tout_Rout),
    FPR     = Tout_Rok  / (Tout_Rok  + Tout_Rout)
    ), by=.(sample_id, method, z_f)]
  sample_perf2   = merge(sample_perf2, props_dt, by=c('sample_id', 'z_f')) %>%
    .[, prop_is_1   := ifelse(prop_z == 1, 'one celltype', 'multiple types') ] %>%
    .[, prop_is_1   := fct_relevel(prop_is_1, 'one celltype') ]
}

.not_used_plots <- function() {
  # ```{r plot_outlier_comparisons}
  perf_list   = c('precision', 'recall')
  g_pr1 = wrap_plots(perf_plot_fn(perf_list), nrow=1)

  g_pr2 = ggplot(pr_dt) +
    aes(y = precision, x = recall, size = N_Tout, fill = min_p_jk) +
    geom_point( shape=21, colour='black' ) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    scale_fill_distiller( palette='Reds', breaks=pretty_breaks() ) + 
    # scale_size( trans='log10' ) +
    facet_wrap( ~ method ) +
    theme_bw() + 
    labs(
      size  = '# true\noutliers\nin sample',
      fill  = 'minimum\ncelltype\npropn.'
      )
  g = g_pr1 / g_pr2
  for (f in formats)
    ggsave(sprintf(comparisons_pr, f, f), g, h=6, w=7)

  perf_list   = c('FOR', 'FPR')
  g = wrap_plots(perf_plot_fn(perf_list), nrow=1)
  for (f in formats)
    ggsave(sprintf(comparisons_f2, f, f), g, h=3, w=8)
  # ```

  # ```{r plot_outliers_by_celltype_1, fig.height=6, fig.width=8}
  perf_melt   = sample_perf %>%
    melt( id=c('sample_id', 'z_f', 'prop_z', 'prop_is_1', 'N_z', 'method', 'N_Tout', 'N_Rout'),
      variable.name='metric', value.name='value')
  plot_dt   = perf_melt %>%
    dcast( sample_id + z_f + prop_z + prop_is_1 + metric ~ method, value.var='value' ) %>%
    setorder(-'prop_z')

  g = ggplot(plot_dt[ metric %in% c('precision', 'recall') ]) +
    aes( y=SampleQC, x=scater, colour=z_f, size=prop_z ) +
    geom_abline( intercept=0, slope=1, colour='black', linetype='dashed' ) +
    geom_point( alpha=0.5 ) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    scale_colour_brewer( palette='Set1' ) +
    # scale_size( limits=c(0,1), breaks=seq(0,1,0.2), range=c(0,6) ) +
    scale_size( limits=c(0,1), breaks=seq(0,1,0.2), range=c(0,4) ) +
    coord_fixed(ratio=1, xlim=c(0,1), ylim=c(0,1)) +
    facet_grid( prop_is_1 ~ metric ) +
    theme_bw() + theme(aspect.ratio=0.9) +
    labs(
      title     = 'Performance metrics by QC celltype', 
      subtitle  = 'task: identifying outliers',
      colour='QC celltype', size='celltype\npropn.'
      )
  print(g)
  for (f in formats)
    ggsave(sprintf(per_celltype_fpr, f, f), g, h=5, w=10)
  # ```

  # ```{r plot_outliers_by_celltype_1B, fig.height=6, fig.width=8}
  perf_melt   = sample_perf2 %>%
    melt( id=c('sample_id', 'z_f', 'prop_z', 'prop_is_1', 'N_z', 'method', 'N_Tok', 'N_Rok'),
      variable.name='metric', value.name='value')
  plot_dt   = perf_melt %>%
    dcast( sample_id + z_f + prop_z + prop_is_1 + metric ~ method, value.var='value' ) %>%
    setorder(-'prop_z')

  g = ggplot(plot_dt[ metric %in% c('precision', 'recall') ]) +
    aes( y=SampleQC, x=scater, colour=z_f, size=prop_z ) +
    geom_abline( intercept=0, slope=1, colour='black', linetype='dashed' ) +
    geom_point( alpha=0.5 ) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    scale_colour_brewer( palette='Set1' ) +
    # scale_size( limits=c(0,1), breaks=seq(0,1,0.2), range=c(0,6) ) +
    scale_size( limits=c(0,1), breaks=seq(0,1,0.2), range=c(0,4) ) +
    coord_cartesian(xlim=c(0,1), ylim=c(0.5,1)) +
    facet_grid( prop_is_1 ~ metric ) +
    theme_bw() + theme(aspect.ratio=0.9) +
    labs(
      title     = 'Performance metrics by QC celltype', 
      subtitle  = 'task: identifying good cells',
      colour='QC celltype', size='celltype\npropn.'
      )
  print(g)
  for (f in formats)
    ggsave(sprintf(per_celltype_fpr_rev, f, f), g, h=5, w=10)
  # ```

  # ```{r plot_outliers_by_celltype_2, fig.height=6, fig.width=7}
  g = ggplot(sample_perf2) +
    aes( y=precision, x=recall, colour=z_f, size=prop_z ) +
    geom_point( alpha=0.5 ) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    scale_colour_brewer( palette='Set1' ) +
    # scale_size( limits=c(0,1), breaks=seq(0,1,0.2), range=c(0,6) ) +
    scale_size( limits=c(0,1), breaks=seq(0,1,0.2), range=c(0,4) ) +
    facet_grid( prop_is_1 ~ method ) +
    theme_bw() + 
    labs( title='Precision and recall by celltype', 
      colour='QC celltype', size='celltype\npropn.' )
  print(g)
  for (f in formats)
    ggsave(sprintf(per_celltype_pr, f, f), g, h=6, w=7)
  # ```

  # ```{r plot_outliers_by_celltype_3, fig.height=6, fig.width=7, eval=FALSE}
  plot_dt   = sample_perf2 %>% 
    dcast( sample_id + z_f + N_z + prop_z + prop_is_1 ~ method, value.var=c('precision', 'recall') )
  g = ggplot(plot_dt) +
    aes(
      y=precision_scater, x=recall_scater, 
      yend=precision_SampleQC, xend=recall_SampleQC, 
      alpha=log10(N_z) ) +
    geom_segment( arrow=arrow( length=unit(0.03, "npc") ) ) +
    scale_x_continuous( breaks=pretty_breaks() ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    scale_colour_brewer( palette='Set1' ) +
    scale_alpha( range=c(0.1, 0.9) ) +
    # scale_size( limits=c(0,1), breaks=seq(0,1,0.2), range=c(0,6) ) +
    # scale_size( limits=c(0,1), breaks=seq(0,1,0.2), range=c(0,4) ) +
    facet_grid( . ~ prop_is_1 ) +
    theme_bw() + theme( aspect.ratio=1 ) +
    labs( title='Effect of changing from scater to SampleQC', 
      x='recall', y='precision',
      size='celltype\npropn.' )
  print(g)
  for (f in formats)
    ggsave(sprintf(per_celltype_effect, f, f), h=4, w=8)
  # ```

  # ```{r plot_outliers_by_celltype_5, fig.height=8, fig.width=8}
  plot_dt   = copy(sample_perf) %>% 
    .[, delta_out := N_Tout_Rout - N_Tout ]
  g = ggplot(plot_dt) +
    aes(x=fct_rev(sample_id), y=delta_out, fill=z_f ) +
    geom_col( position=position_dodge2(preserve='single') ) +
    scale_y_continuous( breaks=pretty_breaks() ) +
    scale_fill_brewer( palette='Set1' ) +
    facet_grid( prop_is_1 ~ method, scales='free_y', space='free' ) +
    coord_flip() +
    theme_bw() +
    labs( title='Bias in outliers reported', 
      x=NULL, y='no. reported - no. true outliers',
      fill='celltype' )
  print(g)
  for (f in formats)
    ggsave(sprintf(per_celltype_bias_out, f, f), g, h=8, w=8)
  # ```
}
