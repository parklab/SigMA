devtools::load_all()

# this script demonstrates different types of NNLS calculations one can do with SigMA
# using simulated data, which allows comparison to the calculation in WGS data
#
# you can do the same calculations with actual panel data you would only miss the part
# where calculations in panel are compared to WGS
#
# please run the following script after test_tune_example_forPfizer.R script run use the following
# which creates the simulated data used below

tumor_type <- 'breast'
catalog_name <- 'cosmic_v3p2_inhouse'
check_msi <- F # set T if you might have msi samples
data <- 'fo'
maf_percent <- 0.001 # replace with NULL for matched normal panels
genome_file <- paste0('matrix_96dim_', data, '_', tumor_type, '_', catalog_name, '_nomsisimulation.csv')
df_genome <- read.csv(genome_file)

## use the lines below if you have not first tried the test_tune_example_forPfizer.R script
outfile <- run(genome_file = genome_file, 
               tumor_type = tumor_type, data = 'msk', check_msi = F,
               catalog_name = catalog_name)
	       
outfile <- run(genome_file,
                      data = paste0('gbm_for_', tumor_type,'_', data, '_', catalog_name, '_maf_percent_', maf_percent),
                      tumor_type = tumor_type,
                      catalog_name = catalog_name,
                      do_mva = T, # set to T for msk
                      do_assign = T,
                      check_msi = check_msi,
                      custom = T,
                      norm96 = weight_3Nfreq$fo)
print('ran')		      
df_sigma_out <- read.csv(outfile)
df_sigma_out <- llh_max_characteristics(df_sigma_out,tumor_type = tumor_type, catalog_name = catalog_name)
df_sigma_lite <- lite_df(df_sigma_out)

df_exps <- get_sig_exps(df_sigma_out, col_sigs = 'sigs_all', col_exps = 'exps_all')
df_exps_cluster <- get_sig_exps(df_sigma_out, col_sigs='cluster_sigs_all', col_exps = 'cluster_exps_all')
colnames(df_exps_cluster) <- paste0('cluster_', colnames(df_exps_cluster))

df_wgs <- get_wgs_data_for_llh(tumor_type = tumor_type, check_msi = check_msi)
df_wgs <- weight_3Nfreq[data] * df_wgs
df_wgs_norm <- t(t(df_wgs)/colSums(df_wgs))
colnames(df_wgs_norm) <- colnames(df_wgs)
df_wgs <- df_wgs_norm

df_likelihoods <- match_to_catalog(genomes = df_genome,
                                   signatures = df_wgs,
               		           method = 'median_catalog',
				   data = data)


df <- llh_max_characteristics_wgs_data(df_likelihoods,
                                       tumor_type = tumor_type,
				       catalog_name = catalog_name,
				       check_msi = check_msi)
df <- df[,colnames(df) %in% c('sig_max_ml','max_ml','is_sig3_sig_max','Signature_3_ml_wgs_data', 'max_ml_sample_sigs_all', 'max_ml_sample_exps_all', 'is_msi_sig_max', 'is_pole_sig_max') | grepl('ave_ml_', colnames(df))]

df_exps_sample_max <- get_sig_exps(df, col_exps = 'max_ml_sample_exps_all', col_sigs = 'max_ml_sample_sigs_all')
colnames(df_exps_sample_max) <- paste0('max_ml_sample_', colnames(df_exps_sample_max))
inds <- match(df_genome$tumor, df_sigma_lite$tumor)
df_out <- cbind(df_genome,
                df_sigma_lite[inds,!(colnames(df_sigma_lite) %in% c('tumor', 'exp_sig3','tumor'))],
		df_exps[inds,],
		df_exps_cluster[inds,],
		df[inds,],
		df_exps_sample_max[inds,])

df_out$rat_sig3 <- df_out$exp_sig3/df_out$total_snvs

# up to here the same script can be used for data the lines following this are comparisons to WGS calculation

# the is_sig3 tag comes from the WGS analysis while the rest of the exposures are calculated from panel simulations
#
#   rat_sig3 is exposure of sig3 calculated with NNLS directly (exp_sig3) divided by total_snvs
#
#   cluster_exp_sig3 is the relative exposure of sig3 in the most likely WGS cluster to which
#   to which the sample is matched to
#
#   max_ml_sample_exp_sig3 is the relative exposure of sig3 in the WGS sample which has the

df_sig3_exp_long <- reshape2::melt(df_out[,c('rat_sig3', 'cluster_exp_sig3', 'max_ml_sample_exp_sig3','is_sig3', 'Signature_3_mva')], id.var = 'is_sig3')
plot <- ggplot2::ggplot(df_sig3_exp_long, ggplot2::aes(x = is_sig3, y = value))
plot <- plot + ggplot2::geom_jitter(size = 0.5, alpha = 0.7, ggplot2::aes(color = is_sig3))
plot <- plot + ggplot2::geom_boxplot(outlier.color = NA, lwd = 0.2, alpha = 0.4, ggplot2::aes(fill = is_sig3))
plot <- plot + ggplot2::facet_wrap(~variable, nrow = 1) 
plot <- plot + ggplot2::theme_bw()
ggplot2::ggsave(plot, file = 'compare_nnls_exposures_calculated_with_different_model.pdf', width =4, height = 3)
write.table(df_out, 'example_exposures_from_WGS.txt', row.names = F, sep = '\t', quote = F)

# The same calculation is done for other signatures as well replace sig3 with any other signature
# e.g. exp_sig17b, cluster_exp_sig17b, max_ml_sample_exp_sig17b,  ave_ml_exp_sig17b, but is_sig3 column is only
# present for sig3 in the built in matrices so the below figure cannot be made. If you want to do the figure below
# for other signatures see the calculation below

df_out <- read.delim('example_exposures_from_WGS.txt')
signatures <- signames_per_tissue_per_catalog[[catalog_name]][[tumor_type]]
if(check_msi)
  signatures <- unique(c(signatures,
                  signames_per_tissue_per_catalog[[catalog_name]][['msi']],
                  signames_per_tissue_per_catalog[[catalog_name]][['msi_extra']],
                  signames_per_tissue_per_catalog[[catalog_name]][['pole']]))

cosmic_catalog <- catalogs[[catalog_name]]

df_wgs <- get_wgs_data_for_llh(tumor_type = tumor_type, check_msi = check_msi)
df_wgs_is_sig3 <- get_Sig3_llh_from_wgs_data(colnames(df_wgs), tumor_type = tumor_type, check_msi = check_msi)

df_exps_wgs <- data.frame(matrix(0, dim(df_wgs)[[2]], length(signatures)))
colnames(df_exps_wgs) <- gsub(signatures, pattern = 'Signature_', replace = 'exp_sig')
for(i in 1:dim(df_wgs)[[2]]){
  signatures_this <- signatures
  spec <- df_wgs[,i]
  if(!df_wgs_is_sig3$is_sig3[df_wgs_is_sig3$tumor == colnames(df_wgs)[[i]]])
    signatures_this <- signatures_this[signatures_this != "Signature_3"]
  if(check_msi){
    if(!df_wgs_is_sig3$is_msi[df_wgs_is_sig3$tumor == colnames(df_wgs)[[i]]])
      signatures_this <- signatures_this[!(signatures_this %in% c(signames_per_tissue_per_catalog[[catalog_name]][['msi']], signames_per_tissue_per_catalog[[catalog_name]][['msi_extra']]))]
    if(!df_wgs_is_sig3$is_pole[df_wgs_is_sig3$tumor == colnames(df_wgs)[[i]]])
      signatures_this <- signatures_this[!(signatures_this %in% c(signames_per_tissue_per_catalog[[catalog_name]][['msi']], signames_per_tissue_per_catalog[[catalog_name]][['msi_extra']], signames_per_tissue_per_catalog[[catalog_name]][['pole']]))]
  }
  
  df_decompose <- match_to_catalog(data.frame(t(df_wgs))[i,], cosmic_catalog[,signatures_this], data = 'wgs', method = 'decompose')
  df_decompose <- cbind(df_decompose, get_sig_exps(df_decompose, col_sigs='sigs_all', col_exps = 'exps_all'))
  for(sig in signatures_this){
    df_decompose[which(df_decompose[,paste0(sig, '_l_rat')] < 0.505), gsub(sig, pattern = 'Signature_', replace = 'exp_sig')]  <- 0
  }
  paste0(gsub(signatures_this, pattern = 'Signature_', replace = 'exp_sig'))
  df_exps_wgs[i, gsub(signatures_this, pattern = 'Signature_', replace = 'exp_sig')] <- df_decompose[,gsub(signatures_this, pattern = 'Signature_', replace = 'exp_sig')]
}
colnames(df_exps_wgs) <- paste0('wgs_', gsub(colnames(df_exps_wgs), pattern = 'Signature_', replace = ''))
df_out <- cbind(df_out, df_exps_wgs[match(df_out$tumor, colnames(df_wgs)),])

# an example comparing wgs exposure in the same sample wgs_exp_sig13 to the one calculated in panel
# note that for clustered mutational processes such as APOBEC (e.g. Sig2, Sig13) a subset of genomic
# regions can have higher exposure than the whole genome if a cluster of mutations are captured
sig_selected <- 'exp_sig17b' # change this to see other signatures of interest

df_out$exp_sig_selected <- df_out[,sig_selected]
df_out$wgs_exp_sig_selected <- df_out[,paste0('wgs_', sig_selected)]
df_out$cluster_exp_sig_selected <- df_out[,paste0('cluster_', sig_selected)]

plot_nnls_panel <- ggplot2::ggplot(df_out, ggplot2::aes(x = wgs_exp_sig_selected, y = exp_sig_selected/total_snvs))
plot_nnls_panel <- plot_nnls_panel + ggplot2::geom_point(size = 0.5, alpha = 0.5)
plot_nnls_panel <- plot_nnls_panel + ggplot2::geom_abline(slope = 1)
plot_nnls_panel <- plot_nnls_panel + ggplot2::geom_smooth(method = 'lm')
plot_nnls_panel <- plot_nnls_panel + ggplot2::theme_bw()
ggplot2::ggsave(plot_nnls_panel, file = paste0(sig_selected, '_scatter_nnls_panel_vs_WGS.pdf'), width = 2, height = 2)

# you can try cleaning out assignments with potential overfitting problem
df_out$selected_l_rat <- df_sigma_out[match(df_out$tumor, df_sigma_out$tumor), paste0(gsub(sig_selected, pattern = 'exp_sig', replace = 'Signature_'), '_l_rat')] 
df_out$exp_sig_selected_cleaned <- df_out$exp_sig_selected
df_out$exp_sig_selected_cleaned[df_out$selected_l_rat < 0.6] <- 0

plot_nnls_panel_clean <- ggplot2::ggplot(df_out, ggplot2::aes(x = wgs_exp_sig_selected, y = exp_sig_selected_cleaned/total_snvs))
plot_nnls_panel_clean <- plot_nnls_panel_clean + ggplot2::geom_point(size = 0.5 , alpha = 0.5)
plot_nnls_panel_clean <- plot_nnls_panel_clean + ggplot2::geom_abline(slope = 1)
plot_nnls_panel_clean <- plot_nnls_panel_clean + ggplot2::geom_smooth(method = 'lm')
plot_nnls_panel_clean <- plot_nnls_panel_clean + ggplot2::theme_bw()
ggplot2::ggsave(plot_nnls_panel_clean, file = paste0(sig_selected, '_scatter_nnls_panel_cleaned_vs_WGS.pdf'), width = 2, height = 2)

# you can also compare to the matched cluster
plot_nnls_cluster <- ggplot2::ggplot(df_out, ggplot2::aes(x = wgs_exp_sig_selected, y = cluster_exp_sig_selected))
plot_nnls_cluster <- plot_nnls_cluster + ggplot2::geom_point(size =0.5, alpha = 0.5)
plot_nnls_cluster <- plot_nnls_cluster + ggplot2::geom_abline(slope = 1)
plot_nnls_cluster <- plot_nnls_cluster + ggplot2::geom_smooth(method = 'lm')
plot_nnls_cluster <- plot_nnls_cluster + ggplot2::theme_bw()
ggplot2::ggsave(plot_nnls_cluster, file = paste0(sig_selected, '_scatter_nnls_matched_cluster_WGS_vs_WGS.pdf'), width = 2, height = 2)

# For the signatures you are interested please compare all the approaches above to decide which gives more accurate estimates from panel