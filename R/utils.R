# gets the clusters with maximum likelihood and decomposes them with NNLS to determine relative exposures of signatures
llh_max_characteristics <- function(df, tumor_type, cosmic_version){
  # calculate the signature exposures of the cluster with maximum likelihood

  if(cosmic_version == "v3")
    cosmic_catalog <- pcawg_catalog
  if(cosmic_version == "v2")
    cosmic_catalog <- cosmic_catalog_v2

  if('sig_max_ml_msi' %in% colnames(df))
    df$max_clust <- df$sig_max_ml_msi
  else if('sig_max_ml' %in% colnames(df))
    df$max_clust <- df$sig_max_ml

  signatures <- cosmic_catalog[,signames_per_tissue[[tumor_type]]]
  if('sig_max_ml_msi' %in% colnames(df))
    signatures_msi_pole <- cbind(signatures, cosmic_catalog[,c(signames_per_tissue[['msi']], signames_per_tissue[['pole']])])
  df$cluster_sigs_all <- '' 
  df$cluster_exps_all  <- ''

  for(max_clust in unique(df$max_clust)){
    if(max_clust %in% colnames(all_catalogs[['msi']])){
      median_clust <- all_catalogs[['msi']][,max_clust]
      signatures_this <- signatures_msi_pole
    }
    else if(max_clust %in% colnames(all_catalogs[['pole']])){
      median_clust <- all_catalogs[['pole']][,max_clust]
      signatures_this <- signatures_msi_pole
    }
    else{
      median_clust <- all_catalogs[[tumor_type]][,max_clust]
      signatures_this <- signatures
    }
    exps <- coefficients(nnls::nnls(as.matrix(signatures_this), median_clust))
  
    df$cluster_sigs_all[df$max_clust == max_clust] <- paste0(colnames(signatures_this)[exps > 0.1], collapse = '.')
    df$cluster_exps_all[df$max_clust == max_clust] <- paste0(exps[exps > 0.1], collapse = '_')
  }
  return(df)
}

# converts exps_all, sigs_all or exps_all_msi, sigs_all_msi columns into an exposure table
get_sig_exps <- function(df, col_sigs, col_exps){
  # returns a data table splitting the exposures and signatures stored in compact columns e.g. sigs_all/exps_all or sigs_all_msi/exps_all_msi

  signames <- unique(unlist(unlist(strsplit(as.character(df[,col_sigs]), split = '\\.'))))
  df_sigs <- data.frame(matrix(0, dim(df)[[1]], length(signames)))
  colnames(df_sigs) <- gsub(signames, pattern = 'Signature_', replace ='exp_sig')

  for(i in 1:dim(df)[[1]]){
    sigs_this <- as.character(unlist(strsplit(as.character(df[i,col_sigs]), split = '\\.')))
    exps_this <- as.numeric(unlist(strsplit(as.character(df[i, col_exps]), split = '\\_')))
    df_sigs[i, match(gsub(sigs_this, pattern = 'Signature_', replace = 'exp_sig'), colnames(df_sigs))] <- exps_this
  }
  return(df_sigs)
}

# returns signature catalogs 
get_catalog <- function(cosmic_version = 'v3', output_file = NULL){
  # calculate the signature exposures of the cluster with maximum likelihood

  if(cosmic_version == "v3")
    cosmic_catalog <- pcawg_catalog
  if(cosmic_version == "v2")
    cosmic_catalog <- cosmic_catalog_v2
  if(!is.null(output_file)){
    message(paste0('catalog is saved to ', output_file))
    write.csv(cosmic_catalog, file = output_file, row.names = F, sep = ',', quote = F)
  }
  else{
    return(cosmic_catalog)
  }
}

