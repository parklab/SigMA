# gets the clusters with maximum likelihood and decomposes them with NNLS to determine relative exposures of signatures
llh_max_characteristics <- function(df, tumor_type, catalog_name = 'cosmic_v2_inhouse'){
  # calculate the signature exposures of the cluster with maximum likelihood
  cosmic_catalog <- catalogs[[catalog_name]]
  
  if('sig_max_ml_msi' %in% colnames(df))
    df$max_clust <- df$sig_max_ml_msi
  else if('sig_max_ml' %in% colnames(df))
    df$max_clust <- df$sig_max_ml

  signatures <- cosmic_catalog[,signames_per_tissue_per_catalog[[catalog_name]][[tumor_type]]]
  if('sig_max_ml_msi' %in% colnames(df))
    signatures_msi_pole <- cbind(signatures, cosmic_catalog[,c(signames_per_tissue_per_catalog[['msi']], signames_per_tissue_per_catalog[['pole']])])
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
    
    if(!grepl('Signature_3_', max_clust))
      signatures_this <- signatures_this[,colnames(signatures_this) != "Signature_3"]

    exps <- coefficients(nnls::nnls(as.matrix(signatures_this), median_clust))
  
    df$cluster_sigs_all[df$max_clust == max_clust] <- paste0(colnames(signatures_this), collapse = '.')
    df$cluster_exps_all[df$max_clust == max_clust] <- paste0(exps, collapse = '_')
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
get_catalog <- function(catalog_name = 'cosmic_v2_inhouse', output_file = NULL){
  # calculate the signature exposures of the cluster with maximum likelihood
  if(!(catalog_name %in% names(catalogs)))
    stop(paste0('Available catalog names are ', paste0(names(catalogs), collapse = ' ')))
    
  cosmic_catalog <- catalogs[[catalog_name]]
  
  if(!is.null(output_file)){
    message(paste0('catalog is saved to ', output_file))
    write.csv(cosmic_catalog, file = output_file, row.names = F, sep = ',', quote = F)
  }
  else{
    return(cosmic_catalog)
  }
}

# the tumor types where the specified signature was discovered
get_tumor_types_for_signature <- function(catalog_name = 'cosmic_v2_inhouse', signature = NULL){
  is_in_tt <- lapply(signames_per_tissue_per_catalog[[catalog_name]], function(x, signature){ signature %in% x }, signature = signature)
  return(unlist(names(is_in_tt))[c(unlist(is_in_tt))])
}

# gets the expected probabilities based on per sample spectra in wgs data
get_wgs_data_for_llh <- function(tumor_type = NULL, check_msi = F){
  data_dir <- system.file("extdata/matrices/matrices_96dim.rda", package="SigMA")
  load(data_dir)
  m_wgs <- matrices_96dim[['matched_normal']][['wgs']][[tumor_type]]
  if(check_msi){
    m_wgs <- rbind(m_wgs,
                   matrices_96dim[['matched_normal']][['wgs']][['msi']],
                   matrices_96dim[['matched_normal']][['wgs']][['pole']])
  }
  m_wgs[,1:96] <- m_wgs[,1:96]/rowSums(m_wgs[,1:96])
  df_wgs <- data.frame(t(m_wgs[,1:96]))
  colnames(df_wgs) <- m_wgs$tumor
  return(df_wgs)
}

# gets the total Sig3, msi, and pole likelihoods based on likelihood of Sig3+, msi, and pole samples in WGS data
# to be used after get_wgs_data_for_llh() function
get_Sig3_llh_from_wgs_data <- function(tumors, tumor_type = NULL, check_msi = F){
  data_dir <- system.file("extdata/matrices/matrices_96dim.rda", package="SigMA")
  load(data_dir)
  m_wgs <- matrices_96dim[['matched_normal']][['wgs']][[tumor_type]]
  df <- data.frame(tumor = tumors)
  df$is_sig3 <- m_wgs$is_sig3[match(df$tumor, m_wgs$tumor)]
  df$is_sig3[is.na(df$is_sig3)] <- F
  if(check_msi){
    m_wgs_msi <- matrices_96dim[['matched_normal']][['wgs']][['msi']]
    df$is_msi <- !is.na(match(df$tumor, m_wgs_msi$tumor))
    m_wgs_pole <- matrices_96dim[['matched_normal']][['wgs']][['pole']]
    df$is_pole <- !is.na(match(df$tumor, m_wgs_pole$tumor))
  }
  return(df)
}

# matches to individual samples and finds the maximum likely sample or calculates an average over likelihoods
llh_max_characteristics_wgs_data <- function(df,
                                             tumor_type,
					     catalog_name = 'cosmic_v2_inhouse',
					     check_msi = F,
					     add_average_over_all_samples = F){
  # calculate the signature exposures of the sample with maximum likelihood
  cosmic_catalog <- catalogs[[catalog_name]]
  signatures <- cosmic_catalog[,signames_per_tissue_per_catalog[[catalog_name]][[tumor_type]]]

  df$max_clust <- gsub(df$sig_max_ml, pattern = '_ml', replace = '')
  
  df_wgs <- get_wgs_data_for_llh(tumor_type = tumor_type, check_msi = check_msi)
  df_wgs_is_sig3 <- get_Sig3_llh_from_wgs_data(colnames(df_wgs), tumor_type = tumor_type, check_msi = check_msi)

  if(add_average_over_all_samples){
    df_exps_wgs <- data.frame(matrix(0, dim(df_wgs)[[2]], length(colnames(signatures))))
    colnames(df_exps_wgs) <- colnames(signatures)
  
    for(i in 1:dim(df_wgs)[[2]]){
      signatures_this <- signatures
      spec <- df_wgs[,i]
      if(!df_wgs_is_sig3$is_sig3[which(df_wgs_is_sig3$tumor == colnames(df_wgs)[[i]])])
        signatures_this <- signatures_this[,colnames(signatures_this) != "Signature_3"]
      if(check_msi){
        if(!df_wgs_is_sig3$is_msi[which(df_wgs_is_sig3$tumor == colnames(df_wgs)[[i]])])
          signatures_this <- cbind(signatures_this,
                                   cosmic_catalog[,c(signames_per_tissue_per_catalog[['msi']], signames_per_tissue_per_catalog[['msi_extra']])])
        if(!df_wgs_is_sig3$is_pole[which(df_wgs_is_sig3$tumor == colnames(df_wgs)[[i]])])
          signatures_this <- cbind(signatures_this,
                                   cosmic_catalog[,c(signames_per_tissue_per_catalog[['msi']], signames_per_tissue_per_catalog[['msi_extra']], signames_per_tissue_per_catalog[['pole']])])
      }
      exps <- coefficients(nnls::nnls(as.matrix(signatures_this), spec))
      df_exps_wgs[i,colnames(signatures_this)] <- exps
    }

    df_ave <- as.matrix(df[,paste0(colnames(df_wgs), '_ml')]) %*% as.matrix(df_exps_wgs)
    colnames(df_ave) <- paste0('ave_ml_', gsub(colnames(signatures), pattern = 'Signature_', replace = 'exp_sig'))
  }
  
  df$max_ml_sample_sigs_all <- ''
  df$max_ml_sample_exps_all  <- ''

  df_is_sig3 <- get_Sig3_llh_from_wgs_data(df$max_clust, tumor_type = tumor_type, check_msi = check_msi)
  
  for(i in 1:dim(df)[[1]]){
    add_sig3 <- F
    add_msi <- F
    add_pole <- F
    
    if(df_is_sig3$is_sig3[i] | sum(df[i, paste0(df_wgs_is_sig3$tumor[df_wgs_is_sig3$is_sig3], '_ml')]) > 0.5)
      add_sig3 <- T
      
    if(check_msi){
      if(df_is_sig3$is_msi[i] | sum(df[i, paste0(df_wgs_is_sig3$tumor[df_wgs_is_sig3$is_msi], '_ml')]) > 0.9) 
        add_msi <- T
      if(df_is_sig3$is_pole[i] | sum(df[i, paste0(df_wgs_is_sig3$tumor[df_wgs_is_sig3$is_pole], '_ml')]) > 0.9)
        add_pole <- T
    }
      
    signatures_this <- signatures
    if(!add_sig3)
      signatures_this <- signatures[,colnames(signatures) != "Signature_3"]
    if(check_msi){
      if(add_msi)
        signatures_this <- cbind(signatures_this,
                                 cosmic_catalog[,c(signames_per_tissue_per_catalog[['msi']], signames_per_tissue_per_catalog[['msi_extra']])])
      if(add_pole)
        signatures_this <- cbind(signatures_this,
                                 cosmic_catalog[,c(signames_per_tissue_per_catalog[['msi']], signames_per_tissue_per_catalog[['msi_extra']], signames_per_tissue_per_catalog[['pole']])])
    }

    median_clust <- c(unlist(df_wgs[, df$max_clust[i]]))    
    exps <- coefficients(nnls::nnls(as.matrix(signatures_this), median_clust))
    
    df$max_ml_sample_sigs_all[i] <- paste0(colnames(signatures_this), collapse = '.')
    df$max_ml_sample_exps_all[i] <- paste0(exps, collapse = '_')
  }
  if(add_average_over_all_samples)
    df <- cbind(df, df_ave)
  return(df)
}
