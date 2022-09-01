#' The quick_simulation() generates a simulated 96-dimensional
#' matrix specific for that tumor_type and based on the total
#' number of mutations in the dataset for which SigMA will be 
#' run on. It is useful for tuning a new model or determining
#' sensitivity and false positive rate for an existing model
#' see SigMA/test/test_determine_threshold.R for an example use case
#' 
#' @param input_file the path to the file that contains the
#' data table that will be used in run.R function
#' @param output_file the path to an output file to save the
#' simulated data 
#' @param tumor_type see list_tumor_types() for more information
#' @param data sequencing platform see list_data_options() for more 
#' information
#' @param remove_msi_pole can be set to TRUE to remove
#' microsattelite instable and POLE-exo domain mutated hypermutated
#' tumors because these bias the mutation count. If the data does
#' not contain hypermutated samples it can be set to FALSE to speed
#' up the calculation
#' @param run_SigMA boolean indicating whether the input file
#' is the output of an initial SigMA calculation (to be set to FALSE)
#' or whether it is an unprocessed file (to be set to TRUE).

quick_simulation <- function(input_file, 
                             output_file = NULL, 
                             tumor_type = NULL,  
                             data = NULL,
			     catalog_name = NULL,
                             remove_msi_pole = T, 
                             run_SigMA = T,
                             input_df = NULL, 
                             return_df = F, 
                             snv_cutoff = 5,
                             below_cutoff = NULL,
                             maf_percent = NULL){

  if(run_SigMA){
    if(is.null(catalog_name))
      stop('please provide a catalog_name argument options: ',  paste0(names(catalogs), sep = ' '))
    
    if(!(catalog_name %in% names(catalogs)))
      stop(paste0('catalog ', catalog_name, ' does not exist provide one of the following: ', paste0(names(catalogs), sep = ' ')))
      
    if(!is.null(input_file)) df_in <- read.csv(input_file)
    else df_in <- input_df 
 
    if(grepl('op_', data)) data_this <- 'op'
    else data_this <- data

    df <- run(input_file,
              tumor_type = tumor_type,
              data = data_this,
	      catalog_name = catalog_name,
              check_msi = remove_msi_pole,
              do_mva = F, #has_model(data, tumor_type),
              do_assign = T, 
              snv_cutoff = snv_cutoff,
              input_df = input_df,
              return_df = T)

    inds <- match(df_in$tumor, df$tumor)
    if(sum(is.na(inds)) != 0) 
      below_cutoff <- df_in[is.na(inds),]
  }
  else{
    if(!is.null(input_file) & is.null(input_df)){
      df <- read.csv(input_file)
    }
    else if(is.null(input_file) & !is.null(input_df)){
      df <- input_df
    }
    else{ 
      stop('only one of the following should be provided: input_df or input file')
    }
  }

  if(remove_msi_pole){ 
    inds <- which((df$Signature_msi_ml > 0.99 | df$Signature_pole_ml > 0.9999) & df$total_snvs > median(df$total_snvs))
    if(length(inds) != 0) df <- df[-inds,]
  }

  # determine approximate Sig3 estimates to scale mutation counts
  if(length(grep(paste0(paste0('Signature_3_c', 1:10, '_ml'), collapse = '|'), colnames(df)) == 1)){
    df$Signature_3_ml <- df$Signature_3_c1_ml
  }
  else{
    df$Signature_3_ml <- rowSums(df[, grep(paste0(paste0('Signature_3_c', 1:10, '_ml'), collapse = '|'), colnames(df))])
  }

  if(length(which(colnames(df) == "Signature_3_mva")) != 0) var <- 'Signature_3_mva'
  else var <- 'Signature_3_ml'
  

  # determine approximately the number of Sig3 pos samples
  stat_this <- stat[stat$tumor_type == tumor_type,]
  Sig3_frac <- mean(stat_this$Sig3_frac)  
 
  if(data == "wgs" | data == "wgs_pancan") sens <- 1
  else if(data == "seqcap" | data == "tcga_mc3" | data == "seqcap_probe") sens <- 0.9
  else if(grepl('op_', data) | data == "msk" | data == "fo") sens <- 0.7
  else sens <- 1

  if(is.null(below_cutoff)){
    n_pos <- round(Sig3_frac * sens * dim(df)[[1]], digit = 0) 
  }
  else{
    n_pos <- round(Sig3_frac * sens * (dim(df)[[1]] + dim(below_cutoff)[[1]]), digit = 0) 
  }
  
  val_cutoff <- sort(df[,var])[dim(df)[[1]] - n_pos]
  df$pass <- df[,var] > val_cutoff

  if(data == "seqcap_probe" | data == "tcga_mc3")
    data_matrix <- 'seqcap'
  else
    data_matrix <- data

  df_simul <- simulate_from_wgs(df, tumor_type = tumor_type, data = data_matrix, below_cutoff, maf_percent)
   
  if(return_df) return(df_simul)
  else{
    if(!is.null(output_file))
      write.table(df_simul, file = output_file, row.names = F, sep = ',', quote = F)
    else{
      output_file <- gsub(input_file, pattern = '.csv', replace = 'simulation.csv')
      write.table(df_simul, file = output_file, row.names = F, sep = ',', quote = F)
    }
    return(output_file)
  }
}

# adjusts the snv count based on the differences in the
# median SNV counts in the downsampled simulation and data
determine_scale <- function(df, m_ref, tumor_type, below_cutoff, scale_separately = F){
  message('determining SBS scale for quick simulation')
  # calculate the scale to match the data to simulations 
  if(sum(grepl(tumor_type, tumor_types_HRD)) > 0 & sum(df$pass, na.rm = T) > 0){
    scale_plus <- NULL
    scale_neg <- NULL
    if(sum(m_ref$is_sig3) > 0){
      total_snvs_plus <- median(df$total_snvs[df$pass], na.rm = T)
      median_pos_built_in <- median(rowSums(m_ref[m_ref$is_sig3, 1:96]), na.rm = T)
      scale_plus <- (total_snvs_plus - median_pos_built_in + 1*as.integer(median_pos_built_in == 0))/(median_pos_built_in + 1*as.integer(median_pos_built_in == 0))
    }
    if(sum(!m_ref$is_sig3) > 0){      
      total_snvs_neg <- median(df$total_snvs[!df$pass], na.rm = T)
      if(!is.null(below_cutoff))
        total_snvs_neg <- median(c(df$total_snvs[!df$pass], rowSums(below_cutoff[,1:96])), na.rm = T)
      median_neg_built_in <- median(rowSums(m_ref[which(!m_ref$is_sig3), 1:96]), na.rm = T)
      scale_neg <- (total_snvs_neg - median_neg_built_in + 1*as.integer(median_neg_built_in == 0))/(median_neg_built_in + 1*as.integer(median_neg_built_in == 0))
    }    
    if(!is.null(scale_plus) & !is.null(scale_neg) & !scale_separately){      
      if( abs(2*(scale_plus - scale_neg)/(scale_plus + scale_neg)) > 0.5 ){
        warning('The scales for positive and negative group disagree')
        match_total_snvs_combined <- T
      }
      else{
        match_total_snvs_combined <- F
        scale <- (scale_plus + scale_neg)/2
      }
    }
    if(!is.null(scale_plus) & is.null(scale_neg)){
      match_total_snvs_combined <- F 
      scale_separately <- F
      scale <- scale_plus
    }
    if(is.null(scale_plus) & !is.null(scale_neg)){
      match_total_snvs_combined <- F
      scale_separately <- F
      scale <- scale_neg
    }
  }
  else{
    match_total_snvs_combined <- T
  }
  if(match_total_snvs_combined){
    total_snvs <- median(rowSums(df[,1:96]), na.rm = T)
    if(!is.null(below_cutoff)){
      total_snvs <- median(c(rowSums(df[,1:96]), rowSums(below_cutoff[,1:96])), na.rm = T)
    }
    median_total_snvs <- median(rowSums(m_ref[,1:96]), na.rm = T) 
    scale <- (total_snvs - median_total_snvs)/median_total_snvs
  }
  m_ref$scale <- scale
  if(scale_separately){
    m_ref$scale[m_ref$is_sig3] <- scale_plus
    m_ref$scale[!m_ref$is_sig3] <- scale_neg
  }
  message('determining indel scale for quick simulation')

  # if the model will include mmej and nhej indel counts these counts should be present
  # in the data input from which simulations are generated
  if('mmej' %in% colnames(df) & 'nhej' %in% colnames(df)){
    inds_add_nonzero_mmej <- integer()
    inds_add_nonzero_nhej <- integer()
  
    if(sum(grepl(tumor_type, tumor_types_HRD)) > 0 & sum(df$pass, na.rm = T) > 0){
      scale_plus_mmej <- NULL
      scale_neg_mmej <- NULL
      scale_plus_nhej <- NULL
      scale_neg_nhej <- NULL
      
      if(sum(m_ref$is_sig3) > 0){
        # below is only used if mmej and nhej counts are provided in the input matrix 
        frac_mmej_nonzero <- sum(df$pass & df$mmej > 0, na.rm = T)/sum(df$pass, na.rm = T)
	frac_nhej_nonzero <- sum(df$pass & df$nhej > 0, na.rm = T)/sum(df$pass, na.rm = T)
	frac_simul_mmej_nonzero <- sum(m_ref$is_sig3 & m_ref$mmej > 0, na.rm = T)/sum(m_ref$is_sig3, na.rm = T)
	frac_simul_nhej_nonzero <- sum(m_ref$is_sig3 & m_ref$nhej > 0, na.rm = T)/sum(m_ref$is_sig3, na.rm = T)
	if(frac_simul_mmej_nonzero < frac_mmej_nonzero ){
	  inds_zero_mmej <- which(m_ref$is_sig3 & m_ref$mmej == 0)
	  inds_add_nonzero_mmej <- c(inds_add_nonzero_mmej,
	                          inds_zero_mmej[unique(round(runif(round(1.1*length(inds_zero_mmej)*(frac_mmej_nonzero - frac_simul_mmej_nonzero)), 0, length(inds_zero_mmej)), digit = 0))])
  	}
	if(frac_simul_nhej_nonzero < frac_nhej_nonzero){
	  inds_zero_nhej <- which(m_ref$is_sig3 & m_ref$nhej == 0)
	  inds_add_nonzero_nhej <- c(inds_add_nonzero_nhej,
	                          inds_zero_nhej[unique(round(runif(round(1.1*length(inds_zero_nhej)*(frac_nhej_nonzero - frac_simul_nhej_nonzero)), 0, length(inds_zero_nhej)), digit = 0))])
	}
      }
      if(sum(!m_ref$is_sig3) > 0){      
        frac_mmej_nonzero <- sum(!df$pass & df$mmej > 0, na.rm = T)/sum(!df$pass, na.rm = T)
	frac_nhej_nonzero <- sum(!df$pass & df$nhej > 0, na.rm = T)/sum(!df$pass, na.rm = T)
	frac_simul_mmej_nonzero <- sum(!m_ref$is_sig3 & m_ref$mmej > 0, na.rm = T)/sum(!m_ref$is_sig3, na.rm = T)
	frac_simul_nhej_nonzero <- sum(!m_ref$is_sig3 & m_ref$nhej > 0, na.rm = T)/sum(!m_ref$is_sig3, na.rm = T)
	if(frac_simul_mmej_nonzero < frac_mmej_nonzero){
	  inds_zero_mmej <- which(!m_ref$is_sig3 & m_ref$mmej == 0)
	  inds_add_nonzero_mmej <- c(inds_add_nonzero_mmej,
	                        inds_zero_mmej[unique(round(runif(round(1.1*length(inds_zero_mmej)*(frac_mmej_nonzero - frac_simul_mmej_nonzero)), 0, length(inds_zero_mmej)), digit = 0))])
	}
	if(frac_simul_nhej_nonzero < frac_nhej_nonzero){
	  inds_zero_nhej <- which(!m_ref$is_sig3 & m_ref$nhej == 0)
	  inds_add_nonzero_nhej <- c(inds_add_nonzero_nhej,
	                        inds_zero_nhej[unique(round(runif(round(1.1*length(inds_zero_nhej)*(frac_nhej_nonzero - frac_simul_nhej_nonzero)), 0, length(inds_zero_nhej)), digit = 0))])
	}        
      }    
      if(sum(!m_ref$is_sig3) > 0 & sum(m_ref$is_sig3) > 0){
        m_ref$mmej[inds_add_nonzero_mmej] <- 1
        m_ref$nhej[inds_add_nonzero_nhej] <- 1
        if(sum(m_ref$is_sig3) > 0){      
          total_mmej_plus <- mean(df$mmej[df$pass & df$mmej > 0], na.rm = T)
          total_nhej_plus <- mean(df$nhej[df$pass & df$nhej > 0], na.rm = T)
          mean_mmej_pos_built_in <- mean(m_ref$mmej[m_ref$is_sig3 & m_ref$mmej > 0], na.rm = T)
  	  mean_nhej_pos_built_in <- mean(m_ref$nhej[m_ref$is_sig3 & m_ref$nhej > 0], na.rm = T)
          scale_plus_mmej <- (total_mmej_plus - mean_mmej_pos_built_in)/mean_mmej_pos_built_in
          scale_plus_nhej <- (total_nhej_plus - mean_nhej_pos_built_in)/mean_nhej_pos_built_in
        }
        if(sum(!m_ref$is_sig3) > 0){      
          total_mmej_neg <- mean(df$mmej[!df$pass & df$mmej > 0], na.rm = T)
          total_nhej_neg <- mean(df$nhej[!df$pass & df$nhej > 0], na.rm = T)
          mean_mmej_neg_built_in <- mean(m_ref$mmej[!m_ref$is_sig3 & m_ref$mmej > 0], na.rm = T)
	  mean_nhej_neg_built_in <- mean(m_ref$nhej[!m_ref$is_sig3 & m_ref$nhej > 0], na.rm = T)
          scale_neg_mmej <- (total_mmej_neg - mean_mmej_neg_built_in)/mean_mmej_neg_built_in
          scale_neg_nhej <- (total_nhej_neg - mean_nhej_neg_built_in)/mean_nhej_neg_built_in
        }
        m_ref$scale_mmej[m_ref$is_sig3] <- scale_plus_mmej
        m_ref$scale_mmej[!m_ref$is_sig3] <- scale_neg_mmej
        m_ref$scale_nhej[m_ref$is_sig3] <- scale_plus_nhej
        m_ref$scale_nhej[!m_ref$is_sig3] <- scale_neg_nhej
      }
      else{
        frac_mmej_nonzero <- sum(df$mmej > 0, na.rm = T)/dim(df)[[1]]
        frac_nhej_nonzero <- sum(df$nhej > 0, na.rm = T)/dim(df)[[1]]
        frac_simul_mmej_nonzero <- sum(m_ref$mmej > 0, na.rm = T)/dim(m_ref)[[1]]
        frac_simul_nhej_nonzero <- sum(m_ref$nhej > 0, na.rm = T)/dim(m_ref)[[1]]
        if(frac_simul_mmej_nonzero - frac_mmej_nonzero < -0.05){
	  inds_zero_mmej <- which(m_ref$mmej == 0)
	  inds_add_nonzero_mmej <- c(inds_add_nonzero_mmej,
    	                            inds_zero_mmej[unique(round(runif(round(1.1*length(inds_zero_mmej)*(frac_mmej_nonzero - frac_simul_mmej_nonzero)), 0, length(inds_zero_mmej)), digit = 0))])
        }
        if(frac_simul_nhej_nonzero - frac_nhej_nonzero < -0.05){
          inds_zero_nhej <- which(m_ref$nhej == 0)
          inds_add_nonzero_nhej <- c(inds_add_nonzero_nhej,
	                             inds_zero_nhej[unique(round(runif(round(1.1*length(inds_zero_nhej)*(frac_nhej_nonzero - frac_simul_nhej_nonzero)), 0, length(inds_zero_nhej)), digit = 0))])
        }
        m_ref$mmej[inds_add_nonzero_mmej] <- 1
        m_ref$nhej[inds_add_nonzero_nhej] <- 1
        total_mmej <- mean(df$mmej[df$mmej > 0], na.rm = T)
        total_nhej <- mean(df$nhej[df$nhej > 0], na.rm = T)
        mean_mmej_built_in <- mean(m_ref$mmej[m_ref$mmej > 0], na.rm = T)
        mean_nhej_built_in <- mean(m_ref$nhej[m_ref$nhej > 0], na.rm = T)
        m_ref$scale_mmej <- (total_mmej - mean_mmej_built_in)/mean_mmej_built_in
        m_ref$scale_nhej <- (total_nhej - mean_nhej_built_in)/mean_nhej_built_in
      }
    }
  }
  return(m_ref)
}
 
# For tumor types where number of Sig3 cases are small
# Sig3 samples with similar signature composition from other
# tumor types are supplemented in the base matrix
# but with adjusted SNV counts to account for the tissue
# specific differences in SNV load in tumors
get_inter_tt_suppl <- function(tumor_type, data, matrices_96dim, maf_percent = NULL){
  
  data_dir <- system.file("extdata/matrices/matrices_96dim.rda", package="SigMA")
  load(data_dir)
  if(is.null(maf_percent)) matrices_96dim <- matrices_96dim[['matched_normal']] 
  else matrices_96dim <- matrices_96dim[[paste0('maf_', maf_percent)]]
  
  additional_tts <- names(inter_tt_suppl[[tumor_type]])
  if(is.null(additional_tts)) return(NULL)
  if(exists('m_comb')) rm(m_comb)
  for(tt in additional_tts){
    tumors <- inter_tt_suppl[[tumor_type]][[tt]] 
    m_this <- matrices_96dim[[data]][[tt]]
    m_this$tumor_type <- tt
    m_this <- m_this[na.omit(match(tumors, m_this$tumor)),]
    if(exists('m_comb')) m_comb <- rbind(m_comb, m_this)
    else m_comb <- m_this
  }
  if(!exists('m_comb')) return(NULL)
  else{
    return(m_comb)
  }
}

# Get the 96-dimensional matrix that will be used in the 
# simulation based on the tumor_type 
get_base_matrix <- function(df, tumor_type, data,  below_cutoff = NULL, main = F, maf_percent = NULL){

  data_dir <- system.file("extdata/matrices/matrices_96dim.rda", package="SigMA")
  load(data_dir)
  if(is.null(maf_percent)) matrices_96dim <- matrices_96dim[['matched_normal']]
  else matrices_96dim <- matrices_96dim[[paste0('maf_', maf_percent)]]
   
  m_out <- matrices_96dim[[data]][[tumor_type]]
  if(main){ 
    m_out <- determine_scale(df, m_out, tumor_type, below_cutoff)
  }
  m_out$main <- T
  m_out_suppl <- get_inter_tt_suppl(tumor_type, data, matrices_96dim, maf_percent) 
  if(!is.null(m_out_suppl)){
    m_out_suppl$main <- F
    if(main){
      if(sum(colnames(m_out_suppl) == "tumor_type") > 0){
        for(tt in unique(m_out_suppl$tumor_type)){
          this <- m_out_suppl[m_out_suppl$tumor_type == tt,]
          this <- determine_scale(df, this, tumor_type, below_cutoff)
          if(exists('m_out_suppl_tmp')){
            m_out_suppl_tmp <- rbind(m_out_suppl_tmp, this)
          }
          else{
            m_out_suppl_tmp <- this
          }
        }
        m_out_suppl <- m_out_suppl_tmp
        rm('m_out_suppl_tmp')
      }
      else{
        m_out_suppl <- determine_scale(df, m_out_suppl, tumor_type, below_cutoff)
      }
    }
    if(sum(colnames(m_out_suppl) == tumor_type) > 0){
      m_out_suppl <- m_out_suppl[, -which(colnames(m_out_suppl) == "tumor_type")]
      m_out <- rbind(m_out, m_out_suppl) 
    }
    else{
      m_out <- rbind(m_out, m_out_suppl[,na.omit(match(colnames(m_out), colnames(m_out_suppl)))]) 
    }
  }
  return(m_out)
}

# The simulations are generated by adjusting the SNV count with
# MC simulations
simulate_from_wgs <- function(df, tumor_type, data, below_cutoff, maf_percent){
  if(data == "wgs"){
    m_main <- get_base_matrix(df = df, tumor_type = tumor_type, data = 'wgs', main = T, maf_percent = maf_percent)
    if(sum(m_main$scale < 0 & !m_main$main) > 0){
      inds <- which(m_main$scale < 0 & !m_main$main)
      m_main[inds,] <- adjust_snvs(m_main[inds,])
    } 
    if(sum(m_main$scale > 0 & !m_main$main) > 0){
      inds <-m_main$scale > 0 & !m_main$main
      m_main[inds, 1:96] <- (1 + m_main$scale[inds]) * m_main[inds, 1:96]
    }
    if(sum(c('mmej','nhej') %in% colnames(df)) == 2){
      if(sum(m_main$scale_mmej < 0 & !m_main$main) > 0){
        inds <- which(m_main$scale_mmej < 0 & !m_main$main)
        m_main[inds,] <- adjust_ids(m_main[inds,], type = 'mmej')
      } 
      if(sum(m_main$scale_nhej < 0 & !m_main$main) > 0){
        inds <- which(m_main$scale_nhej < 0 & !m_main$main)
        m_main[inds,] <- adjust_ids(m_main[inds,], scale = 'nhej')
      } 
      if(sum(m_main$scale_mmej > 0 & !m_main$main) > 0){
        inds <-m_main$scale_mmej > 0 & !m_main$main
        m_main[inds, 1:96] <- (1 + m_main$scale_mmej[inds]) * m_main[inds, 1:96]
      }
      if(sum(m_main$scale_nhej > 0 & !m_main$main) > 0){
        inds <-m_main$scale_nhej > 0 & !m_main$main
        m_main[inds, 1:96] <- (1 + m_main$scale_nhej[inds]) * m_main[inds, 1:96]
      }
    }
    return(m_main)
  }

  # if data is not specified wgs is downsampled
  if(is.null(data)) data = 'wgs'

  m_main <- get_base_matrix(df = df, tumor_type = tumor_type, data = data, main = T, below_cutoff = below_cutoff, maf_percent = maf_percent)

  if(max(m_main$scale) < 0){
    m_main <- adjust_snvs(m_main)
  }
  else{
    if(data == "seqcap" | data == "seqcap_probe"){
      m_extra <- get_base_matrix(df = df, tumor_type = tumor_type, data = 'wgs', below_cutoff = below_cutoff, maf_percent = maf_percent)
    }
    if(data == "msk" | data == "op" | data == "op_test" | data == "fo"){ 
      m_extra <- get_base_matrix(df = df, tumor_type = tumor_type, data = 'seqcap', below_cutoff = below_cutoff, maf_percent = maf_percent)
    }

    m_main <- adjust_snvs(m_main[!is.na(match(m_main$tumor, m_extra$tumor)),], 
                          m_extra[na.omit(match(m_main$tumor, m_extra$tumor)),]) 
  }

  
  if(sum(c('mmej','nhej') %in% colnames(df)) == 2){    
    if(max(m_main$scale_mmej) < 0){
      m_main <- adjust_ids(m_main, type = 'mmej')
    }
    else{
      if(data == "seqcap" | data == "seqcap_probe"){
        m_extra <- get_base_matrix(df = df, tumor_type = tumor_type, data = 'wgs', below_cutoff = below_cutoff, maf_percent = maf_percent)
      }
      if(data == "msk" | data == "op" | data == "op_test" | data == "fo"){ 
        m_extra <- get_base_matrix(df = df, tumor_type = tumor_type, data = 'seqcap', below_cutoff = below_cutoff, maf_percent = maf_percent)
      }
      m_main <- adjust_ids(m_main[!is.na(match(m_main$tumor, m_extra$tumor)),], 
                            m_extra[na.omit(match(m_main$tumor, m_extra$tumor)),],
			    type = 'mmej')
    }

    if(max(m_main$scale_nhej) < 0){
      m_main <- adjust_ids(m_main, type = 'nhej')
    }
    else{
      if(data == "seqcap" | data == "seqcap_probe"){
        m_extra <- get_base_matrix(df = df, tumor_type = tumor_type, data = 'wgs', below_cutoff = below_cutoff, maf_percent = maf_percent)
      }
      if(data == "msk" | data == "op" | data == "op_test" | data == "fo"){ 
        m_extra <- get_base_matrix(df = df, tumor_type = tumor_type, data = 'seqcap', below_cutoff = below_cutoff, maf_percent = maf_percent)
      }
      m_main <- adjust_ids(m_main[!is.na(match(m_main$tumor, m_extra$tumor)),], 
                            m_extra[na.omit(match(m_main$tumor, m_extra$tumor)),],
			    type = 'nhej')
    }

  }

  return(m_main)
}

adjust_snvs <- function(m, m_extra = NULL){
  if(is.null(m_extra)){
    for(i in 1:dim(m)[[1]]){
      spec <- m[i, 1:96]
      scale <- m$scale[[i]]
      if(abs(scale) < 0.5){
        count <- round(sum(spec) * scale, digit = 0)
        m[i, 1:96] <- spec - subsample(spec, abs(count))
      }
      else{
        count <- round(sum(spec) * (1 + scale), digit = 0)
        m[i, 1:96] <- subsample(spec, abs(count))
      }
    }
  }
  else{
    for(i in 1:dim(m)[[1]]){
      spec <- m[i, 1:96]
      spec_extra <- m_extra[i, 1:96]

      count <- round(sum(spec) * m$scale[[i]], digit = 0)
      if(count < 0) m[i, 1:96] <- spec - subsample(spec, abs(count))
      else m[i, 1:96] <- spec + subsample(spec_extra, count)
      m[i,m[i, 1:96] < 0] <- 0
    }
  }
  return(m)
}

subsample <- function(spec, count){
  count <- min(sum(spec), count)
  spec_out <- rep(0, 96)
  max_val <- max(spec) * 1.1
  while(sum(spec_out) < count){
    ind <- round(runif(1, 0.5, 96.5), digit = 0)
    val <- runif(1, 0, max_val)
    if(spec[ind] > val & spec_out[ind] < spec[ind])
      spec_out[ind] <- spec_out[ind] + 1
  }
  return(spec_out)
}

adjust_ids <- function(m, m_extra = NULL, type = NULL){
  if(is.null(m_extra)){
    for(i in 1:dim(m)[[1]]){ 
      id_count <- m[i, type]
      scale <- m[i, paste0('scale_', type)]
#      scale <- rpois(lambda = scale, n = 1)
      if(abs(scale) < 0.5){
        count <- round(id_count * scale, digit = 0)
        m[i, type] <- id_count - subsample_id(id_count, abs(count))
      }
      else{
        count <- round(id_count * (1 + scale), digit = 0)
        m[i, type] <- subsample_id(id_count, abs(count))
      }
    }
  }
  else{
    for(i in 1:dim(m)[[1]]){
      id_count <- m[i, type]
      id_count_extra <- m_extra[i, type]
      scale <- m[i, paste0('scale_', type)]
#      scale <- rpois(lambda = scale, n = 1)
      
      count <- round(id_count * scale, digit = 0)
      if(count < 0) m[i, type] <- id_count - subsample_id(id_count, abs(count))
      else m[i, type] <- id_count + subsample_id(id_count_extra, count)
      if(m[i, type] < 0) m[i,type] <- 0     
    }
  }
  return(m)
}


subsample_id <- function(id_count, count){
  count <- min(id_count, count)
  id_out <- 0
  max_val <- id_count * 1.1
  while(sum(id_out) < count){
    val <- runif(1, 0, max_val)
    if(id_count > val & id_out < id_count)
      id_out <- id_out + 1
  }
  return(id_out)
}

