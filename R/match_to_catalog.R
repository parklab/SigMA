#' Calculates the compatibility of a list of genomes
#' to an input catalog based on likelihood and cosine 
#' similarity
#' 
#' @param genomes a data table or matrix with snv  
#' spectra in the first ntype columns and genomes in each row 
#' @param signatures the input catalog, a data table
#' with signature spectra in each column 
#' @param method can be 'median_catalog', 'weighted_catalog'
#' 'cosine_simil' or 'decompose. 'median_atalog' uses the 
#' signature catalog formed by clustering genome SNV spectra and 
#' using it as a probability distribution. The 'median_catalog'
#' method can be used with any custom signatures data frame 
#' if the user intends to provide their own signature table.
#' @param data sets the type of sequencing platform used, options 
#' are 'msk', 'seqcap', 'wgs'
#'
#' @return A data frame that contains the input genomes
#' and in addition columns associated to each signature in
#' in the catalog with likelihood and cosine simil values

match_to_catalog <- function(genomes, signatures, method = 'median_catalog', data){
  nsig <- dim(signatures)[[2]]
  ntype <- dim(signatures)[[1]]

  #calculates the likelihood
  calc_prob <- function(this_genome, signatures, normalize = T, success = NULL){
    eps <- 0.00001
    probs <- apply(signatures, 2, 
	           function(x, pow){ 
                     if(sum(is.na(x)) > 0) return(0)
                     prob <- 0
                     for(i in 1:length(x)){
                       if(x[[i]] == 0) x[[i]] <- eps
	               if(pow[[i]] > 0){
                         prob <- prob + pow[[i]] * log(x[[i]]) 
                       }
	             }
	             return(prob)
                   },
                   pow = this_genome)
 

    if(normalize){ # normalizes the total prob for all clusters to be 1
      mean_probs <- mean(probs)
      
      q1_probs <- mean(probs[probs < mean(probs)])
      max_probs <- max(probs)
   
      # often having too small values leads to infinities so
      # only the probabilities in the first quartile are kept
      # to be nonzero
      inds_keep <- which(probs >= q1_probs)
      inds_rm <- which(probs < q1_probs)
    
      probs[inds_keep] <- exp(probs[inds_keep] - max_probs)
      probs[inds_keep] <- probs[inds_keep]/sum(probs[inds_keep])

      probs[inds_rm] <- 0
    }

    # the success settings are currently not used it allows to select
    # a subset of signatures from the whole table
    if(!is.null(success)){
      if(sum(success) == 0){
        return(list(probs = probs, ind_max = NA, sig_max = NA, max_val = NA))
      }
      ind_max <- which(max(probs[success]) == probs) 
    }else
      ind_max <- which(max(probs) == probs)

    # signature with highest likelihood
    if(length(ind_max) > 0) ind_max <- ind_max[[1]]
    sig_max <- colnames(signatures)[[ind_max]]
    max_val <- probs[[ind_max]]
 
    return(list(probs = probs, ind_max = ind_max, sig_max = sig_max, max_val = max_val))
  }

  # cosine similarity
  cosine <- function(x, y){
    return(x%*%y / sqrt(x%*%x * y%*%y))
  }

  # calculates cosine similarity for all signature options
  calc_cos <- function(this_genome, signatures, success = NULL){
    simils <- apply(signatures, 2,
                    function(x){ cosine(this_genome, x) })
    
    # the success settings are currently not used it allows to select
    # a subset of signatures from the whole table
    if(!is.null(success)){
      if(sum(success) == 0){
        return(list(simils = simils, ind_max = NA, sig_max = NA, max_val = NA))
      }
      ind_max <- which(max(simils[success]) == simils)
    }else
      ind_max <- which(max(simils) == simils)

    # signature with maximum cosine similarity and its value
    if(length(ind_max) > 0) ind_max = ind_max[[1]]
    sig_max <- colnames(signatures)[[ind_max]]
    max_val <- simils[[ind_max]]

    return(list(simils = simils, 
                ind_max = ind_max, 
                sig_max = sig_max, 
                max_val = max_val))
  }

  signature_names <- colnames(signatures)
  genome_matrix <- genomes[, 1:ntype]

  names_comb <- rep('', dim(genome_matrix)[[1]])
  combined_simil <- rep(0, dim(genome_matrix)[[1]])

  # each genome is decomposed with NNLS into corresponding signatures
  # and likelihood of the decomposition with and without each signature
  # is calculated
  
  if(method == 'decompose'){
    output <- data.frame(sigs_all = character(), 
                                    exps_all = character(),
                                    comb_all_l = double(), 
                                    comb_all_c = double(),
                                    exp_sig3 = double(),
                                    matrix(0, 0, nsig),
                                    matrix(0, 0, nsig))

    pb <- txtProgressBar(min = 0, max = dim(genome_matrix)[[1]], style = 3)

    for(igenome in 1:dim(genome_matrix)[[1]]){
      x <- unlist(genome_matrix[igenome, ])
      # decompose using all signatures
      decomposed_all <- decompose(x, signatures, data)
      sigs <- signatures[, decomposed_all$signatures]
      exps <- decomposed_all$exposures
      if(!is.null(dim(sigs))){
        if(!is.null(dim(sigs))){
          sigs <- sigs[, exps > 0]
          exps <- exps[exps > 0]
        }
      }  
      if(length(na.omit(match('Signature_3', colnames(sigs)))) > 0){
        sig3_exp <- exps[[which(colnames(sigs) == 'Signature_3')]]
      }
      else sig3_exp <- 0

      comb_all <- as.matrix(sigs) %*% exps
      if(length(exps) > 1){
        sigs_all <- paste0(colnames(sigs), collapse = '.')
        exps_all <- paste0(exps, collapse = '_')
      }else{
        sigs_all <- decomposed_all$signatures[[which(decomposed_all$exposures > 0)]]
        exps_all <- as.character(exps)
      }

      # calculate likelihood and similarity of the combination
      probs_comb_all <- calc_prob(x, comb_all/sum(comb_all, normalize = F))$probs
      cos_simil_comb_all <- calc_cos(x, comb_all/sum(comb_all))$simils
      
      probs_comb_without <- rep(0, nsig)
      cos_simil_without <- rep(0, nsig)
      prob_rats <- rep(0, nsig)

      count <- 1
      for(isig in 1:(nsig)){
        # decompose using all except the current signature
        if(nsig > 2){
          decomposed_without <- decompose(x, signatures[, -isig], data)
          sigs <- signatures[, decomposed_without$signatures]
          exps <- decomposed_without$exposures
          if(!is.null(dim(sigs))){
            sigs <- sigs[, exps > 0]
            exps <- exps[exps > 0]
          }
           
          comb_without <- as.matrix(sigs) %*% exps
       
          # calculate_likelihood and similarity of the combination
          probs_comb_without[[isig]] <- calc_prob(x, comb_without/sum(comb_without))$probs
          cos_simil_without[[isig]] <- calc_cos(x, comb_without/sum(comb_without))$simils

          # calculate a ratio of likelihood obtained with and without this
          # signature if it is one of those that were picked out of all signatures
          if(grepl(colnames(signatures)[[isig]], sigs_all)){
            average_logl <- (probs_comb_all + probs_comb_without[[isig]])/2
            rat_prob <- exp(probs_comb_all - average_logl)/(exp(probs_comb_all - average_logl) 
                                                            + exp(probs_comb_without[[isig]] - average_logl))
          }else
            rat_prob <- 0

          prob_rats[[isig]] <- rat_prob
        }  
      }

      output_this <- data.frame(sigs_all = sigs_all,
                            exps_all = exps_all,
                            comb_all_l = probs_comb_all,
                            comb_all_c = cos_simil_comb_all,
                            exp_sig3 = sig3_exp,  
                            matrix(cos_simil_comb_all - cos_simil_without, 1, length(cos_simil_without)),
                            matrix(prob_rats, 1, length(prob_rats))) 
      output <- rbind(output, output_this)
      setTxtProgressBar(pb, igenome)
    }       

    message('\n')
    colnames(output)[6:(5 + nsig)] <- paste0(colnames(signatures), '_c_diff')
    colnames(output)[(6 + nsig):(5 + 2 * nsig)] <- paste0(colnames(signatures), '_l_rat')
  }

  if(method == 'weighted_catalog'){
    ind_bc <- grep('Signature_', colnames(weights_560_bc_cooccur_PCAWG_sig))
    ind_bc_2 <- na.omit(match(colnames(signatures), colnames(weights_560_bc_cooccur_PCAWG_sig)[ind_bc]))
    signatures <- signatures[, colnames(weights_560_bc_cooccur_PCAWG_sig)[ind_bc[ind_bc_2]]]
    weights <- weights_560_bc_cooccur_PCAWG_sig[, ind_bc[ind_bc_2]]
    signatures <- signatures * weights
  }

  if(method == 'median_catalog' | method == "weighted_catalog"){
    probs_all    <- apply(genome_matrix, 1, 
                          function(x, y){ calc_prob(x, y)$probs }, y = signatures)
    max_sigs_all <- apply(genome_matrix, 1, 
                          function(x, y){ calc_prob(x, y)$sig_max }, y = signatures)
    max_val_all  <- apply(genome_matrix, 1, 
                          function(x, y){ calc_prob(x, y)$max_val }, y = signatures)
  }

  if(method == 'cosine_simil'){
    simils_all    <- apply(genome_matrix, 1, 
                          function(x, y){calc_cos(x, y)$simils}, y = signatures)
    max_sigs_cos_all <- apply(genome_matrix, 1, 
                          function(x, y){calc_cos(x, y)$sig_max}, y = signatures)
    max_val_cos_all  <- apply(genome_matrix, 1, 
                          function(x, y){calc_cos(x, y)$max_val}, y = signatures) 
  }
  
  if(method == 'weighted_catalog' | method == 'median_catalog'){
    output <- data.frame(t(probs_all),
                         sig_max = max_sigs_all,
                         max = max_val_all)

    if(method == 'weighted_catalog') mname = 'wl'
    if(method == 'median_catalog') mname = 'ml'
    colnames(output) <- c(paste0(colnames(signatures), '_', mname),
                          paste0('sig_max_', mname),
                          paste0('max_', mname))
  }

  if(method == 'cosine_simil'){ 
    output <- data.frame(t(simils_all), 
                         sig_max_c = max_sigs_cos_all,
                         max_c = max_val_cos_all)
    colnames(output)[1:dim(simils_all)[[1]]] <- paste0(colnames(signatures), '_c')
  }

  return(output)
}