#' Calculates the compatibility of a list of genomes
#' to an input catalog based on likelihood and cosine 
#' similarity
#' 
#' @param genomes a data table or matrix with snv  
#' spectra in the first ntype columns and genomes in each row 
#' @param signatures the input catalog, a data table
#' with signature spectra in each column 
#' @param method can be 'median_catalog', 'weighted_catalog'
#' 'cosine_simil'. 'median_catalog' uses a custom signature
#' catalog formed by clustering genome SNV spectra and 
#' using it as a probability distribution. The 'median_catalog'
#' method can be used with any custom signatures data frame 
#' if the user intends to provide their own signature table.
#'
#' @return A data frame that contains the input genomes
#' and in addition columns associated to each signature in
#' in the catalog with likelihood and cosine simil values

match_to_catalog <- function(genomes, signatures, method = 'median_catalog'){
  nsig <- dim(signatures)[[2]]
  ntype <- dim(signatures)[[1]]

  calc_prob <- function(this_genome, signatures, normalize = T, success = NULL){
    eps <- 0.00001
    probs <- apply(signatures, 2, 
	           function(x, pow){ 
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
 

    if(normalize){
      mean_probs <- mean(probs)
      q1_probs <- mean(probs[probs < mean(probs)])
      max_probs <- max(probs)

      inds_keep <- which(probs >= q1_probs)
      inds_rm <- which(probs < q1_probs)
    
      probs[inds_keep] <- exp(probs[inds_keep] - max_probs)
      probs[inds_keep] <- probs[inds_keep]/sum(probs[inds_keep])

      probs[inds_rm] <- 0
    }

    if(!is.null(success))
      ind_max <- which(max(probs[success]) == probs)
    else
      ind_max <- which(max(probs) == probs)

    if(length(ind_max) > 0) ind_max <- ind_max[[1]]
     
    sig_max <- colnames(signatures)[[ind_max]]
    max_val <- probs[[ind_max]]
    return(list(probs = probs, ind_max = ind_max, sig_max = sig_max, max_val = max_val))
  }

  cosine <- function(x, y){
    return(x%*%y / sqrt(x%*%x * y%*%y))
  }

  calc_cos <- function(this_genome, signatures, success = NULL){
    simils <- apply(signatures, 2,
                    function(x){ cosine(this_genome, x) })
    
    if(!is.null(success))
      ind_max <- which(max(simils[success]) == simils)
    else
      ind_max <- which(max(simils) == simils)

    if(length(ind_max) > 0) ind_max = ind_max[[1]]

    sig_max <- colnames(signatures)[[ind_max]]
    max_val <- simils[[ind_max]]
    return(list(simils = simils, ind_max = ind_max, sig_max = sig_max, max_val = max_val))
  }

  error <- function(this_genome, signatures, exposures){
 
    reco <- (as.matrix(signatures) %*% exposures)
    error_frac <- sqrt(c(this_genome - reco) %*% c(this_genome - reco))/sqrt(this_genome %*% this_genome)
    return(error_frac)
  }

  decompose <- function(this_genome, signatures){
    dim2 <- dim(signatures)[[2]]

    if(dim2 == 2){
      exps <- coef(nnls::nnls(as.matrix(signatures), this_genome))
      inds <- which(exps != 0)
      exps <- exps[inds] 
      sigs <- colnames(signatures)[inds]

      error <- error(this_genome, signatures[, sigs], exps)
      return(list(signatures = sigs,
                  exposures = exps,
                  error = error))
    }

    min_error <- 1
    min_indices <- NULL
    min_exposures <- NULL
    for(i in 1:(dim2 - 2)){
      for(j in (i+1):(dim2 - 1)){
        for(k in (j+1):(dim2)){
          exposures <- coef(nnls::nnls(as.matrix(signatures[,c(i, j, k)]), this_genome))

          error_this <- error(this_genome, as.matrix(signatures[,c(i, j, k)]), exposures)
          if(error_this < min_error){
            min_error <- error_this
            min_indices <- c(i, j, k)
            min_exposures <- exposures
          }
        }
      }
    }
    return(list(signatures = colnames(signatures)[min_indices],
                exposures = min_exposures,
                error = min_error))
  }

  signature_names <- colnames(signatures)
  genome_matrix <- genomes[, 1:ntype]

  names_comb <- rep('', dim(genome_matrix)[[1]])
  combined_simil <- rep(0, dim(genome_matrix)[[1]])

  ##addition
  if(method == 'decompose'){
    output <- data.frame(sigs_all = character(), 
                                    exps_all = character(),
                                    comb_all_l = double(), 
                                    comb_all_c = double(),
                                    matrix(0, 0, nsig),
                                    matrix(0, 0, nsig),
                                    matrix(0, 0, nsig),
                                    matrix(0, 0, nsig * (nsig - 1) / 2),
                                    matrix(0, 0, nsig * (nsig - 1) / 2),
                                    sig_max_dl = character(),
                                    max_dl = double(),
                                    sig_max_dc = character(),
                                    max_dc = double())
    for(igenome in 1:dim(genome_matrix)[[1]]){
      x <- unlist(genome_matrix[igenome, ])

      # decompose using all signatures
      decomposed_all <- decompose(x, signatures)
      sigs <- signatures[, decomposed_all$signatures]
      exps <- decomposed_all$exposures
      sigs <- sigs[, exps > 0]
      exps <- exps[exps > 0]
      
      comb_all <- as.matrix(sigs) %*% exps
      sigs_all <- paste0(colnames(sigs), collapse = '_')
      exps_all <- paste0(exps, collapse = '_')

      # calculate likelihood and similarity of the combination
      probs_comb_all <- calc_prob(x, comb_all/sum(comb_all, normalize = F))$probs
      cos_simil_comb_all <- calc_cos(x, comb_all/sum(comb_all))$simils
      
      signatures_pair <- matrix(0,
                                ntype, 
                                nsig * (nsig - 1) / 2)
      signames_pair <- rep('', nsig * (nsig - 1) / 2)
      pair_success <- rep(T, nsig * (nsig - 1) / 2)
      probs_comb_without <- rep(0, nsig)
      cos_simil_without <- rep(0, nsig)
      prob_rats <- rep(0, nsig)

      count <- 1
      for(isig in 1:(nsig)){
        # decompose using all except the current signature
        
        decomposed_without <- decompose(x, signatures[, -isig])
        sigs <- signatures[, decomposed_without$signatures]
        exps <- decomposed_without$exposures
        sigs <- sigs[, exps > 0]
        exps <- exps[exps > 0]
               
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
  
        if(isig <= nsig - 1){
          for(jsig in (isig + 1):nsig){
            decomposed_pair <- decompose(x, signatures[, c(isig, jsig)])
            sigs <- signatures[, decomposed_pair$signatures]
            exps <- decomposed_pair$exposures
            if(length(exps) > 1){
              sigs <- sigs[, exps > 0]
              exps <- exps[exps > 0]
            }else{
              pair_success[[count]] <- F
            }

            comb_pair <- as.matrix(sigs) %*% exps
            signatures_pair[, count] <- comb_pair/sum(comb_pair)
            signames_pair[[count]] <- sprintf('%s_%s', 
                                       colnames(signatures)[[isig]],
                                       colnames(signatures)[[jsig]])
            
            count <- count + 1
          }
        }
      }

      colnames(signatures_pair) <- signames_pair
      tmp_pair <- calc_prob(x, signatures_pair, success = pair_success)
      probs_pair <- tmp_pair$probs
      max_sigs_pair <- tmp_pair$sig_max
      max_val_pair <- tmp_pair$max_val

      tmp_pair <- calc_cos(x, signatures_pair, success = pair_success)
      simils_pair <- tmp_pair$simils
      max_sigs_cos_pair <- tmp_pair$sig_max
      max_val_cos_pair <- tmp_pair$max_val 

      output_this <- data.frame(sigs_all = sigs_all,
                            exps_all = exps_all,
                            comb_all_l = probs_comb_all,
                            comb_all_c = cos_simil_comb_all,
                            matrix(probs_comb_without, 1, length(probs_comb_without)), 
                            matrix(cos_simil_without, 1, length(cos_simil_without)), 
                            matrix(prob_rats, 1, length(prob_rats)), 
                            matrix(probs_pair, 1, length(probs_pair)),
                            matrix(simils_pair, 1, length(simils_pair)),
                            sig_max_dl = max_sigs_pair,
                            max_dl = max_val_pair,
                            sig_max_dc = max_sigs_cos_pair,
                            max_dc = max_val_cos_pair)
      output <- rbind(output, output_this)
    }       
    colnames(output)[5:(4 + nsig)] <- paste0(colnames(signatures), '_wout_l')
    colnames(output)[(5 + nsig):(4 + 2 * nsig)] <- paste0(colnames(signatures), '_wout_c')
    colnames(output)[(5 + 2 * nsig):(4 + 3 * nsig)] <- paste0(colnames(signatures), '_l_rat')
    colnames(output)[(5 + 3 * nsig):(4 + 3 * nsig + nsig * (nsig - 1) / 2)] <- paste0(signames_pair, '_l')
    colnames(output)[(5 + 3 * nsig + nsig * (nsig - 1) / 2):(4 + 3 * nsig + nsig * (nsig - 1))] <- paste0(signames_pair, '_c')
  }

  ##addition
  if(method == 'weighted_catalog'){
    ind_bc <- grep('Signature_', colnames(weights_560_bc_cooccur_PCAWG_sig))
    signatures <- signatures[, colnames(weights_560_bc_cooccur_PCAWG_sig)[ind_bc]]
    weights <- weights_560_bc_cooccur_PCAWG_sig[, ind_bc]
    signatures <- signatures * weights 
  }

  if(method == 'weighted_catalog' | method == 'median_catalog'){
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