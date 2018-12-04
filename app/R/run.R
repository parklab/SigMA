#' Runs SigMA: (1) calculates likelihood, cosine similarity,
#' NNLS exposures, and likelihood of the decomposition. (2) These
#' features are later used in multivariate analysis. (3) Based
#' on scores a final decision on existence of the signature.
#'
#' @param genome_file a csv file with snv spectra info can be created
#' from vcf file using @make_genome_matrix() function
#' see ?make_genome_matrix
#' @param output_file the output file name, can be NULL in which
#' case input file name is used and appended with "_output"
#' @param data the options are "msk" (for a panel that is similar 
#' size to MSK-Impact panel with 410 genes), "seqcap" (for whole  
#' exome sequencing), or "wgs" (for whole genome sequencing)
#' @param tumor_type the options are "bladder", "bone_other" (Ewing's 
#' sarcoma or Chordoma), "breast", "crc", "eso", "gbm", "lung", 
#' "lymph", "medullo", "osteo", "ovary", "panc_ad", "panc_en",
#' "prost", "stomach", "thy", or "uterus". The exact correspondance
#' of these names can be found in https://github.com/parklab/SigMA
#' @param do_assign boolean for whether a cutoff should be applied 
#' to determine the final decision or just the features should 
#' be returned
#' @param do_mva a boolean for whether multivariate analysis 
#' should be run
#' @param check_msi is a boolean which determines whether the user
#' wants to identify micro-sattelite instable tumors
#'
#' @examples
#' run(genome_file = "input_genomes.csv", 
#'     data = "msk",
#'     tumor_type = "ovary")
#' run(genome_file = "input_genomes.csv", 
#'     data = "seqcap", 
#'     tumor_type = "bone_other")

run <- function(genome_file, 
                output_file = NULL,
                do_assign = T,
                data = "msk",
                tumor_type = "breast",
                do_mva = T,
                check_msi = F, 
                weight_cf = F,
                lite_format = F,
                add_sig3 = F){

  # give a warning if weight_cf is TRUE but the sequencing platform
  # is not a panel
  if(data != 'msk') weight_cf = T
  
  # if a custom output file is not defined 
  # use the input file name
  if(is.null(output_file))
    output_file = gsub(genome_file,
                       pattern = ".csv",
                       replacement = paste0("_output_tumortype_",
                                            tumor_type,
                                            "_platform_",
                                            gsub(platform_names[[data]], 
                                                 pattern = " ", replace = ""),
                                            "_cf", as.integer(weight_cf),
                                            ".csv"))

  if((tumor_type == "medullo" | tumor_type == "bone_other") & data == "msk")
    stop('tumor_type medullo and bone_other cannot be used msk data setting
          use seqcap or wgs data options instead')

  genomes <- read.csv(genome_file)

  # remove genomes with no mutation 
  if(sum(is.na(rowSums(genomes[, 1:96]))) > 0){
    genomes <- genomes[!is.na(rowSums(genomes[, 1:96])), ]
  }
  if(sum(rowSums(genomes[, 1:96]) == 0) > 0){
    genomes <- genomes[which(rowSums(genomes[, 1:96]) > 0), ]
  }

  # lower cutoff on number mutations for SigMA
  if(do_assign | do_mva){
    if(data == "msk"){
      if(tumor_type == "prost")
        genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 4), ]
      else if(tumor_type == "osteo" | tumor_type == "panc_en")
        genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 3), ]
      else
        genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 5), ]      
    }else{ # for exomes and wgs a larger lower cutoff is applied
      if(tumor_type == "bone_other" | tumor_type == "medullo")
        genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 5), ]
      else
        genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 10), ]
    }
  }

  message(paste0("You are running SigMA for ", dim(genomes)[[1]], 
                 " ", tissue_names[[tumor_type]],  "(s) sequenced by ", 
                 platform_names[[data]]))


  # method names to calculate different features
  methods <- c("median_catalog", "cosine_simil", "decompose")

  # signature catalogs to be used for each method
  sig_catalogs <- c("average", "cosmic", "cosmic_tissue")

  # first the samples are assumed to be mss and later msi 
  # calculations are done
  steps <- c("mss", "mss", "mss")

  # signature to be identified from mss samples
  if(tumor_type == "gbm"){
    signames <- rep("Signature_8", 3)
  }
  else{
    signames <- rep("Signature_3", 3) 
  }


  # for breast cancer there is an additional feature that is 
  # the weighted likelihood
  if(tumor_type == "breast"){
    methods <- c(methods, "weighted_catalog")
    sig_catalogs <- c(sig_catalogs, "cosmic_tissue")
    steps <- c(steps, "mss")
    signames <- c(signames, "Signature_3")
  }
    
  if(do_mva){
    methods <- c(methods, "mva")
    sig_catalogs <- c(sig_catalogs, "none")
    signames <- c(signames, "Signature_3")
    steps <- c(steps, "mss")
  }

  if(check_msi){
    methods <- c(methods, "median_catalog", "decompose")
    sig_catalogs <- c(sig_catalogs, "average", "cosmic_tissue")
    steps <- c(steps, "msi", "msi")
    signames <- c(signames, "Signature_msi")
  }

  for(imethod in 1:length(methods)){
    
    sig_catalog <- sig_catalogs[[imethod]]
    method <- methods[[imethod]]
    step <- steps[[imethod]]

    if(method == "median_catalog"){
      average_catalog <- all_catalogs[[tumor_type]]
      if(add_sig3) average_catalog$Signature_3_c1 <- all_catalogs$breast$Signature_3_c1
      message("Calculating likelihoods for each cluster")
    }else if(method == "cosine_simil"){
      message("Calculating cosine similarities")
    }else if(method == "decompose"){
      message("Calculating exposures with NNLS")
    }else if(method == "mva"){
      message("Predicting MVA score")
    }

    if(sig_catalog == "cosmic"){
      signatures <- cosmic_catalog
    }

    if(sig_catalog == "cosmic_tissue"){
      if(step == "mss"){
        signatures <- cosmic_catalog[, signames_per_tissue[[tumor_type]]]
      }else{
        signatures <- cosmic_catalog[, c(signames_per_tissue[[tumor_type]],
                                       signames_per_tissue[["msi"]],
                                       signames_per_tissue[["pole"]])]
      }
      if(add_sig3) signatures$Signature_3 <- cosmic_catalog$Signature_3
    }

    if(sig_catalog == "average"){
      if(step == "mss") signatures <- average_catalog
      else signatures <- cbind(average_catalog, 
                               all_catalogs[["msi"]],
                               all_catalogs[["pole"]])
    }

    # scale for the tri-nucleotide context
    if(data == "seqcap" | data == "seqcap_probe"){
      signatures <- weight_exome*signatures
    }else if(data == "msk"){
      signatures <- weight_msk*signatures
    }else if(data == "fo"){
      signatures <- weight_fo*signatures
    }
    signatures_norm <- t(t(signatures)/colSums(signatures))
    colnames(signatures_norm) <- colnames(signatures)
    signatures <- signatures_norm 
    rm(signatures_norm)

    # set the column with total number of mutations if this column
    # doesn't exist
    if(sum(colnames(genomes) == "total_snvs") == 0)
      genomes$total_snvs <- rowSums(genomes[, 1:96])

    #calculate the likelihood/cos simil/decomposition
    if(method != "mva"){
      if(method == 'median_catalog'){

        # When likelihoods are calculated they are scaled with the number
        # of tumors in the public datasets that are in that specific cluster
        # with respect to which likelihoods are calculated. For WES data
        # this is the default for panel data there is the option of not 
        # having this weighting to make the model more robust 
        if(weight_cf){
          cluster_fractions_this <- cluster_fractions[[tumor_type]]
          if(check_msi & step == 'msi'){
            msi_count <- cluster_fractions$msi_tissue_weights[[tumor_type]]
            pole_count <- cluster_fractions$pole_tissue_weights[[tumor_type]]
            counts <- c(cluster_fractions_this, 
                        msi_count*cluster_fractions$msi/sum(cluster_fractions$msi), 
                        pole_count*cluster_fractions$pole/sum(cluster_fractions$pole))
            names(counts) <- c(names(cluster_fractions_this), 
                               names(cluster_fractions$msi), 
                               names(cluster_fractions$pole))
            cluster_fractions_this <- counts
          }
        }else{
          cluster_fractions_this <- NULL
        }    
        output <- match_to_catalog(genomes, 
                                   signatures,  
                                   method = method, 
                                   data = data, 
                                   cluster_fractions = cluster_fractions_this)
        
      }else{
        output <- match_to_catalog(genomes, 
                                   signatures,  
                                   method = method, 
                                   data = data)
      }

      if(step == "msi"){
        colnames(output) <- paste0(colnames(output), "_msi")        
      }
    }
    else{ # calculate the multivariate analysis score
      if(exists("merged_output")){
        output <- predict_mva(cbind(genomes, merged_output), 
                              signames[[imethod]], 
                              data, 
                              tumor_type, 
                              weight_cf)
      }
      else{ 
        output <- predict_mva(genomes, 
                              signames[[imethod]], 
                              data, 
                              tumor_type, 
                              weight_cf)
      }
    }
    
    # calculates the pass/fail boolean from MVA score or likelihood

    if(do_assign & 
       (method == "median_catalog" | method == "mva")){
      output_comb <- cbind(genomes, output)
      assignments <- assignment(output_comb, 
                                method = method, 
                                signame = signames[[imethod]],
                                data = data, 
                                tumor_type = tumor_type, 
                                weight_cf = weight_cf)
      output <- cbind(output, assignments)
    }


    if(!exists("merged_output")) merged_output <- output
    else merged_output <- cbind(merged_output, output)
  }

  merged_output <- cbind(genomes, merged_output)

  # lite format removes cosine similarities and NNLS + likelihood results for
  # signatures other than Signature 3 it keeps likelihoods with respect to
  # cluster averages, but combines the different clusters in different 
  # categories together
  if(lite_format){   
    lite <- lite_df(merged_output)
    output_file <- gsub(output_file, pattern = '\\.csv', replace = '_lite\\.csv')
    write.table(lite,
                output_file, 
                row.names = F,
                col.names = T,
                quote = F,
                sep = ",") 

    if(file.exists(output_file))
      message(paste0("SigMA output is in: ", output_file))

  }else{
    write.table(merged_output,
                output_file, 
                row.names = F,
                col.names = T,
                quote = F,
                sep = ",")
    if(file.exists(output_file))
      message(paste0("SigMA output is in: ", output_file))
  }

  return(output_file)

}