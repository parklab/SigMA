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
#' @param data determines the type of sequencing platform
#' used for dataset see list_data_options() and find_data_setting()
#' @param tumor_type see list_tumor_types() for available options
#' @param do_mva a boolean for whether multivariate analysis 
#' should be run, has_model() function can be used to check whether
#' there is an available MVA model for that data and tumor_type setting
#' see SigMA/examples/example_find_settings.R
#' @param do_assign when true a cutoff is applied on MVA score 
#' to convert these values into binary classification. The thresholds
#' correspond to positive rates of 0.1 and 0.01 for MVA score, 
#' the columns generated for these are pass_mva and pass_mva_strict,
#' respectively
#' @param check_msi is a boolean which determines whether the user
#' wants to identify micro-sattelite instable tumors which often
#' confound the HRD classification
#' @param weight_cf determines whether the likelihood calculation
#' will take into account the number of tumors in each cluster
#' when it is F the clusters get equal weights and when it's T
#' they are weighted according to the fraction of tumors in each
#' cluster
#' @param lite_format saves the output in a lite format when set
#' to true
#' @param add_sig3 should be set to T when the likelihood of 
#' Signature 3 is calculated for tumor types for which Signature 3 
#' was not discovered by NMF in their WGS data
#' @param norm96 the normalization for custom classifiers to weight
#' signatures taking into account differences in trinucleotide 
#' frequency in whole genome versus the the sequencing platform. 
#' When built in data setting are used nomr96 is not necessary.
#' @param readjust when set TRUE the cutoffs are readjusted to 
#' match the median SNV counts in the data which may differ from
#' the dataset used in tuning the MVA classifier. Set to F by
#' default in order to avoid readjusting for dataset with too 
#' few samples which may bias the prediction.  
#' @param return_df when set TRUE instead of saving the output 
#' to file it returns the data frame
#' @param input_df the input data frame, can be used instead of
#' genome file. Provide one or the other.
#'
#' @examples
#' run(genome_file = "input_genomes.csv", 
#'     data = "msk",
#'     tumor_type = "ovary")
#' run(genome_file = "input_genomes.csv", 
#'     data = "seqcap", 
#'     tumor_type = "bone_other")

run <- function(genome_file = NULL, 
                output_file = NULL,
                data = "msk",
                tumor_type = "breast",
		catalog_name = NULL,
                do_assign = T,
                do_mva = T,
                check_msi = F, 
                weight_cf = F,
                lite_format = F,
                add_sig3 = F,
                norm96 = NULL,
                custom = F,
                readjust = F,
                return_df = F,
                input_df = NULL,
                snv_cutoff = 5){

  if(is.null(catalog_name))
    stop('please provide a catalog_name argument options: ',  paste0(names(catalogs), sep = ' '))

  if(!(catalog_name %in% names(catalogs)))
    stop(paste0('catalog ', catalog_name, ' does not exist provide one of the following: ', paste0(names(catalogs), sep = ' ')))

  cosmic_catalog <- catalogs[[catalog_name]]

  # when readjust is set to FALSE the parameters below are not used

  cut_var = NULL
  limits = NULL
  cutoffs_recalculated = NULL
 
  # only if readjust is set to TRUE the data setting is controlled
  # against other options to pick the best option and overwritten
  best_model <- NULL
  below_cutoff <- NULL

  if(!check_msi & readjust){  
    best_model <- find_data_setting(input_file = genome_file,
                                    remove_msi_pole = F,
                                    tumor_type = tumor_type,
                                    input_df = input_df)

    if(best_model != data){
      warning(paste0('best data setting differs from the data setting 
              used so the setting was changed to: ', best_model))
      data <- best_model
    }
  }
 

  if((data == 'seqcap' | data == "seqcap_probe") &  do_mva & !weight_cf){
    warning('weight_cf was set to FALSE but the seqcap and seqcap_probe requires
    weight_cf to be TRUE, so this setting was change to TRUE')
    weight_cf = T
  }else if(data %in% c('wgs', 'wgs_pancan', 'tcga_mc3')){
    weight_cf = F
  }

  # There are trained MVA models for tumor_type "eso", "osteo", "ovary",
  # "panc_ad", "panc_en", "prost", "stomach", "uterus", "breast", "bladder" 
  # for other tumor types do_mva should be FALSE
  # Signature 3 likelihood and all other variables can still be calculated 
  # by setting add_sig3 = T

  if(do_mva & !custom & sum(tumor_type == names(gbm_models[[data]])) == 0){
    stop('No built-in MVA models for the tumor_type selected for
          targetted gene panels for "medullo" or "ewing" whole
          exome sequencing is available for others set do_mva to FALSE')
  }

  if(!add_sig3 & do_assign & ('Signature_3' %in% signames_per_tissue_per_catalog[[catalog_name]][[tumor_type]]) & sum(grepl('Signature_3_',colnames(all_catalogs[[tumor_type]]))) ==0 ){
     add_sig3 <- T
  }
  
  if(is.null(genome_file) & !is.null(input_df)){
    genomes <- input_df
  }
  else if(!is.null(genome_file) & is.null(input_df)){
    genomes <- read.csv(genome_file)
  }
  else{
    stop('genome_file or input file should point to the input dataset')
  }

  # remove genomes with no mutation 
  if(sum(is.na(rowSums(genomes[, 1:96]))) > 0){
    genomes <- genomes[!is.na(rowSums(genomes[, 1:96])), ]
  }
  if(sum(rowSums(genomes[, 1:96]) == 0) > 0){
    if(readjust)
      below_cutoff <- genomes[is.na(rowSums(genomes[,1:96])),]
    genomes <- genomes[which(rowSums(genomes[, 1:96]) > 0), ]
  }

  # lower cutoff on number mutations for SigMA
  if(do_assign | do_mva){
    if(!is.null(snv_cutoff)){
      if(readjust){
        if(is.null(below_cutoff))
          below_cutoff <- genomes[rowSums(genomes[,1:96]) < snv_cutoff,]
        else{
          below_cutoff <- rbind(below_cutoff, 
                                genomes[rowSums(genomes[,1:96]) < snv_cutoff,])
        }
      }
      genomes <- genomes[rowSums(genomes[,1:96]) >= snv_cutoff,]
    }
    else{
      if(data == "msk" | data == "op" | data == "fo"){
        if(tumor_type == "prost") snv_cutoff <- 4
        else if(tumor_type == "osteo" | tumor_type == "panc_en") snv_cutoff <- 3
        else snv_cutoff <- 5
      }else if(data %in% names(platform_names)){ # for exomes and wgs a larger lower cutoff is applied
        if(tumor_type == "bone_other" | tumor_type == "medullo") snv_cutoff <- 5
        else snv_cutoff <- 10
      }
    }
  }

  if(!custom){
    message(paste0("You are running SigMA for ", dim(genomes)[[1]], 
                   " ", tissue_names[[tumor_type]],  "(s) sequenced by ", 
                   platform_names[[data]]))
  }else{
    message(paste0("You are running SigMA for ", dim(genomes)[[1]], 
                   " ", tissue_names[[tumor_type]]))
  }

  # method names to calculate different features
  methods <- c("median_catalog", "cosine_simil", "decompose")

  # signature catalogs to be used for each method
  sig_catalogs <- c("average", "cosmic_tissue", "cosmic_tissue")

  # first the samples are assumed to be mss and later msi 
  # calculations are done
  steps <- c("mss", "mss", "mss")

  # signature to be identified from mss samples
  signames <- rep("Signature_3", 3) 


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
      if(data == "op_multisig" & tumor_type == "other"){
        average_catalog <- all_catalogs[['other_multisig']]
      }
      else{
        average_catalog <- all_catalogs[[tumor_type]]
      }
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
      if(data == 'op_multisig'){
        if(step == "mss" & method != "cosine_simil"){
          signatures <- cosmic_catalog[, signames_per_tissue_per_catalog[[catalog_name]][[tumor_type]]]
        }else{
          signatures <- cosmic_catalog[, unique(c(signames_per_tissue_per_catalog[[catalog_name]][[tumor_type]],
                                         signames_per_tissue_per_catalog[[catalog_name]][["msi"]],
                                         signames_per_tissue_per_catalog[[catalog_name]][["pole"]]))]
        }						
      }
      else{
        if(step == "mss" & method != "cosine_simil"){
          signatures <- cosmic_catalog[, signames_per_tissue_per_catalog[[catalog_name]][[tumor_type]]]
        }else{
          signatures <- cosmic_catalog[, unique(c(signames_per_tissue_per_catalog[[catalog_name]][[tumor_type]],
                                         signames_per_tissue_per_catalog[[catalog_name]][["msi_extra"]],
                                         signames_per_tissue_per_catalog[[catalog_name]][["msi"]],
                                         signames_per_tissue_per_catalog[[catalog_name]][["pole"]]))]
        }
        if(add_sig3) signatures$Signature_3 <- cosmic_catalog$Signature_3
      }
    }

    if(sig_catalog == "average"){
      signatures_msi <- all_catalogs[["msi"]]

      if(step == "mss") signatures <- average_catalog
      else signatures <- cbind(average_catalog, 
                               signatures_msi,
                               all_catalogs[["pole"]])
      signatures <- unique(signatures)
    }

    # scale for the tri-nucleotide context
    if(data != "wgs" & data != "wgs_pancan"){
      if(!is.null(norm96)){ 
        signatures <- norm96*signatures
      }
      else{
        if(data %in% names(weight_3Nfreq)){
          signatures <- unlist(weight_3Nfreq[data]) * signatures
        }else if(data == "op_multisig"){
	  data_norm <- 'op'
	  signatures <- unlist(weight_3Nfreq[data_norm]) * signatures
	}else{
          stop('choose a built-in model or provide trinucleotide normalization')
        }
      }
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
      output <- predict_mva(cbind(genomes, merged_output), 
                            signames[[imethod]], 
                            data, 
                            tumor_type, 
                            weight_cf,
                            custom)
    }

    # if check_msi is set to FALSE meaning that user has provided only the
    # mismatch repair proficient samples then the adjustment is done initially
    if(method == "mva" & readjust){
      # if readjust is set to TRUE the cutoffs are determined based on average
      # SNV cuts. 1. MSI and POLE-exo mutated sampels are removed. 2. The best model is 
      # determined, 3. cutoffs are recalculated if the average count of this dataset
      # differs from the data used to tune the models 
      combined <- cbind(genomes, merged_output, output)
      if(check_msi){ 
        inds <- which((merged_output$Signature_msi_ml > 0.99 | merged_output$Signature_pole_msi > 0.99) 
                       & merged_output$total_snvs > median(merged_output$total_snvs))
        
        if(length(inds) > 0) df_no_msi_pole <- combined[-inds,]
        else df_no_msi_pole <- combined
      }
      else{
        df_no_msi_pole <- combined
      }
      if(is.null(best_model)){
        best_model <- find_data_setting(input_file = NULL,
                                        input_df = df_no_msi_pole,
                                        remove_msi_pole = F,
                                        tumor_type = tumor_type)
        warning(paste0('Changed to the best data setting is:', best_model))
      }
      if(best_model != data){
        if(exists('merged_output')) rm('merged_output')
        if(exists('output')) rm('output')
        SigMA_output <- run(genome_file = genome_file, 
            output_file = output_file, 
            data = best_model, 
            tumor_type = tumor_type, 
            do_assign = do_assign, 
            do_mva = do_mva,
            check_msi = check_msi, 
            weight_cf = weight_cf,
            lite_format = lite_format,
            add_sig3 = add_sig3, 
            norm96 = norm96,
            custom = custom,
            readjust = readjust,
            return_df = return_df,
            input_df = input_df)
        return(SigMA_output)
      }
      else{
        adjusted_cutoff <- adjust_cutoff(df_no_msi_pole, data, tumor_type, below_cutoff)
        if(!is.null(adjusted_cutoff)){
          custom <- T
          cut_var <- adjusted_cutoff$cut_var
          limits <- adjusted_cutoff$limits
          cutoffs_recalculated <- adjusted_cutoff$cutoffs
        }
        else{
          simul_df <- quick_simulation(input_file = NULL, 
                                       input_df = df_no_msi_pole,
                                       tumor_type = tumor_type,
                                       data = best_model,
                                       remove_msi_pole = F, 
                                       return_df = T,
                                       run_SigMA = F,
                                       below_cutoff = below_cutoff)

          # algorithm is run on simulations
          output_simul_df <- run(genome_file = NULL,
                                input_df =simul_df,
                                data = best_model,
                                tumor_type = tumor_type,
                                do_mva = T,
                                do_assign = T,
                                check_msi = F, 
                                return_df = T)

          cut_var <- 'fpr'
          limits <- c(0.1, 0.01)
       
          thresh <- get_threshold(output_simul_df,
                                limits, var = 'Signature_3_mva', cut_var = cut_var)
       
          cutoffs_recalculated <- thresh$cutoff
          warning('adjusted cutoff was not available because the SNV count 
                  differs from the simulations used to tune the
                  classifier. Results are likely to be subobtimal try retuning see
                  SigMA/examples/test_tune_new.R for a tutorial')
        }
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
                                weight_cf = weight_cf, 
                                custom_model = custom,
                                cut_var = cut_var, 
                                limits = limits,
                                cutoffs_custom = cutoffs_recalculated)
      output <- cbind(output, assignments)
    }

    if(!exists("merged_output")) merged_output <- output
    else merged_output <- cbind(merged_output, output)

    if(!('rat_sig3' %in% colnames(merged_output))){
      if('exp_sig3' %in% colnames(merged_output)){
        merged_output$rat_sig3 <- merged_output$exp_sig3/genomes$total_snvs
      }
    }

    if(!("Signature_3_ml" %in% colnames(merged_output))){
      inds <- na.omit(match(paste0('Signature_3_c', 1:10, '_ml'), colnames(merged_output)))
      if(length(inds) > 1)
        merged_output$Signature_3_ml <- rowSums(merged_output[, inds])
      if(length(inds) == 1)
        merged_output$Signature_3_ml <- merged_output$Signature_3_c1_ml
    }
  }
  
  merged_output <- cbind(genomes, merged_output)
  # lite format removes cosine similarities and NNLS + likelihood results for
  # signatures other than Signature 3 it keeps likelihoods with respect to
  # cluster averages, but combines the different clusters in different 
  # categories together

  if(return_df) return(merged_output)
  else{
    if(!is.null(best_model)) data = best_model
    if(data %in% names(platform_names)){
      platform_name = gsub(platform_names[[data]], 
                         pattern = " ", replace = "")
    }
    else{
      platform_name = paste0('Custom_', data)
    }

    if(is.null(output_file) & !is.null(genome_file)){
      output_file = gsub(genome_file,
                         pattern = ".csv",
                         replacement = paste0("_output_tumortype_",
                                              tumor_type,
                                              "_platform_",
                                              platform_name,
                                              "_cf", as.integer(weight_cf),
                                              ".csv"))
    }else if(is.null(output_file) & is.null(genome_file)){
      output_file = paste0("SigMA_output_tumor_type_",
                           tumor_type,
                           "_platform_",
                           platform_name,
                          "_cf", as.integer(weight_cf), 
                          ".csv")
    }
                     
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
}