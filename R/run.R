#' Loads signatures and genome snv matrices and runs the matching
#'
#' @param genome_file a file with snv spectra info can be created
#' from vcf file using @make_genome_matrix() function:
#' see ?make_genome_matrix
#' @param output_file the output file name
#' @param sig_catalog an array of 'cosmic' 'pcawg' 'default
#' 'custom', 'default option uses the intersection of cosmic
#' and pcawg catalogs adding the subsignatures from pcawg of 
#' existing cosmic signatures, e.g. Signature 7a, 7b, etc, but
#' the new signatures introduced by PCAWG catalog are not used
#' @param custom_sig_file a file with signature distributions
#' @param do_assign produces a boolean of whether the tumor 
#' passes the test based on the measures obtained from the 
#' methods
#' used 
#'
#' @examples
#' run(genome_file = 'input_genomes.csv', 
#'     sig_catalog = c('cosmic', 'pcawg'))
#' run(genome_file = 'input_genomes.csv', 
#'     sig_catalog = 'custom', 
#'     custom_sig_file = 'signatures_me.csv')
#' run(genome_file = 'input_genomes.csv',
#'     sig_catalog = c('pcawg')
#'     rm_sigs = 'Signature_40')

run <- function(genome_file, 
                output_file = NULL,
                method = 'all',
                sig_catalog = 'default', 
                custom_sig_df = NULL,
                custom_tune_file = NULL,
                use_weight = F, 
                do_assign = F,
                exome = F,
                gbm = F){
  
  genomes <- read.csv(genome_file)

  if(sum(rowSums(genomes[, 1:96]) == 0) > 0){
    print('Removing rows with 0 mutations')
    genomes <- genomes[which(rowSums(genomes[, 1:96]) > 0), ]
  }

  if(do_assign){
    genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 3), ]
    print('assign is true samples with more than 3 SNV are used')
  }

  if(method == 'weighted_catalog'){
    methods <- c('weighted_catalog')
    sig_catalogs <- c('cosmic')
    signames <- c('Signature_3')
    use_weight_vec <- c(T)
  }
  else if(method == 'median_catalog'){
    methods <- c('median_catalog')
    sig_catalogs <- c('custom')
    signames <- c('Signature_3_c1')
    use_weight_vec <- c(F)
  }
  else if(method == 'cosine_simil'){
    methods <- c('cosine_simil')
    sig_catalogs <- c('cosmic')
    signames <- c('Signature_3')
    use_weight_vec <- c(F)
  }
  else if(method == 'gbm'){
    methods <- c('gbm')
    sig_catalogs <- c('none')
    signames <- c('Signature_3')
    use_weight_vec <- c(F)
  }
  else if(method == 'all'){
    methods <- c('median_catalog', 'weighted_catalog', 'cosine_simil')
    sig_catalogs <- c('custom', 'cosmic', 'cosmic', 'cosmic_breast')
    signames <- c('Signature_3_c1', 'Signature_3', 'Signature_3', 'Signature_3')
    use_weight_vec <- c(F, T, F)
    if(exome){
      methods <- c(methods, 'decompose')
      sig_catalogs <- c(sig_catalogs, 'cosmic_breast')
      signames <- c(signames, 'Signature_3')
      use_weight_vec <- c(use_weight_vec, F)
    }
    if(gbm){
      methods <- c(methods, 'gbm')
      sig_catalogs <- c(sig_catalogs, 'none')
      signames <- c(signames, 'Signature_3')
      use_weight_vec <- c(use_weight_vec, F)
    }
  }
  else if(method == 'all' & !exome){
    methods <- c('median_catalog', 'weighted_catalog', 'cosine_simil')
    sig_catalogs <- c('custom', 'cosmic', 'cosmic')
    signames <- c('Signature_3_c1', 'Signature_3', 'Signature_3')
    use_weight_vec <- c(F, T, F)
  }
  else if(method == 'custom'){
    methods <- c('custom')
    if(is.null(sig_catalog))
      sig_catalogs <- c('custom')
    else
      sig_catalogs <- c(sig_catalog)
    use_weight_vec <- c(use_weight)
  }
  else
    stop('method can be all, median_catalog, weighted_catalog, 
         cosine_simil or custom')


  for(imethod in 1:length(methods)){
    sig_catalog <- sig_catalogs[[imethod]]
    method <- methods[[imethod]]
    use_weight <- use_weight_vec[[imethod]]

    if(method == 'median_catalog'){
      if(exome) custom_sig_df <- median_catalog_for_exome
      else custom_sig_df <- median_catalog_720bc
    }
    if(method == 'custom'){
      if(is.null(custom_sig_df)) 
        stop('custom signature data frame is empty')
      else custom_sig_df <- custom_sig_df
    }
    if(sig_catalog == "cosmic")
      signatures <- cosmic_catalog
    if(sig_catalog == "cosmic_breast")
      signatures <- cosmic_catalog_breast

    # custom overwrites all the options above and uses the 
    # user defined input file 
    if(sig_catalog == "custom"){
      signatures <- custom_sig_df
    }

    # set the tune based on the method, the tunes are read from
    # sysdata.Rda file if one of the default methods are used if
    # method is custom then a file needs to be provided
    if(method == 'median_catalog')
      tune_df <- tune_median_catalog
    else if(method == 'weighted_catalog')
      tune_df <- tune_weighted_catalog
    else if(method == 'cosine_simil')
      tune_df <- tune_cosine_simil
    else if(method == 'gbm')
      tune_df <- tune_gbm
    else if(method == 'custom'){
      if(is.null(custom_tune_file)) 
        stop('method set to custom but tune is not provided
              set do_assign to F to run without a tune')
      tune_df <- read.csv(custom_tune_file)
    }

    # set the column with total number of mutations if this column
    # doesn't exist
    if(sum(colnames(genomes) == 'total_snvs') == 0)
      genomes$total_snvs <- rowSums(genomes[, 1:96])

    #calculate the likelihood/cos simil/gbm prob
    if(method != 'gbm') 
      output <- match_to_catalog(genomes, 
                                 signatures,  
                                 method = method)

    else{
      if(exists('merged_output'))
        output <- get_gbm_prediction(cbind(genomes, merged_output), signames[[imethod]])
      else 
        output <- get_gbm_prediction(genomes, signames[[imethod]])
    }
    # calculates the pass/fail boolean based on the tune
    if(do_assign){
      if(!is.null(output_file)){
        output_comb <- cbind(genomes, output)
        write.table(output_comb, 
                     sprintf('%s.csv',
                             gsub(output_file, 
                                  pattern = '.csv',
                                  replace = paste0(method, '.csv'))), 
                     sep = ',', 
                     row.names = F, 
                     col.names = T, 
                     quote = F)
        rm(output_comb)
      }else 
        stop('give an output file name')

      assignments <- assignment(sprintf('%s.csv',
                                        gsub(output_file,
                                             pattern = '.csv',
                                             replace = paste0(method, '.csv'))),
                                method = method, 
                                tune = tune_df, 
                                signame = signames[[imethod]])
      output <- cbind(output, assignments)
    }


    if(!exists('merged_output')) merged_output <- output
    else merged_output <- cbind(merged_output, output)

  }

  merged_output <- cbind(genomes, merged_output)


  if(!is.null(output_file)) 
    write.table(merged_output,
                output_file, 
                row.names = F,
                col.names = T,
                quote = F,
                sep = ',')
  return(merged_output)

}