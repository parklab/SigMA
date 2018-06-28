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
                do_assign = F,
                data = "msk",
                tissue = "breast",
                gbm = F){
  
  genomes <- read.csv(genome_file)

  if(sum(rowSums(genomes[, 1:96]) == 0) > 0){
    genomes <- genomes[which(rowSums(genomes[, 1:96]) > 0), ]
  }

  if(do_assign){
    if(data == "msk"){
      if(tissue == "prost")
        genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 4), ]
      else if(tissue == "osteo")
        genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 3), ]
      else
        genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 5), ]
      
    }else{
      genomes <- genomes[which(rowSums(genomes[, 1:96]) >= 10), ]
    }
  }

  if(method == 'all'){
    methods <- c('median_catalog', 'cosine_simil', 'decompose')
    sig_catalogs <- c('average', 'cosmic', 'cosmic_tissue')
    steps <- c('mss', 'mss', 'mss')
    if(tissue == "breast"){
      methods <- c(methods, 'weighted_catalog')
      sig_catalogs <- c(sig_catalogs, 'cosmic_tissue')
      steps <- c(steps, 'mss')
    }
    if(tissue == "crc"){
      signames <- c('Signature_N1', rep('Signature_N1', 2), 
                    'Signature_N1', 'Signature_N1')
    }
    else if(tissue == "gbm"){
      signames <- c('Signature_8', rep('Signature_8', 2), 
                    'Signature_8', 'Signature_8')
    }
    else{
      signames <- c('Signature_3', rep('Signature_3', 2), 
                    'Signature_3', 'Signature_3')
    }
    
    if(tissue == "breast" & gbm){
      methods <- c(methods, 'gbm')
      sig_catalogs <- c(sig_catalogs, 'none')
      signames <- c(signames, 'Signature_3')
      steps <- c(steps, 'mss')
    }

    methods <- c(methods, 'median_catalog', 'decompose')
    sig_catalogs <- c(sig_catalogs, 'average', 'cosmic_tissue')
    steps <- c(steps, 'msi', 'msi')
    signames <- c(signames, 'Signature_msi')

#    if(gbm){
#      methods <- c(methods, "gbm")
#      sig_catalogs <- c(sig_catalogs, 'none')
#      signames <- c(signames, 'Signature_msi')
#      steps <- c(steps, 'msi')
#    }
  }
  else if(method == 'custom'){
    methods <- c('custom')
    if(is.null(sig_catalog))
      sig_catalogs <- c('custom')
    else
      sig_catalogs <- c(sig_catalog)
  }
  else
    stop('method can be all or custom')


  for(imethod in 1:length(methods)){
    sig_catalog <- sig_catalogs[[imethod]]
    method <- methods[[imethod]]
    step <- steps[[imethod]]
    print(method)

    if(method == 'median_catalog'){
      average_catalog <- all_catalogs[[tissue]]
    }
    if(method == 'custom'){
      if(is.null(custom_sig_df)) 
        stop('custom signature data frame is empty')
      else custom_sig_df <- custom_sig_df
    }
    if(sig_catalog == "cosmic"){
      signatures <- cosmic_catalog
    }
    if(sig_catalog == "cosmic_tissue"){
      if(step == "mss"){
        signatures <- cosmic_catalog[, signames_per_tissue[[tissue]]]
      }else{
        signatures <- cosmic_catalog[, c(signames_per_tissue[[tissue]],
                                       signames_per_tissue[['msi']],
                                       signames_per_tissue[['pole']])]
      }
    }

    print(step)
    # custom overwrites all the options above and uses the 
    # user defined input file 
    if(sig_catalog == "average"){
      if(step == "mss") signatures <- average_catalog
      else signatures <- cbind(average_catalog, 
                               all_catalogs[['msi']],
                               all_catalogs[['pole']])
    }

    # scale for the tri-nucleotide context
    if(data == "seqcap"){
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
    if(sum(colnames(genomes) == 'total_snvs') == 0)
      genomes$total_snvs <- rowSums(genomes[, 1:96])

    #calculate the likelihood/cos simil/gbm prob
    if(method != 'gbm'){
      output <- match_to_catalog(genomes, 
                                 signatures,  
                                 method = method, 
                                 data = data)
      if(step == "msi") 
        colnames(output) <- paste0(colnames(output), '_msi')
    }
    else{
      if(exists('merged_output'))
        output <- get_gbm_prediction(cbind(genomes, merged_output), signames[[imethod]], data, step)
      else 
        output <- get_gbm_prediction(genomes, signames[[imethod]], data, step)
    }
    
    # calculates the pass/fail boolean based on the tune
    if(do_assign & 
       (method == "median_catalog" | method == "gbm")){
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
      print('assign')
      assignments <- assignment(output_comb, 
                                method = method, 
                                signame = signames[[imethod]],
                                data = data)
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