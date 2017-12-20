#' Randomly downsample the genomes given in the rows of a
#' dataframe 
#'
#' @param genome_file the file with a data frame that has the 
#' first 96 rows need to be the 96 dimensional SNV spectra
#' @param min_snv is the lower bound of snvs to be generated
#' in the downsampling
#' @param max_snv is the upper bound of snvs to be generated 
#' in the downsampling
#' @param niter is the number of different downsamplings to be
#' generated
#' @param output_dir is the location of outputs to be saved

downsample_random <- function(genome_file, 
                              min_snv,
                              max_snv, 
                              niter,
                              output_dir){

  genomes <- read.csv(genome_file)

  apply(genomes, 1, 
        function(x){
          .generate_from_96spec(as.numeric(x[1:96]), 
                                min_snv = min_snv,
                                max_snv = max_snv,
                                niter = niter, 
                                tumor = x[['tumor']],
                                output_dir = output_dir)
        })
}

.generate_from_96spec <- function(mutations, 
                                 min_snv,
                                 max_snv,
                                 niter, 
                                 tumor, 
                                 output_dir){
  total_mutations <- sum(mutations)
  max_snv <- min(max_snv, total_mutations)  
  max_mutations <- max(mutations)

  max_index <- (max_snv - min_snv + 1)*niter

  mc_snvs <- matrix(0, max_index, 96)
  index <- 1

  for(iter in 1:niter){
    for(isnv in min_snv:max_snv){
      while(sum(mc_snvs[index,]) < isnv){
        location <- as.integer(round(runif(1, 0.5, 96.5), 
                               digit = 0))
        while(mutations[[location]] <= 0){
          location <- as.integer(round(runif(1, 0.5, 96.5), 
                              digit = 0))
        }

        rand_prob <- runif(1, 0, max_mutations)
        if(rand_prob < mutations[[location]]){
          mc_snvs[index:max_index, location] <- mc_snvs[index:max_index, location] + 1
          index <- index + 1
          mutations[[location]] <- mutations[[location]] - 1
          max_mutations <- max(mutations)
        }
      }
    }
  }

  df <- data.frame(mc_snvs)
  colnames_snv <- read.table('snv_colnames.txt')$V1
  colnames(df)[1:96] <- colnames_snv
  write.table(df, 
              sprintf('%s/%s.csv', output_dir, tumor),
              sep = ',', 
              row.names = F, 
              col.names = T, 
              quote = F)
}
