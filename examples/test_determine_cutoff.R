devtools::load_all()
data_dir <- system.file("extdata/matrices/matrices_96dim.rda", package="SigMA")
load(data_dir)

tumor_type <- 'breast'
remove_msi_pole <- F

m <- matrices_96dim[['tcga']][[tumor_type]]

write.table(m, 'tmp.csv', row.names = F, sep = ',', quote = F)
input_file <- 'tmp.csv'

# find the data parameter to be used in the run() function
# so that the median count of mutations in the samples that 
# will be analyzed is closest to the simulations used in the
# tuning of the model
data_val <- find_data_setting(input_file, 
  tumor_type,  
  remove_msi_pole = remove_msi_pole)

# OR if you know which data setting you will be using set it 
# to the preferred value by putting in the line below
# data_val <- 'seqcap' 


# to determine the optimal cutoff for the data setting
# that was calculated above, a new simulation that matches
# the average mutation count of the dataset at hand is 
# performed. This is a quick simulation which adjusts
# a set of already existing simulations that were generated
# specifically for some sequencing platforms. For a more
# reliable mutation simulations should be done by sampling
# the WGS data according to the library of the sequencing
# platform, contact authors for more information.
simul_file <- quick_simulation(input_file = input_file,
  tumor_type = tumor_type,
  data = data_val,
  remove_msi_pole = remove_msi_pole)

# algorithm is run on simulations
output_simul <- run(simul_file,
                    data = data_val,
                    tumor_type = tumor_type,
                    do_mva = T,
                    do_assign = T,
                    check_msi = F)

# using the simulated data the thresholds which correspond to 
# specific false positive rate or sensitivity values can be 
# obtained
cut_var <- 'fpr'
limits <- c(0.1, 0.05, 0.0149)

df <- read.csv(output_simul)
thresh <- get_threshold(df, limits, var = 'Signature_3_mva', cut_var = cut_var)

cutoffs <- thresh$cutoff
sens <- thresh$sen
fprs <- thresh$fpr
fdrs <- thresh$fdr

df_sen_fpr <- data.frame(cutoffs = cutoffs, sens = sens, fprs = fprs, fdrs = fdrs)
write.table(df_sen_fpr, 'df_sen_fpr_example.csv', row.names = F, sep = ',', quote = F)

# the algorithm is run on the actual dataset and 
output_file <- run('tmp.csv', 
                   data = data_val, 
                   do_mva = T, 
                   do_assign = F, 
                   check_msi = F, 
                   weight_cf = T) 

df <- read.csv(output_file)

message('\nassignments with limits')
df <- cbind(df, assignment(df,
                          method = 'mva',
                          cut_var = cut_var,
                          cutoffs = cutoffs,
                          limits = limits))

write.table(df, output_file, row.names = F, sep = ',', quote = F)

message(paste0('See pass_mva_', cut_var, '_X.XX columns for predictions with specific FPR thresholds in ', output_file))
for(i in 1:length(limits)){
  message(paste0('cutoff for ', round(limits[[i]], digit = 2), ' FPR is ', round(fprs[[i]], digit = 3), ' corresponding to sensitivity ', round(sens[[i]], digit = 2)))
}

if(cut_var == 'fpr' | cut_var == 'fdr'){
  strict_limit <- min(limits)
  loose_limit <- max(limits)
}else if(cut_var == 'sen'){
  strict_limit <- max(limits)
  loose_limit <- min(limits)
}else{
  stop('cut_var can be sen, fpr or fdr')
}

colnames(df)[colnames(df) == paste0('pass_mva_', cut_var, '_', round(strict_limit, digit = 2))] <- 'pass_mva_strict'
colnames(df)[colnames(df) == paste0('pass_mva_', cut_var, '_', round(loose_limit, digit = 2))] <- 'pass_mva'

df_lite <- lite_df(df)
write.table(df_lite, file = gsub(output_file, pattern = '.csv', replace = '_lite.csv'), row.names = F, sep = ',', quote = F)
message(paste0('\n lite file written to ', gsub(output_file, pattern = '.csv', replace = '_lite.csv')))