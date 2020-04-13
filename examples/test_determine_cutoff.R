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
data_val <- 'seqcap' 


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
fpr_limits <- c(0.1, 0.05, 0.0149)

df <- read.csv(output_simul)
thresh <- get_threshold(df, fpr_limits, var = 'Signature_3_mva', cut_var = cut_var)

cutoffs <- thresh$cutoff
sens <- thresh$sen
fprs <- thresh$fpr
fdrs <- thresh$fdr

# the algorithm is run on the actual dataset and 
output_file <- run('tmp.csv', 
                   data = data_val, 
                   do_mva = T, 
                   do_assign = F, 
                   check_msi = F, 
                   weight_cf = T) 

df <- read.csv(output_file)
message(paste0('file saved to ', output_file))

message('assignments with limits')
df <- cbind(df, assignment(df,
                          method = 'mva',
                          cut_var = cut_var,
                          cutoffs = cutoffs,
                          limits = fpr_limits))

write.table(df, output_file, row.names = F, sep = ',', quote = F)

message(paste0('See pass_mva_fpr_X.XX columns for predictions with specific FPR thresholds in ', output_file))
for(i in 1:length(fpr_limits)){
  message(paste0('cutoff for ', round(fpr_limits[[i]], digit = 2), ' FPR is ', round(cutoffs[[i]], digit = 3), ' corresponding to sensitivity ', round(sens[[i]], digit = 2)))
}