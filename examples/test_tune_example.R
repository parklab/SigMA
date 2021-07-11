# This example is for users who would like to tune their model 
# for a specific dataset where the average number of mutations
# differ significantly from the datasets for which the built in
# models exists. For those cases it is possible to determine a 
# threshold for an existing model (follow test_determine_cutoff.R
# but a new tune would perform better.


devtools::load_all()

data_file <- system.file("extdata/examples/tcga_mc3_brca.maf", package="SigMA")

m <- make_matrix(data_file, file_type='maf', ref_genome_name = 'hg19')
df <- conv_snv_matrix_to_df(m)

file_name <- 'matrix_96dim_tcga_brca.csv'
tumor_type <- 'breast'

write.table(df, file_name, row.names = F, sep = ',', quote = F)

# First run the built in model
output_file_built_in <- run(file_name,
                            data = 'tcga_mc3', # if it is a panel try msk
                            tumor_type = tumor_type, 
                            do_mva = T,
                            do_assign = T,
                            check_msi = T)



# then assuming that we do not have the model let's try to tune
# a new model for our data and we will then compare to the
# built in calculation

# First determine the closest 'data' setting, which
# indicates the sequencing platform, to use 

# remove MSI samples if you have other sources of hypermutations e.g.
# temozolomide treated patients etc you want to remove them as well
m <- read.csv(output_file_built_in)
lite <- lite_df(m)
m <- m[!(lite$categ %in% c("Signature_msi", "Signature_pole", 'MMRD', 'POLE')),]
file_name_no_msi <- gsub(file_name, pattern = '.csv', replace = '_nomsi.csv')
write.table(m[,1:97], file_name_no_msi, row.names = F, sep  = ',', quote = F)

data <- find_data_setting(input_file = file_name_no_msi, 
                          tumor_type = tumor_type,
                          remove_msi_pole = F) # since we removed these samples above otherwise has to be T


# if the best data option for a panel turns out to be exome sequencing
# you may want to check the germline contamination in your mutation calls
message(paste0('Best data option is ', data))

# Once this is identified it allows us to get an 
# estimate of the total snv counts in Sig3+ and Sig3-
# samples identified at first order. To tune a new model
# that fits the SNV count in our dataset we first
# simulate a new dataset from WGS data for which 
# the sig3 is known from WGS analysis and is more
# reliable
simul_file <- quick_simulation(file_name_no_msi, 
                               tumor_type = tumor_type, 
                               data = data,
                               run_SigMA = T, 
                               remove_msi_pole = F) #since we removed these samples above otherwise has to be T


message('tuning')
message(simul_file)
# Using the imulations tune a new Gradient Boosting Machine
tune_new_gbm(simul_file, 
             tumor_type = tumor_type, 
             data = data,
             run_SigMA = T, 
             rda_file = 'test.rda')

message('tuned')
load('test.rda')

simul_file_output <- gsub(simul_file, pattern = '.csv', replace = '_predictions.csv')
df_predict <- read.csv(simul_file_output)

# Get thresholds for given false positive rate
limits_fpr <- c(0.1, 0.05) 
cut_var <- 'fpr' # you can alternatively use 'sen' to select on specific sensitivity value


thresh <- get_threshold(df_predict, limits_fpr, var = 'prob', cut_var = cut_var, signal = 'is_sig3') # FPR < 0.1
cutoff <- thresh$cutoff[[1]] 
cutoff_strict <- thresh$cutoff[[2]]

                      
# then we add the new gbm model in the system files together with the
# cutoffs we determined so that these can be used in the future for
# new data sequenced with the same sequencing platform
add_gbm_model('example_gbm', 
              tumor_type = tumor_type, 
              gbm_model = gbm_model, 
              cutoff  = cutoff, 
              cutoff_strict = cutoff_strict)
    
# test the following command and replace the bed file with that of your libraries:
bed_file <- system.file("extdata/examples/example_small_bedfile.bed", package="SigMA")
norm96 <- get_trinuc_norm(bed_file)

# then imagine we received a new dataset which is in our example
# the same dataset we analyzed earlier with the built in model
output_new_tune <- run(file_name, 
                      data = 'example_gbm',
                      tumor_type = tumor_type,
                      do_mva = T, 
                      do_assign = T, 
                      check_msi = T,
                      custom = T,  
                      norm96 = norm96) # replace_with norm96 using the bed file that defines the library coverage

