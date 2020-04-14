# This example is for users who would like to tune their model 
# for a specific dataset where the average number of mutations
# differ significantly from the datasets for which the built in
# models exists. For those cases it is possible to determine a 
# threshold for an existing model (follow test_determine_cutoff.R
# but a new tune would perform better.


devtools::load_all()

data_file <- system.file("extdata/examples/tcga_mc3_brca.maf", package="SigMA")

m <- make_matrix(data_file, file_type='maf')
df <- conv_snv_matrix_to_df(m)

file_name <- 'matrix_96dim_tcga_brca.csv'
tumor_type <- 'breast'

write.table(df, file_name, row.names = F, sep = ',', quote = F)

# First run the built in model
output_file_built_in <- run(file_name,
                            data = 'tcga_mc3',
                            tumor_type = tumor_type, 
                            do_mva = T,
                            do_assign = T,
                            check_msi = T)



# then assuming that we do not have the model let's try to tune
# a new model for our data and we will then compare to the
# built in calculation

# First determine the closest 'data' setting, which
# indicates the sequencing platform, to use 

#remove MSI samples
m <- read.csv(output_file_built_in)
lite <- lite_df(m)
m <- m[lite$categ != "Signature_msi" & lite$categ != "Signature_pole",]
write.table(m[,1:97], file_name, row.names = F, sep  = ',', quote = F)


data <- find_data_setting(input_file = file_name, 
                          tumor_type = tumor_type,
                          remove_msi_pole = F) # since we removed these samples above otherwise has to be T


# Once this is identified it allows us to get an 
# estimate of the total snv counts in Sig3+ and Sig3-
# samples identified at first order. To tune a new model
# that fits the SNV count in our dataset we first
# simulate a new dataset from WGS data for which 
# the sig3 is known from WGS analysis and is more
# reliable
simul_file <- quick_simulation(file_name, 
                               tumor_type = 'breast', 
                               data = data,
                               run_SigMA = T, 
                               remove_msi_pole = F) #since we removed these samples above otherwise has to be T

simul_file_output <- run(simul_file, 
                         data = data,
                         tumor_type = tumor_type,
                         do_mva = F,
                         do_assign = F)

message('tuning')
message(simul_file_output)
# Using the imulations tune a new Gradient Boosting Machine
tune_new_gbm(simul_file_output, 
             tumor_type = tumor_type, 
             data = data,
             run_SigMA = T, 
             rda_file = 'test.rda')

message('tuned')
load('test.rda')

df_predict <- read.csv(gsub(simul_file_output, pattern = '.csv', replace = '_predictions.csv'))
tumor_type <- 'breast'

# Get thresholds for given false positive rate
limits_fpr <- c(0.1, 0.05) 
cut_var <- 'fpr' # you can alternatively use 'sen' to select on specific sensitivity value


thresh <- get_threshold(df_predict, limits_fpr[[1]], var = 'prob', cut_var = cut_var) # FPR < 0.1
thresh_strict <- get_threshold(df_predict, limits_fpr[[2]], var = 'prob', cut_var = cut_var) #FPR < 0.05
cutoff <- thresh$cutoff 
cutoff_strict <- thresh_strict$cutoff
sen <- thresh$sen 
fpr <- thresh$fpr
sen_strict <- thresh_strict$sen
fpr_strict <- thresh_strict$fpr

                      
# then we add the new gbm model in the system files together with the
# cutoffs we determined so that these can be used in the future for
# new data sequenced with the same sequencing platform
add_gbm_model('example_gbm', 
              tumor_type = tumor_type, 
              gbm_model = gbm_model, 
              cutoff  = cutoff, 
              cutoff_strict = cutoff_strict)
    
# example bed file:
bed_file <- system.file("extdata/examples/seqcap_capture.bed", package="SigMA")

norm96 <- get_trinuc_norm(bed_file)

# then imagine we received a new dataset which is in our example
# the same dataset we analyzed earlier with the built in model
output_new_tune <- run('matrix_96dim_tcga_brca_msi_removed.csv', 
                      data = 'example_gbm',
                      tumor_type = 'breast',
                      do_mva = T, 
                      do_assign = T, 
                      check_msi = T,
                      custom = T,
                      norm96 = norm96)

df_built_in <- read.csv(output_file_built_in)
df_new_tune <- read.csv(output_new_tune)
df_new_tune <- lite_df(df_new_tune)

df_new_tune$built_in_mva <- df_new_tune$built_in_mva[match(df_new_tune$tumor, df_built_in_mva$tumor)]

plot <- ggplot(df_new_tune, aes(x = Signature_3_mva, y = built_in_mva)) + geom_point(aes(color = categ)) + theme_def
ggsave(plot, 'compare_old_tune_new_tune.pdf', width = 9.8/2.4, height = 9.8/2.4)
