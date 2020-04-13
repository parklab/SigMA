# This example is for users who would like to tune their model 
# for a specific dataset where the average number of mutations
# differ significantly from the datasets for which the built in
# models exists. For those cases it is possible to determine a 
# threshold for an existing model (follow test_determine_cutoff.R
# but a new tune would perform better.

file_name <- 'matrix_96dim_tcga_brca.csv'
tumor_type <- 'breast'

devtools::load_all()
if(0){
m <- make_matrix('tcga_mc3_brca.maf', file_type='maf')
df <- conv_snv_matrix_to_df(m)


write.table(df, file_name, row.names = F, sep = ',', quote = F)

# First run the built in model
output_file_built_in <- run('matrix_96dim_tcga_brca.csv',
                            data = 'tcga_mc3',
                            tumor_type = tumor_type, 
                            do_mva = T,
                            do_assign = T,
                            check_msi = T)


#file_name_out <- 'matrix_96dim_tcga_brca_output_tumortype_breast_platform_TCGAMC3mutationcalls_cf0.csv'
output_file_built_in <- 'matrix_96dim_tcga_brca_output_tumortype_breast_platform_TCGAMC3mutationcalls_cf0.csv'

# then assuming that we do not have the model let's try to tune
# a new model for our data and we will then compare to the
# built in calculation

# First determine the closest 'data' setting, which
# indicates the sequencing platform, to use 

#remove MSI samples
m <- read.csv('matrix_96dim_tcga_brca_output_tumortype_breast_platform_TCGAMC3mutationcalls_cf0.csv')
lite <- lite_df(m)
m <- m[lite$categ != "Signature_msi" & lite$categ != "Signature_pole",]
write.table(m[,1:97], 'matrix_96dim_tcga_brca_msi_removed.csv', row.names = F, sep  = ',', quote = F)


data <- find_data_setting(input_file = 'matrix_96dim_tcga_brca_msi_removed.csv', 
                          tumor_type = tumor_type,
                          remove_msi_pole = F) # since we removed these samples above otherwise has to be T


# Once this is identified it allows us to get an 
# estimate of the total snv counts in Sig3+ and Sig3-
# samples identified at first order. To tune a new model
# that fits the SNV count in our dataset we first
# simulate a new dataset from WGS data for which 
# the sig3 is known from WGS analysis and is more
# reliable
simul_file <- quick_simulation('matrix_96dim_tcga_brca_msi_removed.csv', 
                               tumor_type = 'breast', 
                               data = data,
                               run_SigMA = T, 
                               remove_msi_pole = F) #since we removed these samples above otherwise has to be T

}

data <- 'tcga_mc3'

simul_file <- 'matrix_96dim_tcga_brca_msi_removedsimulation.csv'
simul_file_output <- run('matrix_96dim_tcga_brca_msi_removedsimulation.csv', 
                         data = data,
                         tumor_type = tumor_type,
                         do_mva = F,
                         do_assign = F)

message('tuning')
message(simul_file_output)
# Using the imulations tune a new Gradient Boosting Machine
gbm_model <- tune_new_gbm(simul_file_output, 
                          tumor_type = tumor_type, 
                          data = data,
                          run_SigMA = T, 
                          rda_file = 'test.rda')

message('tuned')
load('test.rda')

message('predicting')
# Using the new model predict the signature 3 scores
df_predict <- predict_prob(file = simul_file_output,
                           gbm_model = gbm_model)


# Get thresholds for given false positive rate
thresh <- get_threshold(df_predict, 0.1, var = 'prob', cut_var = 'fpr') # FPR < 0.1
thresh_strict <- get_threshold(df_predict, 0.05, var = 'prob', cut_var = 'fpr') #FPR < 0.05
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
    
weight_this <- get_trinuc_norm('../../../SeqCap_EZ_Exome/primary_targets/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed')

# then imagine we received a new dataset which is in our example
# the same dataset we analyzed earlier with the built in model
output_new_tune <- run('matrix_96dim_tcga_brca_msi_removed.csv', 
                      data = 'example_gbm',
                      tumor_type = 'breast',
                      do_mva = T, 
                      do_assign = T, 
                      check_msi = T,
                      custom = T,
                      norm96 = weight_this)

df_built_in <- read.csv(output_file_built_in)
df_new_tune <- read.csv(output_new_tune)

df_built_in$Signature_3_mva_new <- df_new_tune$Signature_3_mva[match(df_built_in$tumor, df_new_tune$tumor)]

ggplot(df_built_in, aes(x = Signature_3_mva, y = Signature_3_mva_new)) + geom_point()