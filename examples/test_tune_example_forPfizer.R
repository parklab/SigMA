devtools::load_all()

tumor_type <- 'breast'
data <- 'fo' #replace with msk if you want to use the example dataset below
# cosmic_v2_inhouse is the catalog used by the initial package

catalog_name  <- 'cosmic_v3p2_inhouse'  
maf_percent <- NULL # replace with 0.001 for panels without matched normal

# for using your own dataset remove the 
# file_name <- PATH_TO_CSV_FILE_OF_SIGMA_INPUT

########## remove the lines below if you had defined the file_name above
file_name <- paste0('matrix_96dim_', data , '_',tumor_type,'_', catalog_name,'.csv')
data_file <- system.file("extdata/examples/msk_impact_2017_clinical_data_selected_Oncotree_Codes.tsv", package="SigMA")
df_muts <- read.delim(data_file)
data_file <- paste0('test_msk_Oncotree_Code_', tumor_type, '.maf')
write.table(df_muts[df_muts$Oncotree.Code %in% c("BRCA", "IDC"),], file = data_file, row.names = F, sep = '\t', quote = F)
m <- make_matrix(data_file, file_type='maf', ref_genome_name = 'hg19')
df <- conv_snv_matrix_to_df(m)
write.table(df, file_name, row.names = F, sep = ',', quote = F)
########## remove till here

# First run the built in model
output_file_built_in <- run(file_name,
                            data = data, # if it is a panel try msk
                            tumor_type = tumor_type,
			    catalog_name = catalog_name,
                            do_mva = F, # T when data = 'msk'
                            do_assign = F, # T when data = 'msk'
                            check_msi = T)




# remove MSI samples if you have other sources of hypermutations e.g.
# temozolomide treated patients etc you want to remove them as well
m <- read.csv(output_file_built_in)
lite <- lite_df(m)
m <- m[!(lite$categ %in% c("Signature_msi", "Signature_pole", 'MMRD', 'POLE')),]

m_bef <- read.csv(file_name)
file_name_no_msi <- gsub(file_name, pattern = '.csv', replace = '_nomsi.csv')
write.table(m[,1:dim(m_bef)[[2]]], file_name_no_msi, row.names = F, sep  = ',', quote = F)


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
			       catalog_name = catalog_name,
                               run_SigMA = T,
			       maf_percent = maf_percent,
                               remove_msi_pole = F) #since we removed these samples above otherwise has to be T


message('simulation finished, tuning the model')
message(simul_file)

maf_percent_tag <- 'matched_normal'
if(!is.null(maf_percent)) maf_percent_tag <- paste0('maf_percent_', maf_percent)

# Using the imulations tune a new Gradient Boosting Machine
tune_new_gbm(simul_file, 
             tumor_type = tumor_type, 
             data = data,
             run_SigMA = T,
	     catalog_name = catalog_name,
             rda_file = paste0('test_', tumor_type,'_', data, '_', catalog_name,'_', maf_percent_tag, '.rda'),
	     use_indels = T) # change to T after adding mmej and nhej in your input file see examples/make_snv_id_input_HRD.R

message('tuned')
load(paste0('test_', tumor_type ,'_', data, '_', catalog_name,'_', maf_percent_tag, '.rda'))

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
add_gbm_model(paste0('gbm_for_', tumor_type,'_', data, '_', catalog_name, '_maf_percent_', maf_percent_tag),
              tumor_type = tumor_type, 
              gbm_model = gbm_model, 
              cutoff  = cutoff, 
              cutoff_strict = cutoff_strict)

message('new model added')
# you can calculate a trinucleotide normalization, a 'fo' normalization already exists based on a single patient's panel
# bed_file <- system.file("extdata/examples/example_small_bedfile.bed", package="SigMA")
# norm96 <- get_trinuc_norm(bed_file)

# then imagine we received a new dataset which is in our example
# the same dataset we analyzed earlier with the built in model
output_new_tune <- run(file_name, 
                      data = paste0('gbm_for_', tumor_type, '_', data, '_', catalog_name, '_maf_percent_', maf_percent_tag),
                      tumor_type = tumor_type,
		      catalog_name = catalog_name,
                      do_mva = T, 
                      do_assign = T, 
                      check_msi = T,
                      custom = T, 
                      norm96 = weight_3Nfreq[[data]]) 
#                      norm96 = norm96) # you can replace_with norm96 using the bed file that defines the library coverage

