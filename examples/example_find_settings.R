devtools::load_all()

data_file <- system.file("extdata/examples/test_mutations_50sample.maf", package="SigMA")

genomes_matrix <- make_matrix(data_file, file_type = 'maf')
genomes <- conv_snv_matrix_to_df(genomes_matrix)

genome_file = 'example_maf.csv'

write.table(genomes,
            genome_file,
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)

# Print out the available tumor_type settings
list_tumor_types()

# Print out the available data settings
list_data_options()

# If you have a selection that matches your dataset 
# lets say sequencing platform was MSK-Impact panels
# and the tumor_type is breast cancer, then choose
# tumor_type = 'breast' and data = 'msk'. You
# can check the additional statistics to confirm
# your selection
info_stat(tumor_type = 'breast', data = 'msk')

total_snvs <- rowSums(genomes[,1:96])
median_counts <- median(total_snvs[total_snvs > 0])
sd_counts <- sd(total_snvs[total_snvs > 0])
median_counts <- median(total_snvs[total_snvs < 5*sd_counts])

paste0('median total snvs in data:', median_counts)
paste0('SD total snvs in data:', median_counts)

# If the above is not informative enough you can find the best 
# data setting using find_data_setting() function, change
# remove_msi_pole to T if you think your dataset may contain 
# mismatch repair deficient samples.
best_data_setting <- find_data_setting(input_file = genome_file, 
                  tumor_type = 'breast', 
                  remove_msi_pole = F)

best_data_setting

# you can check whether and MVA model is available for this
# tumor type for that particular data setting
has_model(data = best_data_setting, tumor_type = 'breast')

# then this can be used to run the code
run(genome_file = genome_file,
    data = best_data_settting, 
    readjust = T, # T allows readjustment of the cutoff if the SNV counts are still different
    tumor_type = 'breast')

