devtools::load_all()

data_dir <- system.file("extdata/vcf_examples/", package="SigMA")

m <- make_matrix(data_dir, ref_genome_name = 'hg19')
df <- conv_snv_matrix_to_df(m)

df_msi <- run_overlap_repeats(data_dir, output_file = NULL, file_type = 'vcf', return_df = T)
df <- cbind(df, df_msi)
write.table(df, file = 'input_for_SigMA_classification.csv', row.names = F, sep = ',', quote = F)

df <- run_classifier(df)
get_info(df, df$tumor[[1]])