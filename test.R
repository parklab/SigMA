# short example with 2 simulated panels for testing the tool

# load SigMA
devtools::load_all()

# data_dir can be replaced by the character containing 
# the directory defined by the user
data_dir <- system.file("extdata", package="SigMA")

# prepare the input file from vcf files, maf files can be used
# by changing the file_type to 'maf'
# obtain a data.frame with 96 dimensional spectra
genomes_matrix <- make_matrix(data_dir, file_type = 'vcf')
genomes <- conv_snv_matrix_to_df(genomes_matrix)

# write down the data.frame in a csv file
genome_file = 'example.csv'
write.table(genomes,
            genome_file,
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)
message(paste0('96-dimensional matrix is saved in ', genome_file))

# run SigMA 
message('Running SigMA')
output_file <- run(genome_file, data = "msk", do_assign = T, do_mva = T)

# plot the summary of the output
plot_summary(output_file)
