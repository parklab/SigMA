# short example 2 simulated panels for testing the tool
devtools::load_all()

# data_dir can be replaced by the character containing 
# the directory defined by the user
data_dir <- system.file("extdata", package="SigMA")

genomes_matrix <- make_matrix(data_dir, file_type = 'vcf')
genomes <- conv_snv_matrix_to_df(genomes_matrix)

genome_file = 'example.csv'

write.table(genomes,
            genome_file,
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)


message(paste0('96-dimensional matrix is saved in ', genome_file))
message('Running SigMA')

run(genome_file, 
    data = "msk", 
    do_assign = T, 
    do_mva = T, 
    lite_format = T)

