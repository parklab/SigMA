# Example 50 MSK-IMPACT (Zehir et al. Nat Med (2017)) breast cancer
# panels for testing the tool on a single maf file input

devtools::load_all()

# for data that is a single maf file rather than a directory with
# multiple files, data_dir can be replaced by the character containing
# the path to the file
data_dir <- system.file("extdata/examples/test_mutations_50sample.maf", package="SigMA")

data_file <- 'test_mutations_50sample.maf'

genomes_matrix <- make_matrix(data_file, file_type = 'maf')
genomes <- conv_snv_matrix_to_df(genomes_matrix)

genome_file = 'example_maf.csv'

write.table(genomes,
            genome_file,
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)


message(paste0('96-dimensional matrix is saved in ', genome_file))
message('Running SigMA')

# this dataset contains one MSI sample so we can set check_msi to be TRUE
# this makes the algorithm run extra steps of NNLS and the run time 
# increases. If the MSI tumors are not present in the dataset, as they
# might be identified in other ways, it is not necessary to have this 
# setting TRUE

run(genome_file, 
    data = "msk", 
    do_assign = T, 
    do_mva = T, 
    lite_format = T, 
    check_msi = T)

