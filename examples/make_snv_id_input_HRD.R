devtools::load_all()

# make snv count matrix using a directory that contains vcf files
data_dir <- ''

genomes_matrix <- make_matrix(data_dir, file_type = 'vcf', ref_genome_name = 'hg19')
genomes <- conv_snv_matrix_to_df(genomes_matrix)

genome_file = 'example.csv'


write.table(genomes,
            genome_file,
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)

# below command adds mmej nhej columns to the snv matrix
add_mmej_nhej_id_counts(data_dir, input_matrix_file = genome_file)



