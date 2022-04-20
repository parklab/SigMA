devtools::load_all()

# make snv count matrix using a directory that contains vcf files add the directory path below:
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
mmej_nhej_indels_from_vcfs(data_dir,
                           output_file = 'test_mmej_nhej_output.csv',
                           ref_genome_name = 'hg19',
                           min_size_mh = 2,
                           min_size_del = 5,
                           snv_matrix_file = genome_file)


# output is saved in test_mmej_nhej_output.csv