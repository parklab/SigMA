#after devtools:load_all()

directory <- system.file("extdata/test", package = "LowRStat")
genomes_matrix <- make_matrix(directory, file_type = 'maf')

genomes <- conv_snv_matrix_to_df(genomes_matrix)

write.table(genomes, 
            'genomes_example.csv', 
            sep = ',', 
            row.names = F, 
            col.names = T , 
            quote = F)

#calculates likelihoods and cosine similarity
run(genome_file = 'genomes_example.csv', 
    output_file = 'output_example.csv',
    sig_catalog = 'default')

# plot heatmap for genomes likelihood or cosine similarity values
# only uses the samples that fall within the minimum and maximum 
# of snv_ranges
plot_heatmap('output_example.csv', 
             matching = 'l', 
             snv_ranges = c(5, 10, 50, 100),
             file_name = 'heat_map_example.jpg')