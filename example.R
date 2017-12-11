#after devtools:load_all()

###############Running on maf or vcf files###################
directory <- system.file("extdata/test", package = "LowRStat")
#for vcf replace 'maf' with 'vcf'
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


################Generating MC with single signature#############
# generate Signature 3 samples (more settings are available)
generate_single_sig_mc('Signature_3')
run(genome_file = 'mc_Signature_3_snvs.csv', 
    output_file = 'mc_Signature_3_output.csv')
# generate Signature 5 samples
generate_single_sig_mc('Signature_5')
run(genome_file = 'mc_Signature_5_snvs.csv', 
    output_file = 'mc_Signature_5_output.csv')


###################Sensitivity and false positive rate############
# takes as input two files, one that has the signature of interest
# as the truth, e.g. a MC of signature 3 samples, or BRCA -/- tumors
# and the second one to be used for calculating the false positive 
# rate e.g. a MC of signature 5 samples, or for BRCA -/- samples 
# a highly balanced tumor group which is very unlikely to have HR
# deficiency can be used. Then the signame1 is the signature to be
# used in identifying the positives and signame2 is the truth of 
# the sample which is used to calculate false negatives
plot_sn_fp(file1 = 'mc_Signature_3_output.csv', 
           file2 = 'mc_Signature_5_output.csv',
           signame1 = 'Signature_3',
           signame2 = 'Signature_5')
plot_sn_fp(file1 = 'mc_Signature_3_output.csv', 
           file2 = 'mc_Signature_5_output.csv',
           signame1 = 'Signature_3',
           signame2 = 'Signature_5',
           do_matching = FALSE,
           dependence = 'cutoff')

##########Generating MC based on exposure exposure matrix##########


##########Calculating likelihood of NMF exposures #################
