# See https://github.com/parklab/SigMA/wiki/MMRD-input-file-format

devtools::load_all()

m <- make_matrix('example_muts_mmrd.maf', file_type = 'maf')
df <- conv_snv_matrix_to_df(m)

file_SBS = 'example_mmrd_SBS_matrix.csv'
file_overlap_repeat='counts_indel_repeat.bed.csv'
mut_bed_file='example_muts_mmrd.bed'

write.table(df, file_SBS, row.names = F, sep = ',', quote = F)
combine_SBS_indel(file_SBS, file_overlap_repeat, mut_bed_file)

# Comment out the lines below to add msisensor file too
# msisensor_file='example_msisensor_file.txt'
# combine_SBS_indel(file_SBS = file_SBS, 
#                    file_overlap_repeat = file_overlap_repeat, 
#                    mut_bed_file = mut_bed_file, 
#                    msisensor_file = msisensor_file)

