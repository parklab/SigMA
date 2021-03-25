# This example is for users who would like to tune their model 
# for a specific dataset where the average number of mutations
# differ significantly from the datasets for which the built in
# models exists. For those cases it is possible to determine a 
# threshold for an existing model (follow test_determine_cutoff.R
# but a new tune would perform better.


devtools::load_all() #'~/Dropbox/AstraZeneca/MC/master/update_tune/SigMA')

df_norm <- read.delim('norm96_SNUH.txt')

tumor_type <- 'breast'
maf_percent <- 0.001

file_name <- 'snuh_genome_1e_5.csv'

output_new_tune <- run(file_name,
                      data = 'snuh_0.001',
                      tumor_type = tumor_type,
                      do_mva = T, 
                      do_assign = T, 
                      check_msi = T,
                      custom = T,
                      snv_cutoff = 3,
                      norm96 = df_norm$norm96)

