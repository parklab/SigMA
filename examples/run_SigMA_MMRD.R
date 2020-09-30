# This is an example from TCGA-COAD dataset.

devtools::load_all()

input_file <- 'example_mmrd_SBS_matrix.csv'

output_SigMA <- run(genome_file = 'example_mmrd_SBS_matrix.csv', 
                     data = 'tcga_mc3',
                     tumor_type = 'crc', 
                     do_mva = F, # can also be set to T for Sig3 but not necessary for MMRD detection
                     check_msi = T)
predict_mmrd(input_file = output_SigMA,
             data = 'tcga_mc3')

