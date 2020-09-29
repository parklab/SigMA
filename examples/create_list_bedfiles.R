args = commandArgs(trailingOnly=TRUE)
output_file <- args[[1]]

devtools::load_all()
file_path <- system.file("extdata/repeat_bed_files/",
                         package="SigMA")

files <- list.files(file_path)
files <- paste0(file_path, '/', files)
write.table(files, file = args[[1]], row.names = F, sep = '\t', quote = F, col.names = F)