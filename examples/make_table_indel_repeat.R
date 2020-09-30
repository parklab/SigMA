args = commandArgs(trailingOnly=TRUE)
input <- args[[1]]
output <- args[[2]] 

df <- read.delim(input, header = F)
df <- unique(df[, 4:6])

t_ins <- table(df$V4[df$V5 == "INS"])
t_del <- table(df$V4[df$V5 == "DEL"])

df_nmsi <- data.frame(nmsi_ins = c(unlist(t_ins)),
                       nmsi_del = c(unlist(t_del)),
                       tumor = names(t_ins))

write.table(df_nmsi, file = output,  row.names = F, sep = ',', quote = F)
