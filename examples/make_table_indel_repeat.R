args = commandArgs(trailingOnly=TRUE)
input <- args[[1]]
output <- args[[2]] 

df <- read.delim(input, header = F)
df <- unique(df[, 4:6])

t_ins <- table(df$V4[df$V5 == "INS"])
t_del <- table(df$V4[df$V5 == "DEL"])
tumors = unique(df$V4)

df_msi <- data.frame(nmsi_ins = c(unlist(t_ins[match(tumors, names(t_ins))])),
                     nmsi_del = c(unlist(t_del[match(tumors, names(t_del))])),
                     tumor = tumors)

df_msi$nmsi_ins[is.na(df_msi$nmsi_ins)] <- 0
df_msi$nmsi_del[is.na(df_msi$nmsi_del)] <- 0

write.table(df_msi, file = output,  row.names = F, sep = ',', quote = F)
