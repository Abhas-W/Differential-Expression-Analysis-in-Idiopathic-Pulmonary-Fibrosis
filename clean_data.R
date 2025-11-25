data <- read.table('GSE150910_gene-level_count_file.csv',
                   header = TRUE, sep = ',', row.names = 1)

samples <- colnames(data)
genes <- rownames(data)
classes <- c()

for (sample in samples) {
  tmp <- unlist(strsplit(sample, "_"))
  classes <- append(classes, tmp[1])
}

# strsplit splits into two strings using '_' as a marker and 
# unlist converts it into a vector
# then next line takes the 1st element of that vector
# and append it to 'classes' variable

filtered_data <- data[, classes %in% c("control","ipf")]
filtered_classes <- classes[classes %in% c("control","ipf")]
# classes %in% c("control","ipf") this returns column names 
# 'control' and 'ipf' as TRUE and else as FALSE

saveRDS(filtered_data, file = "rds_objects/filtered_data.RDS")
saveRDS(filtered_classes, file = "rds_objects/filtered_classes.RDS")