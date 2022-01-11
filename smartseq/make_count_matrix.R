# Assemble count matrix from processed Smartseq run
# Author: S Chiou

library("stringr")

dataset_name <- "KIR1"
dataset_path <- "/data/smartseq2/out/KIR1"
output_dir <- "/data/smartseq2/out"

# Split each subdirectory's name by underscores and extract second entry
setwd(dataset_path)
listed_folders <- list.files()
sample_ids <- as.data.frame(str_split_fixed(listed_folders, "_", n = 3))

# Concatenate subdirectory's extracted name to dataset_name
sample_ids <- paste(dataset_name, as.character(sample_ids$V2), sep = "_")

## create data frame from each subdirectory
missing_files <- c()
for (i in seq(1, length(listed_folders))) {
  tmp_file_path <- paste(
    dataset_path,
    listed_folders[i],
    "output_counts_transcript_id.txt",
    sep = "/"
  )

  tmp_file_size <- file.info(tmp_file_path, extra_cols = FALSE)[1]
  tmp_file_size <- as.numeric(tmp_file_size)

  if (tmp_file_size != 0) {
    tmp_df <- read.delim(
      tmp_file_path,
      sep = "\t",
      row.names = 1,
      header = F
    )

    if (i == 1) {
      df <- tmp_df
    } else {
      df <- cbind(df, tmp_df)
    }
  } else {
    # Empty sample
    print(sample_ids[i])
    missing_files <- c(missing_files, sample_ids[i])
  }

  rm(tmp_df, tmp_file_size)
  print(i)
}

colnames(df) <- sample_ids[!(sample_id %in% missing_files)]

## write file
setwd(output_dir)
write.table(
  df,
  file = paste(dataset_name, ".txt", sep = ""),
  sep = "\t",
  col.names = T, row.names = T,
  quote = F
)
