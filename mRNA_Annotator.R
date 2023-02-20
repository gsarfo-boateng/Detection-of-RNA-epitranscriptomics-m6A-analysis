# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
exp_data <- args[1] # file1 field in terminal input
gtf_data <- args[2] # file2 field in terminal input

# Define output file name
output_file <- paste0(tools::file_path_sans_ext(exp_data), ".csv")

# Define acceptable genes
acceptable_genes <- c("3'UTR", "5'UTR", "CDS")

# Open output file for writing
csv_file <- file(output_file, "w")
writer <- csv::writer(csv_file)

# Read experimental data
exp_data_info <- file(exp_data, "r")
headerinfo_old <- strsplit(readLines(exp_data_info, n = 1), "\t")[[1]]
headerinfo_updated <- gsub("\"", "", headerinfo_old[-length(headerinfo_old)])
columns_appended <- c(headerinfo_old[length(headerinfo_old)], "gene_name", "Overlap Type", "Overhang", "Overlap", "Gene Type Assign", "Transcript ID", "Gene Biotype")
headerinfo_updated <- paste(c(headerinfo_updated, columns_appended), collapse = ",")
csv_file <- file(output_file, "w")
writeLines(headerinfo_updated, csv_file)
writeLines("", csv_file) # Write an empty line

exp_data_tsv <- readLines(exp_data_info)[-1]

for (wt_lines in exp_data_tsv) {
  # Define columns
  exp_data_strand <- gsub("\"", "", strsplit(wt_lines, "\t")[[1]][4])
  chromosome <- gsub("\"", "", strsplit(wt_lines, "\t")[[1]][1])
  exp_gene_initial <- as.numeric(gsub("\"", "", strsplit(wt_lines, "\t")[[1]][2]))
  exp_gene_last <- as.numeric(gsub("\"", "", strsplit(wt_lines, "\t")[[1]][3]))
}

# Read GTF data
gtf_data_info <- file(gtf_data, "r")
gtf_data_tsv <- readLines(gtf_data_info)
gtf_data_tsv <- gtf_data_tsv[-(1:4)]

for (lines in gtf_data_tsv) {
  # Define columns
  gtf_data_strand <- strsplit(lines, "\t")[[1]][7]
  gtf_gene_initial <- as.numeric(strsplit(lines, "\t")[[1]][4])
  gtf_gene_last <- as.numeric(strsplit(lines, "\t")[[1]][5])
  relevant_seq <- strsplit(lines, "\t")[[1]][3]
  gtf_chromosome <- strsplit(lines, "\t")[[1]][1]
  transcript_id <- gsub("\"|;", "", strsplit(lines, "\t")[[1]][9])[4]
  gene_symbol_index <- grep("gene_biotype", strsplit(lines, "\t")[[1]][9])
  gene_biotype <- gsub(";|\"|'", "", strsplit(lines, "\t")[[1]][9][gene_symbol_index + 1])
  
  # Define width/depth of genes
  gtf_depth <- gtf_gene_last - gtf_gene_initial
  exp_depth <- exp_gene_last - exp_gene_initial
}

# case A
#   -------------  GTF
#---------         EXP

if (gtf_gene_initial >= exp_gene_initial &
    gtf_gene_last >= exp_gene_last &
    gtf_gene_initial <= exp_gene_last) {
  
  overhang <- as.numeric(gtf_gene_initial - exp_gene_initial) / as.numeric(exp_depth)
  overhang <- round(overhang, 2)
  overlap <- 1 - overhang
  
  cat(gtf_depth, "gtf_depth", exp_depth, "exp_depth\n")
  cat(gtf_gene_initial, gtf_gene_last, "GTF\n")
  cat(exp_gene_initial, exp_gene_last, "EXP\n")
  cat("Type A", overhang, overlap, "\n")
  
  if (overlap <= 0.5) {
    gene_type_assign <- "5 Prime"
  } else {
    gene_type_assign <- "CDS"
  }
  
  overlap_details <- list("Type A", overhang, overlap) #, gene_type_assign]
  
  updated_line <- character()
  gene_name_index <- which(strsplit(strsplit(lines, "\t")[[1]][9], " ")[[1]] == "gene_name")[1]
  gene_name <- strsplit(strsplit(lines, "\t")[[1]][9], " ")[[1]][gene_name_index + 1]
  gene_name <- gsub("\"|;", "", gene_name)
  for (items in strsplit(wt_lines, "\t")[[1]][-length(strsplit(wt_lines, "\t")[[1]])]) {
    items <- gsub("\"", "", items)
    updated_line <- c(updated_line, items)
  }
  gtf_info <- c(substr(strsplit(wt_lines, "\t")[[1]][length(strsplit(wt_lines, "\t")[[1]])], 1, nchar(strsplit(wt_lines, "\t")[[1]][length(strsplit(wt_lines, "\t")[[1]])]) - 1), gene_name)
  updated_line <- c(updated_line, gtf_info)
  updated_line <- paste(updated_line, overlap_details, str(relevant_seq), str(transcript_id), str(gene_biotype), sep = ",")
  updated_line <- gsub("[\\[|\\]|']", "", updated_line)
  cat(updated_line, "\n")
  cat("--------------------------------------------------------------------------------------------\n")
  
  writeLines(updated_line, csv_file)
  writeLines("", csv_file) # to add an empty line after each line
  
}

if ((gtf_gene_initial <= exp_gene_initial) & 
    (gtf_gene_last <= exp_gene_last) &
    (gtf_gene_last >= exp_gene_initial)) {
  
  overhang <- round((exp_gene_last - gtf_gene_last) / exp_depth, 2)
  overlap <- 1 - overhang
  cat(gtf_depth, "gtf_depth", exp_depth, "exp_depth", "\n")
  cat(gtf_gene_initial, gtf_gene_last, "GTF", "\n")
  cat(exp_gene_initial, exp_gene_last, "EXP", "\n")
  cat("Type B", overhang, overlap, "\n")
  
  if (overlap <= 0.5) {
    gene_type_assign <- "3 Prime"
  } else {
    gene_type_assign <- "CDS"
  }
  
  overlap_details <- c("Type B", overhang, overlap)
  
  updated_line <- c()
  gene_name_index <- which(strsplit(strsplit(lines, "\t")[9], " ")[[1]] == "gene_name")
  gene_name <- gsub("\"|;", "", strsplit(strsplit(lines, "\t")[9], " ")[[1]][gene_name_index + 1])
  for (items in strsplit(strsplit(wt_lines, "\t")[[1]][-9], "\"")[[1]][-length(strsplit(strsplit(wt_lines, "\t")[[1]][-9], "\"")[[1]])]) {
    items <- gsub("\"", "", items)
    updated_line <- c(updated_line, items)
  }
  
  gtf_info <- c(strsplit(wt_lines, "\t")[[1]][10], gene_name)
  updated_line <- c(updated_line, gtf_info, overlap_details, relevant_seq, transcript_id, gene_biotype)
  updated_line <- gsub("\\[|\\]|'", "", paste(updated_line, collapse = "\t"))
  cat(updated_line, "\n")
  cat("--------------------------------------------------------------------------------------------", "\n")
  
  writeLines(paste(updated_line, collapse = "\n"), csv_file)
  writeLines("", csv_file)
}

# case C
#------------      GTF
#   ------         EXP
if ((gtf_gene_initial <= exp_gene_initial) & 
    (gtf_gene_last >= exp_gene_last)) {
  overhang = (exp_gene_last - exp_gene_initial) / exp_depth
  overhang = round(overhang, 2)
  overlap = 1 - overhang
  print(gtf_depth, "gtf_depth", exp_depth, "exp_depth")
  print(gtf_gene_initial, gtf_gene_last, "GTF")
  print(exp_gene_initial, exp_gene_last, "EXP")
  print("Type C", overhang, overlap)
  
  gene_type_assign = "CDS"
  
  overlap_details = c("Type C", overhang, overlap) #, gene_type_assign]
  
  updated_line = c()
  gene_name_index = unlist(strsplit(lines.split('\t')[8], ' ')).which(is.element('gene_name', unlist(strsplit(lines.split('\t')[8], ' '))))
  gene_name = strsplit(lines.split('\t')[8], ' ')[gene_name_index+1]
  gene_name = gsub('"', '', gene_name)
  gene_name = gsub(';', '', gene_name)
  for (items in unlist(strsplit(wt_lines, '\t'))[-1]) {
    items = gsub('"', '', items)
    updated_line = c(updated_line, items)
  }
  gtf_info = c(substr(wt_lines, nchar(wt_lines)), gene_name)
  updated_line = c(updated_line, gtf_info, overlap_details, relevant_seq, transcript_id, gene_biotype)
  updated_line = paste(updated_line, collapse = '\t')
  updated_line = gsub('[\\[\\]\']', '', updated_line)
  print(updated_line)
  print("--------------------------------------------------------------------------------------------")
  
  write.csv(updated_line, file = csv_file, append = TRUE, row.names = FALSE, col.names = FALSE)
}


# case D
#   ----           GTF
#----------        EXP
if ((gtf_gene_initial >= exp_gene_initial) & 
    (gtf_gene_last <= exp_gene_last)) {
  print(paste(gtf_depth, "gtf_depth", exp_depth, "exp_depth"))
  print(paste(gtf_gene_initial, gtf_gene_last, "GTF"))
  print(paste(exp_gene_initial, exp_gene_last, "EXP"))
  print("Type D")        
} else {
  # escape cases
  next
  #print("Found None")
}
