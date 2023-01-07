#Annotate bed file reads with 
#NCBI RefSeq table of from UCSC:
#Go to UCSC genome browser > Tools > Table Browser > select correct genome and "NCBI Refseq" in track dropdown menu.
#0-based start, 1-based end.

#NOTE:
##For mRNA, only use transcripts with NM_ prefix in name column
##5_prime_UTR is anything upstream of CDS
##3_prime_UTR is anything downstream of CDS
##intron is anything outside of CDS, 5_prime_UTR and 3_prime_UTR but within txStart and tsEnd.
##In cases where there's more than 1 NM_ prefix transcript, use all but when annotating, give priority
##in this order: CDS>5'UTR,3'UTR,intron.
##e.g. if a site is found in CDS of transcript1 but intron in transcript2 of same gene, denote as CDS.

##For ncRNA, only use transcripts with NR_ prefix in name column
##v3 same as v2 but annotates unannotated sites that are X nucleotides upstream of 5'UTR.


format_gff3 <- function(directory, ref_genome) {
  # Format gff3 file to define 5'UTR, 3'UTR, introns and promoters
  
  # First pass: for each gene, find most confident gene prefix.
  gffs <- list()
  a <- file(paste0(directory, ref_genome, '_NCBI_RefSeq.txt'), 'r')
  line <- readLines(a, n = 1)
  while (length(line) > 0) {
    ar <- read.table(text = line)
    gene <- ar[, 13]
    prefix <- substr(ar[, 2], 1, 2)
    
    if (gene %in% names(gffs)) {
      gffs[[gene]] <- c(gffs[[gene]], prefix)
    } else {
      gffs[[gene]] <- prefix
    }
    line <- readLines(a, n = 1)
  }
  close(a)
  
  for (gene in names(gffs)) {
    # For mRNA
    if ('NM' %in% gffs[[gene]]) {
      gffs[[gene]] <- 'NM'
    }
    # For ncRNA
    else if ('NR' %in% gffs[[gene]]) {
      gffs[[gene]] <- 'NR'
    }
    # Not sure what YP is but it only tags mitochondrial-encoded mRNA.
    else if ('YP' %in% gffs[[gene]]) {
      gffs[[gene]] <- 'YP'
    }
    # For predicted mRNA
    else if ('XM' %in% gffs[[gene]]) {
      gffs[[gene]] <- 'XM'
    }
    # For predicted ncRNA
    else if ('XR' %in% gffs[[gene]]) {
      gffs[[gene]] <- 'XR'
    }
  }
  
  # Second pass, collect all lines that have prefix corresponding to highest confidence gene prefix for each gene, and collect under dictionary "temp".
  temps <- list()
  a <- file(paste0(directory, ref_genome, '_NCBI_RefSeq.txt'), 'r')
  line <- readLines(a, n = 1)
  while (length(line) > 0) {
    ar <- read.table(text = line)
    gene <- ar[, 13]
    prefix <- substr(ar[, 2], 1, 2)
    if (prefix == gffs[[gene]]) {
      if (gene %in% names(temps)) {
        temps[[gene]] <- c(temps[[gene]], line)
      } else {
        temps[[gene]] <- line
      }
    }
    line <- readLines(a, n = 1)
  }
  close(a)
}


# Define all regions and collect under "annotations" dictionary
annotations <- list()
# Generation sub-dictionaries.
for (gene in names(temps)) {
  for (line in temps[[gene]]) {
    ar <- read.table(text = line)
    chro <- ar[, 3]
    strand <- ar[, 4]
    tryCatch({
      annotations[[paste(chro, strand)]] [[gene]] <- list()
    }, error = function(e) {
      annotations[[paste(chro, strand)]] <- list()
      annotations[[paste(chro, strand)]] [[gene]] <- list()
    })
    for (region in c('cds', '3utr', '5utr', 'intron', 'ncRNA', 'intron_ncRNA')) {
      annotations[[paste(chro, strand)]] [[gene]] [[region]] <- list()
    }
  }
}

# Actual defining of regions.
for (gene in names(temps)) {
  for (line in temps[[gene]]) {
    ar <- read.table(text = line)
    chro <- ar[, 3]
    strand <- ar[, 4]
    if (strand == '+') {
      x <- '5utr'
      y <- '3utr'
    } else {
      y <- '5utr'
      x <- '3utr'
    }
    exonstarts <- as.integer(unlist(strsplit(ar[, 10], ','))[-length(strsplit(ar[, 10], ','))])
    exonends <- as.integer(unlist(strsplit(ar[, 11], ','))[-length(strsplit(ar[, 11], ','))])
    # For mRNAs
    if(gffs[gene] %in% c('NM', 'XM')) {
      cdsstart <- as.integer(ar[6])
      cdsend <- as.integer(ar[7])
      # Single exon.
      if(length(exonstarts) == 1) {
        annotations[[chro, strand]][gene]['cds'] <- list(list(cdsstart, cdsend))
        annotations[[chro, strand]][gene][x] <- list(list(txstart, cdsstart))
        annotations[[chro, strand]][gene][y] <- list(list(cdsend, txend))
      } else {
        # Find out if any 5utr or 3utr is interrupted by introns; also denote introns.
        for(i in seq_along(exonstarts)) {
          if(exonstarts[i] <= cdsstart && cdsstart < exonends[i]) {
            leftutr <- i
          }
          if(exonstarts[i] < cdsend && cdsend <= exonends[i]) {
            rightutr <- i
          }
          if(i > 1) {
            annotations[[chro, strand]][gene]['intron'] <- list(list(exonends[i - 1], exonstarts[i]))
          }
        }
      }
    }
    # Denote cds.
    for(i in seq(leftutr, rightutr + 1)) {
      if(i == leftutr) {
        annotations[[chro, strand]][gene]['cds'] <- list(list(cdsstart, exonends[i]))
      } else if(i == rightutr) {
        annotations[[chro, strand]][gene]['cds'] <- list(list(exonstarts[i], cdsend))
      } else {
        annotations[[chro, strand]][gene]['cds'] <- list(list(exonstarts[i], exonends[i]))
      }
    }
    # Denote left-most UTR.
    for(i in seq_len(leftutr + 1)) {
      if(i < leftutr) {
        annotations[[chro, strand]][gene][x] <- list(list(exonstarts[i], exonends[i]))
      } else if(i == leftutr) {
        annotations[[chro, strand]][gene][x] <- list(list(exonstarts[i], cdsstart))
      }
    }
    # Denote right-most UTR.
    for(i in seq(rightutr, length(exonstarts))) {
      if(i == rightutr) {
        annotations[[chro, strand]][gene][y] <- list(list(cdsend, exonends[i]))
      }
      if(i > rightutr) {
        annotations[[chro, strand]][gene][y] <- list(list(exonstarts[i], exonends[i]))
      }
    }
    
    # For ncRNAs.
    if(gffs[gene] %in% c('NR', 'XR')) {
      # Single exon.
      if(length(exonstarts) == 1) {
        annotations[[chro, strand]][gene]['ncRNA'] <- list(list(txstart, txend))
      } else {
        for(i in seq_along(exonstarts)) {
          annotations[[chro, strand]][gene]['ncRNA'] <- list(list(exonstarts[i], exonends[i]))
          if(i > 1) {
            annotations[[chro, strand]][gene]['intron_ncRNA'] <- list(list(exonends[i - 1], exonstarts[i]))
          }
        }
      }
    } else if(gffs[gene] == 'YP') {
      # For YP_ prefixes.
      annotations[[chro, strand]][gene]['cds'] <- list(list(txstart, txend))
    }
    
    rm(gffs)
    return(annotations)

    
## Actual annotations     
ref_genome <- 'hg38'
directory <- paste0('store/annotations/', ref_genome, '/')
priorities <- c('cds', '3utr', '5utr', 'intron', 'ncRNA', 'intron_ncRNA')
annotations <- format_gff3(directory, ref_genome)
filename <- 'Empty-E2+IAV_Pcif1-183B2+IAV'
a <- file('store/quanti/' + filename + '_consolidate_output.txt', 'r')
z <- 3 #For strand index
d <- file('store/quanti/' + filename + '_consolidate_output_annotated.txt', 'w')


# Actual annotating.
sites <- list()
sitecount <- 0
while(length(line <- readLines(a, n = 1)) > 0) {
  sitecount <- sitecount + 1
  sites[[sitecount]] <- read.table(text = line)
}
close(a)

for (site in seq_along(sites)) {
  ar <- sites[[site]]
  chro <- ar[1]
  strand <- ar[4]
  if (!(paste(chro, strand) %in% names(annotations))) {
    sites[[site]] <- c(sites[[site]], 'NA', 'unannotated')
  } else {
    end <- as.integer(ar[3])
    unannotated <- TRUE
    for (priority in priorities) {
      for (gene in names(annotations[[paste(chro, strand)]])) {
        tryCatch({
          for (region in annotations[[paste(chro, strand)]][[gene]][[priority]]) {
            if (region[1] < end <= region[2]) {
              sites[[site]] <- c(sites[[site]], gene, priority)
              unannotated <- FALSE
              break
            }
          }
        }, error = function(e) {})
        if (!unannotated) {
          break
        }
      }
      if (!unannotated) {
        break
      }
    }
  }
}


if (unannotated) {
  for (gene in names(annotations[[paste(chro, strand)]])) {
    tryCatch({
      for (region in annotations[[paste(chro, strand)]][[gene]][['5utr']]) {
        if ((strand == '+' && region[1] - up < end <= region[1]) || (strand == '-' && region[2] < end <= region[2] + up)) {
          sites[[site]] <- c(sites[[site]], gene, '5utr')
          unannotated <- FALSE
          break
        }
      }
    }, error = function(e) {})
    if (!unannotated) {
      break
    }
  }
  
  if (unannotated) {
    sites[[site]] <- c(sites[[site]], 'NA', 'unannotated')
  }
}
rm(annotations)

# To output frequency of each annotation.
represent <- list(unannotated = 0)
for (priority in priorities) {
  represent[[priority]] <- 0
}
for (annotat in sites) {
  if (annotat[[length(annotat)]] == 'unannotated') {
    represent[['unannotated']] <- represent[['unannotated']] + 1
  } else {
    represent[[annotat[[length(annotat)]]]] <- represent[[annotat[[length(annotat)]]]] + 1
  }
  writeLines(paste(annotat, collapse = '\t'), con = d)
}

for (key in names(represent)) {
  cat(paste(key, represent[[key]], sep = '\t'), '\n')
}

# Remove input file to save space.
prog <- 'rm '
command <- paste0(prog, 'store/quanti/', filename, '_consolidate_output.txt')
system(command, intern = TRUE)


