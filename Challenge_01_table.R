# library(dplyr)
# library(docstring)

### Indicate working directory
### this should include a folder called "Challenge samples" and
### a coordination.txt  file
wd <- "D:/Users/Zalva/OneDrive/Amphora/final"
samples_folder <- "Challenge samples/"

### Include new SNPs that are added at the end of the genotype table
### In case you leave this flag as *FALSE* the script will use the SNP
### founds in the first sample as base and new varians found in subsequent
### files will be discarded.
add_new <- FALSE


### if save_table set as TRUE and table_file is a filename
### then a CSV file with the genotype table will be saved
save_table <- TRUE
table_file <- "geno_table.csv"






############################################################
##################### DEFINE FUNCTIONS #####################
############################################################ 

import_geno <- function(file){
  #' Import genotype function
  #'
  #' This function reads a CSV or VCF file from one individual sample
  #' and import the content as a data frame.
  #'
  #' VCF files should include ten columns: "CHROM", "POS", "ID", "REF",
  #' "ALT", "QUAL","FILTER", "INFO", "FORMAT", "individual genotye"
  #' CSV files should include five columns: "CHROM;POS", "REF",
  #' "ALT",	"ALT.1",	"individual genotye"
  #'
  #'@param file input file to be processed. Could be either a CSV or a VCF

  ind <- gsub("\\..*","",basename(file))
  ### check file type and read table
  ### check if file is vcf
  if(grepl(".vcf", file, fixed=TRUE)){
    filetype <- "vcf"
    ### read vcf tables
    t <- read.table(file)
    ### check if vcf file have the expected format
    if(dim(t)[2] != 10){
      cat("WARNING \n
      The number of columns in the VCF file does not match
            what is expected (10).\n
          Please check the format\n
          file", file, "will NOT be included")
    } else {
      ### set col names
      colnames(t) <-  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                        "FILTER", "INFO", "FORMAT", ind)
      t$CHROM <- as.character(t$CHROM)
      t$POS <- as.character(t$POS)
    }
    
    return(t[,c("CHROM", "POS", "ID", "REF", "ALT", ind)])
    ### check if file is csv
  } else if(grepl(".csv", file, fixed=TRUE)){
    filetype <- "csv"
    ### read csv tables
    t <- read.csv(file, head=TRUE)
    ### check if csv file have the expected format
    if(dim(t)[2] != 5){
      cat("WARNING \n
      The number of columns in the CSV file does not match
            what is expected (5).\n
          Please check the format\n
          file", file, "will NOT be included")
    } else {
      ### make new data frame
      t <- data.frame(CHROM=as.character(sub(";.*", "", t[[1]])),
                      POS=as.character(sub(".*;", "", t[[1]])),
                      ID="",
                      REF=t$REF,
                      ALT=t$ALT,
                      ind=sub(",", "|", t[[5]], fixed=TRUE))
      ### set colnames
      colnames(t) <- c("CHROM", "POS", "ID", "REF", "ALT", ind)
      
      
    }
    return(t)
  } else {
    cat("file is different format that expected (VCF/CSV)")
  }
}

geno2num <- function(genotype){
  #' Genotype to numeric function
  #'
  #' Turn a vector of genotype format values (0|0)
  #' into a vector of numeric factors (1,2,3)
  #' 
  #' Values of 0|0 are homozygous for the reference allele and are converted to 0
  #' Values of 0|1 or 1|0 are heterozygous and are converted to 1
  #' Values of 1|1 are homozygous for the alternative allele and are converted to 2
  #'
  #'@param  genotype Vector of genotype values
  #'

  ### turn genotype to numeric factors (0,1,2)
  ### homozygous for reference
  geno_num <- gsub("0|0", "0",  genotype, fixed=TRUE)
  ### heterozygous
  geno_num <- gsub("0|1", "1",  geno_num, fixed=TRUE)
  geno_num <- gsub("1|0", "1",  geno_num, fixed=TRUE)
  ### homozygous for varian
  geno_num <- gsub("1|1", "2",  geno_num, fixed=TRUE)
  return(as.factor(as.numeric(geno_num)))
}


join_geno <- function(base, geno_table, add_new){
  #' Join two genotype tables
  #' 
  #' This function merge two genotype data frames into one new data frame
  #' 
  #' The first data frame is used as base. All the variants and genotypes
  #' present in this data frame will be conserved. 
  #' The second data frame is added to the base. Variants presented in 
  #' both data frames are conserved. If add_new flag is set as TRUE
  #' Variants present in this data frame but absent in the base are
  #' conserved and added at the end of the data. If add_new flag is set as
  #' FALSE then variants present in this data frame but absent in the base are
  #' discarded. All the missing variants of each individual are filled with NA
  #' The first 5 columns of both data frames should be named as follows:
  #' "CHROM", "POS", "ID", "REF", "ALT"
  #' REQUIRES dplyr
  #' 
  #' @param base genotype dataframe used as base
  #' @param geno_table genotype dataframe added to the base
  #' @param add_new If TRUE: new variants are added at the end of the data frame
  #'
  #' 
  
  if(all(colnames(geno_table)[1:5] == c("CHROM", "POS", "ID", "REF", "ALT")) &
     all(colnames(base)[1:5] == c("CHROM", "POS", "ID", "REF", "ALT"))){
    if(add_new == TRUE){
      t <- dplyr::full_join(base, geno_table[,c(1,2,4,5,6)],
                            by=c("CHROM", "POS", "REF", "ALT"))
    } else {
      t <- dplyr::left_join(base, geno_table[,c(1,2,4,5,6)],
                            by=c("CHROM", "POS", "REF", "ALT"))
    }
    return(t)
  } else { 
    cat("column names dont match to expected format\n
    columns should be as follows:\n
    CHROM, POS, ID, REF, ALT genotypes... ")
  }
}


############################################################
########################## START ###########################
############################################################
### set working directory
setwd(wd)

### list samples in samples folder
files <- list.files(samples_folder)

############################################################
####################### read tables ########################
############################################################
### iterate through each file in the folder
for(i in 1:length(files)){
  file <- paste("Challenge samples/", files[i], sep="")
  if(i == 1){
    base_geno <- import_geno(file)
    base_geno[[6]] <- geno2num(base_geno[[6]])
  } else {
    geno_table <- import_geno(file)
    geno_table[[6]] <- geno2num(geno_table[[6]])
    base_geno <- join_geno(base_geno, geno_table, add_new)
  }
}

############################################################
######################## save file #########################
############################################################
if(save_table == TRUE){
  write.csv(base_geno, table_file, row.names=FALSE)
}

