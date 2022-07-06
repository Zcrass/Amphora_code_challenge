# library(scorecard)
# library(factoextra)
# library(class)
# library(dplyr)



### Indicate working directory
### this should include a folder called "Challenge samples" and
### a coordination.txt  file
wd <- "D:/Users/Zalva/OneDrive/Amphora/final"
samples_folder <- "Challenge samples/"
coord_file <- "coordination.txt"

### Indicate the name of the genotype table file
table_file <- "geno_table.csv"

### if you want to omit an individual indicate their UUID
### else, set as FALSE.
### recomended for samples with high number of NAs
omit_ind <- "0abbd799-8737-54f5-b587-7bea7100a61f"




############################################################
##################### DEFINE FUNCTIONS #####################
############################################################ 
load_geno <- function(file, format, omitted = FALSE){
  #' Load genotype tables
  #' 
  #' read a genotype table from a CSV file and turn it to a matrix
  #' 
  #' This function reads a CSV file to load a genotype table and 
  #' return a matrix in different formats:\n
  #' "ind_as_col": each column is an individual genotype
  #' and the rows represent the different loci\n
  #' "snp_as_col": each column is a locus 
  #' and the rows represent each individual genotype\n
  #' "gwas": genotype table required by the GWAS function 
  #' from the rrBLUP package
  #' 
  #' @param file genotype table stored as CSV file
  #' @param format format to return the output data frame (see details)
  #' @omitted character vector of omitted individuals.
  #' 
  #' 
  ### load geno table file
  data <- read.csv(file)
  ### change rownames and colnames
  rownames(data) <- paste(data$CHROM, ";", data$POS, ";",
                          data$REF, ";", data$ALT, sep="")
  colnames(data) <- gsub(".", "-", colnames(data), fixed=T)
  colnames(data) <- gsub("X", "", colnames(data), fixed=T)
  
  ### remove omitted individuals
  if(length(omitted) != FALSE){
    data <- data[,!names(data) %in% omit_ind]
  }
  if(format == "ind_as_col"){
    ### remove SNPs with NAs
    return(t(data[complete.cases(data),6:dim(data)[2]]))
  } else if(format == "snp_as_col"){
    return(data[complete.cases(data),6:dim(data)[2]])
  } else if(format == "gwas"){
    gwas_data <- data.frame(marker = data$ID,
                            chrom =as.numeric(data$CHROM),
                            pos = as.numeric(data$POS))
    gwas_data <- cbind(gwas_data, data[6:ncol(data)])
    return(gwas_data)
  }
}

############################################################
########################## START ###########################
############################################################
### set working directory
setwd(wd)

### set colors for plots
colors <- c("#ff5733", "#1aefdc", "#f0a907", "#7407f0", "#999999", "#046d0e")

### read coordination file
coord <- read.table(coord_file, header = TRUE, sep="\t")

### load genotype table and change rownames
data <- as.data.frame(load_geno(table_file, "ind_as_col", omitted=omit_ind))
colnames(data) <- paste("X", 
                        gsub(";", "_", colnames(data)), sep="")

### add superpopulation code
data$UUID <- rownames(data)
data <- dplyr::left_join(data, coord, by="UUID")
# data$UUID <- NULL

### extract labeled data and unlabeled
data_lab <- data[!is.na(data$Superpopulation.code),]
data_unlab <- data[is.na(data$Superpopulation.code),]
data_unlab$Superpopulation.code <- NULL

### split labeled data in train and test
data_sp <- scorecard::split_df(data_lab, ratio = 0.6, seed = 1234)

train <- data_sp$train; train[,c("UUID", "Superpopulation.code")] <- NULL
test <- data_sp$test; test[,c("UUID", "Superpopulation.code")] <- NULL
country_train <- as.factor(data_sp$train$Superpopulation.code)
country_test <- as.factor(data_sp$test$Superpopulation.code)
pca_vals <- read.csv("pca_vals.csv", row.names = 1)
pca_test <- pca_vals[pca_vals$UUID %in% data_sp$test$UUID, ]

############################################################
##################### Run classifiers ######################
############################################################
### K nearest neighbor classifier
knn <- class::knn(train = train,
                  test = test,
                  k = 1,
                  cl = country_train)
### predicted results for test data
p_knn <- data.frame(UUID = data_sp$test$UUID,
                    p_knn = knn)
### plot clustering
plotly::plot_ly(x=pca_test$PC1, y = pca_test$PC2, z = pca_test$PC3, 
                type="scatter3d", mode = "markers",
                color=p_knn$p_knn,
                colors=colors)

ggplot2::ggplot(pca_test, 
  ggplot2::aes(x=PC1, y=PC2, color=as.factor(p_knn$p_knn))) +
  ggplot2::geom_point() +
  ggplot2::labs(x="PC1", y="PC2") +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                 axis.ticks=ggplot2::element_blank(),
                 axis.text.y=ggplot2::element_blank()) +
  ggplot2::scale_color_manual(values=colors[c(1,2,3,4,6)])

############################################################
### random forest classifier
rf <- randomForest::randomForest(country_train~., 
                                 data=train)
### predicted results for test data
p_rf <- data.frame(UUID = data_sp$test$UUID,
                   p_rf = predict(rf, test))
plotly::plot_ly(x=pca_test$PC1, y = pca_test$PC2, z = pca_test$PC3, 
                type="scatter3d", mode = "markers",
                color=p_rf$p_rf,
                colors=colors)

ggplot2::ggplot(pca_test, 
  ggplot2::aes(x=PC1, y=PC2, color=as.factor(p_rf$p_rf))) +
  ggplot2::geom_point() +
  ggplot2::labs(x="PC1", y="PC2") +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                 axis.ticks=ggplot2::element_blank(),
                 axis.text.y=ggplot2::element_blank()) +
  ggplot2::scale_color_manual(values=colors[c(1,2,3,4,6)])
############################################################
### naive bayes classifier
nb <- naivebayes::naive_bayes(country_train~., data=train)
### predicted results for test data
p_nb <- data.frame(UUID = data_sp$test$UUID,
                   p_nb = predict(nb, test))
plotly::plot_ly(x=pca_test$PC1, y = pca_test$PC2, z = pca_test$PC3, 
                type="scatter3d", mode = "markers",
                color=p_nb$p_nb,
                colors=colors)

ggplot2::ggplot(pca_test, 
  ggplot2::aes(x=PC1, y=PC2, color=as.factor(p_nb$p_nb))) +
  ggplot2::geom_point() +
  ggplot2::labs(x="PC1", y="PC2") +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                 axis.ticks=ggplot2::element_blank(),
                 axis.text.y=ggplot2::element_blank()) +
  ggplot2::scale_color_manual(values=colors[c(1,2,3,4,6)])
############################################################
################# evaluate classification ##################
#############################################################
results <- data.frame(UUID = data_sp$test$UUID,
                      True_val = as.factor(data_sp$test$Superpopulation.code))
results <- dplyr::left_join(results, p_knn, by = "UUID")
results <- dplyr::left_join(results, p_rf, by = "UUID")
results <- dplyr::left_join(results, p_nb, by = "UUID")


str(results)

cm_knn <- caret::confusionMatrix(results$True_val, results$p_knn)
cm_rf <- caret::confusionMatrix(results$True_val, results$p_rf)
cm_nb <- caret::confusionMatrix(results$True_val, results$p_nb)

t(cm_knn$byClass)
t(cm_rf$byClass)
t(cm_nb$byClass)



