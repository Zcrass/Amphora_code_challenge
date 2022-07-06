# library(dplyr)
# library(ggplot2)
# library(plotly)
# library(caret)

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

### indicate if you want to extract the SNPs whit the highest 
### contribution to the PCAs and how many PCAs you want to consider for this.
### recommended to run first the PCA analysis before modify this
identify_SNP <- TRUE
number_PC <- 4

### if you want to extract the SNPs whit the highest 
### contribution to the PCAs and save it as a different genotype
###table set this as TRUE and indicate a filename
save_selected <- TRUE
selected_file <- "geno_selected.csv"



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

### load genotype table
data_pca <- load_geno(table_file, "ind_as_col", omitted=omit_ind)

# ### compute PCA and extract PCA values
pca <- prcomp(data_pca)
pca_vals <- pca$x

### add UUID and Superpopulation code
pca_vals <- cbind(rownames(pca_vals), pca_vals)
colnames(pca_vals)[1] <- "UUID"
pca_vals <- dplyr::left_join(as.data.frame(pca_vals), coord, by="UUID")
pca_vals$Superpopulation.code[is.na(pca_vals$Superpopulation.code)] <- "NA"
### save for future plots
write.csv(pca_vals, "pca_vals.csv")

###################### Plot figures ########################
### plot pca using ggplot
# ggplot(pca_vals, aes(x=PC1, y=PC2, color=Superpopulation.code)) +
#   geom_point() +
#   labs(x="PC1", y="PC2") +
#   theme_classic() +
#   theme(axis.text.x=element_blank(),
#         axis.ticks=element_blank(),
#         axis.text.y=element_blank()) +
#   scale_color_manual(values=colors)
#   

### plot pca using plotly (3D)
plotly::plot_ly(x=pca_vals$PC1, y = pca_vals$PC2, z = pca_vals$PC3, 
        type="scatter3d", mode = "markers",
        color=pca_vals$Superpopulation.code,
        colors=colors)

### plot PCA varianza explained
prop_varianza <- pca$sdev^2 / sum(pca$sdev^2)
PCA_var <- ggplot2::ggplot(data = data.frame(prop_varianza[1:10], pc = 1:10),
  ggplot2::aes(x = pc, y = prop_varianza[1:10])) +
  ggplot2::geom_col(width = 0.3) +
  ggplot2::scale_y_continuous(limits = c(0,0.1)) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Componente principal",
                y = "Prop. de varianza explicada")

### plot PCA cumulative variance
prop_varianza_acum <- cumsum(prop_varianza)
PCA_acum_var <- ggplot2::ggplot(data = data.frame(prop_varianza_acum[1:10], pc = 1:10),
  ggplot2::aes(x = pc, y = prop_varianza_acum[1:10], group = 1)) +
  ggplot2::geom_point() +
  ggplot2::geom_line() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Componente principal",
                y = "Prop. varianza explicada acumulada")

gridExtra::grid.arrange(PCA_var, PCA_acum_var, ncol=2)



############################################################
################# higest contribution SNPs #################
############################################################

if(identify_SNP == TRUE){
  ### if TRUE select the SNPs with the highest contribution to the PCAs
  for(i in 1:number_PC){
    if(i == 1){
      selected <- c()
    }
    selected <- c(selected, 
                  names(sort(abs(pca$rotation[,i]),decreasing = TRUE))[1:100])
  }
  
  ### sort variants selected
  selected <- sort(unique(selected))
  
  ### compute new PCA nly with selected SNPs and extract new values
  pca_selected <- prcomp(data_pca[,selected])
  selected_vals <- pca_selected$x
  
  ### add UUID and Superpopulation code
  selected_vals <- cbind(rownames(selected_vals), selected_vals)
  colnames(selected_vals)[1] <- "UUID"
  selected_vals <- dplyr::full_join(as.data.frame(selected_vals), coord, by="UUID")
  selected_vals$Superpopulation.code[is.na(selected_vals$Superpopulation.code)] <- "NA"
  
  ###################### Plot figures ########################
  ### plot pca using ggplot
  # ggplot(selected_vals, aes(x=PC1, y=PC2, color=Superpopulation.code)) +
  #   geom_point() +
  #   labs(x="PC1", y="PC2") +
  #   theme_classic() +
  #   theme(axis.text.x=element_blank(),
  #         axis.ticks=element_blank(),
  #         axis.text.y=element_blank()) +
  #   scale_color_manual(values=colors)
  
  ### plot pca using plotly (3D)
  plotly::plot_ly(x = selected_vals$PC1, y = selected_vals$PC2, z = selected_vals$PC3, 
          type = "scatter3d", mode = "markers",
          color = selected_vals$Superpopulation.code,
          colors = colors)
  if(save_selected == TRUE){
    write.csv(data_pca[, colnames(data_pca) %in% selected], 
              file = selected_file, 
              row.names = FALSE )
  }
}

############################################################
################# hierarchical clustering ##################
############################################################
### load data for hierarchical clustering
data_clust <- load_geno(table_file, "ind_as_col", omitted=omit_ind)

### compute Nei's distance
data_dis <- poppr::nei.dist(data_clust)

### perform hierarchical clustering and split in 5 clusters
h_clust <- hclust(data_dis, method = "ward.D2")
h_groups <- data.frame(hclust_val = cutree(h_clust, k=5))
h_groups$UUID <- rownames(h_groups)

###################### Plot figures ########################
plotly::plot_ly(x=pca_vals$PC1, y = pca_vals$PC2, z = pca_vals$PC3, 
                type="scatter3d", mode = "markers",
                color=as.factor(h_groups$hclust_val),
                colors=colors[c(1,2,3,4,6)])
### plot pca using ggplot
# ggplot2::ggplot(pca_vals, 
#   ggplot2::aes(x=PC1, y=PC2, color=as.factor(h_groups$hclust_val))) +
#   ggplot2::geom_point() +
#   ggplot2::labs(x="PC1", y="PC2") +
#   ggplot2::theme_classic() +
#   ggplot2::theme(axis.text.x=ggplot2::element_blank(),
#         axis.ticks=ggplot2::element_blank(),
#         axis.text.y=ggplot2::element_blank()) +
#   ggplot2::scale_color_manual(values=colors[c(1,2,3,4,6)])


############################################################
#################### kmeans clustering #####################
############################################################
### load data for kmeans clustering
data_kmeans <- load_geno(table_file, "ind_as_col", omitted=omit_ind)
### perform kmeans clusterin for 5 clusters
kmeans <- kmeans(data_kmeans, centers=5, iter.max = 1000, nstart = 10)
kmeans_group <- data.frame(kmeans_val = kmeans$cluster)
kmeans_group$UUID <- rownames(kmeans_group)

###################### Plot figures ########################
plotly::plot_ly(x=pca_vals$PC1, y = pca_vals$PC2, z = pca_vals$PC3, 
                type="scatter3d", mode = "markers",
                color=as.factor(kmeans_group$kmeans_val),
                colors=colors[c(1,2,3,4,6)])

# ggplot2::ggplot(pca_vals, 
#                 ggplot2::aes(x=PC1, y=PC2, color=as.factor(kmeans_group$kmeans_val))) +
#   ggplot2::geom_point() +
#   ggplot2::labs(x="PC1", y="PC2") +
#   ggplot2::theme_classic() +
#   ggplot2::theme(axis.text.x=ggplot2::element_blank(),
#                  axis.ticks=ggplot2::element_blank(),
#                  axis.text.y=ggplot2::element_blank()) +
#   ggplot2::scale_color_manual(values=colors[c(1,2,3,4,6)])
############################################################
################### evaluate clustering ####################
############################################################
### make data frame for results
results <- data.frame(UUID = coord$UUID,
                      True_val = coord$Superpopulation.code)
results <- dplyr::left_join(results, h_groups, by = "UUID")
results <- dplyr::left_join(results, kmeans_group, by = "UUID")
rownames(results) <- results$UUID

### set groups names based on original data
### for AFR
results$hclust_val[results$hclust_val == results["cfd8c6f1-5e09-5589-92de-67560eb0acb8",3]] <- "AFR"
results$kmeans_val[results$kmeans_val == results["cfd8c6f1-5e09-5589-92de-67560eb0acb8",4]] <- "AFR"

### for AMR
results$hclust_val[results$hclust_val == results["416db4f1-9443-594e-a4b7-bccc5ba59b02",3]] <- "AMR"
results$kmeans_val[results$kmeans_val == results["416db4f1-9443-594e-a4b7-bccc5ba59b02",4]] <- "AMR"

### for EAS
results$hclust_val[results$hclust_val == results["8eb7539a-9b7d-580c-bb32-d940441fcbdc",3]] <- "EAS"
results$kmeans_val[results$kmeans_val == results["8eb7539a-9b7d-580c-bb32-d940441fcbdc",4]] <- "EAS"

### for EUR
results$hclust_val[results$hclust_val == results["75240c09-57b6-5579-bb62-0cc50247dae2",3]] <- "EUR"
results$kmeans_val[results$kmeans_val == results["75240c09-57b6-5579-bb62-0cc50247dae2",4]] <- "EUR"

### for SAS
results$hclust_val[results$hclust_val == results["e07cc175-dea0-5474-a2c9-e1b2cb34baf6",3]] <- "SAS"
results$kmeans_val[results$kmeans_val == results["e07cc175-dea0-5474-a2c9-e1b2cb34baf6",4]] <- "SAS"

### set as factor
results$True_val <- as.factor(results$True_val)
results$hclust_val <- as.factor(results$hclust_val)
results$kmeans_val <- as.factor(results$kmeans_val)

### compute confusion matrix and statistics
cm_hclust <- caret::confusionMatrix(results$True_val, results$hclust_val)
cm_kmeans <- caret::confusionMatrix(results$True_val, results$kmeans_val)
t(cm_hclust$byClass)
t(cm_kmeans$byClass)



