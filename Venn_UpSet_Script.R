#install these packages if needed
#install.packages("readxl")
#install.packages("ggplot2")
#install.packages("data.table")
#install.packages("dplyr")
#install.packages("UpSetR")
library("readxl")
library("data.table")
library("ggplot2")
library("dplyr")
library(eulerr)
set.seed(1)
library(UpSetR)


#example = the Wt alone v. mutants alone def expression data. Has both log2FC and p-values.

#Import excel file with log2foldchange and p-values.
#IMPORTANT: Make sure that all p-values that are NA are just changed to a numeric value higher than p-value threshold (eg. 100)
#IMPORTANT: Make sure that all there are no NAs in the log2foldchange columns in the excel. Change these to 0 or whatever makes sense for your analysis.
#This is required to avoid some class issues.
rnaseq <- read_excel('/Users/dongkyun/Desktop/RNAseq data/MASTER_2017_2020_alone.xlsx')

#dataframe of the log2FC data only from your original excel.
rnaseq_log = rnaseq[, c(3, 5, 7, 9, 11)]

#Make sure all columns are of the double class. Convert any that are not with this line.
rnaseq_log[, 1] <- apply(rnaseq_log[, 1], MAR = 2, as.double)

#dataframe of the p-values only from your original excel; convert p-values from characters to double.
rnaseq_pval = apply(rnaseq[, c(4, 6, 8, 10, 12)], MAR = 2, as.double)

#need to round the p_values because R is dumb 
rnaseq_pval = as.data.frame(round(rnaseq_pval, 3))


#Check to see if it looks good and matches the excel sheets
rnaseq_pval
rnaseq_log

#Extracting each p-value column as a vector. May need to alter column name in the lines below to match what you have.
pval_a = rnaseq_pval$A8_padj
pval_gaf = rnaseq_pval$gaf_padj
pval_b = rnaseq_pval$B8_padj
pval_c = rnaseq_pval$C8_padj
pval_d = rnaseq_pval$D8_padj

#Next few lines: Replacing any p-value in rnaseq_pval that is below the p-value threshold to the character 'notsig'
#This will keep only the values that are below the p-value cut off.
pval_a[pval_a > 0.05] <- 'notsig'
pval_gaf[pval_gaf > 0.05] <- 'notsig'
pval_b[pval_b > 0.05] <- 'notsig'
pval_c[pval_c > 0.05] <- 'notsig'
pval_d[pval_d > 0.05] <- 'notsig'

#Change any NA p-value to notsig. Should have already been done manually in the excel file, but just in case.
pval_a = replace(pval_a, pval_a== 'NA', 'notsig')
pval_gaf = replace(pval_gaf, pval_gaf=='NA', 'notsig')
pval_b = replace(pval_b, pval_b=='NA', 'notsig')
pval_c = replace(pval_c, pval_c=='NA', 'notsig')
pval_d = replace(pval_d, pval_d=='NA', 'notsig')

#make a new data frame with all insignificant p-values replaced with 'notsig'
filtered_pval <- data.frame( pval_a, pval_gaf, pval_b, pval_c, pval_d)


#now need to use filtered_pval to sort through the log2foldchange. 

#Gets rid of anything with a p-value that is too large.
#IF YOU CHANGED THE NAMES OF rnaseq_log AND/OR filtered_pval, YOU ALSO NEED TO CHANGE THAT IN THIS FUNCTION.
check_sig <- function(x_val, p_val){
  x <- rnaseq_log[, x_val]
  p <- filtered_pval[, p_val]
  for (i in 1:nrow(x)) {
    if (p[i] == 'notsig') {
      rnaseq_log[i, x_val] <- 0
    }}
  return(rnaseq_log)}

#for loops take forever to run on R; be patient.
#Basically filtering each column of rnaseq_log with the respective column of p-values from filtered_pval
rnaseq_log <- check_sig(1, 1)
rnaseq_log <- check_sig(2, 2)
rnaseq_log <- check_sig(3, 3)
rnaseq_log <- check_sig(4, 4)
rnaseq_sig <- check_sig(5, 5)

#This dataframe has now been filtered to only retain log2FC data that is significant. 
rnaseq_sig

#Export/save as csv
write.csv(rnaseq_sig, "binary_sig.csv")


#Now only want upregulation. Input is the table with non-signficant data as "0". Outputs table where anything 
#significantly upregulated has a 1
#For all three functions, you can change the number that comes after the logic sign to change your fold-change threshold.
upreg_only <- function(col){
  col[col < 1] <- as.double(0)
  col[col >= 1] <- as.double(1)
  return(col)
}

#downregulation. Input is the table with the non-signficant data as "0". Outputs table where anything
#significantly downregulated has a 1
downreg_only <- function(col){
  col[col > -1] <- as.double(0)
  col[col <= -1] <- as.double(1)
  return(col)
}

#differently expressed. All up/down-regulated genes has value of 1
diff_expr <- function(col) {
  col[(col < 1)&(col > -1)] <- as.double(0)
  col[col >= 1] <- as.double(1)
  col[col <= -1] <- as.double(1)
  return(col)
}

#Run the functions on your rnaseq_sig to get a dataframe with genes that were upregulated, downregulated, etc
all_up = as.data.frame(apply(rnaseq_sig, MAR = 2, upreg_only))
all_down = as.data.frame(apply(rnaseq_sig, MAR = 2, downreg_only))
all_diff = as.data.frame(apply(rnaseq_sig, MAR = 2, diff_expr))

#use these outputed dataframes to make the venn diagram. You can also use to manually check for any mistakes.
all_diff
all_up
all_down

#change colunm names since it's no longer log2FC data.
colnames(all_diff) <- c("CvnA mutant", "GAF mutant","CvnB mutant","CvnC mutant","CvnD mutant")
colnames(all_up) <- c("CvnA mutant", "GAF mutant","CvnB mutant","CvnC mutant","CvnD mutant")
colnames(all_down) <- c("CvnA mutant", "GAF mutant","CvnB mutant","CvnC mutant","CvnD mutant")

#These have the gene locus/info added on.
all_diff_gene = data.frame(rnaseq[, 1], rnaseq[, 2], all_diff[, 1], all_diff[, 2], all_diff[, 3], all_diff[, 4], all_diff[, 5])
all_up_gene = data.frame(rnaseq[, 1], rnaseq[, 2], all_up[, 1], all_up[, 2], all_up[, 3], all_up[, 4], all_up[, 5])
all_down_gene = data.frame(rnaseq[, 1], rnaseq[, 2], all_down[, 1], all_down[, 2], all_down[, 3], all_down[, 4], all_down[, 5])

#Add column names.
colnames(all_diff_int_gene) <- c("Locus", "description", "CvnA Mutant", "GAF Mutant", "CvnB Mutant", "CvnC Mutant", "CvnD Mutant")
colnames(all_up_int_gene) <- c("Locus", "description", "CvnA Mutant", "GAF Mutant", "CvnB Mutant", "CvnC Mutant", "CvnD Mutant")
colnames(all_down_int_gene) <- c("Locus", "description","CvnA Mutant", "GAF Mutant", "CvnB Mutant", "CvnC Mutant", "CvnD Mutant")

#Saved these here.
write.csv(all_diff_gene, "binary_all_diff.csv")
write.csv(all_up_gene, "binary_all_up.csv")
write.csv(all_down_gene, "binary_all_down.csv")


## to make a venn diagram ##

#Run these lines with one dataframe at a time to make Venn/UpSet Diagrams for one set of data at a time.

#Make the Venn diagram object using one of your dataframes.
TheVenn <- euler(all_up, shape="ellipse")

#If you want to do any error-analysis:
TheVenn
error_plot(TheVenn)
coef(TheVenn)

#To check which genes are exclusive to one mutant; can use to double-check your Venns.
#Will need to change A, B, etc to your column names.
nrow(filter(VennDat, A == 1, B == 0, C == 0, D == 0, GAF == 0))

#Plot Venn diagram.
plot(TheVenn, quantities = TRUE)

#Make/plot equi-area Venn diagram (the star-shaped one)
TheVenn2 <- venn(VennDat)
plot(TheVenn2)

##To make UpSet Diagrams
upset(all_up)

#This more complicated UpSet-making line just configurates the diagram differently from the default settings.
upset(all_up, sets = c("GAF mutant", "CvnD mutant", "CvnC mutant", "CvnB mutant", "CvnA mutant"), order.by = 'freq',  sets.x.label = 'Total Genes', mainbar.y.label = "Number of Genes", keep.order = TRUE, empty.intersections = 'on')
