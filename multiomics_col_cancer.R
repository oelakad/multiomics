#Loading packages----
pacman::p_load(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
pacman::p_load(tximport) # package for getting Kallisto results into R
pacman::p_load(biomaRt) # provides access to a wealth of annotation info
pacman::p_load(RColorBrewer) 
pacman::p_load(reshape2) 
pacman::p_load(genefilter)
pacman::p_load(matrixStats)
pacman::p_load(hrbrthemes)
pacman::p_load(reshape2)
pacman::p_load(DT)
pacman::p_load(gt)
pacman::p_load(plotly)
pacman::p_load(limma) #powerful package for differential gene expression using linear modeling
pacman::p_load(edgeR) #another great package for differential gene expression analysis
pacman::p_load(hrbrthemes)
pacman::p_load(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps
pacman::p_load(RColorBrewer) #need colors to make heatmaps
pacman::p_load(limma) #we only use limma in this script for the 'avearrays' function
pacman::p_load(heatmaply) #for making interactive heatmaps using plotly
pacman::p_load(d3heatmap, hrbrthemes) #for making interactive heatmaps using D3
pacman::p_load(knitr, tinytex, rmarkdown) 
devtools::install_github("compgenomr/compGenomRData")#this is where the data used here is
library("compGenomRData")
library("devtools")
#
# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_gex.csv", 
                       package="compGenomRData")
file.exists(csvfile)
typeof(csvfile)
x1 <- read.csv(csvfile, row.names=1)
# Fix the gene names in the data frame
rownames(x1) <- sapply(strsplit(rownames(x1), "\\|"), function(x) x[1])
# Output a table
knitr::kable(head(t(head(x1))), caption="Example gene expression data (head)")
table(head(x1)) %>% as.data.frame() %>% gt()


# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_muts.csv", 
                       package="compGenomRData")
x2 <- read.csv(csvfile, row.names=1)
# Set mutation data to be binary (so if a gene has more than 1 mutation,
# we only count one)
x2[x2>0]=1
x2<- ifelse(x2>0, 1, 0)

# output a table
knitr::kable(head(t(head(x2))), caption="Example mutation data (head)")

# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_cnv.csv", 
                       package="compGenomRData")
x3 <- read.csv(csvfile, row.names=1)
# output a table
knitr::kable(head(t(head(x3))), 
             caption="Example copy number data for CRC samples")


# read in the csv from the companion package as a data frame
csvfile <- system.file("extdata", "multi-omics", "COREAD_CMS13_subtypes.csv",
                       package="compGenomRData")
covariates <- read.csv(csvfile, row.names=1)
# Fix the TCGA identifiers so they match up with the omics data
rownames(covariates) <- gsub(pattern = '-', replacement = '\\.',
                             rownames(covariates))

stringr::str_replace_all(rownames(covariates), "-", ".")

covariates <- covariates[colnames(x1),]
# create a dataframe which will be used to annotate later graphs
anno_col <- data.frame(cms=as.factor(covariates$cms_label))
rownames(anno_col) <- rownames(covariates)