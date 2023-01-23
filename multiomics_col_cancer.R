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
pacman::p_load(heatmaply, cowplot) #for making interactive heatmaps using plotly
pacman::p_load(d3heatmap, hrbrthemes) #for making interactive heatmaps using D3
pacman::p_load(knitr, tinytex, rmarkdown)
pacman::p_load(survminer)
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

rownames(covariates)<- stringr::str_replace_all(rownames(covariates), "-", ".")

covariates <- covariates[colnames(x1),]
# create a dataframe which will be used to annotate later graphs
anno_col <- data.frame(cms=as.factor(covariates$cms_label))
rownames(anno_col) <- rownames(covariates)


pheatmap::pheatmap(x1, show_rownames = F, show_colnames = F, annotation_col = anno_col)
pheatmap::pheatmap(x2, show_rownames = F, show_colnames = F, annotation_col = anno_col)
pheatmap::pheatmap(x3, show_rownames = F, show_colnames = F, annotation_col = anno_col)
pheatmap::pheatmap(dist(t(x1)), annotation_row = anno_col, show_rownames = F, show_colnames = F)
pheatmap::pheatmap(x2, show_rownames = F,
                   show_colnames = F,
                   annotation_col = anno_col,
                   cutree_cols = 3)

covariates %>% 
  dplyr::select(c(age, stage, pt, pn, kras_mut, braf_mut)) %>%
  cor(method = "spearman") %>%
  corrplot::corrplot(addCoef.col = "black")



plot_grid(
  ggplot(covariates, aes(age)) +
  geom_histogram( aes(y = ..density..), bins=30) +
    geom_density(alpha=0.1, fill="red")+
  theme_bw() +
  labs(title="Age")+
    panel_border(),
  
  ggplot(covariates, aes(gender)) +
    geom_bar() +
    geom_density()+
    theme_bw() +
    labs(title="Gender")+
    panel_border(),
  
  ggplot(covariates, aes(cms_label)) +
    geom_bar() +
    geom_density()+
    theme_bw() +
    labs(title="CMS")+
    panel_border(),
  
  ggplot(covariates, aes(cimp)) +
    geom_bar() +
    theme_bw() +
    labs(title="CIMP")+
    panel_border(),
  
  ggplot(covariates, aes(as.factor(stage))) +
    geom_bar() +
    geom_density()+
    theme_bw() +
    labs(title="Stage")+
    panel_border(),
  
  ggplot(covariates, aes(as.factor(kras_mut))) +
    geom_bar() +
    geom_density()+
    theme_bw() +
    labs(title="Kras_mut"),
  
  ggplot(covariates, aes(as.factor(braf_mut))) +
    geom_bar() +
    geom_density()+
    theme_bw() +
    labs(title="Braf_mut")+
    panel_border(),
  
  ggplot(covariates, aes(age, fill=as.factor(stage))) +
    geom_histogram(bins=10) +
    theme_bw() +
    labs(title="age vs stage")+
    panel_border(),
  
  ggplot(covariates, aes(as.factor(stage), age, fill=gender)) +
    geom_boxplot() +
    theme_bw() +
    labs(title="age v stage") +
    panel_border(),
  
  
  ggplot(covariates, aes(as.factor(osStat))) +
    geom_bar() +
    theme_bw() +
    labs(title="age v stage") +
    panel_border()
  
  
  )

splots <- list()
splots[[1]] <-survfit(Surv(osMo, osStat) ~ gender, data = covariates) %>% 
  ggsurvplot(data = covariates, pval = T, ggtheme = theme_bw())

splots[[2]] <-survfit(Surv(osMo, osStat) ~ stage, data = covariates) %>% 
  ggsurvplot(data = covariates, pval = T, ggtheme = theme_bw())

splots[[3]] <-survfit(Surv(osMo, osStat) ~ braf_mut, data = covariates) %>% 
  ggsurvplot(data = covariates, pval = T, ggtheme = theme_bw())

splots[[4]] <-survfit(Surv(osMo, osStat) ~ kras_mut, data = covariates) %>% 
  ggsurvplot(data = covariates, pval = T, ggtheme = theme_bw())

arrange_ggsurvplots(splots, print = TRUE,
                    ncol = 2, nrow = 2, risk.table.height = 0.4)
