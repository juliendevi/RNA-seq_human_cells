##########################################################################################
########### RNA-sequencing Martina human cell samples treated and non-treated ###########
##########################################################################################

# Set working directory ----
setwd("~/Desktop/Current_work/martina_rnaseq/")

### Libraries ----
library(DESeq2)          # Bioconductor v1.36.0 
library(ggplot2)         # CRAN v3.3.6 
library(RColorBrewer)    # CRAN v1.1-3 
library(MetBrewer)       # [github::BlakeRMills/MetBrewer] v0.2.0 
library(pheatmap)        # CRAN v1.0.12 
library(PoiClaClu)       # CRAN v1.0.2.1 
library(tximport)        # Bioconductor v1.24.0 
library(factoextra)      # CRAN v1.0.7 
library(dplyr)           # CRAN v1.0.10 
library(readxl)          # CRAN v1.4.1 
library(VennDiagram)     # CRAN v1.7.3 
library(clusterProfiler) # Bioconductor v4.4.4 
library(org.Hs.eg.db)    # Bioconductor v3.15.0 
library(biomaRt)         # [Annotation::NA/NA] v2.52.0 
library(enrichplot)      # Bioconductor v1.16.2 
library(plotly)          # CRAN v4.10.0  
library(stringr)         # CRAN v1.4.1 
library(cowplot)         # CRAN v1.1.1 



### Import kallisto outputs and experiment information ----

# import experimental scheme
exp = read_excel("~/Desktop/Current_work/martina_rnaseq/experiment.xlsx", col_names = T)
# Make experiment design matrix
treatment = as.factor(exp$treatment)
design = model.matrix(~0 + treatment)
colnames(design)=levels(treatment)
# find file path
path = file.path(exp$sequence, "abundance.tsv")
# check path exist
all(file.exists(path))
# info samples
sample_info <- data.frame(exp$sample, exp$treatment, path)
rownames(sample_info) <- sample_info$exp.sample
colnames(sample_info)=c('sample', 'treatment', 'path')
# import kallisto transcripts
kallisto_out = tximport(path, type = 'kallisto', txOut = T, geneIdCol = T)

### Set up normalisation DESeq2 ----

# Make deseq2 object
dds <- DESeqDataSetFromTximport(kallisto_out,
                                colData = sample_info,
                                design = ~ treatment)

# remove features with low counts
dds = dds[ rowMeans(counts(dds)) > 10, ] 
nrow(dds) #177816 => 70094

# rlog transform counts (by average Tx length and correcting for library size)
rld = rlog(dds, blind=FALSE)




#### Exploration dataset with PCA and clustering ----

### PCA
# Do the PCA
x_pca = prcomp(t(assay(rld)), scale. = T)
# Summary PCA
summary(x_pca)
# Exploring axes explanation
fviz_screeplot(x_pca, addlabels = T) + ggtitle('')
fviz_pca_ind(x_pca, asp=1, pointshape=20, habillage = treatment) + scale_color_discrete() + ggtitle('PCA transcriptomes')
# Do dataframe to visualise all axis at once
pca_df <- x_pca$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sample_info$sample,
             group = sample_info$treatment)
# Pivot the dataframe for plotting
pca_pivot <- pivot_longer(pca_df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)
# Plot one graph by PCA axis
ggplot(pca_pivot) +
  aes(x=sample, y=loadings, fill=group) + 
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA plot for each axis") +
  theme_bw() +
  coord_flip()

### Hierarchical clustering 
# get distance
distance = dist(t(assay(rld)), method = "maximum")
# get cluster
cluster = hclust(distance, method = "ward.D2")
# plot
plot(cluster, labels=exp$sample)




#### Checks DESeq normalisation ----

### two 1st samples plotted against each other to check consistency (for rlog and log2) 
par( mfrow = c( 1, 2 ) )
dds = estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,c(1,2)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
plot(assay(rld)[,c(1,2)],
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
par(mfrow = c(1,1))

### estimate size factors = normalize for dispersion
dds = DESeq2::estimateSizeFactors(dds)
dds = estimateDispersions(dds)
plotDispEsts(dds, xlab= "Mean of Normalised Counts", ylab= "Dispersion", cex=1.0, cex.lab=1.45, cex.axis=1.45)

### check sample distances
# Get distances
sampleDists <- dist(t(assay(rld)))
# Put in a matrix
sampleDistMatrix <- as.matrix(sampleDists)
# Add names
rownames(sampleDistMatrix) <- paste(rld$dex, rld$cell, sep = " - ")
# Choose colors
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
# Visualisation heatmap distances
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_colnames = T,
         fontsize = 18)




#### Differential Gene Expression Analysis ----

# Differential expression (used Benjamin-Hochberg adjustment and a p-val of <0.1)
dds = DESeq2::DESeq(dds, parallel=TRUE)

# Get Human Db
mart=useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# Function get gene name and uniprot ID
getgeneID <- function(x){
  # List transcripts ID
  refseq=as.vector(x$gene)
  # Change transcripts ID for gene ID
  for (i in 1:length(refseq)){
    refseq[i]=unlist(str_split(refseq[i], '\\.')[1])
  }
  gene.df=bitr(refseq, fromType = "REFSEQ", toType = c("ENSEMBL",'UNIPROT', 'ENZYME', 'GENENAME'), OrgDb = org.Hs.eg.db)
  # Match transcript ID to 'gene name' & 'uniprot ID'
  for (i in 1:nrow(x)){
    gene = unlist(str_split(x[i,7], '\\.')[1])[1]
    n = which(gene.df$REFSEQ==gene)
    if(is.null(n)){
      x[i,8]=NA
      x[i,9]=NA
    }
    else{
      x[i,8]=paste(gene.df$UNIPROT[n], collapse = ',')
      x[i,9]=gene.df$GENENAME[n][1]
    }
  }
  # Put the 2 new column names
  colnames(x)=c(colnames(x[1:7]), 'Uniprot', 'Genename')
  return(x)
}




### Explore Test treatment_treated_vs_control: treated=+Log2FC and Control=-Log2FC ----
T_vs_C=results(dds, name="treatment_treated_vs_control")
summary(T_vs_C) # 70094 with nonzero total read count
# Order the results
T_vs_C_sign<- subset(T_vs_C, log2FoldChange > 2 | log2FoldChange < -2) # LogFC threshold
T_vs_C_sign <- subset(T_vs_C_sign, padj < 0.01) # adjusted p-value threshold
T_vs_C_sign_ordered=T_vs_C_sign[order(T_vs_C_sign$padj),] # Order by adjusted p-value
T_vs_C_sign_ordered<-as.data.frame(T_vs_C_sign_ordered) # transform in data frame
T_vs_C_sign_ordered$gene<-rownames(T_vs_C_sign_ordered) # put row names in gene column
rownames(T_vs_C_sign_ordered)<-c() # remove row names
T_vs_C_sign_ordered=getgeneID(T_vs_C_sign_ordered) # Get gene IDs
#write.table(as.data.frame(T_vs_C_sign_ordered), file="diff_exp_output_log2FC_treatment_treated_vs_control.txt",sep="\t", quote = F, col.names = T, row.names = F)

# number of genes biased for each treatment: 2100 DEG, 1452 upregulated in treated, 648 upregulated in control
summary(T_vs_C_sign)



### Explore Test treatment_treated_vs_treated_dexa: treated=+Log2FC and treated_dexa=-Log2FC ----
T_vs_Tdexa=results(dds, contrast = c('treatment', 'treated', 'treated_dexa'))
summary(T_vs_Tdexa) # 70094 with nonzero total read count
T_vs_Tdexa_sign <- subset(T_vs_Tdexa, log2FoldChange > 2 | log2FoldChange < -2)
T_vs_Tdexa_sign <- subset(T_vs_Tdexa_sign, padj < 0.01)
T_vs_Tdexa_sign_ordered=T_vs_Tdexa_sign[order(T_vs_Tdexa_sign$padj),]
T_vs_Tdexa_sign_ordered<-as.data.frame(T_vs_Tdexa_sign_ordered)
T_vs_Tdexa_sign_ordered$gene<-rownames(T_vs_Tdexa_sign_ordered)
rownames(T_vs_Tdexa_sign_ordered)<-c()
T_vs_Tdexa_sign_ordered=getgeneID(T_vs_Tdexa_sign_ordered)
#write.table(as.data.frame(T_vs_Tdexa_sign_ordered), file="diff_exp_output_log2FC_treatment_treated_vs_treated_dexa.txt",sep="\t", quote = F, col.names = T, row.names = F)

# number of genes biased for each treatment
summary(T_vs_Tdexa_sign) # 1161 DEG, 745 upregulated in treated, 416 upregulated in treated_dexa



### Explore Test treatment_treated_vs_treated_IKMO: treated=+Log2FC and treated_IKMO=-Log2FC ----
T_vs_TIKMO=results(dds, contrast = c('treatment', 'treated', 'treated_IKMO'))
summary(T_vs_TIKMO) # 70094 with nonzero total read count
T_vs_TIKMO_sign<- subset(T_vs_TIKMO, log2FoldChange > 2 | log2FoldChange < -2)
T_vs_TIKMO_sign <- subset(T_vs_TIKMO_sign, padj < 0.01)
T_vs_TIKMO_sign_ordered=T_vs_TIKMO_sign[order(T_vs_TIKMO_sign$padj),]
T_vs_TIKMO_sign_ordered<-as.data.frame(T_vs_TIKMO_sign_ordered)
T_vs_TIKMO_sign_ordered$gene<-rownames(T_vs_TIKMO_sign_ordered)
rownames(T_vs_TIKMO_sign_ordered)<-c()
T_vs_TIKMO_sign_ordered=getgeneID(T_vs_TIKMO_sign_ordered)
#write.table(as.data.frame(T_vs_TIKMO_sign_ordered), file="diff_exp_output_log2FC_treatment_treated_vs_treated_IKMO.txt",sep="\t", quote = F, col.names = T, row.names = F)

# number of genes biased for each treatment
summary(T_vs_TIKMO_sign) # 416 DEG, 165 upregulated in treated, 251 upregulated in treated_IKMO




#### Plots for visualization ----

### Venn Diagram ----
# Get all genes DEG 
list_T_vs_C=rownames(T_vs_C_sign)
list_T_vs_Tdexa=rownames(T_vs_Tdexa_sign)
list_T_vs_TIKMO=rownames(T_vs_TIKMO_sign)
# Choose colour pal
a=met.brewer('Egypt', n=3, 'discrete')
# Generate the plot 
venn.diagram(x=list(list_T_vs_C, list_T_vs_Tdexa, list_T_vs_TIKMO), 
             category.names = c('T_vs_C', 'T_vs_TD', 'T_vs_TK'),
             fill=a,
             col=a,
             cat.fontface = "bold",
             fontface = "bold",
             filename = 'VennDiagram')


### Volcano plot ----

## T vs. C ----
# transform for df
T_vs_C.df <- as.data.frame(T_vs_C)
T_vs_C.df$gene <- row.names(T_vs_C.df)
# transform for df
T_vs_C_sign.df <- as.data.frame(T_vs_C_sign)
T_vs_C_sign.df$gene <- row.names(T_vs_C_sign.df)
# put significatif or not
T_vs_C.df$sig <- "no"
T_vs_C.df$sig[T_vs_C.df$gene %in% T_vs_C_sign.df$gene] <- "yes"

#v=ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig, text = paste("Symbol:", gene))) + 
T_vs_C_plot=ggplot(T_vs_C.df, aes(x = log2FoldChange, y = -log10(padj), colour = sig, text = paste("Symbol:", gene))) + 
  geom_point(size=1.5)+
  scale_colour_manual("", breaks=c("no","yes"),
                      values = c("black","red"))+
  xlim(-10,10)+
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 2, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -2, linetype="longdash", colour="#2C467A", size=1) +
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = "none")+
  coord_flip()
#ggplotly(v)

## T vs. Tdexa ----
# transform for df
T_vs_Tdexa.df <- as.data.frame(T_vs_Tdexa)
T_vs_Tdexa.df$gene <- row.names(T_vs_Tdexa.df)
# transform for df
T_vs_Tdexa_sign.df <- as.data.frame(T_vs_Tdexa_sign)
T_vs_Tdexa_sign.df$gene <- row.names(T_vs_Tdexa_sign.df)
# put significatif or not
T_vs_Tdexa.df$sig <- "no"
T_vs_Tdexa.df$sig[T_vs_Tdexa.df$gene %in% T_vs_Tdexa_sign.df$gene] <- "yes"
# Plot 
T_vs_Tdexa_plot=ggplot(T_vs_Tdexa.df, aes(x = log2FoldChange, y = -log10(padj), colour = sig, text = paste("Symbol:", gene))) + 
  geom_point(size=1.5)+
  scale_colour_manual("", breaks=c("no","yes"),
                      values = c("black","red"))+
  xlim(-10,10)+
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 2, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -2, linetype="longdash", colour="#2C467A", size=1) +
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = "none")+
  coord_flip()


## T vs. TIKMO ----
# transform for df
T_vs_TIKMO.df <- as.data.frame(T_vs_TIKMO)
T_vs_TIKMO.df$gene <- row.names(T_vs_TIKMO.df)
# transform for df
T_vs_TIKMO_sign.df <- as.data.frame(T_vs_TIKMO_sign)
T_vs_TIKMO_sign.df$gene <- row.names(T_vs_TIKMO_sign.df)
# put significatif or not
T_vs_TIKMO.df$sig <- "no"
T_vs_TIKMO.df$sig[T_vs_TIKMO.df$gene %in% T_vs_TIKMO_sign.df$gene] <- "yes"
# Plot
T_vs_TIKMO_plot=ggplot(T_vs_TIKMO.df, aes(x = log2FoldChange, y = -log10(padj), colour = sig, text = paste("Symbol:", gene))) + 
  geom_point(size=1.5)+
  scale_colour_manual("", breaks=c("no","yes"),
                      values = c("black","red"))+
  xlim(-10,10)+
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 2, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -2, linetype="longdash", colour="#2C467A", size=1) +
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = "none")+
  coord_flip()


## Plot the 3 volcano plots ----
plot_grid(T_vs_Tdexa_plot,T_vs_C_plot, T_vs_TIKMO_plot, labels=c('T vs T dexa', 'C vs T', 'T vs T IKMO'),vjust=1.2, rel_widths = c(1,2), scale = c(0.9, 0.9, 0.9))

















