##########################################################################################
########### RNA-sequencing Martina human cell samples treated and non-treated ###########
##########################################################################################

# Set working directory ----
setwd("~/Desktop/Current_work/martina_rnaseq/")

### Libraries ----
library(DESeq2) # DESeqDataSetFromTximport counts rlog summary estimateSizeFactors estimateDispersions plotDispEsts DESeq resultsNames results plotMA plotCounts
library(ggplot2) # ggtitle scale_color_discrete ggplot aes geom_bar facet_wrap labs theme_bw coord_flip geom_boxplot position_dodge xlab ylab scale_y_log10 theme element_blank element_text scale_fill_manual geom_point scale_colour_manual xlim geom_hline geom_vline
library(RColorBrewer) # brewer.pal
library(MetBrewer) # met.brewer
library(pheatmap) # pheatmap
library(PoiClaClu) # PoissonDistance
library(tximport) # tximport 
library(factoextra) # fviz_screeplot fviz_pca_ind
library(dplyr) # %>% as_tibble
library(readxl) # read_excel
library(VennDiagram) # venn.diagram
library(clusterProfiler) # %>% bitr enrichGO goplot cnetplot emapplot
library(org.Hs.eg.db) # No used functions found
library(biomaRt) # useDataset useMart
library(enrichplot) # ggtitle goplot cnetplot pairwise_termsim emapplot
library(plotly)
library(stringr)



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
### Same thing with the poisson method
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste(rld$dex, rld$cell, sep=" - ")
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         show_colnames = T,
         fontsize = 18)



#### Differential Gene Expression Analysis ----

# Differential expression (used Benjamin-Hochberg adjustment and a p-val of <0.1)
dds = DESeq2::DESeq(dds, parallel=TRUE)
# Show name of the results (all comparisons)
resultsNames(dds)


### Explore Test treatment_treated_vs_control: treated=+Log2FC and Control=-Log2FC
res=results(dds, name="treatment_treated_vs_control")
summary(res)
#distribution of coefficents of the model
plotMA(res, ylim=c(-5,5),cex=1.0, cex.lab=1.45, 
       cex.axis=1.45, xlab='Mean of Normalised Counts',
       ylab='Log2 Fold Change')
# plot of p-vals excluding genes with very small counts
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
     cex.axis=1.45, cex.lab=1.45, main="")
# Order the results by fold change to see the biggest changing genes
res_ordered=res[order(res$log2FoldChange),]
head(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_ordered$gene<-rownames(res_ordered)
rownames(res_ordered)<-c()
#write.table(as.data.frame(res_ordered), file="diff_exp_output_log2FC_treatment_treated_vs_control.txt",sep="\t", quote = F, col.names = T, row.names = F)

#out of 70094 with nonzero total read count
res_significant<- subset(res, log2FoldChange > 2 | log2FoldChange < -2)
res_significant <- subset(res_significant, padj < 0.01)
nrow(res_significant) #2100 total differentially expressed genes
nrow(res_significant[res_significant$log2FoldChange > 2,]) #1452 upregulated in treated
nrow(res_significant[res_significant$log2FoldChange < -2,]) #648 upregulated in control



### Explore Test treatment_ttt2_vs_control: ttt1=+Log2FC and Control=-Log2FC
#res2=results(dds, name="treatment_ttt2_vs_control")
#summary(res2)
#distribution of coefficents of the model
#plotMA(res2, ylim=c(-5,5),cex=1.0, cex.lab=1.45, 
#       cex.axis=1.45, xlab='Mean of Normalised Counts',
#       ylab='Log2 Fold Change')
# plot of p-vals excluding genes with very small counts
#hist(res2$pvalue[res$baseMean > 1], breaks = 0:20/20,
#     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
#     cex.axis=1.45, cex.lab=1.45, main="")
# Order the results by fold change to see the biggest changing genes
#res2_ordered=res2[order(res2$log2FoldChange),]
#head(res2_ordered)
#res2_ordered<-as.data.frame(res2_ordered)
#res2_ordered$gene<-rownames(res2_ordered)
#rownames(res2_ordered)<-c()
#write.table(as.data.frame(res_ordered), file="diff_exp_output_log2FC_treatment_ttt1_vs_control.txt",sep="\t", quote = F, col.names = T, row.names = F)

#out of 70094 with nonzero total read count
#res2_significant<- subset(res2, log2FoldChange > 2 | log2FoldChange < -2)
#res2_significant <- subset(res2_significant, padj < 0.01)
#nrow(res2_significant) #1603 total differentially expressed genes
#nrow(res2_significant[res2_significant$log2FoldChange > 2,]) #868 upregulated in ttt2
#nrow(res2_significant[res2_significant$log2FoldChange < -2,]) #735 upregulated in control



### Explore Test treatment_ttt1_vs_control: ttt1=+Log2FC and Control=-Log2FC
#res3=results(dds, name="treatment_ttt3_vs_control")
#summary(res3)
#distribution of coefficents of the model
#plotMA(res3, ylim=c(-5,5),cex=1.0, cex.lab=1.45, 
#       cex.axis=1.45, xlab='Mean of Normalised Counts',
#       ylab='Log2 Fold Change')
# plot of p-vals excluding genes with very small counts
#hist(res3$pvalue[res$baseMean > 1], breaks = 0:20/20,
#     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
#     cex.axis=1.45, cex.lab=1.45, main="")
# Order the results by fold change to see the biggest changing genes
#res3_ordered=res[order(res3$log2FoldChange),]
#head(res3_ordered)
#res3_ordered<-as.data.frame(res3_ordered)
#res3_ordered$gene<-rownames(res3_ordered)
#rownames(res_ordered3)<-c()
#write.table(as.data.frame(res_ordered), file="diff_exp_output_log2FC_treatment_ttt1_vs_control.txt",sep="\t", quote = F, col.names = T, row.names = F)

#out of 70094 with nonzero total read count
#res3_significant<- subset(res3, log2FoldChange > 2 | log2FoldChange < -2)
#res3_significant <- subset(res3_significant, padj < 0.01)
#nrow(res3_significant) #1526 total differentially expressed genes
#nrow(res3_significant[res3_significant$log2FoldChange > 2,]) #1197 upregulated in ttt2
#nrow(res3_significant[res3_significant$log2FoldChange < -2,]) #329 upregulated in control






#### treated vs. treated_dexa
res1vs2=results(dds, contrast = c('treatment', 'treated', 'treated_dexa'))
summary(res1vs2)
#out of 70094 with nonzero total read count
res1vs2_significant<- subset(res1vs2, log2FoldChange > 2 | log2FoldChange < -2)
res1vs2_significant <- subset(res1vs2_significant, padj < 0.01)
nrow(res1vs2_significant) #1161 total differentially expressed genes
nrow(res1vs2_significant[res1vs2_significant$log2FoldChange > 2,]) #745 upregulated in treated_dexa
nrow(res1vs2_significant[res1vs2_significant$log2FoldChange < -2,]) #416 upregulated in treated



#### treated vs. treated_IKMO
res1vs3=results(dds, contrast = c('treatment', 'treated', 'treated_IKMO'))
summary(res1vs3)
res1vs3_significant<- subset(res1vs3, log2FoldChange > 2 | log2FoldChange < -2)
res1vs3_significant <- subset(res1vs3_significant, padj < 0.01)
nrow(res1vs3_significant) #416 total differentially expressed genes
nrow(res1vs3_significant[res1vs3_significant$log2FoldChange > 2,]) #165 upregulated in treated_IKMO
nrow(res1vs3_significant[res1vs3_significant$log2FoldChange < -2,]) #251 upregulated in treated


#### Plots for visualisation ----

### Venn Diagram 
# Get all genes DEG ttt vs. control
list_ttt1=rownames(res_significant)
list_ttt1vs2=rownames(res1vs2_significant)
list_ttt1vs3=rownames(res1vs3_significant)

a=met.brewer('Egypt', n=3, 'discrete')
# Generate the plot 
venn.diagram(x=list(list_ttt1, list_ttt1vs2, list_ttt1vs3), 
             category.names = c('U_vs_T', 'T_vs_TD', 'T_vs_TK'),
             fill=a,
             col=a,
             cat.fontface = "bold",
             fontface = "bold",
             filename = 'VennDiagram')


### Heatmap

###

#ATTENTION: Here I just use the most differentially express between control and treated (T)

###

# set number of lines to show
n=100
# take the best DEG 
topdiff = head(c(1:nrow(res))[order(res$padj)],n)
# set color palette
a=met.brewer('Egypt', n=4, 'discrete')
my_colors = list(treatment = c(control = a[3] , treated =a[1], treated_dexa=a[2], treated_IKMO=a[4]))
# prepare dataset for plotting
mat = assay(rld)[topdiff, ]
mat = mat - rowMeans(mat)
df2 = as.data.frame(colData(rld)[,c("treatment"),drop=FALSE])
colnames(df2)<-c("treatment")
# Plot the heatmap
pheatmap(mat, annotation_col=df2,
         show_rownames = F,annotation_colors = my_colors,
         col=met.brewer("Manet",n=100, type="continuous"),
         legend_breaks = c(-0.2, 0, 0.2,0.4, max(mat)),
         legend_labels = c("-0.2", "0", "0.2","0.4", "z-score"),
         fontsize = 16, annotation_legend = F, angle_col = 0, annotation_names_col = F, border_color = "white")

### Boxplot best differentially expressed genes

###

#ATTENTION: Here I just use the most diffrentially express between control and ttt1 (T)

###

n=20
selGenes = head(rownames(res)[order(res$padj)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene,plotCounts(dds, gene=gene, intgroup=c("treatment"), returnData=TRUE))))
ggplot(data, aes(x=treatment, y=count, fill=treatment)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) + facet_wrap(~gene) +
  xlab("Treatment") + ylab("Normalised read count") + 
  scale_y_log10() + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22),
        legend.position = 'none')+
  scale_fill_manual("", breaks=c("control","treated", "treated_dexa", "treated_IKMO"),values = c(a[3],a[1], a[2], a[4]))

### Volcano plot

###

#ATTENTION: Here I just use the most diffrentially express between control and ttt1 (T)

###

res_df <- as.data.frame(res)
res_df$gene <- row.names(res_df)

res_significant_df <- as.data.frame(res_significant)
res_significant_df$gene <- row.names(res_significant_df)

res_df$sig <- "no"
res_df$sig[res_df$gene %in% res_significant_df$gene] <- "yes"

#v=ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig, text = paste("Symbol:", gene))) + 
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig, text = paste("Symbol:", gene))) + 
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
        legend.position = "none")
#ggplotly(v)





#### GO term enrichment analysis ----

# Get Db and correspondence with transcriptomes ID
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# list genes treated vs. control
refseq <- as.vector(rownames(res_significant))
# change ID
for (i in 1:length(refseq)){
  refseq[i]=unlist(str_split(refseq[i], '\\.')[1])
}
# Get ENSEMBL IDs
gene.df <- bitr(refseq, fromType = "REFSEQ", toType = c("ENSEMBL",'UNIPROT'), OrgDb = org.Hs.eg.db)

# Do GO enrichment analysis 
enrich_analysis=enrichGO(gene=gene.df$ENSEMBL, 
                         OrgDb = org.Hs.eg.db, 
                         keyType = 'ENSEMBL',
                         ont='MF', 
                         pAdjustMethod = 'BH', 
                         pvalueCutoff = 0.01, 
                         qvalueCutoff = 0.05, 
                         readable = T)

# Basic plot GO analysis
goplot(enrich_analysis)
# Other plot 
cnetplot(enrich_analysis, showCategory = 20, node_label="gene")
# Best representation (Own preference)
edo <- pairwise_termsim(enrich_analysis)
emapplot(edo, node_label="category", cex_label_category=0.7)



#### KEGG analysis ----
gene.df.bis = bitr_kegg(gene.df$UNIPROT, fromType = "uniprot", toType = "ncbi-geneid", organism     = 'hsa')
kk <- enrichKEGG(gene         = gene.df.bis$`ncbi-geneid`,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
kp=which(kk@result$ID=='hsa00380')
kk@result[kp,]

kk2 <- gseKEGG(geneList     = gene.df.bis$`ncbi-geneid`,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

mkk2 <- gseMKEGG(geneList = gene.df.bis$`ncbi-geneid`,
                 organism = 'hsa',
                 pvalueCutoff = 1)

mkk <- enrichMKEGG(gene = gene.df.bis$`ncbi-geneid`,
                   organism = 'hsa',
                   pAdjustMethod = "fdr",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 1)
?enrichMKEGG
head(mkk)
kp=which(mkk@result$ID=='M00038')
mkk@result[kp,]
library("pathview")
hsa00380 <- pathview(gene.data  = gene.df.bis$`ncbi-geneid`,
                     pathway.id = "hsa00380",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1))
?enrichKEGG

head(kk)

