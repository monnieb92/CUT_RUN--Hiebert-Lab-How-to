```{r, echo=FALSE}
library(tidyverse)
library(viridis)
library(ggplot2)
ibrary(ggplot2)
library(ComplexHeatmap)
library(circlize)
library("DESeq2")
library(tidyverse)
library(readr)
```


```{r}
theme_monica <- function(){
  theme_classic() %+replace%    #replace elements we want to change
    #font <- "Times", 
    theme(
      
  plot.title = element_text(             #title
                  family = "sans",              #set font family
                   size = 20,                #set font size
                   face = 'bold',            #bold typeface
                   hjust = 0.5,                #left align
                   vjust = 2),               #raise slightly
      
plot.subtitle = element_text(          #subtitle
                   family = "sans",            #font family
                   size = 14),               #font size
      
plot.caption = element_text(           #caption
                   family = "sans",            #font family
                   size = 10,                 #font size
                   hjust = 1),               #right align
      legend.text = element_text(             #axis titles
                   family =  "sans",            #font family
                   face = 'bold',
                   size = 12),
      axis.title = element_text(             #axis titles
                   family = "sans",            #font family
                   face = 'bold',
                   size = 16),               #font size
      
      axis.text = element_text(              #axis text
                   family = "sans",            #axis famuly
                   face = 'bold',
                   size = 16),                #font size
      
      #axis.text.x = element_text(            #margin for axis text
                    #margin=margin(10, b = 12))
      
    )}
```
```{r}
##GB CHANGES
PS3exps_KasT6h_4repsgb <- PS3exps_KasT6h_4reps %>% 
  mutate(FoldChange_gb = ((2^(abs(log2FoldChange_gb)))*sign(log2FoldChange_gb)))

PS3exps_KasT6h_4repsgb <- transform(PS3exps_KasT6h_4repsgb, category = ifelse((padj_gb <= 0.05 &  FoldChange_gb <= -2), "Down", ifelse((padj_gb<= 0.05 &  FoldChange_gb >= 2), "Up",  "NS")))


PS3exps_KasT6h_4repsgbd <- PS3exps_KasT6h_4repsgb %>% mutate(avgGBD = (gbd_3191.MB.3.hg19.sorted.F4q10.BLfiltered+gbd_3191.MB.4.hg19.sorted.F4q10.BLfiltered+gbd_3568.MB.1.hg19.sorted.F4q10.BLfiltered+gbd_3568.MB.2.hg19.sorted.F4q10.BLfiltered+gbd_5641.MB.1.hg19.F4q10.sorted.BLfiltered+gbd_5641.MB.3.hg19.F4q10.sorted.BLfiltered+gbd_5641.MB.4.hg19.F4q10.sorted.BLfiltered+gbd_3568.MB.5.hg19.sorted.F4q10.BLfiltered+gbd_3568.MB.6.hg19.sorted.F4q10.BLfiltered)/9)

PS3exps_KasT6h_4repsgbd <- PS3exps_KasT6h_4repsgbd %>% dplyr::select(log2FoldChange_gb,padj_gb,FoldChange_gb,category,avgGBD,pvalue_gb)


PS3exps_KasT6h_4repsgbd <- na.omit(PS3exps_KasT6h_4repsgbd)

PS3exps_KasT6h_4repsgbdMAplot <- ggplot(PS3exps_KasT6h_4repsgbd, aes(y = log2FoldChange_gb,x = avgGBD, color = category)) +   geom_point( size = 0.35) + theme_monica() + 
  theme(panel.border = element_rect(fill=NA, colour = "black", size=1)) +
  scale_color_manual(breaks = c("Up", "NS", "Down"), values=c('#ca0020','#bababa','#0571b0')) +
  geom_hline(yintercept=0, linetype="solid", color = "black", size=1) +
  scale_y_continuous(breaks=seq(-4, 4, 1), limits=c(-4,4)) +
  #scale_x_continuous(breaks=seq(0, 650, 50), limits=c(0,655)) 
  scale_x_log10(breaks = c(0.1,1,10,100)) +
  expand_limits(x=c(0,100))

ggsave("/Volumes/Bomber4TB/SMARCA5_PRO-seq_test/3191_5641_3568_T0hvT6hdT.cgbd0.004/3191-3568-5641-T0h.vs.T6hdT4reps.MAplot.GBD.log10xaxis.tiff", plot = PS3exps_KasT6h_4repsgbdMAplot,width = 4.25, height = 3, dpi = 300)

PS3exps_KasT6h_4repsgbd <- as.data.frame(PS3exps_KasT6h_4repsgbd)
dplyr::count(PS3exps_KasT6h_4repsgbd, category)

print(PS3exps_KasT6h_4repsgbdMAplot)
```



## Heatmap FPKM gb up gb down RNAseq
```{r}
set.seed(123)
Mtg16Nullrestable.sigchanges <- as_data_frame(Mtg16Nullrestable.sigchanges)
AA=Mtg16Nullrestable.sigchanges
AA
rnames=Mtg16Nullrestable.sigchanges$genenames...1.

mat_data=data.matrix(AA[,25:40])
mat_data
glimpse(mat_data)
rownames(mat_data)=rnames

head(mat_data)

#dend2h = hclust(dist(mat_data), method = "ward.D2")
#rownames(mat_scaled)=rnames

my_palette2 = circlize::colorRamp2(c(-2,-1.5,-1 ,-0.5, 0, 0.5,1,1.5,2), c("#2166ac", "#4393c3","#92c5de","#d1e5f0","#f7f7f7","#fddbc7","#f4a582","#d6604d","#b2182b"))

heatmap_sigppdn_log2foldclust <- Heatmap(mat_data, col = my_palette2, column_title = "Mtg16 Changed Genes", cluster_columns = FALSE, clustering_method_rows = "ward.D2",  cluster_row_slices = FALSE, show_row_names = FALSE, show_column_names = TRUE,split=Mtg16Nullrestable.sigchanges$Category,row_names_gp = gpar(fontsize = 5), width = unit(50, "mm"))
heatmap_sigppdn_log2foldclust

pdf("/Volumes/Bomber4TB/RNAseq8127-Mouse/8127-AllNulls_andhomo-R/Mtg16RNA-seq_log2.FPKMoveravgWT.wardD2clustering.pdf")
heatmap_sigppdn_log2foldclust
dev.off()


heatmap_sigppdn_log2fold <- Heatmap(mat_data, col = my_palette2, column_title = "Mtg16 Changed Genes", cluster_columns = FALSE, cluster_rows = FALSE, cluster_row_slices = FALSE, show_row_names = FALSE, show_column_names = TRUE,split=Mtg16Nullrestable.sigchanges$Category,row_names_gp = gpar(fontsize = 5), width = unit(50, "mm"))
heatmap_sigppdn_log2fold

pdf("/Volumes/Bomber4TB/RNAseq8127-Mouse/8127-AllNulls_andhomo-R/Mtg16RNA-seq_log2.FPKMoveravgWT.noclustering.pdf")
heatmap_sigppdn_log2fold
dev.off()

```
#PCA plot

```{r}
result_DE <- results(dds_DE)
vsd_DE <- varianceStabilizingTransformation(dds_DE)

PCA <- plotPCA(vsd_DE, returnData = TRUE)
glimpse(PCA)
```

```{r, fig.align="center", out.width = "80%"}
percentVar <- round(100 * attr(PCA, "percentVar"))

virdispalette_6 <- viridis(10)
  
PCAplot <- ggplot(PCA, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size=3) +
  geom_text_repel(color = "black") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c(virdispalette_6)) +
  ggtitle("SMARCA5_PCA_allreps")+
  theme_classic()
PCAplot
```
