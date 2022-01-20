# 08.0.2020 TM #

library(tidyverse)
library(limma)
library(plotly)
library(pheatmap)
library(wesanderson)
library(pcaMethods)
library(car)
library(impute)

##### load file ################################################################

RawD = read_delim("RES_wSites.csv", delim=",")


##### Function definitions ####################################################

# FC cal 
fc_cal = function(row_before = NULL, row_after = NULL) {
  
  before = unlist(row_before)
  after = unlist(row_after)
  
  FC = "ERROR"
  
  if( is.numeric(before) & is.numeric(after)  ) {
    
    if( sum(is.na(after)) == length(after) & sum(is.na(before)) < length(before)  ) { FC = 0.01 }
    if( sum(is.na(after)) < length(after) & sum(is.na(before)) == length(before)  ) { FC = 100 }
    if( sum(is.na(after)) == length(after) & sum(is.na(before)) == length(before)  ) { FC = NA }
    
    if( sum(is.na(after)) < length(after) & sum(is.na(before)) < length(before)  ) { FC = median(after, na.rm=T)/median(before, na.rm=T) }
    
  }
  
  return(FC)
  
}

# t-test
t_test = function(row_before=NULL, row_after=NULL) {
  
  before = unlist(row_before)
  after = unlist(row_after)
  
  pvalue = 6666666666
  
  if( !is.character(before) & !is.character(after)  ) {
    
    if( (sum(is.na(after)) < length(after) &
         sum(is.na(before)) < length(before) ) &
        ( sum( is.na( c(before) ) ) < 0.6*length(before) | sum( is.na( c(after) ) ) < 0.6*length(after) ) ) 
    {ttest = t.test(x=after, y=before, var.equal = TRUE, alternative = "two.sided")
    pvalue = ttest$p.value}

    if( (sum(is.na(after)) < length(after) & 
         sum(is.na(before)) < length(before) ) &
        ( sum( is.na( c(before) ) ) >= 0.6*length(before) | sum( is.na( c(after) ) ) >= 0.6*length(after) )  ) 
    { pvalue = 1 }
    
    if( (sum(is.na(after)) == length(after) | 
         sum(is.na(before)) == length(before) ) &
        ( sum( is.na( c(before) ) ) >= 0.6*length(before) | sum( is.na( c(after) ) ) >= 0.6*length(after) )  ) 
    { pvalue = 1 }
        
    if( ( sum(is.na(after)) == length(after) &
          sum(is.na(before)) < length(before) ) &
        ( sum( is.na( c(before) ) ) < 0.6*length(before) | sum( is.na( c(after) ) ) < 0.6*length(after) ) )  
    { pvalue = 0.00001 }
    
    if( ( sum(is.na(after)) < length(after) &
          sum(is.na(before)) == length(before) ) &
        ( sum( is.na( c(before) ) ) < 0.6*length(before) | sum( is.na( c(after) ) ) < 0.6*length(after) ) ) 
    { pvalue = 0.00001 }
     
  }
  
  return(pvalue)
  
}




##### make results table #######################################################


# normalize data

RawD_norm = RawD
RawD_norm[6:53] = normalizeVSN(RawD[6:53])
RawD_norm[6:53] = apply(RawD_norm[6:53], c(1,2), function(x) 2^x )



# check intensity distribution before and after normalization

raw_int = RawD[6:53] %>% as_vector() %>% log10()
norm_int = RawD_norm[6:53] %>% as_vector() %>% log10()


ggplot() +
  geom_histogram( aes(x=raw_int), alpha=0.75, fill="darkred", bins = 100 ) +
  geom_histogram( aes(x=norm_int), alpha=0.75, fill="navy", bins = 100 ) +
  labs(title = "log10 of intensities",
       subtitle = "before (red) and after (blue) variance stability normalization") +
  theme_minimal()

ggsave("00_intensity_VSnormalization.png", width=6, height=6.5)


# calculate AKT stim vs WT
RES_AKT = RawD_norm %>%
  na_if(0) %>%
  mutate(
    
    FC = map2_dbl(split(RawD_norm[18:29], 1:nrow(RawD_norm)), split(RawD_norm[6:17], 1:nrow(RawD_norm)), fc_cal ),
    pval = map2_dbl(split(RawD_norm[18:29], 1:nrow(RawD_norm)), split(RawD_norm[6:17], 1:nrow(RawD_norm)), t_test ),
    qval =  p.adjust(pval, method="BH")
    
  )



# calculate CD19 stim vs WT

RES_CD19 = RawD_norm %>%
  na_if(0) %>%
  mutate(
    
    FC = map2_dbl(split(RawD_norm[42:53], 1:nrow(RawD_norm)), split(RawD_norm[30:41], 1:nrow(RawD_norm)), fc_cal ),
    pval = map2_dbl(split(RawD_norm[42:53], 1:nrow(RawD_norm)), split(RawD_norm[30:41], 1:nrow(RawD_norm)), t_test ),
    qval =  p.adjust(pval, method="BH")
    
  )


# calculate AKT vs WT stim
RES_AKT_WT_stim = RawD_norm %>%
  na_if(0) %>%
  mutate(
    
    FC = map2_dbl(split(RawD_norm[30:41], 1:nrow(RawD_norm)), split(RawD_norm[6:17], 1:nrow(RawD_norm)), fc_cal ),
    pval = map2_dbl(split(RawD_norm[30:41], 1:nrow(RawD_norm)), split(RawD_norm[6:17], 1:nrow(RawD_norm)), t_test ),
    qval =  p.adjust(pval, method="BH")
    
  )


# calculate AKT vs CD19 WT
RES_AKT_CD19_wT = RawD_norm %>%
  na_if(0) %>%
  mutate(
    
    FC = map2_dbl(split(RawD_norm[42:53], 1:nrow(RawD_norm)), split(RawD_norm[18:29], 1:nrow(RawD_norm)), fc_cal ),
    pval = map2_dbl(split(RawD_norm[42:53], 1:nrow(RawD_norm)), split(RawD_norm[18:29], 1:nrow(RawD_norm)), t_test ),
    qval =  p.adjust(pval, method="BH")
    
  )





##### RESULTS for each comparison seaprate #####################################


RES = RES_AKT_WT_stim


# check adjustment

q1 = RES %>% filter(pval < 0.05 & qval > 0.05) %>% nrow()
q2 = RES %>% filter(pval > 0.05 & qval > 0.05) %>% nrow()
q3 = RES %>% filter(pval > 0.05 & qval < 0.05) %>% nrow()
q4 = RES %>% filter(pval < 0.05 & qval < 0.05) %>% nrow()

ggplot(RES) +
  geom_point( aes(x=pval, y=qval), size =0.1 ) +
  labs(title = "p-value adjustment",
       subtitle = "lines indicate p-value of 0.05") +
  geom_hline(yintercept = 0.05) +
  geom_vline(xintercept = 0.05) +
  geom_text(label=q1, aes(x=0,y=0.1) ) +
  geom_text(label=q2, aes(x=0.5,y=0.5) ) +
  geom_text(label=q4, aes(x=0,y=-0.05) ) +
  theme_minimal()


ggsave("03_adjustment_check_AKT_CD19_wt.png", width=6, height=6.5)



# export results

write_csv(RES, "result_table_AKT_WT_stim.csv")

RES = read_csv("result_table.csv") # recover results is needed



######## volcano plot ##########################################################

PvalCo <- 0.05       # 5%
FcCo <- 2            # 2-fold OR 0.5-fold

curr_log2 = log2( RES$FC )
curr_pval = -log10(RES$qval)

linetype_costum = "dashed"
color_costum = "darkgrey"


dot_col = function(pval, FC) {
  
  ret = "not_significant"
  
  if(pval<=0.05 & FC>=2) { ret = "up" }
  if(pval<=0.05 & FC<=0.5) { ret = "down" }
  
  return(ret)
  
}

RES["dot_col"] = map2_chr( RES$qval, RES$FC, dot_col )

no_up = RES$dot_col %>% str_detect("up") %>% sum(na.rm=TRUE)
no_down = RES$dot_col %>% str_detect("down") %>% sum(na.rm=TRUE)

dot_col_map = c(
  "up" = "red",
  "not_significant" = "#696969",
  "down" = "navy"
)


tooltip_costum = paste(
  paste("Uniprot = ", RES$prot.names, sep=""),
  paste("FC = ", round(RES$FC, digits=2), sep=""),
  paste("qval = ", round(RES$qval, digits=4), sep=""),
  sep="\n"
) 

p1 = ggplot(RES) +
  labs(title = "Phosphopeptides AKT vs WT",
       subtitle = "stimulated",
       caption = paste("\n The lines indicate a p-value of ", PvalCo, " and a ", FcCo, "-fold change", sep=""), 
       x = "log2(FC)", 
       y = "-log10(q-value)") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_point( aes(x=curr_log2, y=curr_pval, color=dot_col, text=tooltip_costum  ), size=0.5) +
  scale_color_manual(values=dot_col_map) +
  scale_x_continuous( breaks = seq(-10,10,1), limits = c(-10,10) ) +
  scale_y_continuous( breaks = seq(0,20,1) ) +
  geom_segment(  aes( x= -log2(FcCo), xend= -log2(FcCo), y= -log10(PvalCo), yend=  Inf ), linetype=linetype_costum ,colour =color_costum) +
  geom_segment(  aes( x= log2(FcCo), xend= log2(FcCo), y= -log10(PvalCo), yend=  Inf ), linetype=linetype_costum ,colour =color_costum) +
  geom_segment(  aes( x= -Inf, xend= -log2(FcCo), y= -log10(PvalCo), yend=  -log10(PvalCo) ), linetype=linetype_costum ,colour =color_costum) +
  geom_segment(  aes( x= log2(FcCo), xend= Inf, y= -log10(PvalCo), yend=  -log10(PvalCo) ), linetype=linetype_costum ,colour =color_costum) +
  geom_text( aes(x=-8, y=6, label=paste(no_down, "downregulated", sep=" ")), color="navy"  ) +
  geom_text( aes(x=8, y=6, label=paste(no_up, "upregulated", sep=" ")), color="red"  )

ggsave("04_volcano_AKT_WT_stim.png", width=5.5, height=5.5)

ggplotly(p1, tooltip = "tooltip_costum")




##### heatmap ##################################################################


cal_z_score <- function(x){
  (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
}


# check % of missing values

mv = RES[6:23] %>% apply(2, function(x) sum( is.na(x) )/ nrow(RES) )

mv_tib = tibble(
  sample = names(mv),
  mv = mv
)

ggplot(mv_tib) +
  geom_col( aes( x=sample, y=mv*100 ) ) +
  labs(title = "Missing Values %") +
  theme_minimal() +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,5) ) +
  theme( axis.text.x = element_text(angle=90, vjust=0.5) )

ggsave("01_missingValues.png", width=8, height=5)


# missing value imputation #

impute_with_knn= function(data, k=10){
  knn_out= impute.knn(as.matrix(data),k, colmax = 0.9, rowmax = 0.9)
  matrix(as.integer(knn_out$data),ncol=ncol(data))
}

RES_imputed = RES
RES_imputed[6:23] = impute_with_knn( RES_imputed[6:23] )
# RES_imputed[9:38] = impute.MinProb(RES_imputed[9:38] , q = 0.01, tune.sigma = 1)


# prepare data for heatmap #

hmapdat = RES_imputed %>% 
  filter(qval<=0.05 & abs( log2(FC) ) >= 1 ) %>%
  arrange(desc(FC)) %>%
  select(6:23) %>% 
  as.matrix() %>% 
  apply(c(1,2), log2) %>% 
  apply(2, cal_z_score)

rownames(hmapdat) = RES_imputed %>% 
  filter(qval<=0.05 & abs( log2(FC) ) >= 1 ) %>%
  arrange(desc(FC)) %>%
  select(prot.names) %>%
  as.matrix()

hmap = pheatmap(hmapdat,
                color=wes_palette("Zissou1", 10, type = "continuous"),
                breaks=seq(-5,5, 1),
                cluster_cols = F,
                cluster_rows = T,
                show_rownames = F,
                border_color = NA,
                na_col="lightgrey") 




pdf("05_heatmap_significant_only.pdf", height = 0.15*nrow(hmapdat) , width = 7)
hmap
dev.off()





#######  PCA ####################################################################


# create data matrix for pca

curr_data = c(6:23)

pca_data <- as.matrix( RES[curr_data] )
rownames(pca_data) <- as.matrix( RES$prot.names )
colnames(pca_data) <- as.matrix( names(RES[curr_data]) )
t_pca_data <- t(pca_data)

t_pca_data = t_pca_data[, colSums(is.na(t_pca_data)) != nrow(t_pca_data)]

# percentage missing values
sum(is.na(t_pca_data)) / sum(!is.na(t_pca_data))


# Run this code if data should be scaled:
t_pca_data = apply(t_pca_data,c(1,2), log10)


# perform PCA with svdImpute method and centering

pca_result <- pca(t_pca_data, method="svdImpute", center=T, nPcs=5)

pca_result_tibble = as.data.frame(scores(pca_result) ) %>% as_tibble()

pca_result_tibble["sample"] = names(RES[curr_data])
pca_result_tibble["condition"] = c( rep("control",9), rep("DMAT",9) ) %>% as.factor()
pca_result_tibble["replicate"] = rep(1:3,2, each=3) %>% as.factor()
pca_result_tibble["run"] = rep(1:3,6)  %>% as.factor()




png("06_PCA_r2.png", height = 7 , width = 7, units="in", res=300 )
plotPcs(pca_result, pcs = 1:3 )
dev.off()

png("06_PCA_loadings.png", width=12, height=7, units="in", res=300 )
slplot(pca_result)
dev.off()

ggplot(pca_result_tibble) +
  labs(title = "PCA",
       subtitle = "DMAT vs. control \nMissing value imputation: SVD Impute | Centering: yes | Data Transformation: log10 | Missing values: 35% \n\n",
       x = "PC1 (89%)", 
       y = "PC2 (9%)") +
  scale_color_manual( values=c("#335584", "darkgreen") ) +
  theme_minimal() +
  geom_point( aes(x=PC1, y=PC2, color=condition, shape=replicate ), size = 3 )

ggsave("06_PCA_plot.png", width=8, height=6, dpi=300)


scatter3d(x = pca_result_tibble$PC1, 
          y = pca_result_tibble$PC2, 
          z = pca_result_tibble$PC3, 
          groups = pca_result_tibble$condition,
          grid = FALSE, 
          surface = FALSE,
          ellipsoid = TRUE,
          surface.col = c("#2a548a", "#8a602a") )


# write PC loadings into a data frame

loading_tib = tibble(
  entry = RES$prot.names
)

loadings_raw = pca_result@loadings
loadings_raw = cbind(rownames(loadings_raw), data.frame(loadings_raw[,c(1:3)], row.names=NULL) )
colnames(loadings_raw) = c("entry","PC1","PC2","PC3")

loading_tib = full_join(loading_tib, loadings_raw, by="entry")

write_csv(loading_tib, "loadings_table.csv")

