##Part 6
#replicate figure 1D 
Ad_1 <- read.csv("/project/bf528/project_2/data/samples/Ad_1/genes.fpkm_tracking",sep="\t",header=T)
Ad_2 <- read.csv("/project/bf528/project_2/data/samples/Ad_2/genes.fpkm_tracking",sep="\t",header=T)
P0_1 <- read.csv("/projectnb2/bf528/users/dreadlocks/project-2-dreadlocks/outputs/genes.fpkm_tracking",sep="\t",header=T)
P0_2 <- read.csv("/project/bf528/project_2/data/samples/P0_2/genes.fpkm_tracking",sep="\t",header=T)
P4_1 <- read.csv("/project/bf528/project_2/data/samples/P4_1/genes.fpkm_tracking",sep="\t",header=T)
P4_2 <- read.csv("/project/bf528/project_2/data/samples/P4_2/genes.fpkm_tracking",sep="\t",header=T)
P7_1 <- read.csv("/project/bf528/project_2/data/samples/P7_1/genes.fpkm_tracking",sep="\t",header=T)
P7_2 <- read.csv("/project/bf528/project_2/data/samples/P7_2/genes.fpkm_tracking",sep="\t",header=T)

#subset data
Ad_1 <- Ad_1[,c(1,5,10)]
Ad_2 <- Ad_2[,c(1,5,10)]
P0_1 <- P0_1[,c(1,5,10)]
P0_2 <- P0_2[,c(1,5,10)]
P4_1 <- P4_1[,c(1,5,10)]
P4_2 <- P4_2[,c(1,5,10)]
P7_1 <- P7_1[,c(1,5,10)]
P7_2 <- P7_2[,c(1,5,10)]

library(dplyr)
P0_1_1 <- P0_1 %>% group_by(gene_short_name) %>% summarise(P0_1=mean(FPKM))
P0_2_1 <- P0_2 %>% group_by(gene_short_name) %>% summarise(P0_2=mean(FPKM))
Ad_1_1 <- Ad_1 %>% group_by(gene_short_name) %>% summarise(Ad_1=mean(FPKM))
Ad_2_1 <- Ad_2 %>% group_by(gene_short_name) %>% summarise(Ad_2=mean(FPKM))
P4_1_1 <- P4_1 %>% group_by(gene_short_name) %>% summarise(P4_1=mean(FPKM))
P4_2_1 <- P4_2 %>% group_by(gene_short_name) %>% summarise(P4_2=mean(FPKM))
P7_1_1 <- P7_1 %>% group_by(gene_short_name) %>% summarise(P7_1=mean(FPKM))
P7_2_1 <- P7_2 %>% group_by(gene_short_name) %>% summarise(P7_2=mean(FPKM))

#merge all subset
merge_df1<-merge(P0_1_1,P0_2_1,by="gene_short_name")
merge_df2<-merge(merge_df1,P4_1_1,by="gene_short_name")
merge_df3<-merge(merge_df2,P4_2_1,by="gene_short_name")
merge_df4<-merge(merge_df3,P7_1_1,by="gene_short_name")
merge_df5<-merge(merge_df4,P7_2_1,by="gene_short_name")
merge_df6<-merge(merge_df5,Ad_1_1,by="gene_short_name")
final_df<-merge(merge_df6,Ad_2_1,by="gene_short_name")

final_df$P0<-(final_df$P0_1+final_df$P0_2)/2
final_df$P4<-(final_df$P4_1+final_df$P4_2)/2
final_df$P7<-(final_df$P7_1+final_df$P7_2)/2
final_df$Ad<-(final_df$Ad_1+final_df$Ad_2)/2

fpkm_df<-final_df[,c(1,10:13)]
sarcomere <- c('Pdlim5', 'Pygm', 'Myoz2', 'Des', 'Csrp3', 'Tcap', 'Cryab')
mitochondria <- c("Mpc1","Prdx3","Acat1","Echs1","Slc25a11","Phyh")
cell_cycle <- c("Cdc7","E2f8","Cdk7","Cdc26","Cdc6","Cdc27","E2f1","Cdc45","Rad51","Aurkb","Cdc23")

sarcomere_fpkm<-fpkm_df[fpkm_df$gene_short_name%in%sarcomere,]
mitochondria_fpkm<-fpkm_df[fpkm_df$gene_short_name%in%mitochondria,]
cell_cycle_fpkm<-fpkm_df[fpkm_df$gene_short_name%in%cell_cycle,]

#plot
library(ggplot2)
ggplot(sarcomere_fpkm) 
library(reshape)
mdata_sar <- melt(sarcomere_fpkm, id=c("gene_short_name"))
ggplot(mdata_sar, aes(x=factor(variable), y=value,fill=gene_short_name)) +
  geom_line(aes(color = gene_short_name, group = gene_short_name))+
  geom_point()+
  xlab("Stage") + ylab("FPKM") + # Set axis labels
  ggtitle("Sarcomere") +     # Set title
  theme_bw()

mdata_mito <- melt(mitochondria_fpkm, id=c("gene_short_name"))
ggplot(mdata_mito, aes(x=factor(variable), y=value,fill=gene_short_name)) +
  geom_line(aes(color = gene_short_name, group = gene_short_name))+
  geom_point()+
  xlab("Stage") + ylab("FPKM") + # Set axis labels
  ggtitle("Mitochondria") +     # Set title
  theme_bw()

mdata_cell <- melt(cell_cycle_fpkm, id=c("gene_short_name"))
ggplot(mdata_cell, aes(x=factor(variable), y=value,fill=gene_short_name)) +
  geom_line(aes(color = gene_short_name, group = gene_short_name))+
  geom_point()+
  xlab("Stage") + ylab("FPKM") + # Set axis labels
  ggtitle("Cell cycle") +     # Set title
  theme_bw()

#hierarchical clustering 
heat_data<-final_df[,c(1:9)]
gene_exp_diff <- read.table("/projectnb2/bf528/users/dreadlocks/project-2-dreadlocks/outputs/gene_exp.diff", 
                            header = TRUE)
gene_sig<-gene_exp_diff[gene_exp_diff$significant=="yes",]
gene_top_1000<-gene_sig[order(gene_sig$p_value),][1:1000,]

heat_data_final<-heat_data[heat_data$gene_short_name%in%gene_top_1000$gene,]

rownames(heat_data_final)<-heat_data_final[,1]
heat_data_final<-heat_data_final[,-1]

library(gplots)
png(file="Heatmap_ind.png")
heatmap(as.matrix(heat_data_final),
        main="Hierarchical clustering for differentially expressed genes",
        scale = "row",
        labRow=T, Colv = NA,
        xlab="Sample", ylab="Gene")
dev.off()


