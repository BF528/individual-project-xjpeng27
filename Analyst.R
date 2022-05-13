df<-read.table(file.choose(),header = TRUE)
#Part5
##1.
df<-df[order(df$q_value),]
top10<-df[1:10,c("gene","value_1","value_2","log2.fold_change.","p_value","q_value")]
write.csv(top10,"/Desktop/BF528/individual project/top10.csv")

##2.
hist(df$log2.fold_change.,
     breaks = 60,
     xlab = "log fold change",
     main = "histogram of log fold change for all genes")

##3.
df.sub<-subset(df,df$significant=="yes")

##4.
hist(df.sub$log2.fold_change.,
     breaks = 60,
     xlab = "log fold change",
     main = "histogram of log fold change for significant genes")

##5.
upgenes<-df.sub[which(df.sub$log2.fold_change.>0),c("gene")]
length(upgenes)#2830
downgenes<-df.sub[which(df.sub$log2.fold_change.<0),c("gene")]
length(downgenes)#2597

##6.
write(upgenes,"Desktop/BF528/individual project/upgenes.txt")
write(downgenes,"Desktop/BF528/individual project/downgenes.txt")
