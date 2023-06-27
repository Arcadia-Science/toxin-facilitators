## This script was designed to carry out the analysis of the accessory sequences clusters
# it uses the mmseqs cluster output file reformated as a csv file to determine and vizualise the size of each accessory sequences cluster
# it uses the Pfam annotations obtained for each of the 232 representative sequences that got an Pfam annotation and highlights which clusters are associated with a toxin-assocaited Pfam.
# it conducted per cluster species and original toxin diversity analysis

# Environment setup
library(ggplot2)
library(plotly)
library(tidyverse)
library(dplyr)
library('extrafont')
loadfonts(device='win',quiet=T)
library(htmlwidgets)

#Required datasets
Data_C=read.csv('results/Results_R/Analysis/Clusters_Accessory_UpdatedJune2023.csv') # cluster list output from the mmseq cluster module (identifies the representative sequence of a cluster (column 1) and all the other sequences contained in the cluster (column 2))
Data_C_Pfam=read.csv('metadata/R_scripts_data/Cluster_Acc_Ref_PFAM_UpdatedJune2023.csv') # list of the representative sequences that obtained a Pfam annotation (column 6) and whether or not this annoation is a toxin-associated Pfam (column 8)
Data_seq_o=read.csv('metadata/R_scripts_data/Origin_data_accessory_seq_clusters_UpdatedJune2023.csv') # table that summarize for each accessory sequence cluster, the associated accessory sequences
# the outlier sequence associated with the accessory sequence, its rep number (one outlier can have multiple rep)
# the reference toxin of the toxin cluster the outlier sequence sequence was identified from


### Step 1: Summarizing the number of accessory sequences in each cluster ###

## Creates a dataframe that provides for each cluster (identified by its representative sequence)

Clust_list=c(unique(Data_C$Clust_Ref)) #creates a vector identifying each cluster by its representative sequence ID

n_clust=length(Clust_list)
mat=matrix(nrow=n_clust,ncol=2)
colnames(mat)=c('Clust_ref','nb_hits')

for (i in 1:n_clust){
  Acc=Clust_list[i]
  Dat=subset(Data_C, Data_C$Clust_Ref==Acc)
  mat[i,1]=Acc
  mat[i,2]=dim(Dat)[1]
  
}

Clust_size=as.data.frame(mat)
Clust_size$nb_hits=as.numeric(Clust_size$nb_hits)


### Step 2: Cluster Pfam annotation analysis and visualization - Figure 4 ###
## Add information concerning Pfam annotation to the previous cluster size table

Data_C_Pfam_filtered=subset(Data_C_Pfam, Data_C_Pfam$e_value<=1.0e-5)
Data_all=left_join(Clust_size,Data_C_Pfam_filtered,by='Clust_ref')
Data_all$Toxin.PFAM=Data_all$Toxin.PFAM%>%
  replace_na('NA_PFAM') #add a specific label mentioning that no Pfam annotation was found for any representative without Pfam annotation

## Visualization of cluster with and without Pfam annotation based on cluster length distribution.
hpt=ggplot(Data_all, aes(x=nb_hits, fill=Toxin.PFAM))+
  geom_histogram(binwidth = 1,col='black',size=0.1)+ theme_light()+ theme(text=element_text(family="Arial"),legend.position='none')+
  xlab('Cluster size') + ylab('Number of clusters') + ggtitle('Distribution of cluster size of accessory sequences')+
  scale_fill_manual(values=c('#5088C5','#F28360','#3B9886'))




### Step 3: Filtering out toxin-associated clusters  ###
## Add information concerning Pfam annotation to the previous cluster size table

Data_Clust_NoT=subset(Data_all, Data_all$Toxin.PFAM!='YES' )

write.csv(Data_Clust_NoT,'results/Results_R_Analysis/Accessory_clusters_nonToxin_UpdatedJune2023.csv')


### Step 4: Identifying and counting outlier sequences from accessory sequence in each cluster ###
## The same outlier sequences may have generated multiple accessory sequences each identified by a _rep number. 
## The clusters contain accessory sequences, some may come from the same outlier and 1 cluster may actually contain a single outlier while having multiple accessory sequences
## For the rest of the analysis, we want to work with the outlier sequences per cluster

Clust_291=c(Data_Clust_NoT$Clust_ref) # vector of all representative sequences 
Data_seq=subset(Data_seq_o, Data_seq_o$Clust_Ref%in%Clust_291) 

Data_291_short=Data_Clust_NoT[,1:2]

Data_per_outlier=data.frame() # this table will show all outlier sequences in each accessory sequence cluster (identified by its representative sequence) and how many rep of that outlier were found

for (i in 1:length(Clust_291) ){
  Clust_id=Clust_291[i]
  Dat=subset(Data_seq, Data_seq$Clust_Ref==Clust_id)
  Seqs_or=c(unique(Dat$Seq_o))
  n=length(Seqs_or)
  
  Dat_clust=subset(Data_291_short, Data_291_short$Clust_ref==Clust_id)
  Clust_s=Dat_clust[1,2]
  
  if (n==1){
    Seq=Seqs_or
    rep=dim(Dat)[1]
    spe=unique(Dat$Species)
    tox=unique(Dat$Clust_o)
    
    Dat_tmp=data.frame('Clust_ref'=Clust_id,'Clust_size'=Clust_s,
                       'Outlier_seq'=Seq,'Nb_acc'=rep,'Species'=spe,'Origin_C'=tox)
    
  }
  
  if (n>1){
    Dat_tmp=data.frame()
    for (j in 1:n){
      Seq=Seqs_or[j]
      Dat2=subset(Dat, Dat$Seq_o==Seq)
      rep=dim(Dat2)[1]
      spe=unique(Dat2$Species)
      tox=unique(Dat2$Clust_o)
      
      Dat_tmptmp=data.frame('Clust_ref'=Clust_id,'Clust_size'=Clust_s,
                            'Outlier_seq'=Seq,'Nb_acc'=rep,'Species'=spe,'Origin_C'=tox)
      
      Dat_tmp=rbind(Dat_tmp, Dat_tmptmp)
      
    }
  }
  
  Data_per_outlier=rbind(Data_per_outlier, Dat_tmp)
}


### Step 5: Summarizing outlier, species and toxin information per cluster ###
## For each accessory cluster we summarize the number of present unique outliers, the number of different species, and the number of different toxins the outliers are associated with

Data_summary_acc_cluster=data.frame()

for (i in 1:length(Clust_291)){
  Clustid=Clust_291[i]
  Dat=subset(Data_per_outlier, Data_per_outlier$Clust_ref==Clustid)
  
  Dat_clust=subset(Data_291_short, Data_291_short$Clust_ref==Clustid)
  Clust_s=Dat_clust[1,2]
  
  n_outs=length(c(unique(Dat$Outlier_seq)))
  n_spe=length(c(unique(Dat$Species))) 
  n_tox=length(c(unique(Dat$Origin_C)))
  
  Data_tmp=data.frame('Clust_ref'=Clustid,"Clust_size"=Clust_s,
                      "nb_outlier"=n_outs,"nb_spe"=n_spe,'nb_OCT'=n_tox)
  
  Data_summary_acc_cluster=rbind(Data_summary_acc_cluster, Data_tmp)
}


### Step 5: Filtering out cluster with single outlier sequence ###

Data_summary_acc_cluster_no1=subset(Data_summary_acc_cluster, Data_summary_acc_cluster$nb_outlier>1)
write.csv(Data_summary_acc_cluster_no1, 
          "results/Results_R_Analysis/Data_summary_acc_cluster_no1_UpdatedJune2023.csv")

### Step 6: Cluster data summary visualization ###
# Histograms 


# H2: Species per cluster
hp_spe=ggplot(Data_summary_acc_cluster_no1, aes(x=nb_spe))+
  geom_histogram(binwidth = 1,fill='#8A99AD',col='black')+ theme_light()+ theme(text=element_text(family="Arial"))+
  xlab('# of different species/cluster') + ylab('Number of clusters') + 
  scale_x_continuous(limits = c(0, 11), breaks = seq(0,10,1))

pdf('figures/Fig5A.pdf')
hp_spe
dev.off()

# H3: Toxins per cluster
hp_tox=ggplot(Data_summary_acc_cluster_no1, aes(x=nb_OCT))+
  geom_histogram(binwidth = 1,fill='#8A99AD',col='black')+ theme_light()+ theme(text=element_text(family="Arial"))+
  xlab('# of different original toxin clusters/cluster') + ylab('Number of clusters') +
  scale_x_continuous(limits = c(0, 3), breaks = seq(0,3,1)) 

pdf('figures/Fig5B.pdf')
hp_tox
dev.off()

### Step 7: Data export - Creation of Table 2 ##

Data_FHM=Data_summary_acc_cluster_no1

for (i in 1:dim(Data_FHM[1])){
  Data_FHM$id[i]=paste('Cluster',i) 
}

write.csv(Data_FHM, "results/Results_R_Analysis/ASClusters_metrics_summary_UpdatedJune2023.csv")
