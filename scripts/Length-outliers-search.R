## This R script carries out the length-outlier search within toxin-clusters 

#Environment setup
library('tidyverse')
library('ggpubr')
library('plotly')
library('extrafont')
loadfonts(device='win',quiet=T)

# Required datasets
Data_clust=read_tsv('results/2023-06-15-updated-uniprot-vs-venom-tick-proteins-summaries.tsv') ##  (this corresponds to the summary generated after the clustering)
Ref_size=read.csv('metadata/R_scripts_data/Uniprot_Cluster_Info.csv') ## (this corresponds to the metadata of Uniprot reference toxin assocaited with each toxin cluster)


#### Step 1: Data reduction and cluster filtering ####

#Data reduction: keeping only the highest match for proteins with multiple hits
  ##Some proteins had multiple hits and would be redundant in multiple clusters. Thus we only want to keep the best hit for these proteins
  ## We use a for loop and nested if condition that goes through all venom proteins that got a hit to look whether multiple hits are found to only keep the best hit (lowest e-value). 
Prot_ID=c(unique(Data_clust$protein_id))
Data_B=data.frame()

for (i in 1:dim(Data_clust)[1]){
  Prot=Prot_ID[i]
  Table=subset(Data_clust, Data_clust$protein_id==Prot)
  nb=dim(Table)[1]
  if (nb==1){
    Data_B=rbind(Data_B,Table)
  }
  if (nb>1){
    Mn=min(Table$evalue)
    TableG=subset(Table, Table$evalue==Mn)
    Data_B=rbind(Data_B,TableG)
  }
  
}

Data_BH=Data_B

  
#Plotting the distribution of cluster size 
  ## We create a dataframe Count_Hits, that provides for each cluster the accession number of its reference toxin and the number of venom proteins
Acc_hits=c(unique(Data_BH$accession_hit))
n_hits=length(Acc_hits)

mat=matrix(nrow=n_hits,ncol=2)
colnames(mat)=c('Acc_hits','nb_hits')

for (i in 1:n_hits){
  Acc=Acc_hits[i]
  Dat=subset(Data_BH, Data_BH$accession_hit==Acc)
  mat[i,1]=Acc
  mat[i,2]=dim(Dat)[1]
}

Counts_Hits=as.data.frame(mat)
Counts_Hits$nb_hits=as.numeric(Counts_Hits$nb_hits)

 ## Plot for visualization (Figure not included in Pub)
hp=ggplot(Counts_Hits, aes(x=nb_hits))+
  geom_histogram(binwidth = 1,fill='#8A99AD')+ theme_light()+
  xlab('Cluster size') + ylab('Number of clusters') + ggtitle('Distribution of cluster size')

ggplotly(hp)


#Filtering out clusters with les than 5 sequences
Clust5=subset(Counts_Hits, Counts_Hits$nb_hits>=5)
Acc_5=c(Clust5$Acc_hits)
Data_C5=subset(Data_BH, Data_BH$accession_hit%in%Acc_5)


#### Step 2: Protein length ratio calculation ####

# Calculating protein size ratio between clustered proteins and cluster reference toxin within each cluster independently

Data_C5$Ref_L=0 # adds an empty column were associated reference toxin length will be stored
Data_C5$Ratio=0 # adds an empty column were protein length ratio will be stored

Data_build=data.frame()

for (i in 1:dim(Data_C5)[1]){
  Acc=Acc_5[i]
  DatC=subset(Data_C5, Data_C5$accession_hit==Acc)
  Ref=subset(Ref_size, Ref_size$RefSeq==Acc)
  Ref_L=as.numeric(Ref$RefSeq_L)
  
  DatC[,'Ref_L']=Ref_L
  DatC[,'Ratio']=DatC$length/DatC$Ref_L
  
  Data_build=rbind(Data_build,DatC)
}

Data_HBR=Data_build # dataframe that contains for any venom protein in a cluster > 5 sequences different metadata, including the cluster it belongs too (identified as it reference toxin accession - accession_hit) and the size ration (Ratio) 



### Step 3: Identification of outliers ###

# Data organization/reformatting for easier visualization and outliers search
  ## We want to organize sequences by cluster and then organize the clusters based on the reference toxin length. 
  ## For further visualization, each cluster is also associated with a Plot number as we will visualize clusters 10 a time

  #1- Assigning each cluster to a given plot
Ref_length_5=subset(Ref_size, Ref_size$RefSeq%in%Clust5$Acc_hits)

Ref_L5=Ref_length_5[order(Ref_length_5$RefSeq_L),] # sort reference toxin based on their length
chunk <- 10
n <- nrow(Ref_L5)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
Split_T <- split(Ref_L5,r) # creates a list of elements that contains each 10 reference toxins sorted based on length

Data_B=data.frame() # empty dataframe that will be filled in the for loop to assign the right plot number for each cluster.
for (i in 1:length(Split_T)){
  Tab=Split_T[[i]]
  Tab$Plot=paste('Plot',i)
  Data_B=rbind(Data_B,Tab)
}

  #2- Reformating/Simplifying data 
Cluster_plotnum=Data_B
Cluster_plotnum=Cluster_plotnum[,c(-6,-4,-3,-1)]
colnames(Cluster_plotnum)=c('accession_hit', 'Ref_length','Plot')

  #3- Data preparation for outlier search and plotting
Data_plotnum=left_join(Data_HBR,Cluster_plotnum, by='accession_hit')
colnames(Clust5)=c('accession_hit','nb_hits')
Data_plotnum=left_join(Data_plotnum,Clust5,by='accession_hit')


# Identification of outliers in each cluster individually
## Outliers are labelled in an extra column as 'Yes' or 'No' based on the following outlier criteria:
  #- ratio > Q3+1.5*IQ ( Q3: 3rd quartile of the cluster ratio distribution; IQ=interquartile range)
  #- ratio > 1.5

data_out_cat=data.frame()
Data_plotnum=Data_plotnum[order(Data_plotnum$Ref_length),] # sort the table according to toxin reference length
CID=c(unique(Data_plotnum$accession_hit)) # get the vector with unqiue cluster ID

for (i in 1:length(CID)){
  Clust=CID[i]
  TempD=subset(Data_plotnum, Data_plotnum$accession_hit==Clust)
  quant=quantile(TempD$Ratio, probs=c(0.25,0.75))
  Q1=quant[1]
  Q3=quant[2]
  IQ=(Q3-Q1)
  UB_out=Q3+1.5*IQ
  
  TempD$Outlier=ifelse(TempD$Ratio>=UB_out & TempD$Ratio>=1.5 ,"Yes","No")
  data_out_cat=rbind(data_out_cat,TempD)
}

# Plot: distribution of length ratio per cluster with highlighted length-outlier proteins(10 cluster at a time) 
  ## The for loop generates and print 1 plot every 10 clusters
Plot_num=c(unique(Data_plotnum$Plot))

pdf('figures/Outlier_visualization.pdf')

for(i in 1:length(Plot_num)){
  PlN=Plot_num[i]
  DP=subset(data_out_cat, data_out_cat$Plot==PlN)
  DP=DP[order(DP$Ref_L),]
  DP$accession_hit=as.factor(DP$accession_hit)
  
  p=ggplot(DP, aes(x=fct_inorder(accession_hit),y=Ratio, col=Outlier))+
    geom_jitter(position=position_jitter(0.2), cex=0.5) + theme_light()+
    theme(axis.text.x=element_text(size=7,angle=90))+ theme(text=element_text(family="Arial"))+
    geom_hline(yintercept = 1, col='gray63') + xlab('ClusterID') + ggtitle (paste(PlN))+
    scale_color_manual(values=c('#8A99AD','#FFB984'))
  
  print(p)
  
  
}
dev.off()

### Step 4: Extracting outliers information ###

Outliers_seq=subset(data_out_cat, data_out_cat$Outlier=='Yes')
write.csv(Outliers_seq,'results/Results_R_Analysis/Venomproteins_ticks_toxins_outliers_June2023.csv')

### Step 5: Preparation of Figure 3 ###

# Figure 3A
i=16 # select 1 random plot among all the plots generated in the Outlier_visualization.pdf document to illustrate how outliers are identified.
PlN=Plot_num[i]
DP=subset(data_out_cat, data_out_cat$Plot==PlN)
DP=DP[order(DP$Ref_L),]
DP$accession_hit=as.factor(DP$accession_hit)

p=ggplot(DP, aes(x=fct_inorder(accession_hit),y=Ratio, col=Outlier, 
                 text=paste("protein_seq_id:",protein_id)))+
  geom_jitter(position=position_jitter(0.15), cex=0.5,show.legend = FALSE) + theme_light()+
  theme(axis.text.x=element_text(size=7,angle=90), text=element_text(family="Arial"))+
  geom_hline(yintercept = 1, col='gray63') + xlab('\nClusterID') +
  scale_color_manual(values=c('#8A99AD','#FFB984'))

Fig3A=ggplotly(p)

htmlwidgets::saveWidget(
  widget = Fig3A, 
  file = "figures/Fig3A.html",
  selfcontained = TRUE 
)


