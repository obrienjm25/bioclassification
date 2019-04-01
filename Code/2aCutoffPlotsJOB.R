########################################
# Create dendrogram and obtain species groupings for benthic assemblages 
# Exploratory plots to choose cutoff size
# Katie Gale (katie.gale@dfo-mpo.gc.ca) _ Created 2015
# Formatted and uploaded 2017 
#
# As published in 
# Rubidge, E., Gale, K.S.P., Curtis, J.M.R., McClelland, E., Feyrer, L., Bodtker, K., and Robb, C.
# 2016. Methodology of the Pacific Marine Ecological Classification System and its Application
# to the Northern and Southern Shelf Bioregions. DFO Can. Sci. Advis. Sec. Res. Doc.
# 2016/035. xi + 124 p.
#
## Adapted code for Biological classification analysis - Eastern Canadian groundfish - Stanley, Heaslip, Jeffery, O'Brien
#########################################

pkgs <- list('vegan', 'simba', 'maptools', 'pvclust', 'dendroextras', 'dendextend', 'reshape', 'reshape2', 'dplyr', 'tidyr', 'NbClust','ggplot2')
invisible(lapply(pkgs, library, character.only = T))

benthtree<-readRDS("Data/benthtree4km.rds")
SiteXSpecies<-read.csv("Data/ClusterData4km.csv", stringsAsFactors = F, row.names = 1)
grid.pa.simp<-sim(SiteXSpecies,  method='simpson')

#Examine internal cluster validity index as a function of cluster size
#calculate various CVI's for partititions with k=2-20 

cvi <- c('ch', 'cindex', 'ptbiserial','db','silhouette')
nclust <- c(2:20)
vals <- list()

for (i in cvi) {
  
  vals[[i]] <- NbClust(data = SiteXSpecies, diss = grid.pa.simp, distance = NULL, 
                       method = 'average',index = i, min.nc = 2, max.nc = 20)$All.index
  
}

cvi.df <- data.frame(NC = rep(nclust, length(cvi)), Value = unlist(vals))
cvi.df$Index <- gsub('[.][0-9]+','',row.names(cvi.df))
p.cvi <- ggplot(cvi.df) +
  geom_point(aes(x = NC, y = Value)) +
  labs(x = 'Number of clusters', y = 'Index value') +
  facet_wrap(~ Index, scales = 'free_y')+
  theme_bw(); p.cvi

#2, 4, & 9 clusters all possible

ggsave("Output/CVI_Mar.tiff",p.cvi)

#Pinpoint cut-off value by examining number of sites assigned to top clusters and evenness among clusters
#as a function of dissimilarity

####################################################
### Figures to help pick Cutoff
###  create graphs and, if wanted, choose automatic cutoff for dendrogram
#STEPS
# cut tree at h
# make table
# count top 10 clusters
# table of height vs n cells in top 10 clusters

ncells1<-list()
ncells2<-list()
ncells3<-list()
ncells4<-list()
ncells5<-list()
ncells6<-list()
ncells7<-list()
ncells8<-list()
ncells9<-list()
ncells10<-list()
nClusters <- list()

s<-seq(30,100, by=0.1)
for (i in 1:length(s)){
  clus<-dendroextras::slice(benthtree, h=(s[i]/100))
  cluscount<-as.data.frame(table(clus))
  cluscount2<-cluscount[order(-cluscount$Freq),]
  nClusters[i] <- (length(cluscount2$Freq))
  ncells1[i]<-sum(cluscount2$Freq[1:1])
  ncells2[i]<-sum(cluscount2$Freq[2:2])
  ncells3[i]<-sum(cluscount2$Freq[3:3])
  ncells4[i]<-sum(cluscount2$Freq[4:4])
  ncells5[i]<-sum(cluscount2$Freq[5:5])
  ncells6[i]<-sum(cluscount2$Freq[6:6])
  ncells7[i]<-sum(cluscount2$Freq[7:7])
  ncells8[i]<-sum(cluscount2$Freq[8:8])
  ncells9[i]<-sum(cluscount2$Freq[9:9])
  ncells10[i]<-sum(cluscount2$Freq[10:10])
}
table<-cbind.data.frame(s/100, unlist(ncells1),unlist(ncells2),unlist(ncells3),unlist(ncells4),unlist(ncells5),unlist(ncells6),unlist(ncells7),unlist(ncells8),unlist(ncells9),unlist(ncells10),unlist(nClusters))
colnames(table)<-c("distance_Bsim", "ncells_top_1","ncells_top_2","ncells_top_3","ncells_top_4","ncells_top_5","ncells_top_6","ncells_top_7","ncells_top_8","ncells_top_9","ncells_top_10","nClusters")

head(table) #number of sites in top cluster, 2nd top cluster, 3rd top cluster etc
table$distance_Bsim[table$nClusters>=2 & table$nClusters<=9] #Look for cut-off between 0.587-0.945

tablea<-table[,c(-1,-12)]
tablea$SD <- apply(tablea,1, sd, na.rm = TRUE)
tablea$rowsums<-rowSums(tablea[,c(1:10)], na.rm=T)
tablea$nunclassified<-nrow(SiteXSpecies)-tablea$rowsums
tablea$ratio<-tablea$rowsums/tablea$SD
tablea$Bsim<-table$distance_Bsim
tablea$Bsim100<-tablea$Bsim*100

#Exploratory plots

cutoff1<-0.587
cutoff2<-0.691
cutoff3<-0.724

#Plot that shows ratio of assigned sites to variation (SD) in cluster size
plot(tablea$Bsim,tablea$ratio, type="l", xlab="Bsim", ylab="Ratio of sites in top 10 clusters to cluster size SD")
minRatio<-subset(tablea, ratio==min(tablea$ratio, na.rm=T))
points(minRatio$Bsim, minRatio$ratio,col="red")
abline(v = c(cutoff1, cutoff2, cutoff3))

tableb<-tablea[complete.cases(tablea),]
plot(tableb$Bsim,tableb$SD, xlab = "Bsim", yab = "SD of top 10 cluster size") #Plot StDev of cluster size by Bsim cutoff
plot(tableb$Bsim,tableb$rowsums, type="l", col="red", xlab = "Bsim", ylab = "Number of sites in top 10 clusters") #Plot Number of sites included in top 10 by Bsim cutoff


#######
#Figures for paper
######
table2<-reshape2::melt(table[,-12], id="distance_Bsim")
table2
table2$value[is.na(table2$value)] <- 0

library(grid)
library(ggplot2)

####
#Number of sites in each of the top 10 clusters for dendrograms cut at increasing Bsim values
#Area plot
#####
tiff("Output/CutOff_Bsim_byNSites_Maritimes.tiff",
  res=1200, width=180, height=84, units="mm", pointsize=8)
ggplot(table2) +
  geom_area(aes(y = value,x = distance_Bsim,fill=variable), color = 'black', position = position_stack(vjust = 0.5, reverse = T), stat="identity", na.rm = T) +  theme_classic()+ 
  scale_fill_grey(start=0.95, end=0.3,labels=c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10"))+
  scale_x_continuous(expression(paste(beta["sim"]," cut-off")),expand = c(0, 0))+
  scale_y_continuous("Number of sites in each cluster", limits=c(0,2600), expand=c(0,0))+labs(fill="") +
  geom_segment(aes(x = cutoff1-0.005, y = max(table2$value, na.rm=T), xend = cutoff1-0.005, yend = 0), linetype=2)+
  geom_segment(aes(x = cutoff1+0.005, y = max(table2$value, na.rm=T), xend = cutoff1+0.005, yend = 0), linetype=2)+
  geom_segment(aes(x = cutoff1-0.005, y = max(table2$value, na.rm=T), xend = cutoff1+0.005, yend = max(table2$value, na.rm=T)), linetype=2)
dev.off()

####
#Variation (SD) of the number of sites in the top 10 clusters vs. number of sites 
# retained in those clusters for increasing Bsim cut-off values.
#####
#
tiff("Output/CutOff_NSitesBySD_Maritimes.tiff",
    res=1200, width=84, height=84, units="mm", pointsize=8)
plot(tablea$rowsums,tablea$SD, type="n",xlab="Number of sites in top 10 clusters", ylab="SD of number of sites in top 10 clusters")
points(tablea$rowsums[tablea$Bsim100>=30&tablea$Bsim100<=39],tablea$SD[tablea$Bsim100>=30&tablea$Bsim100<=39], type="l",  lty=4, lwd=2)
points(tablea$rowsums[tablea$Bsim100>38&tablea$Bsim100<=49],tablea$SD[tablea$Bsim100>38&tablea$Bsim100<=49], type="l",   lty=1,lwd=1)
points(tablea$rowsums[tablea$Bsim100>48&tablea$Bsim100<=59],tablea$SD[tablea$Bsim100>48&tablea$Bsim100<=59], type="l",  lty=2, lwd=1.5)
points(tablea$rowsums[tablea$Bsim100>58&tablea$Bsim100<=69],tablea$SD[tablea$Bsim100>58&tablea$Bsim100<=69], type="l",  lty=1, lwd=2)
points(tablea$rowsums[tablea$Bsim100>68&tablea$Bsim100<=78],tablea$SD[tablea$Bsim100>68&tablea$Bsim100<=78], type="l",   lty=3, lwd=2)
legend("topleft",legend=c((expression(paste(beta["sim"]," cut-off"))),"0.30 to 0.39","0.40 to 0.49","0.50 to 0.59", "0.60 to 0.69","0.70 to 0.78"),lty=c(0,4,1,5,1,3),lwd=c(0,2,1,1.5,2,2), bty="n")
points(2527,332.3078,col="black", pch=22, cex=3)
# points(tablea$rowsums[tablea$Bsim==cutoff1],tablea$SD[tablea$Bsim==cutoff1],col="black", pch=21, cex=3)
text(2527-160,332.3078,col="black", labels=c("0.58"))
# text(tablea$rowsums[tablea$Bsim==cutoff1]-130,tablea$SD[tablea$Bsim==cutoff1]+15,col="black", labels=c("0.37"))
dev.off()


