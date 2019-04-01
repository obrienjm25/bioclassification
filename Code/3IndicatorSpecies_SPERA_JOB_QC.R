########################################
# Indicator Species Analysis
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

library(labdsv)
library(dplyr)

####Indicator species analysis on top 3 clusters####

#Import Site by Species matrix with cluster assignments
colorGridData<-read.csv("Data/Coloredgrid4km_cluster_assignment_QC.csv")
head(colorGridData)
names(colorGridData)
colorGridData<-colorGridData[,-1] #remove superfluous index variable
top3 <- names(sort(table(colorGridData$cl), decreasing = T))[1:3] #list of top 3 clusters

colorGridData_top3<-colorGridData[colorGridData$cl %in% top3,] #create subset df of top 3 clusters
PA_data<- dplyr::select(colorGridData_top3,-id,-cl,-assigned) #remove columns to create strictly sitexspecies PA matrix

colorGridData_top3$cl2[colorGridData_top3$cl==top3[1]]<-1 #assign top 3 clusters new numeric cluster membership
colorGridData_top3$cl2[colorGridData_top3$cl==top3[2]]<-2 #between 1 and 3 in order of largest to smallest
colorGridData_top3$cl2[colorGridData_top3$cl==top3[3]]<-3


#Provide a name to each cluster associated with broad geomorphic features
colorGridData_top3$ecoregion[colorGridData_top3$cl==top3[1]]<- "Deep Channels" 
colorGridData_top3$ecoregion[colorGridData_top3$cl==top3[2]]<- "Shallow Banks & Straits"
colorGridData_top3$ecoregion[colorGridData_top3$cl==top3[3]]<- "Channel Heads & Slopes"

colorGridData_top3$cl2<-as.factor(colorGridData_top3$cl2)
colorGridData_top3$cl2<-droplevels(colorGridData_top3$cl2)
summary(colorGridData_top3$cl2)
dim(colorGridData_top3)
dim(PA_data)

indicators2<-indval(PA_data, as.numeric(colorGridData_top3$cl2)) #indicator analysis as proposed by Dufrene & Legendre 1997

summary.indval(indicators2, type="long") #species x cluster table with indicator values
summary.indval(indicators2) #indicator value and probability for significant indicator species grouped by cluster

relfreq<-indicators2$relfrq #get relative frequency for that species in that cluster
relabund<-indicators2$relabu #get relative abundace for that species in that cluster
maxcls<-indicators2$maxcls #get cluster with max indval for that species
pval<-indicators2$pval #get pvalue
indval<-indicators2$indval #get indicator species value
indmax <- indicators2$indcls #get indVal for species in its maximum class
maxfreq <- pmax(relfreq[1], relfreq[2],relfreq[3]) # gives relative frequency for species in its max class
names(maxfreq) <- "maxfreq"

names1<-c(paste("freq_", c("1","2",
                           "3"), sep=""))
names2<-c(paste("abund_", c("1","2",
                            "3"), sep=""))
names3<-c(paste("indval_", c("1","2",
                             "3"), sep=""))
names(relfreq)<-names1
names(relabund)<-names2
names(indval)<-names3

#Combine relative frequency, abundance, indicator value, cluster with max indVal and pval into 1 df
combine<-as.data.frame(cbind(maxcls,pval,relfreq, relabund,indval, indmax, maxfreq))
combine$species<-rownames(combine) #Create species name column
head(combine)

#filter out observations with p values > 0.05 and IndVal < 0.25
combine <- combine %>% 
  dplyr::filter(.,pval<=0.05) %>% 
  dplyr::filter(.,indmax >= 0.25) 
combine<-combine[order(combine$maxcls, -combine$indmax),] #group by cluster
combine$species<-gsub("[.]"," ",combine$species)

#Get common names for indicator species

QC <- readRDS("Data/QC_invertsAdded.rds")
SpeciesIndex <- QC %>% dplyr::select(species, common_name) %>% distinct()

for (i in 1:length(combine$species)){
  combine$Common_name[i] <- SpeciesIndex[agrep(combine$species[i],SpeciesIndex$species),2]
}

IndicatorTable <- combine %>% dplyr::select(maxcls,species,Common_name,maxfreq,indmax)
IndicatorTable$indmax <- round(IndicatorTable$indmax,3)
IndicatorTable$maxfreq <- round((IndicatorTable$maxfreq*100),1)
colnames(IndicatorTable) <- c("Ecoregion","Species","Common name", "Relative Frequency (%)","IndVal")
write.table(IndicatorTable, file = "Output/IndicatorSpecies_QC.txt",
            row.names = F, sep = ",", qmethod = "escape")

#On to Random Forest Modeling