setwd("E:/Biological Classification")

#load packages
pkgs <- list("rgdal", "sp", "ggplot2", "grid", "gridExtra", "gtable")
invisible(lapply(pkgs, library, character.only = T))

#load grid cell file and wide format survey data
load("Data/MaritimeGridsWholeCell.RData")
Mardata <- readOGR(dsn = "Data", layer = "WideMaritimes_wgt_ptshp")

#create empty vectors to fill in loop below
CellSize <- vector(mode = 'numeric', 0)
NumCell <- vector(mode = 'numeric', 0)
Populated <- vector(mode = 'numeric', 0)
Mean <- vector(mode = "numeric",0)
SD <- vector(mode = "numeric", 0)
Max <- vector(mode = 'numeric',0)
Median <- vector(mode = 'numeric',0)
q25 <- vector(mode = 'numeric',0)
q75 <- vector(mode = 'numeric',0)
q10 <- vector(mode = "numeric", 0)
q90 <- vector(mode = "numeric", 0)

#Create data frame with variables characterizing sampling effort and coverage over range of grid sizes
for (i in 1:length(Griddata)){
  
  print(paste("working on grid ~ ",(sqrt(Griddata[[i]]@polygons[[1]]@area))/1000,"km resolution",sep=" "))
  
  CellSize <- c(CellSize, (sqrt(Griddata[[i]]@polygons[[1]]@area))/1000)
  MaritimeGrid <- Griddata[[i]] 
  NumCell <- c(NumCell, length(MaritimeGrid))
  
  GridID<- over(Mardata,sp::geometry(MaritimeGrid))
  joinedgrid<-cbind(Mardata@data,GridID )
  populated<-(unique(joinedgrid$GridID))[complete.cases(unique(joinedgrid$GridID))]
  Populated <- c(Populated, length(populated))
  
  goodco <- joinedgrid[!is.na(joinedgrid$GridID),] 
  print(paste(nrow(goodco), 'sets in study area'))
  
  NumGrid <- data.frame(table(goodco$GridID))
  colnames(NumGrid) <- c("Grid","Frequency")
  NumGrid$Grid <- as.numeric(as.character(NumGrid$Grid))
  Mean <- c(Mean, mean(NumGrid$Frequency))
  Median <- c(Median, median(NumGrid$Frequency))
  SD <- c(SD, sd(NumGrid$Frequency))
  Max <- c(Max, max(NumGrid$Frequency))
  q25 <- c(q25,stats::quantile(NumGrid$Frequency, probs = 0.25))
  q75 <- c(q75, stats::quantile(NumGrid$Frequency, probs = 0.75))
  q10 <- c(q10, stats::quantile(NumGrid$Frequency, probs = 0.1))
  q90 <- c(q90, stats::quantile(NumGrid$Frequency, probs = 0.9))}

SamEffort <- data.frame(CellSize, NumCell, Populated, Mean, SD, Median, Max, q25, q75, q10, q90)

#Add column for average Simpson dissiliarity in site X site matrix using data
#superimposed on grids of various sizes
SamEffort$Mean.diss <- c(NA,NA,NA,NA,NA,0.328,NA,0.369,NA,0.416,NA,
                         0.477,0.511,0.566,0.579,0.591,0.605,NA,0.620,0.632,0.637)

theme_update(panel.border = element_blank(), panel.background = element_blank(), 
             panel.grid = element_blank(), plot.background = element_blank(), 
             axis.line = element_line(colour = "black", linetype = 1), 
             axis.text = element_text(size = 12, colour = "black"), 
             axis.title = element_text(size = 16, colour = "black"), 
             plot.title = element_text(size = 16, colour = "black", hjust = 0.5))

p1 <- ggplot(SamEffort[8:21,]) +
  geom_linerange(aes(x = CellSize, ymin = Mean - SD, ymax = Mean + SD)) +
  geom_crossbar(aes(x = CellSize, y = Median, ymin = q25, 
                    ymax = q75, fill = (Populated/NumCell)*100), fatten = 0.5) +
  geom_point(aes(x = CellSize, y = Mean)) + 
  guides(fill = guide_colorbar(title = 'Cells Populated (%)', title.position = "top")) +
  labs(x = "Grid cell size (km)", y = "Survey sets per grid cell") + 
  theme(legend.position = c(0.5,0.85), legend.direction = "horizontal")

p2 <- ggplot(SamEffort[8:21,]) +
  geom_point(aes(x = CellSize, y = Mean.diss), size = 2, col = 'black') +
  geom_smooth(aes(x = CellSize, y = Mean.diss), method = 'loess', se = F) +
  labs(x = "Grid cell size (km)", y = "Simpson distance")
  

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g <- rbind(g1, g2, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths)

tiff('E:/Biological Classification/Output/GridSelection_Maritimes.tiff', height = 664, width = 780, units = 'px')
grid.newpage()
grid.draw(g)
dev.off()



