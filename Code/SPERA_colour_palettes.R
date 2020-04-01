#### Colour Palettes for SPERA Bioclassification Project ####

# library(colorspace)
# library(wesanderson)

# Cluster colours in 4 regions are selected from following palettes #
# Rush <- as.character(wes_palette('Rushmore1',6,'continuous'))
# Rush2 <- as.character(wes_palette('Rushmore1',10,'continuous'))
# Rush3 <- darken(Rush1, amount = 0.5)
# Rush4 <- darken(Rush2, amount = 0.5)
# Royal <- as.character(wes_palette('Royal1'))
# Zissou <- as.character(wes_palette('Zissou1')) %>% 
#   desaturate(., amount = 0.2)
  

MAR.palette <- data.frame(
  assigned = c("#9FA682", # light olive green
               "#E1BD6D", # light goldenrod
               "#649373", # dark sea green
               "#1B5656", # teal
               "#F2300F", # orangered2
               "#5A283E", # plum
               "#e0e0e0"), #gray88 
  cl = factor(c(9,8,4,2,5,1,NA)),
  name = c("WSS/Outer BoF",
         "WSS: Banks/Inner BoF",
         "ESS",
         "ESS: Banks",
         "Laurentian Channel/Shelf Break",
         "Slope",
         "Unclassified"), 
  stringsAsFactors = F
  )

QC.palette <- data.frame(
  assigned = c("#C93312", # orangered3
               "#FAEFD1", # papayawhip
               "#DC863B", # peru
               "#e0e0e0"), # gray88 
  cl = factor(c(6,1,7,NA)),
  name = c("Deep Channels", 
         "Shallow Banks & Straits",
         "Channel Heads & Slopes",
         "Unclassified"),
  stringsAsFactors = F
  )

NL.palette <- data.frame(
  assigned = c("#85B5C0", # light sky blue
               "#3F7F92", # cadet blue
               "#E1AF00", # sub yellow
               # "#EBCC2A", # yellow
               # "#DAB051", # tan
               "#3F0A25", # merlot
              "#E23B35", # reddish brown
              "#e0e0e0"), # gray88 
  cl = factor(c(7,8,6,1,4,NA)),
  name = c("Inner Shelf",
         "Outer Shelf",
         "Grand Banks",
         "Slope", 
         "Laurentian Channel/Shelf Break",
         "Unclassified"),
  stringsAsFactors = F
)

GULF.palette <- data.frame(
  assigned = c("#899DA4", # gray61
               "#002F2F", # dark slate gray
               "#7F1100", # dark red
               "#725900", # orange4
               "#e0e0e0"), # gray88
  cl = factor(c(7,6,9,3,NA)),
  name = c("Magdalen Shallows" ,
         "Inshore/Magdalen Is.",
         "Laurentian Channel",
         "Northumberland Strait/St. George's Bay",
         "Unclassified"),
  stringsAsFactors = F
)
