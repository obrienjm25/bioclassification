#Function to get basic taxonomic information

Classify <- function(x,newsciname=T){
  
  output <- taxize::classification(x,db="itis",return_id=FALSE)
  output2 <- as.data.frame(output[[1]],stringsAsfactors=F)
  output3 <- as.data.frame(t(output2$name),stringsAsFactors = F)
  names(output3) <- output2$rank
  
  PhyloNames <- c("phylum","subphylum","class","order","family","genus","species") #Info we want
  
  PhyloDiff <- intersect(PhyloNames,names(output3)) # the columns which are available
  
  #If some information is not available, fill it with NAs
  if(length(PhyloDiff)<length(PhyloNames)){
    
    ouput3 <- output3[,PhyloDiff]
    
    #what columns are missing
    MissingCols <- setdiff(PhyloNames,PhyloDiff)
  
    #create a dummy dataframe with NAs and the corresponding missing columns
    MissingDF <- data.frame(t(rep(NA,length(MissingCols))))
    colnames(MissingDF) <- MissingCols
    
    #append dummy 'NA' dataframe to the data we do have
    output3 <- cbind(output3,MissingDF)
    
    #ensure the output is in the right order
    output3 <- output3[,PhyloNames]
  
  
  } else {output3 <- output3[,PhyloNames]} # if we actually had all the columns
  
 
  if(newsciname){output3$newsciname=x} #if we want to add the database entry name for merging later (default = T)
  
  return(output3) #return the dataframe row
  
}


