# Calculate the species richness in a sample
Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  count=sum(x[x>0.5]^0)
  return(count)
}

# Calculate the Effective species richness in each individual sample
Eff.Species.richness <- function(x)
{
  # Count only the OTUs that are present more than the set proportion
  total=sum(x)
  count=sum(x[x/total>eff.cutoff]^0)
  return(count)
}

# Calculate the Normalized species richness in each individual sample
Norm.Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  # Given a fixed Normalization reads depth
  total=sum(x)
  count=sum(x[norm.cutoff*x/total>0.5]^0)
  return(count)
}


# Calculate the Shannon diversity index
Shannon.entropy <- function(x)
{
  total=sum(x)
  se=-sum(x[x>0]/total*log(x[x>0]/total))
  return(se)
}

# Calculate the effective number of species for Shannon
Shannon.effective <- function(x)
{
  total=sum(x)
  se=round(exp(-sum(x[x>0]/total*log(x[x>0]/total))),digits =2)
  return(se)
}

# Calculate the Simpson diversity index
Simpson.concentration <- function(x)
{
  total=sum(x)
  si=sum((x[x>0]/total)^2)
  return(si)
}

# Calculate the effective number of species for Simpson
Simpson.effective <- function(x)
{
  total=sum(x)
  si=round(1/sum((x[x>0]/total)^2),digits =2)
  return(si)
}

writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, 
                   as.character(paste(">", data[rowNum,"ASV"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"Sequence"]))
  }
  
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}