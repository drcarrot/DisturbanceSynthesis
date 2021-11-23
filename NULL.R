
suppressPackageStartupMessages(lapply(c('plyr','picante','ggfortify', 'psych','viridis', 'ggplot2','RColorBrewer', 'reshape2', 'ape', 'phyloseq', 'vegan', 'cowplot', 'dplyr', 'gplots', 'dada2', 'phangorn', 'data.table', 'hillR', 'betapart', 'tidyverse', 'Carrotworks'), require,character.only=TRUE)) #add as necessary


#Function from metagMisc https://github.com/vmikk/metagMisc/
phyloseq_sep_variable <- function(physeq, variable, drop_zeroes = T){
   ## Check the input
  if(is.null(phyloseq::sample_data(physeq, errorIfNULL = F))){
    stop("Sample data is missing in the phyloseq-object.\n")
  }
  
  ## Extract samle meta-data
  mtd <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")
  
  if(!variable %in% colnames(mtd)){
    stop("Grouping variable is missing from the sample data of phyloseq-object.\n")
  }
  
  if(class(mtd[, variable]) %in% c("integer", "numeric") ){
    if( length( unique(mtd[, variable]) ) > 5){
      stop("Groupping variable is numeric and it has too many levels. Consider transforming it to factor.\n")
    } else {
      warning("Groupping variable is numeric and it was coerced to factor.\n")
      mtd[, variable] <- factor(mtd[, variable])
    }
  }
  
  if(length(table(mtd[, variable])) == 1){
    cat("Warning: there is only one group of samples in the resulting list.\n")
  }
  
  ## Add sample IDs to the meta-data
  smp <- data.frame(
    SID = phyloseq::sample_names(physeq),
    mtd,
    stringsAsFactors = F)
  
  ## Exatract sample names by the specified variable
  svv <- plyr::dlply(.data = smp, .variables = variable, .fun = function(z){ z$SID })
  
  ## Extract samples by groupping variable
  res <- plyr::llply(.data = svv, .fun = function(z){ phyloseq::prune_samples(z, x = physeq) })
  
  ## Remove taxa with zero abundance
  if(drop_zeroes == TRUE){
    res <- plyr::llply(.data = res, .fun = function(x){ phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) })
  }
  
  return(res)
}


#from main notebook
merged.pruned.rp=readRDS("merged.pruned.rp.rds")
merged.pruned.rp.=phyloseq_sep_variable(merged.pruned.rp, "Time_series", drop_zeroes = FALSE)


niter=1000 #iterations 

bb=lapply(merged.pruned.rp., dist2listt) #convert to dataframe
for(i in 1:length(bb)){
  tmp=bb[[i]]
  if (i>1){
    ee.=rbind(ee., tmp)
  } else {
    ee.=tmp}
}


for(i in 1:niter){
  merged.pruned.rp.0=vector("list", length = length(merged.pruned.rp.) )
  for (j in 1: length (merged.pruned.rp.) ){
    #reshuffle 
    merged.pruned.rp.0[[j]]= phyloseq(otu_table(
      picante::randomizeMatrix(merged.pruned.rp.[[j]]@otu_table@.Data, "richness", 1000), taxa_are_rows = FALSE), 
      tax_table(merged.pruned.rp.[[j]]@tax_table), 
      sample_data(merged.pruned.rp.[[j]]@sam_data))
  }
    bb0=lapply(merged.pruned.rp.0, dist2listt)
    for(k in 1:length(bb0)){
      tmp=bb0[[k]]
      if (k>1){
        cc=rbind(cc, tmp)
      } else {
        cc=tmp
        }
    }
    ee.[[paste0("iter",niter)]]=cc$value
}

  
ee.$Disturbance=merged.pruned.rp@sam_data$Disturbance[match(ee.$col,row.names(merged.pruned.rp@sam_data))]

ee.$Time_series=merged.pruned.rp@sam_data$Time_series[match(ee.$col,row.names(merged.pruned.rp@sam_data))]
ee.$Time_series.=merged.pruned.rp@sam_data$Time_series[match(ee.$row,row.names(merged.pruned.rp@sam_data))]

ee.$Environment=merged.pruned.rp@sam_data$Environment[match(ee.$col,row.names(merged.pruned.rp@sam_data))]

ee.$Study=merged.pruned.rp@sam_data$Study[match(ee.$col,row.names(merged.pruned.rp@sam_data))]

ee.$Experimental_unit=merged.pruned.rp@sam_data$Experimental_Unit_ID[match(ee.$col,row.names(merged.pruned.rp@sam_data))]

ee.$Experimental_unit.=merged.pruned.rp@sam_data$Experimental_Unit_ID[match(ee.$row,row.names(merged.pruned.rp@sam_data))]

ee.$Selective_Factor=merged.pruned.rp@sam_data$Selective_Factor[match(ee.$row,row.names(merged.pruned.rp@sam_data))]

ee.$Time_series_Time=merged.pruned.rp@sam_data$Time_series_Time[match(ee.$col,row.names(merged.pruned.rp@sam_data))]

#add rank data to the main file
merged.pruned.rp@sam_data$Rank=ee.$Rank[match(merged.pruned.rp@sam_data$Time_series_Time,ee.$Time_series_Time,)]

Dispersions=ee.[which(ee.$Time_since_dist==ee.$Time_since_dist. &
                        (ee.$Experimental_unit!=ee.$Experimental_unit.|
                           is.na(ee.$Experimental_unit))&
                        ee.$Time_series==ee.$Time_series.),]
#Dispersions for Shane
write.table(Dispersions, "dispersions.txt")
