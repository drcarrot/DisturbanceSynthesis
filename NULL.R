
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


#from https://rdrr.io/github/vmikk/metagMisc/src/R/phyloseq_randomize.R
phyloseq_randomize <- function(physeq, null_model = "phylogeny.pool", verbose = T, ...){
  # c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")
  
  # require(phyloseq)
  # require(picante)
  # require(vegan)
  
  ## Picante models
  pm <- c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")
  
  ## Vegan models
  # Binary null models
  v1 <- c("r00", "r0", "r0_old", "r1", "r2", "c0", "swap", "tswap", "curveball", "quasiswap", "backtracking", "backtrack")
  
  # Quantitative Models for Counts with Fixed Marginal Sums
  v2 <- c("r2dtable", "quasiswap_count")
  
  # Quantitative swap models
  v3 <- c("swap_count", "abuswap_r", "abuswap_c")
  
  # Quantitative Swap and Shuffle
  v4 <- c("swsh_samp", "swsh_both", "swsh_samp_r", "swsh_samp_c", "swsh_both_r", "swsh_both_c")
  
  # Quantitative Shuffle Methods
  v5 <- c("r00_ind", "r0_ind", "c0_ind", "r00_samp", "r0_samp", "c0_samp", "r00_both", "r0_both", "c0_both")
  
  # All vegan models
  vv <- c(v1, v2, v3, v4, v5)
  
  ## Print implementation details
  if(verbose == TRUE){
    if(null_model %in% pm){
      cat("Randomization null model is based on implementation from the 'picante' package.\n")
    }
    if(null_model %in% vv){
      cat("Randomization null model is based on implementation from the 'vegan' package.\n")
    }
    cat("Please cite it in the publications.\n")
  }
  
  ## TO DO - vegan models
  if(null_model %in% vv){
    stop("Vegan null models are currently not yet implemented.\n")
  }
  
  ### Extract phyloseq components
  ## OTU table
  comm <- as.data.frame(phyloseq::otu_table(physeq))
  if(!phyloseq::taxa_are_rows(physeq)){
    comm <- t(comm)
  }
  
  ## Phylo tree
  tree_null <- is.null(phyloseq::phy_tree(physeq, errorIfNULL=F))
  if(!tree_null){
    phy <- phyloseq::phy_tree(physeq)
  }
  
  
  if(null_model == "taxa.labels"){
    if(tree_null){ # No phylogeny -> shuffle taxa names
      rownames(comm) <- sample(phyloseq::taxa_names(physeq))
    } else {       # Phylogeny present -> shuffle tip names
      phy <- picante::tipShuffle(phy)
    }
  } # end of "taxa.labels"
  
  ## Randomize community data matrices with picante::randomizeMatrix
  if(null_model %in% c("frequency", "richness", "independentswap", "trialswap")){
    
    if(null_model == "sample.pool"){ null_model <- "richness" }  # this is the same models?
    # https://github.com/skembel/picante/blob/649edc7938b878429914c617e22c67198a8c189a/R/phylodiversity.R#L179
    
    comm <- t( picante::randomizeMatrix(t(comm), null.model = null_model, ...) )
  }
  if(null_model == "phylogeny.pool"){
    if(tree_null){ stop("Error: phylogeny tree is not available; therefore 'phylogeny.pool' model is not applicable.\n") }
    comm <- t( picante::randomizeMatrix(t(comm), null.model = "richness", ...) )
    phy <- picante::tipShuffle(phy)
  }
  
  ## Replace phyloseq slots with the randomized ones
  if(!phyloseq::taxa_are_rows(physeq)){  # transpose OTU table back
    comm <- t(comm)
  }
  phyloseq::otu_table(physeq) <- phyloseq::otu_table(comm, taxa_are_rows = TRUE)
  
  if(!tree_null){
    phyloseq::phy_tree(physeq) <- phy
  }
  
  return(physeq)
}


dist2listt <-function(physeq){
  dist=dist2list(vegdist(decostand(otu_table(physeq),method = "hellinger")), tri=TRUE)
  dist$Study=physeq@sam_data$Study[1]
  dist$Time_series=physeq@sam_data$Time_series[1]
  dist$Time_since_dist=physeq@sam_data$Time_since_dist[match(dist$col, sample_names(physeq))]
  dist$Time_since_dist.=physeq@sam_data$Time_since_dist[match(dist$row, sample_names(physeq))]
  return(dist)
}


merged.pruned.rp=readRDS("merged.pruned.rp.rds")


merged.pruned.rp.=phyloseq_sep_variable(merged.pruned.rp, "Time_series", drop_zeroes = FALSE)


niter=1 #iterations 

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