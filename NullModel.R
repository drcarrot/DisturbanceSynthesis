
suppressPackageStartupMessages(lapply(c('plyr','picante','ggfortify', 'psych','viridis', 'ggplot2','RColorBrewer', 'reshape2', 'ape', 'phyloseq', 'vegan', 'cowplot', 'dplyr', 'gplots', 'dada2', 'phangorn', 'data.table', 'hillR', 'betapart', 'tidyverse', 'matrixStats'), require,character.only=TRUE)) #add as necessary


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



dist2list <- function (dist, tri=TRUE) {
  if (!class(dist) == "dist") { stop("Error: The input data must be a dist object.\n") }

  dat <- as.data.frame(as.matrix(dist))
  if (is.null(names(dat))) {
    rownames(dat) <- paste(1:nrow(dat))
  }
  value <- stack(dat)$values
  rnames <- rownames(dat)
  namecol <- expand.grid(rnames, rnames)
  colnames(namecol) <- c("col", "row")
  res <- data.frame(namecol, value)

  if(tri == TRUE){    # return only lower triangular part of dist
    res <- res[-which(upper.tri(as.matrix(dist), diag = T)), ]
  }

  return(res)
}


dist2listt <-function(physeq){
  dist=dist2list(vegdist(decostand(otu_table(physeq),method = "hellinger")), tri=TRUE)
  dist$Study=physeq@sam_data$Study[1]
  dist$Time_series=physeq@sam_data$Time_series[1]
  dist$Time_since_dist=physeq@sam_data$Time_since_dist[match(dist$col, sample_names(physeq))]
  dist$Time_since_dist.=physeq@sam_data$Time_since_dist[match(dist$row, sample_names(physeq))]
  return(dist)
}

merged.pruned.rp=readRDS("merged.pruned.rp.rerarefied.rds")

merged.pruned.rp.=phyloseq_sep_variable(merged.pruned.rp, "Time_series", drop_zeroes = FALSE)

#Actual data
bb=lapply(merged.pruned.rp., dist2listt) #a list of linearized distance matrices

key=list()
for(i in 1:length(bb)){
  key[[i]]=data.frame(levels=as.numeric(levels(as.factor(merged.pruned.rp.[[i]]@sam_data$Time_since_dist))),
                 rank=rank(as.numeric(levels(as.factor(merged.pruned.rp.[[i]]@sam_data$Time_since_dist)))))
  bb[[i]]$Rank=key[[i]]$rank[match(bb[[i]]$Time_since_dist, key[[i]]$levels)]#rank of first Time column
  bb[[i]]$Rank.=key[[i]]$rank[match(bb[[i]]$Time_since_dist., key[[i]]$levels)]#rank of second Time column
}

#Convert bb to data frame 
for(i in 1:length(bb)){
  tmp=bb[[i]]
  if (i>1){
    ee.=rbind(ee., tmp)
  } else {
    ee.=tmp}
}

niter=3 #iterations 

for(i in 1:niter){
  merged.pruned.rp.0=vector("list", length = length(merged.pruned.rp.) )
  for (j in 1: length (merged.pruned.rp.) ){
    #reshuffle 
    merged.pruned.rp.0[[j]]= phyloseq(otu_table(
      picante::randomizeMatrix(merged.pruned.rp.[[j]]@otu_table@.Data, "richness", niter), taxa_are_rows = FALSE), 
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
    ee.[ , ncol(ee.) + 1] <- cc$value            
    colnames(ee.)[ncol(ee.)] <- paste0("iter",i)

}


ee.$Disturbance=merged.pruned.rp@sam_data$Disturbance[match(ee.$col,row.names(merged.pruned.rp@sam_data))]
ee.$Environment=merged.pruned.rp@sam_data$Environment[match(ee.$col,row.names(merged.pruned.rp@sam_data))]
ee.$Study=merged.pruned.rp@sam_data$Study[match(ee.$col,row.names(merged.pruned.rp@sam_data))]
ee.$Experimental_unit=merged.pruned.rp@sam_data$Experimental_Unit_ID[match(ee.$col,row.names(merged.pruned.rp@sam_data))]
ee.$Experimental_unit.=merged.pruned.rp@sam_data$Experimental_Unit_ID[match(ee.$row,row.names(merged.pruned.rp@sam_data))]
ee.$Selective_Factor=merged.pruned.rp@sam_data$Selective_Factor[match(ee.$row,row.names(merged.pruned.rp@sam_data))]
ee.$Time_series_Time=merged.pruned.rp@sam_data$Time_series_Time[match(ee.$col,row.names(merged.pruned.rp@sam_data))]

# add a homogenized time
for(i in 1:dim(ee.)[1]){
if (ee.$Rank.[i]>ee.$Rank[i]){
ee.$Time[i]=ee.$Time_since_dist.[i]}
else {ee.$Time[i]=ee.$Time_since_dist[i]}
}

ee.$nmodel.means=rowMeans(ee.%>% select(starts_with('iter')))
ee.$nmodel.sds=rowSds(as.matrix(ee.%>% select(starts_with('iter'))))

ee.$nmodel.sds[ee.$nmodel.sds==0]<-min(ee.$nmodel.sds[ee.$nmodel.sds>0])

ee.$nmodel.Zscores=(ee.$value-ee.$nmodel.means)/ee.$nmodel.sds


write.table(ee., "allnullmodeloutputs.txt")

Resilience=ee.[which(ee.$Time_since_dist!=ee.$Time_since_dist. &
                       (ee.$Rank==1|ee.$Rank.==1)&
                       (ee.$Experimental_unit!=ee.$Experimental_unit.|
                          is.na(ee.$Experimental_unit))),]

Resilience <- Resilience %>% select(-contains("iter"))
write.table(Resilience, "Resilience.zscores.txt")

Dispersions=ee.[which(ee.$Time_since_dist==ee.$Time_since_dist. &
                        (ee.$Experimental_unit!=ee.$Experimental_unit.|
                           is.na(ee.$Experimental_unit))),]
Dispersions <- Dispersions %>% select(-contains("iter"))
write.table(Dispersions, "dispersions.zscores.txt")
