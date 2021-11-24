
suppressPackageStartupMessages(lapply(c('plyr','picante','ggfortify', 'psych','viridis', 'ggplot2','RColorBrewer', 'reshape2', 'ape', 'phyloseq', 'vegan', 'cowplot', 'dplyr', 'gplots', 'dada2', 'phangorn', 'data.table', 'hillR', 'betapart', 'tidyverse', 'matrixStats'), require,character.only=TRUE)) #add as necessary


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



files <- list.files(path=".",pattern=".RDS")

for (f in 1:length(files)) {
  dat <- readRDS(paste0("iter",f,".RDS"))
 ee.[ , ncol(ee.) + 1]<- dat
  colnames(ee.)[ncol(ee.)] = paste0("iter",f)

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
