ee.=read.table("nmodelWithMeans.txt")
merged.pruned.rp=read_rds("merged.pruned.rp.rds")

pairwise=read.table ("pairwise.table.txt")
colnames(pairwise)
ee.$Rank=pairwise$Rank[match(ee.$col,pairwise$col)]
ee.$Rank.=pairwise$Rank.[match(ee.$row,pairwise$row)]
# add a homogenized time
for(i in 1:dim(ee.)[1]){
if (ee.$Rank.[i]>ee.$Rank[i]){
ee.$Time[i]=ee.$Time_since_dist.[i]}
else {ee.$Time[i]=ee.$Time_since_dist[i]}
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
ee.$Rank=pairwise$Rank[match(ee.$col,pairwise$col)]
ee.$Rank.=pairwise$Rank.[match(ee.$row,pairwise$row)]
Dispersions=ee.[which(ee.$Time_since_dist==ee.$Time_since_dist. &
(ee.$Experimental_unit!=ee.$Experimental_unit.|
is.na(ee.$Experimental_unit))&
ee.$Time_series==ee.$Time_series.),]

ee.$nmodel.sds[ee.$nmodel.sds==0]<-min(ee.$nmodel.sds[ee.$nmodel.sds>0])
ee.$nmodel.Zscores=(ee.$value-ee.$nmodel.means)/ee.$nmodel.sds


Resilience=ee.[which(ee.$Time_since_dist!=ee.$Time_since_dist. &
(ee.$Rank==1|ee.$Rank.==1)&
(ee.$Experimental_unit!=ee.$Experimental_unit.|
is.na(ee.$Experimental_unit))&
ee.$Time_series==ee.$Time_series.),]
write.table(Resilience, "Resilience.zscores.txt")

Dispersions=ee.[which(ee.$Time_since_dist==ee.$Time_since_dist. &
(ee.$Experimental_unit!=ee.$Experimental_unit.|
is.na(ee.$Experimental_unit))&
ee.$Time_series==ee.$Time_series.),]
write.table(Dispersions, "dispersions.zscores.txt")

