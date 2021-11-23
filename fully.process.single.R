################################################
## Processing and analyses for Miseq-Short run
################################################
## Should be run in an R session, with working
## directory containing oligo-trimmed fastq files
## and the TRAININGSET file. Available for download @
## https://benjjneb.github.io/dada2/training.html

################################################
## Parameters specific to project.
## Modify as necessary.
################################################


################################################
## Load required packages
################################################
lapply(c('ggfortify', 'psych','msa', 'ggplot2','RColorBrewer', 'reshape2', 'ape', 'phyloseq', 'vegan', 'cowplot', 'dplyr', 'gplots', 'dada2', 'phangorn'), require,character.only=TRUE) #add as necessary

fully.process<-function(forward_read_names=".fastq.gz", seed=1, displayx=10, filetype=".jpg",trimleft=10, studyname="study", trainingset="silva_nr_v132_train_set.fa.gz", truncLenF=240, tree=TRUE,multithread=TRUE){
  print("A log of processing is saved as log.txt")
  sink(file="log.txt")
  paste0("forward_read_names=", forward_read_names, "seed=",seed,"displayx=", displayx,"trimleft=",trimleft, "trainingset=", trainingset, "truncLenF=", truncLenF,  "tree=", tree, "multithread=", multithread)

################################################
##Analyze sequence quality for cutoffs
################################################

print("A log of processing is saved as log.txt")
path= getwd()
fnFs= sort(list.files(path, pattern= forward_read_names))
print(paste(length (fnFs), "samples were identified."))
sample.names = sapply(strsplit(fnFs, ".fastq.gz"), `[`, 1)  #extract sample names
fnFs = file.path(path, fnFs) #specify global path to the sequence files
filt_path= file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs=file.path(filt_path, paste0(sample.names, forward_read_names))#Create a subdirectory and file names for the filtered files

#Plot quality profile for forward reads for the first sample...
Forwardpqp=plotQualityProfile(fnFs[1])
ggsave(paste0("Qualityprof",studyname, filetype),Forwardpqp)

################################################
##Begin Processing
################################################

print(paste("A plot of the sequence quality of the first sample has been plotted and saved as Qualityprof", filetype))
#Filter and trim sequences according to above plots. The current parameters are standard according to tutorial.
#NOTE: this takes time
print("Sequences are being filtered and trimmed. This takes time.")
out=filterAndTrim(fnFs, filtFs,  truncLen=truncLenF, maxN=0, maxEE=2, truncQ=2, trimLeft=trimleft,rm.phix=TRUE, compress=TRUE, multithread=multithread)
print(paste("the number of files filtered is", length(list.files("filtered"))))

#Calculate error rates.
#NOTE: this takes time
print("Error rates are being calculated. This takes time.")
errF = learnErrors(filtFs, multithread=multithread)

#Dereplicate
derepFs = derepFastq(filtFs, verbose=TRUE)
print("Samples are dereplicated")
names(derepFs) = sample.names # Name the derep-class objects by the sample names

#Infer sequence variants
print("Sequence variants are being inferred. This takes time.")
dadaFs = dada(derepFs, err=errF, multithread=multithread)
print("Some insight into the processing of the first Forward sample")
print(dadaFs[[1]])

#Construct a sequence table (OTU table)
seqtab = makeSequenceTable(dadaFs)
print("The sequence table has been constructed.")

#Remove chimeras
print("Chimeras are being removed. This may take time")
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=multithread, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.nochim.rds")
dim(seqtab.nochim)

print(paste0("The percentage of data which remained after chimera removal is", sum(seqtab.nochim)/sum(seqtab)))

#Make tracking table
getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) = c("Original Reads", "Filtered Reads", "Denoised Reads",  "tabled", "Non-chimeric Reads")
rownames(track) = sample.names
track1= track[, -6]
write.table(track1, "tracking_table.txt", sep="\t")
print("This is what your tracking table looks like.")
print(head(track1))

#Assign taxonomy
print("Assigning taxonomy, this may take a while.")
taxa = assignTaxonomy(seqtab.nochim, trainingset, multithread=multithread, tryRC=TRUE)
unname(head(taxa))
taxa=as.data.frame(taxa)
#Modify taxa names so that they carry over their deepest classification down the taxonomy.
taxa[] <- lapply(taxa, as.character)
taxa$Phylum[is.na(taxa$Phylum)]<-taxa$Kingdom[is.na(taxa$Phylum)]
taxa$Class[is.na(taxa$Class)]<-taxa$Phylum[is.na(taxa$Class)]
taxa$Order[is.na(taxa$Order)]<-taxa$Class[is.na(taxa$Order)]
taxa$Family[is.na(taxa$Family)]<-taxa$Order[is.na(taxa$Family)]
taxa$Genus[is.na(taxa$Genus)]<-taxa$Family[is.na(taxa$Genus)]
taxa[] <- lapply(taxa, as.factor)
taxa=as.matrix(taxa)

if (tree=="TRUE"){
  print("Take a seat, this can take several hours.")
  seqs <- getSequences(seqtab.nochim)
  names(seqs) <- seqs # This propagates to the tip labels of the tree
  mult <- msa(seqs, method="ClustalW", type="dna", order="input")
  phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
  dm <- dist.ml(phang.align)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data=phang.align)
  ## negative edges length changed to 0!
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
  detach("package:phangorn", unload=TRUE)
  saveRDS(fitGTR$tree, "tree.RDS")
  merged = phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),  tax_table(taxa), phy_tree(fitGTR$tree))
  merged@sam_data=sample_data(as.data.frame(sample_names(merged)))
  sample_names(sample_data(merged))=sample_names(merged)
  sample_data(merged)[ , 2] <- sample_data(merged)[ ,1]
  colnames(merged@sam_data)=c("names", "dummy")
  merged = subset_taxa(merged, Kingdom== "Bacteria")
  saveRDS(merged, "merged.rds")
} else{
  merged = phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),  tax_table(taxa))
  merged@sam_data=sample_data(as.data.frame(sample_names(merged)))
  sample_names(sample_data(merged))=sample_names(merged)
  sample_data(merged)[ , 2] <- sample_data(merged)[ ,1]
  colnames(merged@sam_data)=c("names", "dummy")

  }
print(paste(c("The number of reads  assigned to Bacteria is", sum(taxa_sums(subset_taxa(merged, Kingdom=="Bacteria"))))))
merged = subset_taxa(merged, Kingdom== "Bacteria")
saveRDS(merged, "merged.rds")
print("Your phyloseq object has been created and saved as merged.RDS. All reads unclassified at the Kingdom level have been removed")

################################################
##Make Basic figures
################################################
#reads per OTU and per sample figure
reads_per_OTU=data.frame(nreads = sort(taxa_sums(merged), TRUE), sorted = 1:ntaxa(merged), type = "OTUs")
rpo=ggplot(data=reads_per_OTU, aes(x = sorted, y = nreads,fill=TRUE))+ geom_area(stat="identity")+scale_x_log10()+ scale_fill_manual(values= "#69BE28")+ theme(legend.position="none")+labs(x = "(log) Unique OTUs", y="Reads per OTU", title="Reads per OTU")
reads_per_sample=data.frame(nreads = sort(sample_sums(merged),TRUE), sorted = 1:nsamples(merged), type = "Samples")
rps=ggplot(data=reads_per_sample, aes(x = sorted, y = nreads,fill=TRUE))+ geom_bar(stat="identity")+ scale_fill_manual(values= "#69BE28")+ theme(legend.position="none")+labs(x = "Samples", y="Reads per sample", title="Reads per sample")
readinfo=plot_grid(rpo, rps, labels=c("A", "B"), ncol = 2, nrow = 1)
ggsave(paste0("readinfo",studyname,filetype), plot=readinfo)


tracking=read.table("tracking_table.txt", header=TRUE, row.names=1, sep="\t")
tracking$pseq=sample_sums(merged)
tracking$study=studyname
write.table(tracking, "tracking_table.txt", sep="\t")


Goods=1-(sum(taxa_sums(merged)==1)/sum(taxa_sums(merged)))
print(paste("Good's coverage index is", Goods))

#Prevalence vs. abundance plot
prev0 = apply(X = otu_table(merged), MARGIN = ifelse(taxa_are_rows(merged), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0, TotalAbundance = taxa_sums(merged), tax_table(merged))
prevalence.plot=ggplot(prevdf, aes(TotalAbundance, Prevalence, color = Phylum)) + geom_point(size = 2, alpha = 0.7, show.legend=FALSE) + scale_y_log10() + scale_x_log10() + xlab("Total Abundance") +facet_wrap(~Phylum)
ggsave(plot=prevalence.plot, paste0("prevalence_plot_Phyla",studyname,filetype))

#Richness plot
Richness=plot_richness(merged, measures="Observed")+ggsave(paste("Unrarefied_Observed_Richness",studyname,filetype))

#Bray-curtis ordination.
merged_ord_PCoA=ordinate(merged, "PCoA", "bray")
merged_ord_plot_PCoA=plot_ordination(merged, merged_ord_PCoA)+ geom_point(size=4)+geom_text(aes(label=names),hjust=-0.3, vjust=-0.5, size=2)+ggsave(paste("Unrarefied_Bray_PCOA",studyname,filetype))

print ("These are the general abundances in your entire experiment")
melted=psmelt(merged)
print(summary(melted))

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

sink()
savehistory(file = "history.Rhistory")
}
