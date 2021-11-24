# DisturbanceSynthesis
Code for the bioinformatics process and null models used in  Disturbance and recovery: a synthesis of microbial community assembly following disturbance across realms

The repository includes a dada2 wrapper (fully.process.single.R), sequence quality data (Tracking_table), the main phyloseq sequence file used (merged.ecological.rds), and code for the calculation of richness, dispersion, turnover (Ecological_metaanalysis.Sept2021.Rmd), and null model-based Z-scores (NullModel.R and ExtractNullToZScores.r). 

For the calculation of Z-Scores, NullModel.R is called iteratively within a high performance computing environment, and the output files are read by ExtractNullToZScores.r to calculate z-scores. 
