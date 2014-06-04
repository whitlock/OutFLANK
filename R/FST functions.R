###### Fst functions for OutFLANK

WC_FST_FiniteSample_Haploids_2AllelesB_MCW<-function(AllCounts){
  #Input a matrix of the counts of each allele (columns) in each population (rows)
  #returns vector instead of list of Fst values, according to Weir
  
  n_pops<-dim(AllCounts)[1]
  r<-n_pops
  counts1 <- AllCounts[,1]
  sample_sizes <- rowSums(AllCounts)
  n_ave <- mean(as.numeric(sample_sizes))
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)  
  p_freqs = counts1/sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  #note: this differs slightly from mean(p_freqs) in R
  He <- 2*p_ave*(1-p_ave)
  
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  #note: this differs slightly from var(p_freqs) in R
  
  T1 <- s2 - 1/(n_ave-1)*(p_ave*(1-p_ave) -(s2*(r-1)/r))
  T2 <- (n_c - 1)*(p_ave*(1-p_ave))/(n_ave-1) + (1 + (r-1)*(n_ave-n_c)/(n_ave-1))*s2/r
  
  FST <- T1/T2 
  
  return(c(He,p_ave, FST, T1, T2))
  
}

#####The following function does not correct for finite local sample size, and
#####therefore creates a biased estimate of Fst that does not go negative.
WC_FST_FiniteSample_Haploids_2AllelesB_NoSamplingCorrection<-function(AllCounts){
  #Input a matrix of the counts of each allele (columns) in each population (rows)
  
  n_pops<-dim(AllCounts)[1]
  r<-n_pops
  counts1 <- AllCounts[,1]
  sample_sizes <- rowSums(AllCounts)
  n_ave <- mean(as.numeric(sample_sizes))
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)  
  p_freqs = counts1/sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  #note: this differs slightly from mean(p_freqs) in R
  He <- 2*p_ave*(1-p_ave)
  
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  #note: this differs slightly from var(p_freqs) in R
  
  T1NoCorr <- s2 
  T2NoCorr <- s2/r+(p_ave*(1 - p_ave))
  
  FSTNoCorr <- T1NoCorr/T2NoCorr 
  
  return(c(HeNoCorr=He,p_aveNoCorr=p_ave, FSTNoCorr=FSTNoCorr, T1NoCorr=T1NoCorr, T2NoCorr=T2NoCorr))
  
}

fstBarCalculatorNoCorr=function(DataList){
  #Calculates mean FstNoCorr from the dataframe, using sum(T1NoCorr) / sum(T2NoCorr) as the estimate of mean Fst.
  #Uses only data for which qvalues > qthreshold (i.e. $OutlierFlag==FALSE)
  #Does not internally screen for low MAF or low He values (but that can be added by only sending the
  #  high MAF rows to this function)
  sum(DataList$T1NoCorr[which(!DataList$OutlierFlag)])/sum(DataList$T2NoCorr[which(!DataList$OutlierFlag)])
}

fstBarCalculator=function(DataList){
  #Calculates mean Fst from the dataframe, using sum(T1) / sum(T2) as the estimate of mean Fst.
  #Uses only data for which qvalues > qthreshold (i.e. $OutlierFlag==FALSE)
  #Does not internally screen for low MAF or low He values (but that can be added by only sending the
  #  high MAF rows to this function)
  sum(DataList$T1[which(!DataList$OutlierFlag)])/sum(DataList$T2[which(!DataList$OutlierFlag)])
}

#####From FODR-- Diploid Fst from Weir and Cockerham 1984##########
##########################################
## WC FST for infinite sample of diploid allele freqs; without a correction for loacl sample size
###########################################  


#' 
#' Calculates FST without correct for local sample sizes. This is necessary for using OutFLANK, which depends on these uncorrected values for reliable function. (Otherwise, sampling corrections can someties cause negative estiamtes of FST.)
#' 
#' @title FSTNoCorr calcaution for biallelic diploid data
#'
#' @param Sample_Mat This is an array with a row for each population, and three values per row: Number of Homozygotes of one type, Nuber of heterozygotes, number of homozygotes of other type.
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   He:  the expected heterozygosity of the locus
#'  \item 	FST:  Fst (without sample size correction)
#'  \item 	T1: The numerator of the uncorrected sample size correction (similar to Weir and Cockerham 1984)
#'  \item   T2: The denominator of the uncorrected sample size correction
#'  }
#'  @export
#'  
WC_FST_FiniteSample_Diploids_2Alleles_NoCorr<-function(Sample_Mat){
  
  #Sample Mat has three columns (homo_p,m heterozygotes, and homo_q) and a row for each population
  
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = (Sample_Mat[,1] + Sample_Mat[,2]/2) /sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  #note: this differs slightly from mean(p_freqs) in R
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  #note: this differs slightly from var(p_freqs) in R
  if(s2==0){return(1); break}  
  
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  #note: this differs slightly from mean(h_freqs) in R  
  
  a <- n_ave/n_c*(s2)
  
  b <- (p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave)/(4*n_ave)*h_ave)
  
  c <- 1/2*h_ave
  
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  
  FST <- a/(a+b+c) 
  return(list(He=He,FST=FST, T1=a, T2=(a+b+c)))
}

#############FSt for diploids with local sample size corrections###############

##########################################
## WC FST for infinite sample of diploid allele freqs
###########################################  

WC_FST_FiniteSample_Diploids_2Alleles<-function(Sample_Mat){
  
  #Sample Mat has three columns (homo_p,m heterozygotes, and homo_q) and a row for each population
  
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = (Sample_Mat[,1] + Sample_Mat[,2]/2) /sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  #note: this differs slightly from mean(p_freqs) in R
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  #note: this differs slightly from var(p_freqs) in R
  if(s2==0){return(1); break}	
  
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  #note: this differs slightly from mean(h_freqs) in R	
  
  a <- n_ave/n_c*(s2 - 1/(n_ave-1)*(p_ave*(1-p_ave)-((r-1)/r)*s2-(1/4)*h_ave))
  
  b <- n_ave/(n_ave-1)*(p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave - 1)/(4*n_ave)*h_ave)
  
  c <- 1/2*h_ave
  
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  
  FST <- a/(a+b+c) 
  return(list(He=He,FST=FST, T1=a, T2=(a+b+c)))
}

