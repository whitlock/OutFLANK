###### Fst functions for OutFLANK

####### Fst for haploids
#' 
#' Calculates FST with correction for local sample sizes, for haploid biallelic data. Based on Weir (1996 - Genetic Data Analysis II)
#' 
#' @title FST calculation for biallelic haploid data
#'
#' @param AllCounts This is an array with a row for each population, and two values per row: Number of alleles in the sample of one type,  number of alleles of other type.
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   He:  the expected heterozygosity of the locus
#'  \item   p_ave: the average allele frequency
#'  \item   FST:  Fst (with sample size correction)
#'  \item   T1: The numerator of the Fst calculation
#'  \item   T2: The denominator of the Fst calculation
#'  }
#'  @export
#'  
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

  He <- 2*p_ave*(1-p_ave)
  
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)

  
  T1 <- s2 - 1/(n_ave-1)*(p_ave*(1-p_ave) -(s2*(r-1)/r))
  T2 <- (n_c - 1)*(p_ave*(1-p_ave))/(n_ave-1) + (1 + (r-1)*(n_ave-n_c)/(n_ave-1))*s2/r
  
  FST <- T1/T2 
  
  return(c(He,p_ave, FST, T1, T2))
  
}

#####The following function does not correct for finite local sample size, and
#####therefore creates a biased estimate of Fst that does not go negative.

#' 
#' Calculates FST without correction for local sample sizes, for haploid biallelic data. This is necessary for using OutFLANK, which depends on these uncorrected values for reliable function. (Otherwise, sampling corrections can sometimes cause negative estimates of FST.)
#' 
#' @title FSTNoCorr calculation for biallelic haploid data
#'
#' @param AllCounts This is an array with a row for each population, and two values per row: Number of alleles in the sample of one type,  number of alleles of other type.
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   HeNoCorr:  the expected heterozygosity of the locus
#'  \item   p_aveNoCorr: the average allele frequency
#'  \item   FSTNoCorr:  Fst (without sample size correction)
#'  \item 	T1NoCorr: The numerator of the FSTNoCOrr calculation 
#'  \item   T2NoCorr: The denominator of the FSTNoCOrr calculation
#'  }
#'  @export
#'  
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

  He <- 2*p_ave*(1-p_ave)
  
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)

  
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
#' Calculates FST without correction for local sample sizes, for diploid biallelic data. This is necessary for using OutFLANK, which depends on these uncorrected values for reliable function. (Otherwise, sampling corrections can sometimes cause negative estimates of FST.)
#' 
#' @title FSTNoCorr calculation for biallelic diploid data
#'
#' @param Sample_Mat This is an array with a row for each population, and three values per row: Number of Homozygotes of one type, number of heterozygotes, number of homozygotes of other type.
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   He:  the expected heterozygosity of the locus
#'  \item 	FSTNoCorr:  Fst (without sample size correction)
#'  \item 	T1NoCorr: The numerator of the uncorrected sample size correction (similar to Weir and Cockerham 1984)
#'  \item   T2NoCorr: The denominator of the uncorrected sample size correction
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
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)

  if(s2==0){return(1); break}  
  
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  
  a <- n_ave/n_c*(s2)
  
  b <- (p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave)/(4*n_ave)*h_ave)
  
  c <- 1/2*h_ave
  
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  
  FST <- a/(a+b+c) 
  return(list(He=He,FSTNoCorr=FST, T1NoCorr=a, T2NoCorr=(a+b+c)))
}

#############FSt for diploids with local sample size corrections###############

##########################################
## WC FST for infinite sample of diploid allele freqs
###########################################  
#' 
#' Calculates FST with correction for local sample sizes, for diploid biallelic data. 
#' 
#' @title FST calculation for biallelic diploid data
#'
#' @param Sample_Mat This is an array with a row for each population, and three values per row: Number of Homozygotes of one type, number of heterozygotes, number of homozygotes of other type.
#' 
#' @return Returns a list of values related to FST:
#'  \itemize{
#'  \item   He:  the expected heterozygosity of the locus
#'  \item   FST:  Fst (with sample size correction)
#'  \item 	T1: The numerator of the Fst calculation (a from Weir and Cockerham 1984)
#'  \item   T2NoCorr: The denominator of the Fst calculation (a+b+c from Weir and Cockerham 1984)
#'  }
#'  @export
#' 
WC_FST_FiniteSample_Diploids_2Alleles<-function(Sample_Mat){
  
  #Sample Mat has three columns (homo_p,m heterozygotes, and homo_q) and a row for each population
  
  sample_sizes = rowSums(Sample_Mat)
  n_ave = mean(sample_sizes)
  n_pops = nrow(Sample_Mat) #r
  r = n_pops
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)
  p_freqs = (Sample_Mat[,1] + Sample_Mat[,2]/2) /sample_sizes
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)

  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  if(s2==0){return(1); break}	
  
  h_freqs = Sample_Mat[,2]/sample_sizes
  h_ave = sum(sample_sizes*h_freqs)/(n_ave*n_pops)
  
  a <- n_ave/n_c*(s2 - 1/(n_ave-1)*(p_ave*(1-p_ave)-((r-1)/r)*s2-(1/4)*h_ave))
  
  b <- n_ave/(n_ave-1)*(p_ave*(1-p_ave) - (r-1)/r*s2 - (2*n_ave - 1)/(4*n_ave)*h_ave)
  
  c <- 1/2*h_ave
  
  He <- 1-sum(p_ave^2, (1-p_ave)^2)
  
  FST <- a/(a+b+c) 
  return(list(He=He,FST=FST, T1=a, T2=(a+b+c)))
}

