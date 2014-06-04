
#OutFLANK:  An Fst outlier approach by Mike Whitlock and Katie Lotterhos, University of British Columbia.
#Development supported by AdapTree, Genome Canada, Genome BC, and an NSERC Discovery Grant to MCW.

########################How to use OutFLANK##################

#This method looks for Fst outliers from a list of Fst's for different loci. It
#assumes that each locus has been genotyped in all populations.

#OutFLANK estimates the distribution of Fst based on a trimmed sample of Fst. It
#assumes that the majority of loci in the center of the distribution are
#neutral, and infers the shape of the distribution of neutral Fst using a trimmed set of
#loci. Loci with the highest and lowest Fst's are trimmed from the data set
#before this inference, and the distribution of Fst df/(mean Fst)  is assumed to
#follow a chi-square distribution. Based on this inferred distribution, each
#locus is given a q-value based on its quantile in the inferred null
#distribution.

#The main procedure is called OutFLANK -- see comments in that
#function immediately below for input and output formats. The other functions
#here are necessary and must be uploaded, but are not necessarily needed by the
#user directly.

#Steps:
# 1. Make sure you have the biocLite package on your computer.  Code for getting it is commented in the next section.
#
# 2. Load all the functions in this script.
#
# 3. Create a file that has a row for each locus in your data set, with the following columns:
#     $LocusName: a character string that uniquely names each locus. 
#     $FST: Fst calculated for this locus. (Kept here to report the unbiased Fst of the results) 
#     $T1: The numerator of the estimator for Fst (necessary, with $T2, to calculate mean Fst) 
#     $T2: The denominator of the estimator of Fst 
#     $FSTNoCorr: Fst calculated for this locus without sample size correction. (Used to find outliers) 
#     $T1NoCorr: The numerator of the estimator for Fst without sample size correction (necessary, 
#                with $T2NoCorr, to calculate mean Fst) 
#     $T2NoCorr: The denominator of the estimator of Fst without sample size correction 
#     $He: The heterozygosity of the locus (used to screen out low heterozygosity loci that have 
#                a different distribution) 

#FstNoCorr, T1NoCorr, and T2NoCorr can be calculated from function given below: 
#WC_FST_FiniteSample_Haploids_2AllelesB_MCW  for the haploid case or
#WC_FST_FiniteSample_Diploids_2Alleles_NoCorr for diploids.



#The procedure will return a list with the following elements:
#     $FSTbar:  the mean Fst of the data from loci with high enough heterozygosity
#     $dfInferred:     the effective number of populations in the data (equals df + 1)
#     $numberLowFstOutliers:  the number of loci flagged by the OutFLANK procedure as having significantly low Fst (i.e. with a q-value less than 0.05)
#     $numberHighFstOutliers:  the number of loci flagged by the OutFLANK procedure as having significantly high Fst (i.e. with a q-value less than 0.05)
#     $results:   a data frame with information about each locus


# This results dataframe includes all of the input data, plus the following three columns:
#     $indexOrder: integer index giving the original order of rows in the input file
#     $OutlierFlag: TRUE if locus is an outlier; FALSE otherwise 
#     $qvalues: q-value for locus against null hypothesis of neutrality


#############LOAD NECESSARY PACAKAGES############# 

#Download the biocLite package at first use. On subsequent uses, run library(qvalue) before 
#using functions in the rest of this file.

#source("http://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)
source()

#############FUNCTIONS###################################
 #'
 #'Takes Fst data for a list of loci to find outliers, using a trimmed likelihood approach.
 #'
 #'This function should take in a dataframe ("FstDataFrame") that 
#'has columns for $LocusName,$Fst,$T1,$T2,$FstNoCorr, $T1NoCorr, $T2NoCorr,$H. It should return a dataframe 
#'with those same columns but also new columns for $LowOutlierFlag, $HighOutlierFlag,and $q.
#'
#'This function requires Fst's calculated without sample size correction. These
#'can be calculated, for example, with WC_FST_FiniteSample_Haploids_2AllelesB_NoSamplingCorrection in this package.
#'
#'This use of the biased FSTs is necessary for the trimming outlier approach 
#'with small samples, because the debiasing sometimes creates negtive Fsts 
#'which do not fit into the chi-square distribution.

#'This will use FST's calculated without sample size correction for outlier tests.
#'Such FSTs will be biased upwards, but as long as the sample size is similar for
#'all loci, the resulting measures ought to be give similar results.

#'This use of the biased FSTs is necessary for the trimming outlier approach with
#'small samples, because the debiasing sometimes creates negtive Fsts which do
#'not fit into the chi-square distribution.
#'
#'@title Fst outliers with trimming
#'
#'@param FstDataFrame A data frame that includes a row for each locus, with columns as follows: 
#'                    $LocusName: a character string that uniquely names each locus. 
#'                    $FST: Fst calculated for this locus. (Kept here to report the unbased Fst of the results) 
#'                    $T1: The numerator of the estimator for Fst (necessary, with $T2, to calculate mean Fst) 
#'                    $T2: The denominator of the estimator of Fst 
#'                    $FSTNoCorr: Fst calculated for this locus without sample
#'                    size correction. (Used to find outliers) 
#'                    $T1NoCorr: The numerator of the estimator for Fst without sample size correction (necessary, with $T2, to 
#'                    calculate mean Fst) 
#'                    $T2NoCorr: The denominator of the estimator of Fst 
#'                    without sample size correction 
#'                    $He: The heterozygosity of the locus (used to screen out low heterozygosity loci that have a different distribution) 
#'                    $indexOrder: integer index giving the original order of rows in the input file.
#'                    
#' @param LeftTrimFraction The proportion of loci that are trimmed from the lower end of the range of Fst before the likelihood funciton is applied.
#' 
#' @param RightTrimFraction The proportion of loci that are trimmed from the upper end of the range of Fst before the likelihood funciton is applied.
#' 
#' @param Hmin The minimum heterozygosity required before including calculations from a locus.
#' 
#' @param NumberOfSamples The number of spatial locations included in the data set.
#' 
#' @param qthreshold The desired false discovery rate threshold for calculating q-values.
#' 
#' @return
#' 
#' The function returns a list with seven elements:
#' \itemize{
#'  \item   FSTbar: the mean FST inferred from loci not marked as outliers 
#'  \item 	FSTNoCorrbar: the mean FST (not corrected for sample size -gives an upwardly biased estimate of FST)
#'  \item 	dfInferred: the inferred number of degrees of freedom for the chi-square distribution of neutral FST
#'   \item  numberLowFstOutliers: Number of loci flagged as having a signficantly low FST (not reliable)
#'   \item  numberHighFstOutliers: Number of loci identified as haivng significantly high FST
#'   \item  results:  a data frame with rows for each locus. This data frame includes all the original columns in the 
#'                    data set, and four new ones: 
#'                    \itemize{
#'              \item $indexOrder (the original order of the input data set),
#'              \item $GoodH (Boolean variable which is TRUE if the expected heterozygosity is greater than the Hemin set by input),
#'              \item $OutlierFlag (TRUE if the method identifies the locus as an outlier, FALSE otherwise), and 
#'              \item $q (the q-value for the test of neutrality for the locus)
#'              }
#'  }
#'  @export


OutFLANK=function(FstDataFrame, LeftTrimFraction=0.05, RightTrimFraction=0.05, Hmin=0.1, NumberOfSamples, qthreshold=0.05){
  
  #
  #
  #Setting up necessary columns in dataframe
  Fstdata= outputDFStarterNoCorr(FstDataFrame,Hmin)

  
  #making working dataframe with real Fst (no NAs), storing NAs to add back later
  workingDataFrame=Fstdata[which(!is.na(Fstdata$FSTNoCorr)),]
  storedDataFrameNA=Fstdata[which(is.na(Fstdata$FSTNoCorr)),]
  
  
  #Finding upper and lower bounds for trimming (eliminating NAs, but not negative FSTs)
  sortedDataFrame=workingDataFrame[order(workingDataFrame$FSTNoCorr),]
  
  NLociTotal=length(sortedDataFrame$FSTNoCorr)
  SmallestKeeper=ceiling(NLociTotal*LeftTrimFraction)
  LargestKeeper=floor(NLociTotal*(1-RightTrimFraction))
  LowTrimPoint=sortedDataFrame$FSTNoCorr[[SmallestKeeper]]
  HighTrimPoint=sortedDataFrame$FSTNoCorr[[LargestKeeper]]
    
  
  if(LowTrimPoint<0) {writeLines("ERROR: The smallest FST in the trimmed set must be > 0. Please use a larger LeftTrimFraction."); return()}
  if(HighTrimPoint>=1) {writeLines("ERROR: The largest FST in the trimmed set must be < 1. Please use a larger RightTrimFraction."); return()}
  
  #finding dfInferred and Fstbar iteratively  
  putativeNeutralListTemp=ifelse(workingDataFrame$FSTNoCorr>0,TRUE,FALSE)
  
  oldOutlierFlag=rep(FALSE,NLociTotal)
  
  
  #Note: All negative FST loci are maked as putative outliers, which will need
  #to be tested with the coalescent model later. In the meantime, they are
  #removed so as to not confuse the likelihood function.
  
  keepGoing=TRUE
  count = 0
  #writeLines(paste(mean(workingDataFrame$FSTNoCorr[putativeNeutralListTemp])))
  
  while(keepGoing){
    count=count+1
    if(count>19) {
      keepGoing=FALSE  
      writeLines("Exceeded iteration maximum.") ###Try with increased maximum value for count two lines above.
    }
    
    FstbarNoCorrTemp=fstBarCalculatorNoCorr(workingDataFrame[putativeNeutralListTemp,])  

    dfInferredTemp=EffectiveNumberSamplesMLE(workingDataFrame$FSTNoCorr[putativeNeutralListTemp],FstbarNoCorrTemp,NumberOfSamples,LowTrimPoint,HighTrimPoint)
    workingDataFrame=pOutlierFinderChiSqNoCorr(workingDataFrame,FstbarNoCorrTemp,dfInferredTemp,qthreshold)

    #### mark all negative FSTs as outliers if lowest nonneg FST is outlier
    #### (because negative Fst estimates can't be evaluated through the
    #### chi-square approach on their own)
    if(any(workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr<LowTrimPoint])) workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr<0]=TRUE
    
    ####Any loci previously marked as $OutlierFlag=TRUE remain so, even if the new iteration doesn''t flag them as outliers
    #     workingDataFrame$OutlierFlag=!as.logical((!workingDataFrame$OutlierFlag)*(!oldOutlierFlag))
    
    #Resetting neutral list, and checking whether the outlier list has stabilized
    putativeNeutralListTemp=ifelse((!workingDataFrame$OutlierFlag),TRUE,FALSE)
    if(sum(putativeNeutralListTemp)==0) {writeLines("No loci in neutral list..."); return("FAIL")}
    
    if(identical(oldOutlierFlag,workingDataFrame$OutlierFlag)) keepGoing=FALSE
    
    ######if all in trimmed get IDed as outlier - return to user with warning
    if(all(workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr<LowTrimPoint])){
      writeLines("All loci with Fst below the lower (lefthand) trim point were marked as outliers. Re-run with larger LeftTrimFraction or smaller qthreshold.")
      return(0)
    }
    
    if(all(workingDataFrame$OutlierFlag[workingDataFrame$FSTNoCorr>HighTrimPoint])){
      writeLines("All loci with Fst above the upper (righthand) trim point were marked as outliers. Re-run with smaller RightTrimFraction or smaller qthreshold.")
      return(0)
    }
    
    oldOutlierFlag=workingDataFrame$OutlierFlag
    
    #writeLines(paste(as.character(count),"   ",as.character(sum(putativeNeutralListTemp))))
  }
  
  if(count>19) writeLines("Loop iteration limit exceeded.")
  
  numberLowFstOutliers=sum(workingDataFrame$OutlierFlag[(workingDataFrame$FSTNoCorr<LowTrimPoint)])
  numberHighFstOutliers=sum(workingDataFrame$OutlierFlag[(workingDataFrame$FSTNoCorr>HighTrimPoint)])
  
  FSTbar=fstBarCalculator(workingDataFrame[putativeNeutralListTemp,])  
  
  
  #merge NA list back to working list, and sort by original order	
  resultsDataFrame=rbind(workingDataFrame,storedDataFrameNA)
  resultsDataFrame=resultsDataFrame[order(resultsDataFrame$indexOrder),]
  #return new dataframe
  list(FSTbar=FSTbar,FSTNoCorrbar=FstbarNoCorrTemp,dfInferred=dfInferredTemp,numberLowFstOutliers=numberLowFstOutliers,numberHighFstOutliers=numberHighFstOutliers,results=resultsDataFrame)
}


outputDFStarterNoCorr=function(FstDataFrame,Hmin=0.1) {
  #This will take a given dataframe with $LocusName, $FST,$He, $T1,  $T2, etc. and 
  #    initialize $indexOrder,$GoodH,$OutlierFlag (to 0), and $q (to 1).
  
  #Hmin is the smallest allowable He for which a locus should be included in 
  #the initial calculations. By default this requires that a locus have 
  #heterozygosity equal to 10% or more.
  
  len=length(FstDataFrame$FSTNoCorr)
  indexOrder=seq(1,len)
  GoodH=ifelse(FstDataFrame$He<Hmin,"lowH","goodH")
  OutlierFlag=ifelse(is.na(FstDataFrame$FSTNoCorr),NA,FALSE)
  qvalues=rep(NA,len)
  cbind(FstDataFrame, indexOrder, GoodH, qvalues,OutlierFlag )
  
}

EffectiveNumberSamplesMLE=function(FstVect, Fstbar, NumberOfSamples, SmallestFstInTrimmedList, LargestFstInTrimmedList){
  #This function should find the maximum likelihood value 
  #of the effective number of samples, for a given list of
  #Fst values.
  
  #The FstVect should already have been purged of NaN values and of loci with 
  #too low heterozygosity or MAF. 
  
  sortedFst=FstVect[order(FstVect)]	
  
  #The Minimum Fst considered in the trimmed data is the larger of the amount
  #specified by the user or the mean FSt over 100. This is to prevent extremely
  #small Fsts from causing estiamtion errors (Especially when R interprets a
  #small Fst as FSt=0.)
  LowTrimPoint=max(Fstbar/100,SmallestFstInTrimmedList)
  
  trimmedFstVect =FstVect[which((FstVect>=LowTrimPoint)&(FstVect<=LargestFstInTrimmedList))]
  
  trimmedFstArray=as.array(trimmedFstVect)
  
  localNLLAllData=function(dfInferred){
    localNLLOneLocus=function(Fst){
      negLLdfFstTrim(Fst,dfInferred,Fstbar,LowTrimPoint,LargestFstInTrimmedList)
    }
    sum(localNLLOneLocus(trimmedFstVect))
  }
  
  optim(NumberOfSamples, localNLLAllData, lower=2, method="L-BFGS-B")$par
}





pOutlierFinderChiSqNoCorr=function(DataList, Fstbar, dfInferred, qthreshold=0.05){
  #Finds outliers based on chi-squared distribution
  #Takes given values of dfInferred and Fstbar, and returns a list of q-values for all loci based on chi-square.
  #Assumes that the DataList input has a column called $FST and tat other columns are set up via the
  # "outputDFStarter" function.
  
  #Divide DataList into 3 lists:  DataListGood has $FST>0; DataListNeg has cases where $FST <=0; and
  #   DataListNA has cases where $FST is NA.
  #DataListNeg is necessary to keep separate here because these cases do not have meaningful  results with the chi-square aprach;
  #   however, they do carry information.
  
  DataListGood=DataList[which(DataList$FSTNoCorr>0),]
  DataListNonPosFst=DataList[which(DataList$FSTNoCorr<=0),]
  DataListNA=DataList[which(is.na(DataList$FSTNoCorr)),]
  
  
  
  pList=pTwoSidedFromChiSq(DataListGood$FSTNoCorr*(dfInferred)/Fstbar,dfInferred)
  
  qtemp=qvalue(pList,fdr.level=qthreshold,pi0.method="bootstrap")
  #Note:  Using the bootstrap method here seems OK, but if this causes problems remove the pi0.method="bootstrap" in the previous line to revert to the default.
  
  DataListGood$qvalues=qtemp$qvalues
  DataListGood$OutlierFlag=qtemp$significant
  rbind(DataListGood,DataListNonPosFst,DataListNA) 
  
}

pTwoSidedFromChiSq=function(x,df){
  #Takes a value x, finds the two-sided p-value for comparison to a chi-square distribution with df degrees of freedom.
  pOneSided=pchisq(x,df)
  ifelse(pOneSided>.5,(1-pOneSided)*2,pOneSided*2)
}




IncompleteGammaFunction=function(a, z) {
  #equivalence to Mathematica Gamma[a,z] according to 
  #   http://r.789695.n4.nabble.com/Incomplete-Gamma-function-td833545.html
  pgamma(z,a,lower=FALSE)*gamma(a)
}

negLLdfFstTrim=function(Fst, dfInferred, Fstbar, LowTrimPoint, HighTrimPoint){
  #Fst is the Fst from a locus, and dfInferred is the candidate value for the
  #degrees of freedom for the chi-squared distribution of neutral Fst, and
  #Fstbar is the mean Fst of all neutral loci (sequentially inferred from
  #non-outlier loci) LowTrimPoint and HighTrimPoint are the values of the lowest
  #and highest Fst values allowed to be included in the Fst list.
  #
  #Finds contribution to the negative log likelihood of a given locus' Fst for a
  #given dfInferred #CHECKED AGAINST MATHEMATICA DERIVATION##
  
  df=dfInferred
  
  1/(2*Fstbar)*(df * Fst +df * Fstbar * log(2) - df * Fstbar *log(df)-(df-2)*Fstbar * log(Fst)+df * Fstbar * log(Fstbar) + 2*Fstbar * log(-IncompleteGammaFunction(df/2,df*HighTrimPoint/(2*Fstbar))+IncompleteGammaFunction(df/2,df*LowTrimPoint/(2*Fstbar))))
}


