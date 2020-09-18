#include <Rcpp.h> 
#include   <stdlib.h>
#include <math.h>  
using namespace Rcpp; 
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <iterator>
#include<vector>


// [[Rcpp::export]]

NumericVector  UXKindex(std::vector< std::string > seq1){
  
  int m=seq1.size();
  IntegerVector  Kindex(m,0);
  
  const char charC1[]="K";
  const char charC2[]="U";
  const char charC3[]="X";
  
  for(int i=0;i<m;i++)
  {
    std::string pepsequen1=seq1[i];
  
    
    for (std::string::iterator it=pepsequen1.begin(); it!=pepsequen1.end()-1; ++it)
    {
     if(*it==*charC1){
       Kindex[i]=1;
      }
    }
    
       for (std::string::iterator it=pepsequen1.begin(); it!=pepsequen1.end(); ++it)
      {
        if((*it==*charC2)|(*it==*charC3)){
          Kindex[i]=0;
         }
      }
   }
  return wrap(Kindex);
}


/*** R

setwd("D:\\database")
library(data.table)
library(Rcpp)
library(inline)
library(ggplot2)
library(e1071)
  
peptideID<-fread("uniprot-proteomehuman_digested_Mass600to4000normal.txt")

indexUXK<-UXKindex(peptideID$Sequence)
index<-which(indexUXK==1)
peptideID1<-peptideID[index,]
write.csv(peptideID1,file="peptidesequence.csv")
*/
