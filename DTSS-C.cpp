#include <Rcpp.h> 
#include   <stdlib.h>
#include <math.h>  
using namespace Rcpp; 
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <map>
#include <algorithm>
#include <iterator>
#include<vector>


// [[Rcpp::export]]

std::vector<int> getcount(std::vector<int> x){
 
  std::vector<int> w;
  int count=1;
  for(int i=1;i<x.size();i++){
    if(x[i]==x[i-1]){
      count++;
      }
    else if(x[i]!=x[i-1]){
      if(count>(6))
        {
        w.push_back(x[i-1]);
        } 
      count=1;
      }
  }
  
  return w;
}





// [[Rcpp::export]]

std::vector<int>  getindex(std::vector<int > x,std::vector<int > y,std::vector<int > z){
int i=0;
int j=0;
int m=x.size();
int n=y.size();
std::vector<int> w;
w.push_back(0);
while(i<m && j<n)
{
  
  if (x[i]<y[j])
  {
    i=i+1;
    }
  else if(x[i]==y[j]){
    w.push_back(z[j]);
    j++;
    } 
  else if(x[i]>y[j]) j=j+1;
}

std::sort(w.begin(),w.end());

return w;
}


// [[Rcpp::export]]

int getcharindex(std::string x,std::vector< std::string > seqy){
  
  int m=seqy.size();
  std::vector<int> w;
  int t=0;

  for(int i=0;i<m;i++){
    if(x==seqy[i]) {
      t=i;
      break;}
    else continue;
  }
    
   return t;
}



// [[Rcpp::export]]

int getCIDpos(IntegerVector x, int t,int pos){
  
  
  int m=x.size();
  int leng=0;
  
  for(int i=pos;i<m;i++)
  {
    if ( x[i]<t) 
    {
      continue;
    }
    
    else if ( x[i]==t) 
    {
      leng++;
    }
    else
      break;
  }
  
  return leng;
  
}

// [[Rcpp::export]]

IntegerVector getmgfindex(IntegerVector scanindex,int t){
  
  IntegerVector z;
  int m=scanindex.size();
  
  
  for(int i=0;i<m;i++)
  {
    
    int scani=scanindex[i];
    
    
    if ( scani<(t-40)) 
    {
      continue;
    }
    
    else if (( scani>=(t-40)) &&( scani<(t+40))) 
    {
      z.push_back(i);
    }
    else if(scani>=40) 
      break;
  }
  
  return wrap(z);
  
}



// [[Rcpp::export]]

NumericVector getspectra(IntegerVector x, NumericVector y,int t){
  
  NumericVector z;
  int m=x.size();
  
  
  for(int i=0;i<m;i++)
  {
    if ( x[i]<t) 
    {
      continue;
    }
    
    else if ( x[i]==t) 
    {
      z.push_back(y[i]);
    }
    else
      break;
  }
  
  return wrap(z);
  
}

// [[Rcpp::export]]

NumericVector posgetspectra(IntegerVector x, NumericVector y,int t,int pos){
  
  NumericVector z;
  int m=x.size();
  int tempos=std::max(0,pos-2);
  
  for(int i=tempos;i<m;i++)
  {
    if ( x[i]<t) 
    {
      continue;
    }
    
    else if ( x[i]==t) 
    {
      z.push_back(y[i]);
    }
    else
      break;
  }
  
  return wrap(z);
  
}


// [[Rcpp::export]]

DataFrame getmzseqz(std::vector<int> scoremgfindex,std::vector<int> scorepepseqindex,std::vector<double> pepexpmz1,
                    std::vector<int> pepexpz1,std::vector<int> mgfpos1,std::vector<int> mgfleng1,std::vector<std::string> phossequence,
                    std::vector<int> Theomzpos1,std::vector<int> Theomzleng1){
  
  std::vector<double> expmz1;
  std::vector<int> pepseqindex;
  std::vector<int> expz1;
  std::vector<int> pos1;
  std::vector<int> leng1;
  std::vector<std::string> sequence1;
  std::vector<int> Tmzpos1;
  std::vector<int> Tmzleng1;
  
  for(int i=0;i<scoremgfindex.size()-1;i++)
  {
    int tempmgfindex=scoremgfindex[i]-1;
    int tempseqindex=scorepepseqindex[i]-1;
    
    double tempmz=pepexpmz1[tempmgfindex];
    double tempz=pepexpz1[tempmgfindex];
    int tempos=mgfpos1[tempmgfindex];
    int temleng=mgfleng1[tempmgfindex];
    std::string tempseq=phossequence[tempseqindex];
    
    expmz1.push_back(tempmz);
    expz1.push_back(tempz);
    pos1.push_back(tempos);
    leng1.push_back(temleng);
    sequence1.push_back(tempseq);
    pepseqindex.push_back(tempseqindex+1);
    
    Tmzpos1.push_back(Theomzpos1[tempseqindex]);
    Tmzleng1.push_back(Theomzleng1[tempseqindex]);
  }
  
  return Rcpp::DataFrame::create(Named("pepexpmz")=expmz1,Named("pepexpz")=expz1,Named("mgfpos")=pos1,Named("mgfleng")=leng1,
                                       Named("pepseqindex")=pepseqindex,Named("peptideseq1")=sequence1,Named("Theomzpos")=Tmzpos1,Named("Theomzleng")=Tmzleng1); 
}


// [[Rcpp::export]]

NumericVector deisotope(NumericVector x, NumericVector y,NumericVector z){
  int m=x.size();
  NumericVector iso(m);
  std::fill(iso.begin(), iso.end(), 0);
  
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>350)&&((fabs(y[j]-y[i]-1.0034) <0.01) |(fabs(y[j]-y[i]-0.50017)<0.01)|(fabs(y[j]-y[i]-0.3344)<0.01)|(fabs(y[j]-y[i]-0.2509)<0.006)|(fabs(y[j]-y[i]-0.20068)<0.006))&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
        iso[j]=1;  
    }
  }
  
  return wrap(iso); 
  
}



// [[Rcpp::export]]

NumericVector calz(NumericVector x, NumericVector y,NumericVector z){
  
  
  int m=x.size();
  NumericVector iso(m);
  std::fill(iso.begin(), iso.end(), 1);
  
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>250)&&((fabs(y[j]-y[i]-0.50017)<0.01))&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
        iso[i]=2;  
    }
  }
  
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>250)&&((fabs(y[j]-y[i]-0.3344) <0.01) )&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
        iso[i]=3;  
    }
  }
  
  for(int i=0;i<m;i++)
  {
    int minsize=std::min(i+5,m-1);
    for(int j=i+1;j<minsize;j++)
    {
      if((y[i]>250)&&((fabs(y[j]-y[i]-0.25085) <0.01) )&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
        iso[i]=4;  
    }
  }
    for(int i=0;i<m;i++)
    {
      int minsize=std::min(i+5,m-1);
      for(int j=i+1;j<minsize;j++)
      {
        if((y[i]>250)&&((fabs(y[j]-y[i]-0.20068) <0.01) )&&(x[i]==x[j])&&(z[i]>2000)&&(z[j]/z[i]>y[i]/8000)&&(z[j]/z[i]<y[i]/300))
          iso[i]=5;  
      }
  }
  
  return wrap(iso);
  
}





// [[Rcpp::export]]

NumericVector  calmass(std::vector< std::string > seq1){
  
  int m=seq1.size();
  std::vector<double> pepmass;
  
  std::map<char,double> AAmass;
  AAmass.insert ( std::pair<char,double>('A',71.03711) );
  AAmass.insert ( std::pair<char,double>('R',156.10111) );
  AAmass.insert ( std::pair<char,double>('N',114.04293) );
  AAmass.insert ( std::pair<char,double>('D',115.02694) );
  AAmass.insert ( std::pair<char,double>('C',160.0338) );
  AAmass.insert ( std::pair<char,double>('E',129.04259) );
  AAmass.insert ( std::pair<char,double>('Q',128.05858) );
  AAmass.insert ( std::pair<char,double>('G',57.02146) );
  AAmass.insert ( std::pair<char,double>('H',137.05891) );
  AAmass.insert ( std::pair<char,double>('I',113.08406) );
  AAmass.insert ( std::pair<char,double>('L',113.08406) );
  AAmass.insert ( std::pair<char,double>('K',128.09496) );
  AAmass.insert ( std::pair<char,double>('M',131.04049) );
  AAmass.insert ( std::pair<char,double>('F',147.06841) );
  AAmass.insert ( std::pair<char,double>('P',97.05276) );
  AAmass.insert ( std::pair<char,double>('S',87.03203) );
  AAmass.insert ( std::pair<char,double>('T',101.04768) );
  AAmass.insert ( std::pair<char,double>('W',186.07931) );
  AAmass.insert ( std::pair<char,double>('Y',163.06333) );
  AAmass.insert ( std::pair<char,double>('V',99.06841) );
  
  
  for(int i=0;i<m;i++)
  {
    std::string pepsequen1=seq1[i];
    double sumass=0;
    
    for (std::string::iterator it=pepsequen1.begin(); it!=pepsequen1.end(); ++it)
    {
      double residulmass= AAmass.find(*it)->second;
      sumass=sumass+residulmass;
    }
    
    sumass=sumass+1.0078+18.0106;
    pepmass.push_back(sumass);
    
    
  }
  
  return wrap(pepmass);
}

// [[Rcpp::export]]

int mass300pos(std::vector<double> x,int cutmass){
  int m=x.size();
  int t=0;
  for(int i=0;i<m;i++){if(x[i]<cutmass) t++;
  else break;}
  
  return(t);
}


// [[Rcpp::export]]

DataFrame  spectragen(std::vector< std::string > seq1){
  
  
  int m=seq1.size();
  std::vector<char> Residue;
  std::vector<int> Index;
  std::vector<std::string> Type;
  std::vector<std::string> Linkseq;
  std::vector< double > mz;
  std::vector< double > pepmass;
  
  std::map<char,double> AAmass;
  AAmass.insert ( std::pair<char,double>('A',71.03711));
  AAmass.insert ( std::pair<char,double>('R',156.10111));
  AAmass.insert ( std::pair<char,double>('N',114.04293));
  AAmass.insert ( std::pair<char,double>('D',115.02694));
  AAmass.insert ( std::pair<char,double>('C',160.0338));
  AAmass.insert ( std::pair<char,double>('E',129.04259));
  AAmass.insert ( std::pair<char,double>('Q',128.05858));
  AAmass.insert ( std::pair<char,double>('G',57.02146));
  AAmass.insert ( std::pair<char,double>('H',137.05891));
  AAmass.insert ( std::pair<char,double>('I',113.08406));
  AAmass.insert ( std::pair<char,double>('L',113.08406));
  AAmass.insert ( std::pair<char,double>('K',128.09496));
  AAmass.insert ( std::pair<char,double>('M',131.04049));
  AAmass.insert ( std::pair<char,double>('F',147.06841));
  AAmass.insert ( std::pair<char,double>('P',97.05276));
  AAmass.insert ( std::pair<char,double>('S',87.03203));
  AAmass.insert ( std::pair<char,double>('T',101.04768));
  AAmass.insert ( std::pair<char,double>('W',186.07931));
  AAmass.insert ( std::pair<char,double>('Y',163.06333));
  AAmass.insert ( std::pair<char,double>('V',99.06841));
  
  for(int i=0;i<m;i++)
  {
  
    std::string pepsequen1=seq1[i];
  
    double sumass=1.0078;
   
    for (std::string::iterator it=pepsequen1.begin(); it!=pepsequen1.end()-1; ++it)
    {
      
      double residulmass= AAmass.find(*it)->second;
      sumass=sumass+residulmass;
      Residue.push_back(*it);
      Index.push_back(i+1);
      Type.push_back("b");
      Linkseq.push_back(pepsequen1);
      mz.push_back(sumass);
    }
    
    
    sumass=19.0184;
    for (std::string::iterator it=pepsequen1.end()-1; it!=pepsequen1.begin(); --it)
    {
      
      double residulmass= AAmass.find(*it)->second;
      sumass=sumass+residulmass;
      
      Residue.push_back(*it);
      Index.push_back(i+1);
      Type.push_back("y");
      Linkseq.push_back(pepsequen1);
      mz.push_back(sumass);
    }
  }
  
  
  
  return Rcpp::DataFrame::create(Named("index")=Index,Named("type")=Type,Named("Residue")=Residue,
                                       Named("seq")=Linkseq,Named("mz")=mz);
}

// [[Rcpp::export]]

NumericVector  singlespectragenmod(std::string  seq1,double pmass){
  
  
  std::vector< double > mz;
  std::map<char,double> AAmass;
  AAmass.insert ( std::pair<char,double>('A',71.03711) );
  AAmass.insert ( std::pair<char,double>('R',156.10111) );
  AAmass.insert ( std::pair<char,double>('N',114.04293) );
  AAmass.insert ( std::pair<char,double>('D',115.02694) );
  AAmass.insert ( std::pair<char,double>('C',160.0338) );
  AAmass.insert ( std::pair<char,double>('E',129.04259) );
  AAmass.insert ( std::pair<char,double>('Q',128.05858) );
  AAmass.insert ( std::pair<char,double>('G',57.02146) );
  AAmass.insert ( std::pair<char,double>('H',137.05891) );
  AAmass.insert ( std::pair<char,double>('I',113.08406) );
  AAmass.insert ( std::pair<char,double>('L',113.08406) );
  AAmass.insert ( std::pair<char,double>('K',128.09496) );
  AAmass.insert ( std::pair<char,double>('M',131.04049) );
  AAmass.insert ( std::pair<char,double>('F',147.06841) );
  AAmass.insert ( std::pair<char,double>('P',97.05276) );
  AAmass.insert ( std::pair<char,double>('S',87.03203) );
  AAmass.insert ( std::pair<char,double>('T',101.04768) );
  AAmass.insert ( std::pair<char,double>('W',186.07931) );
  AAmass.insert ( std::pair<char,double>('Y',163.06333) );
  AAmass.insert ( std::pair<char,double>('V',99.06841) );
  
  const char charC1[]="K";
  const char charC2[]="K";
  
  std::string pepsequen1=seq1;
  double sumass=0;
  IntegerVector Cposb;
  IntegerVector Cposy;
  int j=1;
  
  std::string::iterator itstart=pepsequen1.begin();
  std::string::iterator itend=pepsequen1.end()-1;
  
  for (std::string::iterator it=pepsequen1.begin(); it!=pepsequen1.end()-1; ++it)
  {
    
    if((*it==*charC1)|(*it==*charC2)){Cposb.push_back(j);}
    j++;
  }  
  
  
  j=1;
  
  for (std::string::iterator it=pepsequen1.end()-1; it!=pepsequen1.begin()+1; --it)
  {
    
    if((*it==*charC1)|(*it==*charC2)){Cposy.push_back(j);}
    j++;
  }
  
  
  sumass=1.0078;
  j=1;
  
  for (std::string::iterator it=pepsequen1.begin(); it!=pepsequen1.end()-1; ++it)
  {
    
    double residulmass= AAmass.find(*it)->second;
    sumass=sumass+residulmass;
    
    if(j==*Cposb.begin()) {
      sumass=sumass+pmass;
      mz.push_back(sumass);}
    else if(j>*Cposb.begin()){
      mz.push_back(sumass);
    }
    j++;
  }
  
  
  sumass=19.0184;
  j=1;
  for (std::string::iterator it=pepsequen1.end()-1; it!=pepsequen1.begin(); --it)
  {
    double residulmass= AAmass.find(*it)->second;
    sumass=sumass+residulmass;
    if(j==*Cposy.begin()) {
      sumass=sumass+pmass;
      mz.push_back(sumass);}
    else if(j>*Cposy.begin()){
      mz.push_back(sumass);
    }
    j++;
  }
  std::sort(mz.begin(),mz.end());  
  return wrap(mz);
}



// [[Rcpp::export]]

NumericVector  singlespectragen(std::string pepsequen1){
  
  std::vector< double > mz;
  std::map<char,double> AAmass;
  AAmass.insert ( std::pair<char,double>('A',71.03711) );
  AAmass.insert ( std::pair<char,double>('R',156.10111) );
  AAmass.insert ( std::pair<char,double>('N',114.04293) );
  AAmass.insert ( std::pair<char,double>('D',115.02694) );
  AAmass.insert ( std::pair<char,double>('C',160.0338) );
  AAmass.insert ( std::pair<char,double>('E',129.04259) );
  AAmass.insert ( std::pair<char,double>('Q',128.05858) );
  AAmass.insert ( std::pair<char,double>('G',57.02146) );
  AAmass.insert ( std::pair<char,double>('H',137.05891) );
  AAmass.insert ( std::pair<char,double>('I',113.08406) );
  AAmass.insert ( std::pair<char,double>('L',113.08406) );
  AAmass.insert ( std::pair<char,double>('K',128.09496) );
  AAmass.insert ( std::pair<char,double>('M',131.04049) );
  AAmass.insert ( std::pair<char,double>('F',147.06841) );
  AAmass.insert ( std::pair<char,double>('P',97.05276) );
  AAmass.insert ( std::pair<char,double>('S',87.03203) );
  AAmass.insert ( std::pair<char,double>('T',101.04768) );
  AAmass.insert ( std::pair<char,double>('W',186.07931) );
  AAmass.insert ( std::pair<char,double>('Y',163.06333) );
  AAmass.insert ( std::pair<char,double>('V',99.06841) );
  
  
  double sumass=1.0078;
  for (std::string::iterator it=pepsequen1.begin(); it!=pepsequen1.end()-1; ++it)
  {
    double residulmass= AAmass.find(*it)->second;
    sumass=sumass+residulmass;
    mz.push_back(sumass);
  }
  
  
  sumass=19.0184;
  for (std::string::iterator it=pepsequen1.end()-1; it!=pepsequen1.begin(); --it)
  {
    double residulmass= AAmass.find(*it)->second;
    sumass=sumass+residulmass;
    mz.push_back(sumass);
  }
  
  std::sort(mz.begin(),mz.end());
  
  int m=mz.size();
  int n=std::min(6,m);
  
  
  for(std::vector<double>::iterator it = mz.begin();it!= mz.begin()+n;it++){
    if((fabs(*it-147.1128)<0.001)|(fabs(*it-175.1190)<0.001)){
      mz.erase(it);}
  }
  
  
  
  return wrap(mz);
  
}


// [[Rcpp::export]]

int spectracomp(NumericVector x, NumericVector y,double delta,double eppm){
  
  int m = x.size();
  int n = y.size(); 
  int i=0;
  int j=0;
  std::vector<double> deltamass;
  
  while(i<m && j<n)
  {
    
    if ((x[i]-y[j])<delta&&(x[i]-y[j])>(-delta))
    {
      double xtemp=x[i];
      double ytemp=y[j];
      double error=fabs(xtemp-ytemp);
      error=error/xtemp;
      deltamass.push_back(error);
      i=i+1;
      j=j+1;
    }
    else if(x[i]<=y[j]-delta) i=i+1;
    else if(x[i]>=y[j]+delta) j=j+1;
  }
  
  int t=0;
  int indexsize=deltamass.size();
  for(int w=0;w<indexsize;w++){
    
    if((deltamass[w]*1e6)<eppm){
      t++;}
    else continue;
  }
  return t;
}


// [[Rcpp::export]]

int spectracomp2(std::vector<double> x, NumericVector y,double delta,double eppm){
  
  int m = x.size();
  int n = y.size(); 
 
  int i=0;
  int j=0;
  std::vector<double> deltamass;

  
  while(i<m && j<n)
  {
    
    if ((x[i]-y[j])<delta&&(x[i]-y[j])>(-delta))
    {
      double xtemp=x[i];
      double ytemp=y[j];
      double error=fabs(xtemp-ytemp);
      error=error/xtemp;
      deltamass.push_back(error);
      i=i+1;
      j=j+1;
     }
    else if(x[i]<=y[j]-delta) i=i+1;
    else if(x[i]>=y[j]+delta) j=j+1;
  }
  
  int t=0;
  int indexsize=deltamass.size();
  for(int w=0;w<indexsize;w++){
      
      if((deltamass[w]*1e6)<eppm){
        t++;}
      else continue;
      }
 
  return t;
}


// [[Rcpp::export]]

DataFrame Indexsearch(IntegerVector mgfindex,IntegerVector mgfpos,NumericVector pepmz,IntegerVector mgfindex11,NumericVector mz1,
                          IntegerVector seqindex,NumericVector seqTheorymz,std::vector<int > seqindex11,std::vector<int > index11){
  
  
  int xx=mgfindex.size();
  
  int n = seqTheorymz.size();
  
  std::vector<int > Theoryindex;
  std::vector<int > mgfi;

 
     for(int x=0;x<xx;x++){
      
        NumericVector spectraexpHCD=posgetspectra(mgfindex11,mz1,mgfindex[x],mgfpos[x]);
        int m = spectraexpHCD.size();
        
        std::vector<int > Theoryseqindextemp;
        int i=0;
        int j=0;
        
        while(i<m&&j<n)
        {
          if ((spectraexpHCD[i]-seqTheorymz[j])<0.05&&(spectraexpHCD[i]-seqTheorymz[j])>(-0.05))
          {
            Theoryseqindextemp.push_back(seqindex[j]);
            j++;
           }
          else if(spectraexpHCD[i]<=seqTheorymz[j]-0.05) i=i+1;
          else if(spectraexpHCD[i]>=seqTheorymz[j]+0.05) j=j+1;
        }
        std::sort(Theoryseqindextemp.begin(),Theoryseqindextemp.end());
        std::vector<int> tempindex=getindex(Theoryseqindextemp,seqindex11,index11);
     
        if(tempindex.size()>3){
        std::vector<int> tempcountindex=getcount(tempindex);
        int tempsize=tempcountindex.size();
        if(tempsize>0){
          Theoryindex.insert(Theoryindex.end(), tempcountindex.begin(), tempcountindex.end());
           mgfi.insert(mgfi.end(), tempsize, mgfindex[x]);
          
        }}}
     
  return Rcpp::DataFrame::create(Named("mgfindexexp")=mgfi,Named("pepseqindex")=Theoryindex);
}

// [[Rcpp::export]]

DataFrame Openmodsearch(IntegerVector mgfindexexp,IntegerVector mgfpos,IntegerVector mgfleng,NumericVector pepmz,NumericVector pepmass,IntegerVector mgfindex11,NumericVector mz1,
                        std::vector< int > seqindex,std::vector< std::string > phossequence,NumericVector phosmass,NumericVector Theomz1,IntegerVector fmzpos,IntegerVector fmzleng){
  
  
  int m=pepmz.size();
  std::vector<int> pepseqindex;
  std::vector<int> mgfi;
  std::vector<int> HCDscore;
  std::vector<int> HCDscore2;
  std::vector<int> CIDscore;
  std::vector<int> sizeHCD;
  std::vector<int> sizeCID;
  std::vector<std::string> pepseq;
  std::vector<int> mgfpos1;
  std::vector<int> mgfleng1;
  std::vector<double> mz11(mz1.begin(),mz1.end());
  std::vector<int> mgfindex12(mgfindex11.begin(),mgfindex11.end());
  
  
  for(int i=0;i<m;i++){
    std::vector<double> spectraexpHCD(mz11.begin()+mgfpos[i],mz11.begin()+mgfpos[i]+mgfleng[i]);
    double phosmassi=phosmass[i];
    double modpepmassi=pepmass[i]-phosmassi;
   
   
   std::string phoseqi=phossequence[i];
   NumericVector phosspectra(Theomz1.begin()+fmzpos[i],Theomz1.begin()+fmzpos[i]+fmzleng[i]);;
   NumericVector CIDspectrai=singlespectragenmod(phoseqi,modpepmassi);
   int t=spectracomp2(spectraexpHCD,phosspectra,0.06,20)+spectracomp2(spectraexpHCD,CIDspectrai,0.06,20);
   if(t>5){
     t=t+spectracomp2(spectraexpHCD,phosspectra-18.01056,0.06,20)+spectracomp2(spectraexpHCD,phosspectra-17.02655,0.06,20);
     
      if(modpepmassi>800){
        mgfi.push_back(mgfindexexp[i]);
        pepseq.push_back(phossequence[i]);
        HCDscore.push_back(t);
      
        sizeHCD.push_back(spectraexpHCD.size());
        mgfpos1.push_back(mgfpos[i]);
        mgfleng1.push_back(mgfleng[i]);
        pepseqindex.push_back(seqindex[i]);
        }}}  
  
   return Rcpp::DataFrame::create(Named("mgfindexexp")=mgfi,Named("pepseqindex")=pepseqindex,Named("phossequence")=pepseq,Named("pep1HCDscore")=HCDscore,
                                Named("sizeHCD")=sizeHCD,Named("mgfpos")=mgfpos1,Named("mgfleng")=mgfleng1);
}


// [[Rcpp::export]]
DataFrame Research(NumericVector pepmass,IntegerVector pepindex,std::vector< std::string > peptide1seq,IntegerVector mgfindexexp,NumericVector pepexpmz,NumericVector Theomasspep1,IntegerVector scanindex,
                   IntegerVector scoreHCD,IntegerVector sizeHCD,NumericVector pepmz1,IntegerVector mgfpos,IntegerVector mgfindex11,NumericVector mz1,
                   NumericVector nonphosmass,std::vector< std::string > nonphossequence){
  
  
  std::vector<int> mgfindex;
  std::vector<int> mgfindex2;
  std::vector<int> phosrev;
  std::vector<int> score1HCD;
  std::vector<int> score2HCD;

  std::vector<int> score1CID;
  std::vector<int> score2CID;

  std::vector<double> deltamass;
  std::vector<double> pepmz;
  std::vector<int> scanindexexp;
  std::vector< std::string > pep1seq;
  std::vector< std::string > pep2seq;
  std::vector<int> sizeHCD1;
  std::vector<int> sizeCID1;
  std::vector<int> indicator;
  std::vector<int> scoreprecursor;
  std::vector<int> pepindex1;
  std::vector<int> pepindex2;
  
  std::vector<int> posmgf;
  
  
  int m=scanindex.size();
  int n=nonphossequence.size();

  for(int i=0;i<m;i++)
  {
  
  double p1mass=Theomasspep1[i]-1.0078;
  double pepmassi=pepmass[i]-1.0078;
  std::string phoseqi=peptide1seq[i];
  
  int x=mgfindexexp[i]-1;
  NumericVector expspectra=posgetspectra(mgfindex11,mz1,x+1,mgfpos[i]);

    
  for(int j=i+1;j<m;j++){
   if(mgfindexexp[i]==mgfindexexp[j]){
     std::string phoseqj=peptide1seq[j];
     double p2mass=Theomasspep1[j]-1.0078;
     double residualmass=fabs(p1mass+p2mass-pepmassi+209.97181);
     double abstdeltamass=fabs(p1mass+p2mass-pepmassi+209.97181);
     double tresmass=abstdeltamass-(int(abstdeltamass+0.5))*1.0033;
     tresmass=fabs(tresmass);
     double error=round(1e6*tresmass/(pepmassi+209));
      if((abstdeltamass<1.2)&&(error<10)){
          
          std::string phoseqj=peptide1seq[j];
          mgfindex.push_back(mgfindexexp[i]);
          deltamass.push_back(residualmass);
          pep1seq.push_back(phoseqi);
          pep2seq.push_back(phoseqj);
          scanindexexp.push_back(scanindex[i]);
          pepmz.push_back(pepexpmz[i]);
          score1HCD.push_back(scoreHCD[i]);
          score2HCD.push_back(scoreHCD[j]);
          sizeHCD1.push_back(sizeHCD[i]);
          indicator.push_back(1);
          pepindex1.push_back(pepindex[i]);
          pepindex2.push_back(pepindex[j]);
          posmgf.push_back(mgfpos[i]);
        }}
      else if(mgfindexexp[j]>mgfindexexp[i]){break;}
   }
  
  
  /*for(int j=0;j<n;j++){
  std::string nonphoseqj=nonphossequence[j];
  double p2mass=nonphosmass[j]-1.0078;
  double residualmass=fabs(p1mass+p2mass-pepmassi+209.97181);
  double abstdeltamass=fabs(p1mass+p2mass-pepmassi+209.97181);
  double tresmass=abstdeltamass-(int(abstdeltamass+0.5))*1.0033;
  tresmass=fabs(tresmass);
  double error=round(1e6*tresmass/(pepmassi+209));
  if((abstdeltamass<1.2)&&(error<10)){
   

  NumericVector CIDspectraj=singlespectragenmod(nonphoseqj,pepmassi-p2mass);
  
  NumericVector HCDspectraj=singlespectragen(nonphoseqj);
  

  int scoretemp2=spectracomp(expspectra,HCDspectraj,0.06,20)+spectracomp(expspectra,CIDspectraj,0.06,20);
  if(scoretemp2>5){
  scoretemp2=scoretemp2+spectracomp(expspectra,HCDspectraj-18.01056,0.06,20)+spectracomp(expspectra,HCDspectraj-17.02655,0.06,20);
  
  mgfindex.push_back(mgfindexexp[i]);
  deltamass.push_back(residualmass);
  pep1seq.push_back(phoseqi);
  pep2seq.push_back(nonphoseqj);
  scanindexexp.push_back(scanindex[i]);
  pepmz.push_back(pepexpmz[i]);
  score1HCD.push_back(scoreHCD[i]);
  score2HCD.push_back(scoretemp2);
  sizeHCD1.push_back(sizeHCD[i]);
  indicator.push_back(0);
  pepindex1.push_back(pepindex[i]);
  pepindex2.push_back(j+1);
  }}}*/
  
  }
  return Rcpp::DataFrame::create(Named("mgfindexexp")=mgfindex,Named("scanindex")=scanindexexp,Named("pepexpmz")=pepmz,
                                 Named("Alphapep")=pep1seq,Named("Betapep")=pep2seq,Named("Alphaindex")=pepindex1,Named("Betaindex")=pepindex2,
    Named("residualmass")=deltamass,Named("scoreAlpha")=score1HCD, Named("scoreBeta")=score2HCD,
          Named("expsizeHCD")=sizeHCD1,Named("mgfpos")=posmgf,Named("phosIndicator")=indicator);}



// [[Rcpp::export]]

DataFrame monoisosearch(IntegerVector mgfindex,NumericVector phosmass1,NumericVector phosmass2,IntegerVector mgfpos,IntegerVector mgfindex11,NumericVector mz1,
                        std::vector< std::string > phossequence1,std::vector< std::string > phossequence2){
  
  
  int m=mgfindex.size();
  std::vector<int> score1;
  std::vector<int> score2;
  
  for(int i=0;i<m;i++){
    NumericVector expspectra=posgetspectra(mgfindex11,mz1,mgfindex[i],mgfpos[i]);
    double modpepmass1=phosmass1[i]+209.97181-1.0078;
    double modpepmass2=phosmass2[i]+209.97181-1.0078;
    std::string phoseqi=phossequence1[i];
    NumericVector phosspectrai=singlespectragen(phoseqi);
    NumericVector CIDspectrai=singlespectragenmod(phoseqi,modpepmass2);
    
    
    std::string phoseqj=phossequence2[i];
    NumericVector phosspectraj=singlespectragen(phoseqj);
    NumericVector CIDspectraj=singlespectragenmod(phoseqj,modpepmass1);
    
    int t1=spectracomp(expspectra,phosspectrai,0.06,20)+spectracomp(expspectra,CIDspectrai,0.06,20);
    if(t1>5){
      t1=t1+int(0.5*(spectracomp(expspectra,phosspectrai-18.01056,0.06,10)+spectracomp(expspectra,phosspectrai-17.02655,0.06,10)));}
    
    int t2=spectracomp(expspectra,phosspectraj,0.06,20)+spectracomp(expspectra,CIDspectraj,0.06,20);
    if(t2>5){
      t2=t2+int(0.5*(spectracomp(expspectra,phosspectraj-18.01056,0.06,10)+spectracomp(expspectra,phosspectraj-17.02655,0.06,10)));}
    
    score1.push_back(t1);
    score2.push_back(t2);
  } 
  
  return Rcpp::DataFrame::create(Named("mscoreAlpha")= score1,Named("mscoreBeta")= score2);
}




// [[Rcpp::export]]

DataFrame getproteiname(std::vector< std::string > Proname,std::vector< int > alphaindex,std::vector< int > betaindex){
  
  

  std::vector< std::string > name1;
  std::vector< std::string > name2;
  int m=alphaindex.size();
  
  for(int i=0;i<m;i++){
    int t1=alphaindex[i];
    int t2=betaindex[i];
    name1.push_back(Proname[t1-1]);
    name2.push_back(Proname[t2-1]);
  } 
  
  return Rcpp::DataFrame::create(Named("Proteina")=name1,Named("Proteinb")=name2);
}

// [[Rcpp::export]]

DataFrame getproteiname2(std::vector< std::string > Proname,std::vector< std::string > Proname2,std::vector< int > alphaindex,std::vector< int > betaindex){
  
  
  std::vector< std::string > name1;
  std::vector< std::string > name2;
  int m=alphaindex.size();
  
  for(int i=0;i<m;i++){
    int t1=alphaindex[i];
    int t2=betaindex[i];
    name1.push_back(Proname[t1-1]);
    name2.push_back(Proname2[t2-1]);
  } 
  
  return Rcpp::DataFrame::create(Named("Proteina")=name1,Named("Proteinb")=name2);
}


/*** R

setwd("D:\\Rfiles\\Indexsearch")
library(data.table)
  library(Rcpp)
  library(inline)
  library(ggplot2)
  
  
mgfread<-function(df)
{
  sax<-df[,1:2]
  colnames(sax)<-c("mz","intensity")
  pepmassindex<-which(sax$mz=="TITLE")
  pepmass<-sax$intensity[pepmassindex+1]
  pepmass<-strsplit(pepmass," ")
  pepmass<-sapply(pepmass,function(x) x[1])
  charge<-sax$intensity[pepmassindex+2]
  charge<-strsplit(charge,NULL)
  charge<-sapply(charge,function(x) x[1])
  Rtime<-sax$intensity[pepmassindex+3]
  scanindex<-sax$intensity[pepmassindex+4]
  mgfindex=seq(1:length(pepmass))
  parention<-data.frame(mgfindex,pepmass,charge,Rtime,scanindex)
  remove_index<-c(pepmassindex,pepmassindex+1,pepmassindex+2,pepmassindex+3,pepmassindex+4,pepmassindex+5)
  fragmention<-sax$mz[-remove_index]
  index<-which((fragmention=="BEGIN IONS")|(fragmention=="\nBEGIN IONS"))
  cumdex<-rep(0,times=length(fragmention))
  cumdex[index]<-1
  cumdex[1]<-1
  mgfindex<-cumsum(cumdex)
  fragmention<-as.character(fragmention)
  fragmention<-strsplit(fragmention," ")
  mz<-sapply(fragmention,function(x) x[1])
  intensity<-sapply(fragmention,function(x) x[2])
  fragmention<-data.frame(mgfindex,mz,intensity)
  index<-which(fragmention$intensity=="IONS")
  fragmention<-fragmention[-index,]
  rm(sax)
  return (list("parention"=parention,"fragmention"=fragmention))
}# readmgf files and convert the mgf files to a list of mz values including charge and scanindex et al.
  

  fragmentprocessor<-function(spectralist){
    spectralist1<-as.data.table(spectralist)
    spectralist1$mz<-round(as.numeric(spectralist1$mz),4)
    spectralist1<-spectralist1[,.SD[1],by=.(mgfindex,mz)]
    setkey(spectralist1,mgfindex,intensity)
    spectralist1<-spectralist1[intensity>2000,]
    spectralist1<-spectralist1[,tail(.SD,300),by=mgfindex]
    mgfindex<-as.numeric(spectralist1$mgfindex)
    mz<-as.numeric(spectralist1$mz)
    spectralist1<-spectralist1[order(mgfindex,mz),]
    spectralist1mgf<-as.numeric(spectralist1$mgfindex)
    spectralist1mz<-as.numeric(spectralist1$mz)
    spectralist1intensity<-as.numeric(spectralist1$intensity)
    isotope<-deisotope(spectralist1mgf,spectralist1mz,spectralist1intensity)
    fragz<-calz(spectralist1mgf,spectralist1mz,spectralist1intensity)
    spectralisttemp<-data.frame(spectralist1,isotope,fragz)
    spectralisttemp$mz<-(spectralisttemp$mz-1.0078)*spectralisttemp$fragz+1.0078
    index<-which(spectralisttemp$isotope==0)
    spectralist1<-spectralisttemp[index,]
    spectralist1<-spectralist1[,1:6]
    spectralist1<-as.data.table(spectralist1)
    spectralist1<-spectralist1[,.SD[1],by=.(mgfindex,mz)]
    setkey(spectralist1,mgfindex,mz)
    mgfindex11<-as.numeric(spectralist1$mgfindex)
    mz1<-as.numeric(spectralist1$mz)
    intensity1<-as.numeric(spectralist1$intensity)
    spectrapos<-spectralist1[,.N,by=mgfindex]
    mgfleng<-spectrapos$N
    mgfpos<-cumsum(mgfleng)
    mgfpos<-c(0,mgfpos)
    leng<-length(mgfpos)
    mgfpos<-mgfpos[-leng]
    mgfpos<-as.integer(mgfpos)
    mgfindex<-spectrapos$mgfindex
    return (list("mgfindex11"=mgfindex11,"mz1"=mz1,"intensity1"=intensity1,"mgfpos"=mgfpos,"mgfleng"=mgfleng,"mgfindex"=mgfindex))
  }  #Generate fragmentation mz and intensity information for each scan

  
  Indexgen<-function(phossequence){
    Theoryspectra<-spectragen(phossequence) 
    Theoryspectra<-as.data.table(Theoryspectra)
    Theoryspectra<-Theoryspectra[,start:=mz-0.000001]
    Theoryspectra<-Theoryspectra[,end:=mz+0.000001]
    setkey(Theoryspectra,start,end)
    
    overlapmz<-seq(120,2000,0.02)
    leng<-length(overlapmz)
    seqindex<-seq(1:leng)
    indexmz<-data.frame(seqindex,overlapmz)
    indexmz<-as.data.table(indexmz)
    indexmz<-indexmz[,start:=overlapmz-0.01]
    indexmz<-indexmz[,end:=overlapmz+0.01]
    setkey(indexmz,start,end)
    
    overlapspectra<-foverlaps(Theoryspectra,indexmz)
    overlapspectra<-na.omit(overlapspectra)
    #rm(Theoryspectra)
    #rm(indexmz)
    gc()
    overlapspectra<-as.data.table(overlapspectra)
    setkey(overlapspectra,seqindex)
    overlapspectra<-overlapspectra[,constantoverlap:=1] 
    overlapspectra<-overlapspectra[,commonions:=sum(constantoverlap),by=seqindex]
    overlapspectra<-overlapspectra[commonions>0,] 
    overlapspectra<-overlapspectra[,c(1,2,5,13)]
    overlapspectra2<-overlapspectra[,c(1,3)]
    concisemz<-overlapspectra[,.SD[1],by=seqindex]
    concisemz<-as.data.table(concisemz)
    setkey(concisemz,seqindex,overlapmz)
    seqindex<-concisemz$seqindex
    seqTheomz<-concisemz$overlapmz
    #rm(overlapspectra)
    gc()
    overlapspectra2<-as.data.table(overlapspectra2)
    setkey(overlapspectra2,seqindex,index)
    seqindex11<-overlapspectra2$seqindex
    index11<-overlapspectra2$index
    return (list("seqindex"=seqindex,"seqTheomz"=seqTheomz,"seqindex11"=seqindex11,"index11"=index11))
  }   #Generate fragment ions index information for the following indexsearch
  
  

  sax<-fread("SAMPLE.mgf",sep="=", fill=TRUE,blank.lines.skip=TRUE,data.table=TRUE)
  sax<-data.frame(sax)
  spectra<-mgfread(sax)
  parention1<-spectra$parention
  spectralist1<-spectra$fragmention
  write.csv(spectralist1,file="spectralist1.csv")
  write.csv(parention1,file="parention1.csv")
  rm(spectra,spectralist1,parention1)

  parention1<-read.csv("parention1.csv")
  spectralist1<-fread("spectralist1.csv")
  spectralist<-fragmentprocessor(spectralist1)
  mgfindex11<- spectralist$mgfindex11
  mz1<-spectralist$mz1
  intensity1<-spectralist$intensity1
  mgfpos<-spectralist$mgfpos
  mgfleng<-spectralist$mgfleng
  mgfindex<-spectralist$mgfindex
  parention2<-data.frame(mgfindex,mgfpos,mgfleng)
  parention1<-data.table(parention1)
  setkey(parention1,mgfindex)
  parention2<-data.table(parention2)
  setkey(parention2,mgfindex)
  parention<-merge(parention1,parention2,by="mgfindex",all=TRUE)
  write.csv(parention,file="parention.csv")
  rm(parention)
  parention<-read.csv("parention.csv")
  parention<-as.data.table(parention)
  setkey(parention,mgfindex,scanindex)
  mgfindex1<-as.integer(parention$mgfindex)
  pepmz1<-parention$pepmass
  Rtime1<-parention$Rtime
  charge1<-as.integer(parention$charge)
  scanindex1<-as.integer(parention$scanindex)
  pepmass1<-(pepmz1-1.0078)*charge1+1.0078
  mgfpos<-as.integer(parention$mgfpos)
  mgfleng<-as.integer(parention$mgfleng)
  index<-which(is.na(mgfpos))
  mgfpos[index]<-length(mz1)-1
  mgfleng[index]<-0


  

  peptideID<-read.csv("Type0pep.csv")
  phossequence<-as.character(peptideID$Sequence)
  phosAcession<-as.character(peptideID$Protein_Name)
  phospepmass<-calmass( phossequence)
  phosindex<-seq(1:length(phossequence))
  
  
  
  sequenceindex<-Indexgen(phossequence)
  seqindex<-sequenceindex$seqindex
  seqTheomz<-sequenceindex$seqTheomz
  seqindex11<-sequenceindex$seqindex11
  index11<-sequenceindex$index11
  
  rm(spectralist,parention1,parention2,sequenceindex)
  gc()
  
  t1<-Sys.time()
  print(t1)
  
  score<-Indexsearch(mgfindex1,mgfpos,round(pepmz1,1),mgfindex11,mz1,seqindex,seqTheomz,seqindex11,index11)
  #The first and fast index search, for each scan,only those peptide sequences that have more than three matched ions are kept.
  #Only monolinked peptides are searched and only HCD spectra are searched.

  t2<-Sys.time()
  print(t2)
  gc()
  #write.csv(score,file="score-TP.csv")
  #score<-fread("score-TP.csv")
  
  
  
  
  NTspectra<-spectragen(phossequence)
  NTspectra<-as.data.table(NTspectra)
  setkey(NTspectra,index,mz)
  Theomzleng<-NTspectra[,.N,by=index]$N
  Theomzpos<-cumsum(Theomzleng)
  Theomzpos<-c(0,Theomzpos)
  leng<-length(Theomzpos)
  Theomzpos<-Theomzpos[-leng]
  Theomzpos<-as.integer(Theomzpos)
  Theomz1<-NTspectra$mz
  
  
  scoremgfindex<-as.integer(score$mgfindexexp)
  scorepepseqindex<-as.integer(score$pepseqindex)
  mz_seq_z<-getmzseqz(scoremgfindex,scorepepseqindex,pepmz1,charge1,mgfpos,mgfleng,phossequence,Theomzpos,Theomzleng)
  leng<-length(score$mgfindexexp)-1
  score<-score[1:leng,]
  
  
  
  score<-data.frame(score,mz_seq_z) 
  score$peptideseq1<-as.character(score$peptideseq1)
  phosmass1<-calmass(score$peptideseq1)
  pepmass<-(score$pepexpmz-1.0078)*score$pepexpz+1.0078
  score<-data.frame(score,pepmass,phosmass1)
  rm(phosmass1,mz_seq_z,scoremgfindex,scorepepseqindex,pepmass)
 
  t3<-Sys.time()
  print(t3)
  
  newscorephos<-Openmodsearch(score$mgfindexexp,score$mgfpos,score$mgfleng,score$pepexpmz,score$pepmass,mgfindex11,mz1,score$pepseqindex,score$peptideseq1,score$phosmass1,Theomz1,score$Theomzpos,score$Theomzleng)
  #The second and open modfication search, both HCD and CID information are employed. A list is generated.
  
  #write.csv(newscorephos,file="newscorephos.csv")
  #newscorephos<-read.csv("newscorephos.csv")
  
  
  t4<-Sys.time()
  print(t4)
  rm(score,spectralist1)
  gc()
  
  newscorephos<-as.data.table(newscorephos)

    
  leng<-length(newscorephos$mgfindexexp)
  pepexpmz<-rep(NA,times=leng)
  pepexpz<-rep(NA,times=leng)
  scanindexexp<-rep(NA,times=leng)
    
  for(i in 1:leng){
      index<-newscorephos$mgfindexexp[i]
      pepexpmz[i]<-pepmz1[index]
      pepexpz[i]<-charge1[index]
      scanindexexp[i]<-scanindex1[index]
      }
    
    mgfindexexp<-newscorephos$mgfindexexp
    mgfpos<-newscorephos$mgfpos
    pepmass<-(pepexpmz-1.0078)*pepexpz+1.0078
    peptide1seq<-as.character(newscorephos$phossequence)
    Theomasspep1<-calmass(peptide1seq) 
    scoreHCD<-newscorephos$pep1HCDscore
 
    sizeHCD<-newscorephos$sizeHCD
   pepseqindex<-newscorephos$pepseqindex
      
    t5<-Sys.time()
    print(t5)
    
    peptideID2<-read.csv("Type0pep.csv")
    nonphossequence<-as.character(peptideID2$Sequence)
    nonphosAcession<-as.character(peptideID2$Protein_Name)
    nonphosmass<-calmass( nonphossequence)
    
    newscore<-Research(pepmass,pepseqindex,peptide1seq,mgfindexexp,pepexpmz,Theomasspep1,scanindexexp,scoreHCD,sizeHCD,pepmz1,mgfpos,mgfindex11,mz1,nonphosmass,nonphossequence)
    #Crosslinking search. Each spectra is searched against a combination of spectra.
    #The precursor mass need to be matched.  
    
    t6<-Sys.time()
    print(t6)
    

    
    Alphamass<-calmass(as.character(newscore$Alphapep))
    Betamass<-calmass(as.character(newscore$Betapep))
    
    monoscore<-monoisosearch(as.integer(newscore$mgfindexexp),Alphamass,Betamass,as.integer(newscore$mgfpos),mgfindex11,mz1,as.character(newscore$Alphapep),as.character(newscore$Betapep))
    newscore<-data.frame(newscore,monoscore)
    
    
    
    
    
    Alphasize<-nchar(as.character(newscore$Alphapep))
    Betasize<-nchar(as.character(newscore$Betapep))
    
    
    newscore<-data.frame(newscore,Alphasize,Betasize)
    
    
    
    
    Alphascore<-1-ppois(newscore$mscoreAlpha,(newscore$expsizeHCD*newscore$Alphasize*2/2000))
    Alphascore<-(-7.2)*log10(Alphascore)
    Alphascore<-round(Alphascore,2)
    
    Betascore<-1-ppois(newscore$mscoreBeta,(newscore$expsizeHCD*newscore$Betasize*2/2000))
    Betascore<-(-7.2)*log10(Betascore)
    Betascore<-round(Betascore,2)
    
    
    newscore<-data.frame(newscore,Alphascore,Betascore)
    
    
    p3<-ggplot()+geom_point(data=newscore,aes(x=Alphascore,y=Betascore),colour="Red")+
      geom_text(data=newscore,aes(x=Alphascore,y=Betascore,label=scanindex,size=0.1))
    p4<-ggplot()+geom_point(data=newscore,aes(x=Alphascore,y=Betascore),colour="Red")+
      geom_text(data=newscore,aes(x=Alphascore,y=Betascore,label=phosIndicator,size=0.1))
    
    
    
    leng<-length(newscore$scanindex)
    peptideID$Sequence<-as.character(peptideID$Sequence)
    newscore$Alphapep<-as.character(newscore$Alphapep)
    newscore$Betapep<-as.character(newscore$Betapep)
    peptideID$Protein_Name<-as.character(peptideID$Protein_Name)
    
    index<-which(newscore$phosIndicator==0)
    newscore0<-newscore[index,]
    protein0<-getproteiname2(as.character(peptideID$Protein_Name),as.character(peptideID2$Protein_Name),newscore0$Alphaindex,newscore0$Betaindex)
    newscore0<-data.frame(protein0,newscore0)
    index<-which(newscore$phosIndicator==1)
    newscore1<-newscore[index,]
    protein1<-getproteiname(as.character(peptideID$Protein_Name),newscore1$Alphaindex,newscore1$Betaindex)
    newscore1<-data.frame(protein1,newscore1)

    newscore<-rbind(newscore0,newscore1)
    
    

    psm<-read.csv("type0psms.csv")
    psmscan<-as.integer(psm$FirstScan)
    leng<-length(newscore$scanindex)
    monoscan<-rep(1,times=leng)
    for(i in 1:leng)
    {
      tempscan<-newscore$scanindex[i]
      if(any(psmscan==tempscan)){monoscan[i]<-0}
        
    }
    index<-which(monoscan==1)
      newscore<-newscore[index,]
    
  write.csv(newscore,file="m-mDB.csv")
    
    
    
  */
