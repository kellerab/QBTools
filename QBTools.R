################################################################################
#
# Quantitative Biodiversity Functions Source Code; Spring 2015 Course, taught by
#     Dr. Jay Lennon
#
# Written by: Adrienne Keller
#
# Last update: 2015/06/03
#
# Notes: This file contains functions from alpha diversity week 1 exercise
#
################################################################################

require("vegan")||install.packages("vegan"); require("vegan")

#Observed species richness
S.obs<-function(x=""){
  rowSums(x > 0) * 1
}

#Good's average
C<- function(x=""){
  1- (sum(x == 1)/rowSums(x))
}

#Chao 1
S.chao1<-function(x = ""){
  S.obs(x) + (sum(x == 1)^2)/(2 * sum(x == 2))
}

#Chao 2
S.chao2<- function(site = "", SbyS = ""){
  SbyS = as.data.frame(SbyS)
  x = SbyS[site,]
  SbyS.pa<-(SbyS > 0) *1
  Q1 = sum(colSums(SbyS.pa) == 1)
  Q2 = sum(colSums(SbyS.pa) == 2)
  S.chao2 = S.obs(x) + (Q1^2)/(2 * Q2)
  return(S.chao2)
}

#Rank abundance curve
RAC<-function(x = ""){
  x = as.vector(x)
  x.ab = x[x>0]
  x.ab.ranked = x.ab[order(x.ab, decreasing =T)]
  return(x.ab.ranked)
}

#Simpson's evenness
SimpE<- function(x = ""){
  x = as.data.frame(x)
  D<-diversity(x,"inv")
  S<-S.obs(x)
  E<-(D)/S
  return(E)
}

#Smith and Wilson's evenness index
Evar<-function(x){
  x<-as.vector(x[x>0])
  1-(2/pi)*atan(var(log(x)))
}

#Shannon's diversity
H<-function(x = ""){
  H=0
  for (n_i in x){
    p = n_i/sum(x)
    H=H-p*log(p)
  }
  return(H)
}

#Simpson's diversity
D<-function(x=""){
  D=0
  N=sum(x)
  for (n_i in x){
    D=D+(n_i^2)/(N^2)
  }
  return(D)
}