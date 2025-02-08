########### R codes for Application Two in my manuscript ############
rm(list=ls());
library(readr);  ## for function read.csv
library(pracma); 
library(ggplot2); ## for plotting histgram purpose
source(file.choose()); ##Please choose the local file named Non_self-referential GSP.R
                       ##the functions to be called are NsrGSP_Vk(U,kk) and NsrGSP_COEFFk(U,kk)
ProjFun=function(vec1,vec2)##project vector 1 onto vector 2
{
num=vec1%*%vec2;
den=vec2%*%vec2;
return(num/den)
}


##### real data analysis for classifying 202 countries into 6 categories
##### from Year 2013 to 2023
data=read_csv(file.choose()); ## choose the data from the file data
dataMatrix=as.matrix(data[1:1212,4:14]);
dim(dataMatrix);    ## the dimension is 1212-by-11, i.e., 202 countries 
##against 6 inidcators (ID1-ID6) for 11 years (from Year 2013 to 2023)


## China is a 6-by-11 matrix, i.e., ID1-ID6 for 11 years.
## The same statements also hold for other representative countries
China=dataMatrix[223:228,]; 
France=dataMatrix[373:378,];
GER=dataMatrix[403:408,];
India=dataMatrix[487:492,];
Russia=dataMatrix[901:906,];
USA=dataMatrix[1153:1158,];

##Rep is a list with length 11 (yrs) - each of its element is a
## 6-by-6 matrix (the row is associated with ID1-ID6, and the column is 
## associated with the 6 representative countries. 
Rep=list();       
for(t in 1:11)
{
Rep[[t]]=cbind(China[,t],France[,t],GER[,t],India[,t],Russia[,t],USA[,t]);
}

##Country is a 11-by-202 list: Country[[t,]] records ID1-ID6 for country j
## at time t.
Country=as.list(numeric(11*202));dim(Country)=c(11,202);
for(t in 1:11)
{
Year_t=dataMatrix[,t];
  for(j in 1:202)
  {
  id=(6*(j-1)+1):(6*j);
  Country[[t,j]]=Year_t[id];
  }
}


## DD is a list with length 11 - DD[[t]] is a 202-by-6 matrix which records
## all the projection values onto the 6 representative vectors at time t 
DD=list(); 
for(t in 1:11)
{
Country_t=Country[t,];
Yt=t(matrix(unlist(Country_t),6,202));
Ut=Rep[[t]];
Ut=Ut/colSums(Ut*Ut);
Dt=Yt%*%Ut;
DD[[t]]=Dt;
}

##----- find the Vk converted from the 6 representative vectors
VRep=list();
for(t in 1:11)
{
Ut=Rep[[t]];
V_t=NsrGSP_Vk(Ut,6);
VRep[[t]]=V_t;
}

##----- find the coefficients attached to the 6 representative vectors for
## the orthogonal features (vectors)
Coff=list();
for(t in 1:11)
{
Ut=Rep[[t]];
cof=NsrGSP_COEFFk(Ut,6);
Coff[[t]]=cof;
}

##----- compute the similarities between each country and the orthogonal
## features

CosSim=list();
for(t in 1:11)
{
Dt=DD[[t]];
COFt=Coff[[t]];
cos_t=Dt%*%COFt/colSums(COFt*COFt);
CosSim[[t]]=cos_t;
}

##-----------------------plotting for the results----------------------

cat2013=max.col(CosSim[[1]]);
cat2014=max.col(CosSim[[2]]);
cat2015=max.col(CosSim[[3]]);
cat2016=max.col(CosSim[[4]]);
cat2017=max.col(CosSim[[5]]);
cat2018=max.col(CosSim[[6]]);
cat2019=max.col(CosSim[[7]]);
cat2020=max.col(CosSim[[8]]);
cat2021=max.col(CosSim[[9]]);
cat2022=max.col(CosSim[[10]]);
cat2023=max.col(CosSim[[11]]);

catYear2013to2023=list();
for(t in 1:11)
{
catYear2013to2023[[t]]=max.col(CosSim[[t]]);
}

hist2013=hist(cat2013);
hist2014=hist(cat2014);
hist2015=hist(cat2015);
hist2016=hist(cat2016);
hist2017=hist(cat2017);
hist2018=hist(cat2018);
hist2019=hist(cat2019);
hist2020=hist(cat2020);
hist2021=hist(cat2021);
hist2022=hist(cat2022);
hist2023=hist(cat2023);

cat2013;
cat2014;
cat2015;
cat2016;
cat2017;
cat2018;
cat2019;
cat2020;
cat2021;
cat2022;
cat2023;

par(mfrow=c(4,3));
hist2013=hist(cat2013);
hist2014=hist(cat2014);
hist2015=hist(cat2015);
hist2016=hist(cat2016);
hist2017=hist(cat2017);
hist2018=hist(cat2018);
hist2019=hist(cat2019);
hist2020=hist(cat2020);
hist2021=hist(cat2021);
hist2022=hist(cat2022);
hist2023=hist(cat2023);

##--11 years, 6 representatives, categories assigned
Cat6=matrix(0,11,6);
for(t in 1:11)
{
cat_t=catYear2013to2023[[t]];
Cat6[t,]=cat_t[c(38,63,68,82,151,193)]
}


##---compute the cosine similarity between country and othogonal representatives





