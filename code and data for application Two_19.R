########### R codes for Application Two in my manuscript ############
### some extra programs or data should be downloaded from github shown 
### in the manuscript before running this program
rm(list=ls());
library(readr);  ## for function read.csv
library(pracma); 
library(ggplot2); ## for plotting histgram purpose
print("Please select your local file named Non_self-referential GSP.R");
source(file.choose()); ##Please choose the local file named Non_self-referential GSP.R
                       ##the functions to be called are NsrGSP_Vk(U,kk) and NsrGSP_COEFFk(U,kk)
print("Good! Non_self-referential GSP.R is selected");
ProjFun=function(vec1,vec2)##project vector 1 onto vector 2
{
num=vec1%*%vec2;
den=vec2%*%vec2;
return(num/den)
}


##### real data analysis for classifying 202 countries into 6 categories
##### from Year 2013 to 2023
print("Please select your local file named data");
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
normDt=as.matrix(sqrt(rowSums(Dt*Dt)),default=0);
normCOFt=sqrt(rowSums(COFt*COFt));
cos_t=(Dt%*%t(COFt))*(1/(normDt%*%normCOFt));
CosSim[[t]]=cos_t;
}

##----categorization of all the 202 ocuntries on a yearly basis---------

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

##--11 years, 6 representatives, categories assigned
Cat6=matrix(0,11,6);
for(t in 1:11)
{
cat_t=catYear2013to2023[[t]];
Cat6[t,]=cat_t[c(38,63,68,82,151,193)]
}

##----------Presentations of results -----------
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

##-------for comparison (without Non-self-referntial GSP)
cost=function(vec1,vec2) ##function defining cosine value of vector 1,2
{
norm1=sqrt(vec1%*%vec1);
norm2=sqrt(vec2%*%vec2);
cosine=(vec1%*%vec2)/(norm1*norm2);
return(cosine);
}

## calculate the cosine values
COS=list();    ## a list to collect all the cosine values for each year
ind=c(38,63,68,82,151,193); ## the country codes for the t representatives
for(t in 1:11)
{
Country_t=Country[t,];       
Yt=t(matrix(unlist(Country_t),6,202));
Rept=Yt[ind,];     ## the representative vectors (raw, without GSP)

COS_t=matrix(0,202,6);  ## the cosine values at Year t
   for(i in 1:202)
   {
    veci=Yt[i,];   ## the IDs for country i
      for(j in 1:6)
      {
      vecj=Yt[ind[j],];  ## the IDs for the representative country j
      COS_t[i,j]=cost(veci,vecj);
      }
   }
COS[[t]]=COS_t;
}

## identify the optimal categories
cat2013_raw=max.col(COS[[1]]);
cat2014_raw=max.col(COS[[2]]);
cat2015_raw=max.col(COS[[3]]);
cat2016_raw=max.col(COS[[4]]);
cat2017_raw=max.col(COS[[5]]);
cat2018_raw=max.col(COS[[6]]);
cat2019_raw=max.col(COS[[7]]);
cat2020_raw=max.col(COS[[8]]);
cat2021_raw=max.col(COS[[9]]);
cat2022_raw=max.col(COS[[10]]);
cat2023_raw=max.col(COS[[11]]);


## presentations of the comparative results

cat2013_raw
cat2014_raw
cat2015_raw
cat2016_raw
cat2017_raw
cat2018_raw
cat2019_raw
cat2020_raw
cat2021_raw
cat2022_raw
cat2023_raw

par(mfrow=c(4,3));
hist2013_r=hist(cat2013_raw);
hist2014_r=hist(cat2014_raw);
hist2015_r=hist(cat2015_raw);
hist2016_r=hist(cat2016_raw);
hist2017_r=hist(cat2017_raw);
hist2018_r=hist(cat2018_raw);
hist2019_r=hist(cat2019_raw);
hist2020_r=hist(cat2020_raw);
hist2021_r=hist(cat2021_raw);
hist2022_r=hist(cat2022_raw);
hist2023_r=hist(cat2023_raw);

inde=c(38,63,68,82,151,193);
cat2013_raw[inde];
cat2014_raw[inde];
cat2015_raw[inde];
cat2016_raw[inde];
cat2017_raw[inde];
cat2018_raw[inde];
cat2019_raw[inde];
cat2020_raw[inde];
cat2021_raw[inde];
cat2022_raw[inde];
cat2023_raw[inde];

###-------- for comparasion 2 ------------
ProjtoVspace=list();
for(t in 1:11)
{
Country_t=Country[t,];       
Yt=t(matrix(unlist(Country_t),6,202));
Prj_t=Yt%*%Rep[[t]];

COS2_t=matrix(0,202,6);  ## the cosine values at Year t
   for(i in 1:202)
   {
    veci=Prj_t[i,];   ## the IDs for country i
      for(j in 1:6)
      {
      vecj=Prj_t[ind[j],];  ## the IDs for the representative country j
      COS2_t[i,j]=cost(veci,vecj);
      }
   }
ProjtoVspace[[t]]=COS2_t;
}

cat2013_Vspace=max.col(ProjtoVspace[[1]]);
cat2014_Vspace=max.col(ProjtoVspace[[2]]);
cat2015_Vspace=max.col(ProjtoVspace[[3]]);
cat2016_Vspace=max.col(ProjtoVspace[[4]]);
cat2017_Vspace=max.col(ProjtoVspace[[5]]);
cat2018_Vspace=max.col(ProjtoVspace[[6]]);
cat2019_Vspace=max.col(ProjtoVspace[[7]]);
cat2020_Vspace=max.col(ProjtoVspace[[8]]);
cat2021_Vspace=max.col(ProjtoVspace[[9]]);
cat2022_Vspace=max.col(ProjtoVspace[[10]]);
cat2023_Vspace=max.col(ProjtoVspace[[11]]);

##----presentations of the results

cat2013_Vspace;
cat2014_Vspace;
cat2015_Vspace;
cat2016_Vspace;
cat2017_Vspace;
cat2018_Vspace;
cat2019_Vspace;
cat2020_Vspace;
cat2021_Vspace;
cat2022_Vspace;
cat2023_Vspace;

par(mfrow=c(4,3));
hist2013_V=hist(cat2013_Vspace);
hist2014_V=hist(cat2014_Vspace);
hist2015_V=hist(cat2015_Vspace);
hist2016_V=hist(cat2016_Vspace);
hist2017_V=hist(cat2017_Vspace);
hist2018_V=hist(cat2018_Vspace);
hist2019_V=hist(cat2019_Vspace);
hist2020_V=hist(cat2020_Vspace);
hist2021_V=hist(cat2021_Vspace);
hist2022_V=hist(cat2022_Vspace);
hist2023_V=hist(cat2023_Vspace);



