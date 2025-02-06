########## R codes for Application 1 in the manuscript#########
rm(list=ls());
library(gtools);
library(combinat);
library(pracma);

sqNorm=function(A,B)
{
C=A-B;
return(sum(C*C));
}

### data analysis for application 1
U=matrix(c(-2,7,10,2,4,5,-1,0,5,3,1,-1,-3,1,11,0,2,2,7,0,4,-4,1,-1,6),5,5);
U=t(U);U;

U_ast=matrix(0,5,5);  ## unit U matrix
for(j in 1:5)
{
U_ast[,j]=U[,j]/sqrt(sum(U[,j]*U[,j]));
}
U_ast;

P=permn(1:5);  ##the permuted group/indexes of {1,2,3,4,5}
matbbU_ast=list();
matbbV_ast=list();
for(i in 1:120)
{
matbbU_ast[[i]]=U_ast[,unlist(P[i])];
matbbV_ast[[i]]=gramSchmidt(U_ast[,unlist(P[i])])$Q;
}

sqNormDiff=matrix(0,120,120);
for(i in 1:120)
{
V_i=unlist(matbbV_ast[[i]]);
  for(j in 1:120)
  {
  U_j=unlist(matbbU_ast[[j]]);
  sqNormDiff[i,j]=sqNorm(V_i,U_j)  
  }
}
rS=rowSums(sqNormDiff);
min=min(rS);
ind=which(rS==min);
V_ast=matbbV_ast[[ind]];
U_ast=matbbU_ast[[ind]];
V_ast;U_ast;








