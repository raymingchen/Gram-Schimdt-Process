########### Non-self-referential Gram-Schimdt Process #############
### to fully understand the notations and ideas behind these codes, please refer to the msnuscript.
rm(list=ls()); 
##----------------
U_dot_k=function(U,k) ##U is a matrix with at least k+1 columns
{
U_k_plus1=U[,k+1];
U_1tok=U[,1:k];
u_dot_k=t(U_k_plus1)%*%U_1tok;
return(u_dot_k);
}
##-----------------
normMatrixfromV=function(V_kminus1)
{
n=dim(V_kminus1)[2];
normMatV=matrix(0,n,n);
diag(normMatV)=-1/(colSums(V_kminus1*V_kminus1));
return(normMatV);
}

##-----------------
V_k=function(b_plus_Delta_k,U,k)
{
Vk=U[,1:k]%*%t(b_plus_Delta_k);
return(Vk);
}

##------------------
B_plus_Delta_k=function(b_plus_Delta_k,U,k,Vk)
{
m=dim(b_plus_Delta_k)[1];
n=dim(b_plus_Delta_k)[2];
b_plus_Delta_k2=matrix(0,m+1,n+1);
b_plus_Delta_k2[1:m,1:n]=b_plus_Delta_k;
b_plus_Delta_k2[1,n+1]=0;
b_plus_Delta_k2[m+1,n+1]=1;
u_dot_k=U_dot_k(U,k);   
normMatV=normMatrixfromV(Vk);
a1=u_dot_k%*%t(b_plus_Delta_k);
a2=a1%*%normMatV;
a3=a2%*%b_plus_Delta_k;
b_plus_Delta_k2[m+1,1:n]=a3;
b_plus_Delta_k=b_plus_Delta_k2;
return(b_plus_Delta_k);
}

##-(main programs)---Non-self-referential Gram-Schimdt Process----
##U is a m-by-n matrix consisting column vectors, where n<=m
##kk is the number of generated orthogonal vectors, where k<=n.
##to get output orthogonal vectors Vk from the input vectors U 
NsrGSP_Vk=function(U,kk)
{
b_plus_Delta_k=matrix(c(1),1,1); ##the initial B^+_{Delta_1}=[1];
Vk=as.matrix(U[,1],default=0);  ##the initial Vk=U[,1];
 for(k in 1:(kk-1))
 {
 b_plus_Delta_k=B_plus_Delta_k(b_plus_Delta_k,U,k,Vk);
 Vk=V_k(b_plus_Delta_k,U,k+1);  ##indeed this is V_{k+1}
 }
return(Vk);
}

## to get the coefficients associated with input vectors U
NsrGSP_COEFFk=function(U,kk)
{
b_plus_Delta_k=matrix(c(1),1,1); ##the initial B^+_{Delta_1}=[1];
Vk=as.matrix(U[,1],default=0);  ##the initial Vk=U[,1];
 for(k in 1:(kk-1))
 {
 b_plus_Delta_k=B_plus_Delta_k(b_plus_Delta_k,U,k,Vk);
 Vk=V_k(b_plus_Delta_k,U,k+1);  ##indeed this is V_{k+1}
 }
return(b_plus_Delta_k);
}


######======= Experiment 1 (definite input U)
U=matrix(c(-2,5,1,0,4,7,-1,-1,2,-4,10,0,-3,2,1,2,5,1,7,-1,4,3,11,0,6),5,5);
kk=dim(U)[2];
Vk=NsrGSP_Vk(U,kk);
check=t(Vk)%*%Vk;
COEFFk=NsrGSP_COEFFk(U,kk);
print("InputU=");U;
print("Output V=NsrGSP(U)=");Vk;
print("check the orthogonality of Vk via t(Vk)%*%Vk=");t(Vk)%*%Vk;
print("the coefficients associated with U are COEFF=");COEFFk;
print("the correctness of the assigned coefficients is verified by U%*%t(COEFFk)==Vk??");
U%*%t(COEFFk)==Vk;
print("the specific coefficients C_i associated with U to generate v_i is COEFFk[i,]");
print("for example, C_4="); as.matrix(COEFFk[4,],default=0);
print("check v_4==Ut(C_4)??"); Vk[,4]==U%*%COEFFk[4,];


######======== Experiment 2 (randomly generated input U)
U=matrix(runif(15*9),15,9);
kk=dim(U)[2];
Vk=NsrGSP_Vk(U,kk);
check=t(Vk)%*%Vk;
COEFFk=NsrGSP_COEFFk(U,kk);
print("InputU=");U;
print("Output V=NsrGSP(U)=");Vk;
print("check the orthogonality of Vk via t(Vk)%*%Vk=");t(Vk)%*%Vk;
print("the coefficients associated with U are COEFF=");COEFFk;
print("the correctness of the assigned coefficients is verified by U%*%t(COEFFk)==Vk??");
U%*%t(COEFFk)==Vk;
print("the specific coefficients C_i associated with U to generate v_i is COEFFk[i,]");
print("for example, C_8="); COEFFk[8,];
print("check v_8==Ut(C_8)??"); Vk[,8]==U%*%COEFFk[8,];




